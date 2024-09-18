!  jadeR_Fortran.f90 
!
!  FUNCTIONS:
!  jadeR_Fortran - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: jadeR_Fortran
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

program jadeR_Fortran
    implicit none

    ! Variables
    ! Body of jadeR_Fortran
    
        
    real*8,allocatable::X(:,:),B(:,:)
    integer m, X_size1, X_size2,i
        
    m = 2
    X_size1 = 3
    X_size2 = 4
    allocate(X(X_size1, X_size2),B(m,X_size1))
    X(1,1) = 1
    X(1,2) = 3
    X(1,3) = 2
    X(1,4) = 4
    X(2,1) = 9
    X(2,2) = 7
    X(2,3) = 8
    X(2,4) = 5
    X(3,1) = 9
    X(3,2) = 15
    X(3,3) = 18
    X(3,4) = 13
    ! print *, X(1,1),X(1,2),X(1,3)
    ! print *, X(2,1),X(2,2),X(2,3)
    call jadeR(X, m, X_size1, X_size2)


    
    
    

end program jadeR_Fortran
    
    
    
    
subroutine jadeR(X, m, X_size1, X_size2) 
    use blas95
    use f95_precision
    implicit none
    ! X is the matrix     m is the number of the modes   
    ! X_size1 is the rows , X_size2 is the cols

    integer i,j,i_,j_  ! 循环变量
    integer m, n, T, X_size1, X_size2
    real*8 X(X_size1, X_size2)
    real*8 sum_E
    real*8 , allocatable :: U(:,:), D(:)
    real*8 , allocatable :: Ds(:)
    integer, allocatable :: k(:)
    real*8 , allocatable :: C(:,:),E(:,:)
    integer*8 , allocatable :: PCs(:)  
    real*8 , allocatable :: B(:,:),scales(:),PC(:,:) 
    real*8 , allocatable :: inv_PC(:,:),EOF_(:,:)
    real*8 , allocatable :: X_(:,:),X_tran(:,:)    
    integer dimsymm,nbcm
    real*8 , allocatable :: CM(:,:),R(:,:),Qij(:,:),Xim(:,:),Xijm(:,:),R_vec_i(:,:),R_vec_j(:,:)
    integer, allocatable :: Uns(:),range_(:)
    real*8 , allocatable :: V(:,:),Diag(:,:)
    real*8 On,Off
    integer encore,sweep,updates,upds
    real*8 seuil,c_,s_,ton,toff,theta,Gain
    real*8 , allocatable :: g(:,:),gg(:,:),G_(:,:)
    integer, allocatable :: Ip(:),Iq(:),pair(:)
    real*8 , allocatable :: A(:,:)
    real*8 , allocatable ::b_(:)

        
    n = X_size1
    T = X_size2

    
    if (X_size1 < m) then
        write (*,*) "jade -> Do not ask more sources than sensors here!!!"
        return
    else
        write (*,"(' ',A,' ',I2.2,' ',A)") "jade -> Looking for", m , "sources"
    end if
    

! =======================================================================================================    
    ! Mean removal
    ! 对数据矩阵的行进行去平均操作 
    write (*,*) "jade -> Removing the mean value"    
    call juping(X, n, T)
    ! =======================================================================================================
    

    

! =======================================================================================================
    ! whitening & projection onto signal subspace
    ! 信号子空间的白化与投影
    write (*,*) "jade -> Whitening the data"
    
    ! An eigen basis for the sample covariance matrix
    ! 计算特征值
    ! U为特征向量 按列排列
    ! D为特征值与U的每一列对应   
    ! --------------------------------------------
    allocate(U(X_size1,X_size1))
    allocate(D(X_size1))
    call eig(matmul(X,transpose(X))/T, n, U, D)

    ! Sort by increasing variances
    ! 将特征值从小到大排列
    ! Ds从小到大排序后的特征值序列
    ! k为Ds中每个对应特征值在原始序列D的位置
    ! --------------------------------------------
    allocate(Ds(X_size1))
    allocate(k(X_size1))
    call sort_(Ds,k,D,X_size1)

    ! The m most significant princip. comp. by decreasing variance
    ! 选择最后m个最显著的特征值
    ! --------------------------------------------
    allocate(PCs(m))
    j = 1
    do i = n , n-m+1,-1
        PCs(j) = i
        j = j + 1
    end do

    allocate(C(n,n))
    allocate(E(n,1))
    C = matmul(X,transpose(X))/T
    do i = 1,n
        E(n-i+1,1) = D(i)
    end do
    sum_E = sum(E)
    E     = E *100.0 /sum_E
    
    ! ---  PCA  ----------------------------------------
    allocate(B(m,n))
    B = transpose(U(:,k(PCs)))   ! At this stage, B does the PCA on m components


    ! ---  Scaling  -------------------------------------
    allocate(scales(m))
    allocate(PC(T,m))
    allocate(inv_PC(m,m))
    allocate(EOF_(n,m))
    scales = SQRT(Ds(PCs))  ! The scales of the principal components .
    do i = 1,m 
        B(i,:) = 1.0/scales(i) * B(i,:)  ! Now, B does PCA followed by a rescaling = sphering
    end do    
    PC   = TRANSPOSE(MATMUL(B , X))
    call max_inv(MATMUL(TRANSPOSE(PC),PC),m,inv_PC)    
    EOF_ = MATMUL(X , MATMUL(PC , inv_PC))
    
    ! --- sphering -------------------------------------
    allocate(X_(m,T))
    X_    = MATMUL(B,X)
    
    ! --- release variable memory-----------------------
    deallocate(U,D,Ds,k,PCs,scales,inv_PC,C)
    ! =======================================================================================================


    

! =======================================================================================================
    ! Reshaping of the data, hoping to speed up things a little bit
    ! 改变数据矩阵形状，以加速
    write (*,*) "jade -> Estimating cumulant matrices"
    allocate(X_tran(T,m))
    X_tran = TRANSPOSE(X_)

    dimsymm    = (m*(m+1))/2  ! Dim. of the space of real symm matrices
    nbcm       = dimsymm      ! number of cumulant matrices

    allocate(CM(m,m*nbcm))
    allocate(R(m,m),R_vec_i(m,1),R_vec_j(m,1))
    allocate(Qij(m,m),Xim(T,1),Xijm(T,1),Uns(m),range_(m))
    R          = 0.0
    do i = 1,m
        R(i,i) = 1.0
    end do
    Uns        = 1.0
    do i = 1,m
        range_(i) = i
    end do
    
    
    do i = 1,m
        Xim(:,1)      = X_tran(:,i)
        Xijm          = Xim * Xim
        R_vec_i(:,1)  = R(:,i)  ! 若不单独另开数组，仅仅是引用某个2维数组的单列的话，会被系统认为是一维数组无法进行转置
        Qij           = matmul(transpose(Xijm(:,Uns) * X_tran) , X_tran) / T - R - 2*matmul(R_vec_i,transpose(R_vec_i))       
        CM(:,range_)  = Qij
        range_        = range_  + m 

        do j = 1,i-1
            Xijm(:,1)     = X_tran(:,j)
            Xijm          = Xim * Xijm
            R_vec_j(:,1)  = R(:,j)
            Qij           = sqrt(2.0) * (matmul(transpose(Xijm(:,Uns) * X_tran) , X_tran) / T - matmul(R_vec_i,transpose(R_vec_j)) - matmul(R_vec_j,transpose(R_vec_i)))   
            CM(:,range_)  = Qij
            range_        = range_  + m 
        end do
    end do
    ! Now we have nbcm = m(m+1)/2 cumulants matrices stored in a big m x m*nbcm array.  
    ! =======================================================================================================



    
! =======================================================================================================
    ! joint diagonalization of the cumulant matrices
    ! 累积量矩阵的联合对角化

    
    ! The dont-try-to-be-smart init
    ! --------------------------------------------
    allocate(V(m,m),Diag(m,1))
    do i = 1,m
        V(i,i) = 1.0  ! la rotation initiale
    end do
    Diag       = 0.0
    On         = 0.0
    do i = 1,m
        range_(i) = i
    end do

    do i = 1,nbcm
        do j = 1,m 
            Diag(j,1)   = CM(j,range_(j)) 
        end do
        On      = On + sum(matmul(transpose(Diag),Diag))
        range_  = range_  + m 
    end do
    Off   = sum(CM * CM) - On

    seuil	= 1.0e-6 / sqrt(dble(T))  ! A statistically scaled threshold on `small' angles
    encore  = 1
    sweep   = 0  ! sweep number
    updates = 0  ! Total number of rotations
    upds    = 0  ! Number of rotations in a given seep
    allocate(g(2,nbcm),gg(2,2),G_(2,2),pair(2))
    g       = 0.0
    gg      = 0.0
    G_      = 0.0
    c_      = 0.0
    s_      = 0.0
    ton     = 0.0
    toff    = 0.0
    theta   = 0.0
    Gain    = 0.0


    ! Joint diagonalization proper
    ! 联合对角化
    ! --------------------------------------------
    write (*,*) "jade -> Contrast optimization by joint diagonalization"

    do while (encore)
        encore = 0
        upds   = 0

        do i = 1,m-1
            do j = i+1,m
                allocate(Ip((m*nbcm-i)/2+1))
                allocate(Iq((m*nbcm-j)/2+1))

                j_ = 1
                do i_ = i,m*nbcm,m
                    Ip(j_) = i_
                    j_     = j_ + 1
                end do

                j_ = 1
                do i_ = j,m*nbcm,m
                    Iq(j_) = i_
                    j_     = j_ + 1
                end do

                ! computation of Givens angle
                g(1,:) = CM(i,Ip) - CM(j,Iq) 
                g(2,:) = CM(i,Iq) + CM(j,Ip)
                gg     = matmul(g,transpose(g))
                ton    = gg(1,1) - gg(2,2)
                toff   = gg(1,2) + gg(2,1) 
                theta  = 0.5 * atan2(toff,ton + sqrt(ton*ton+toff*toff))
                Gain   = (sqrt(ton*ton+toff*toff) - ton) / 4


                ! Givens update
                if (abs(theta) > seuil) then
                    encore        = 1
                    upds          = upds + 1
                    c_            = cos(theta)
                    s_            = sin(theta)
                    G_            = reshape([c_,s_,-s_,c_],[2,2])
                    pair          = [i,j]
                    V(:,pair)     = matmul(V(:,pair) , G_)
                    CM(pair,:)    = matmul(transpose(G_),CM(pair,:))
                    CM(:,[Ip,Iq]) = reshape([c_*CM(:,Ip)+s_*CM(:,Iq),-s_*CM(:,Ip)+c_*CM(:,Iq)],[m,m*nbcm]) 

                    On   = On  + Gain
                    Off  = Off - Gain

                end if
                deallocate(Ip,Iq)
            end do ! of the loop on j
        end do ! of the loop on i        
        write (*,"(' ',A,I3.3,A,I3.3,A)") "jade -> Sweep ", sweep," completed in ", upds, " rotations"
        sweep   = sweep + 1 
        updates = updates + upds
    end do ! of the while loop

    write (*,"(' ',A,I3.3,A)") "jade -> Total of ",updates , " Givens rotations"
    ! =======================================================================================================




! =======================================================================================================
    ! A separating matrix
    ! 分离矩阵
    B = matmul(transpose(V),B)
    ! =======================================================================================================
    

! =======================================================================================================
    ! Permute the rows of the separating matrix B to get the most energetic components first.
    ! Here the **signals** are normalized to unit variance.  Therefore, the sort is
    ! according to the norm of the columns of A = pinv(B)
    ! 首先改变分离矩阵B的行次序，以获得能量最高的分量。
    ! 这里，**信号**被归一化为单位方差。因此，排序是根据A=pinv（B）列的范数进行的
    write (*,*) "jade -> Sorting the components"
    allocate(A(n,m),Ds(m),k(m),D(m))
    call max_pinv(B,m,n,A)
    A  = A * A
    do i = 1,m
        D(i) = sum(A(:,i))
    end do  
    call sort_(Ds,k,D,m)
    B  = B(k,:)
    B  = B(m:1:-1,:)

    ! Signs are fixed by forcing the first column of B to have non-negative entries.
    ! 通过强制B的第一列包含非负数来固定标志。
    ! --------------------------------------------
    write (*,*) "jade -> Fixing the signs"
    allocate(b_(m))
    b_  = B(:,1)
    b_  = sign(1.0, sign(1.0,b_) + 0.1)  ! just a trick to deal with sign=0  具体的sign用法 参考百度
    
    do i = 1,m
        B(i,:) = B(i,:) * b_(i)
    end do   
    ! =======================================================================================================


    do i = 1,size(B,1)
        write(*,'(*(f10.4))')  B(i,:)
    end do  
    100 format(' ',A,' ',I2.2,' ',A)
end subroutine

