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
    
        
    real*8,allocatable::X(:,:)
    integer m, X_size1, X_size2
        
    m = 2
    X_size1 = 3
    X_size2 = 4
    allocate(X(X_size1, X_size2))
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
    
    
    ! print *, X(1,1),X(1,2),X(1,3)
    ! print *, X(2,1),X(2,2),X(2,3)
end program jadeR_Fortran
    
    
    
    
subroutine jadeR(X, m, X_size1, X_size2) 
    implicit none
    ! X is the matrix     m is the number of the modes   
    ! X_size1 is the rows , X_size2 is the cols

    integer i,j  ! 循环变量
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
        
    n = X_size1
    T = X_size2

    
    if (X_size1 < m) then
        write (*,*) "jade -> Do not ask more sources than sensors here!!!"
        return
    else
        write (*,100) "jade -> Looking for", m , "sources"
    end if
    
    ! Mean removal
    ! 对数据矩阵的行进行去平均操作 
    ! =======================================================================================================    
    write (*,*) "jade -> Removing the mean value"    
    call juping(X, n, T)
    ! =======================================================================================================
    

    
    ! whitening & projection onto signal subspace
    ! 信号子空间的白化与投影
    ! =======================================================================================================
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
    deallocate(U,D,Ds,PCs,scales)
    
    ! =======================================================================================================

    
    ! Reshaping of the data, hoping to speed up things a little bit
    ! 改变数据矩阵形状，以加速
    ! =======================================================================================================
    allocate(X_tran(T,m))
    X_tran = TRANSPOSE(X_)
    ! =======================================================================================================
    do i = 1,T
        write(*,'(*(f10.4))') X_tran(i,:)
    end do    
    
    


    
    


100 format(' ',A,' ',I2.2,' ',A)
end subroutine

