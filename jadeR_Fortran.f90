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
    X_size1 = 2
    X_size2 = 3
    allocate(X(X_size1, X_size2))
    X(1,1) = 1
    X(1,2) = 2
    X(1,3) = 3
    X(2,1) = 4
    X(2,2) = 5
    X(2,3) = 6
    print *, X(1,1),X(1,2),X(1,3)
    print *, X(2,1),X(2,2),X(2,3)
    call jadeR(X, m, X_size1, X_size2)
    
    
    print *, X(1,1),X(1,2),X(1,3)
    print *, X(2,1),X(2,2),X(2,3)
end program jadeR_Fortran
    
    
    
    
subroutine jadeR(X, m, X_size1, X_size2) 
    implicit none
    ! X is the matrix     m is the number of the modes   
    ! X_size1 is the rows , X_size2 is the cols
    
    integer m, n, T, X_size1, X_size2
    real*8 X(X_size1, X_size2)
    real*8 , allocatable :: U(:,:), D(:)
    real*8 , allocatable :: Ds(:),k(:)

    ! real*8 a(3,3)
    ! a=reshape([1,1,0,0,0,2,0,0,-1],shape(a))
    ! a=transpose(a)
    
        
    n = X_size1
    T = X_size2

    
    if (X_size1 < m) then
        write (*,*) "jade -> Do not ask more sources than sensors here!!!"
        return
    else
        write (*,100) "jade -> Looking for", m , "sources"
    end if
    
    ! 对数据矩阵的行进行去平均操作 Mean removal
    ! =======================================
    write (*,*) "jade -> Removing the mean value"    
    call juping(X, n, T)


    ! whitening & projection onto signal subspace
    ! 信号子空间的白化与投影
    ! =========================================
    allocate(U(X_size1,X_size1))
    allocate(D(X_size1))
    call eig(matmul(X,transpose(X))/T, n, U, D)

    ! Sort by increasing variances
    ! 将特征值从小到大排列
    ! Ds从小到大排序后的特征值序列
    ! k为Ds中每个对应特征值在原始序列D的位置
    ! ============================================
    allocate(Ds(X_size1))
    allocate(k(X_size1))
    call sort_(Ds,k,D)

    print *, U
    print *, D


100 format(' ',A,' ',I2.2,' ',A)
end subroutine

