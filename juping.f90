    subroutine juping(X, X_size1, X_size2)
            
    use blas95
    implicit none
    integer X_size1, X_size2
    real*8 X(X_size1, X_size2)
    real*8 mean_X(X_size1, 1)
    real*8 Y(X_size1, X_size2)
    real*8 l(X_size2, 1) 
    real*8 alpha

    alpha = 1.0 / X_size2
    l     = 1.0
    
    call gemm(X, l, mean_X, 'N', 'N', alpha)
    call gemm(mean_X, l, Y, 'N', 'T')

    X = X - Y
    end subroutine
    