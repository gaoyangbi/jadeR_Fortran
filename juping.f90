subroutine juping(X, X_size1, X_size2)

    ! 对输入数据矩阵进行去平均操作，根据jadeR.m文件
    ! 对每一行进行去平均
    ! X	= X - mean(X')' * ones(1,T); 
    ! T为数据矩阵的col列数
    use blas95
    use f95_precision
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
    