subroutine max_inv(x,x_size,inv_x)
    ! 对方阵进行求逆
    use lapack95
    use f95_precision  
    integer x_size,i
    real*8 x(x_size,x_size)
    real*8 inv_x(x_size,x_size)
    integer ipiv(x_size)
    
    call getrf(x,ipiv)
    call getri(x,ipiv)

    inv_x = x

    
end subroutine

