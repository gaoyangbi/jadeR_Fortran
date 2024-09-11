subroutine eig(x, x_size, U, D)
    ! U为特征向量 按列排列
    ! D为特征值与U的每一列对应   
    use lapack95
    use f95_precision  
    integer x_size
    real*8 x(x_size, x_size), U(x_size,x_size), D(x_size) 

    real*8 wr(x_size), wi(x_size), vr(x_size,x_size), vl(x_size,x_size)
    integer info,i

    call geev(x, wr, wi, vl, vr, info)

    U = vr
    D = wr

    write (*,*) '特征值为:'
    write (*,'(*(f10.4))') wr

    write (*,*) '特征向量为:'

    do i = 1,x_size
        write(*,'(*(f10.4))') vr(i,:)
    end do
    
end subroutine

