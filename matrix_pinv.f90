subroutine max_pinv(x,x_size1,x_size2,pinv_x)
    ! 对矩阵进行求广义逆
    use lapack95
    use f95_precision  
    integer x_size1,x_size2,i
    real*8 x(x_size1,x_size2),U(x_size1,x_size1),S(x_size1,x_size2),x0(x_size1,x_size2)
    real*8 V(x_size2,x_size2),SS(min(x_size1,x_size2),min(x_size1,x_size2))
    real*8 s1(min(x_size1,x_size2))
    real*8 pinv_x(x_size2,x_size1)

    ! write(*,*) "x is ..."
    ! do i = 1,size(x,1)
    !     write(*,'(*(f10.4))')  x(i,:)
    ! end do

    x0 = x
    call gesvd(x0,s1,U,V)  ! //此处返回的V是V的转置Vt  ! SVD会把X改变，注意！！

    ! write(*,*) "U is ..."
    ! do i = 1,size(U,1)
    !     write(*,'(*(f10.4))')  U(i,:)
    ! end do

    ! write(*,*) "S is ..."
    S  = 0.0
    forall (i = 1 : min(x_size1,x_size2)) S(i,i) = s1(i)
    ! do i = 1,size(S,1)
    !     write(*,'(*(f10.4))')  S(i,:)
    ! end do

    ! write(*,*) "Vt is ..."
    ! do i = 1,size(V,1)
    !     write(*,'(*(f10.4))')  V(i,:)
    ! end do

    ! ================================================
    ! write(*,*) "checking, x is ..."
    ! x = matmul(matmul(U,S),V)  ! // x = U*S*Vt
    ! ================================================
    
    ! do i = 1,size(x,1)
    !     write(*,'(*(f10.4))')  x(i,:)
    ! end do

    V = transpose(V)
    forall (i = 1 : min(x_size1,x_size2)) S(i,i) = 1.0 / S(i,i)
    if (x_size1 >= x_size2) then
        SS     = S(1:x_size2,1:x_size2)
        pinv_x = matmul( matmul(V,SS),transpose(U(:,1:x_size2)) )
    else
        SS     = S(1:x_size1,1:x_size1) 
        pinv_x = matmul( matmul(V(:,1:x_size1), SS), transpose(U) )  ! // pinv(X) = V*S*Ut
    end if

    ! write(*,*) "pinv(x) is ..."
    ! do i = 1,size(pinv_x,1)
    !     write(*,'(*(f10.4))') pinv_x(i,:)
    ! end do    
end subroutine

