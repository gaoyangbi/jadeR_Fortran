subroutine sort_(Ds,k,D,x_size)
    ! 将向量D中的元素从小到大排列
    ! Ds是从小到大排序后的特征值序列
    ! k为Ds中每个对应特征值在原始序列D的位置
    integer x_size,i,j
    real*8 Ds(x_size), D(x_size) 
    integer k(x_size), index
    real*8 string

    Ds = D
    i  = 1
    do while(.true.)
        do j = 1,x_size
            if (Ds(i) >= Ds(j)) then
                index = j
                ! print *, index
            end if
        end do 

        string = Ds(i)
        Ds(i) = Ds(index)
        Ds(index) = string

        if (i == index) then 
            i = i + 1
        end if

        if (i == x_size) then
            exit
        end if
        
    end do

    do i = 1,x_size
        do j = 1,x_size
            if (abs(Ds(i)-D(j)) <= 1e-10) then
                k(i) = j
            end if
        end do 
    end do

    write (*,*) '原始特征值为:'
    write (*,'(*(f10.4))') D

    write (*,*) '从小到大排列后的特征值为:'
    write (*,'(*(f10.4))') Ds

    write (*,*) '排列后的特征值在原始特征值序列的位置为:'
    write (*,'(*(i10.3))') k
    
end subroutine

