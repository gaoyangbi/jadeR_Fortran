program example
    implicit none

    ! Variables
    ! Body of jadeR_Fortran
    
        
    real*8,allocatable::X(:,:),B_TWS(:,:),E_TWS(:,:),C_TWS(:,:),S1_TWS(:,:),T1_TWS(:,:)
    integer m, X_size1, X_size2,i
        
    m = 2
    X_size1 = 3
    X_size2 = 4
    allocate(X(X_size1, X_size2),B_TWS(m,X_size1),E_TWS(X_size1,1),C_TWS(X_size1,X_size1))
    allocate(S1_TWS(X_size1,m),T1_TWS(X_size2,m))
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
    call jadeR(X, m, X_size1, X_size2,B_TWS,E_TWS,C_TWS,S1_TWS,T1_TWS)

    write(*,*) "B_TWS:"
    do i = 1,size(B_TWS,1)
        write(*,'(*(f10.4))')  B_TWS(i,:)
    end do  

    write(*,*) "E_TWS:"
    do i = 1,size(E_TWS,1)
        write(*,'(*(f10.4))')  E_TWS(i,:)
    end do  

    write(*,*) "C_TWS:"
    do i = 1,size(C_TWS,1)
        write(*,'(*(f10.4))')  C_TWS(i,:)
    end do  

    write(*,*) "S1_TWS:"
    do i = 1,size(S1_TWS,1)
        write(*,'(*(f10.4))')  S1_TWS(i,:)
    end do  

    write(*,*) "T1_TWS:"
    do i = 1,size(T1_TWS,1)
        write(*,'(*(f10.4))')  T1_TWS(i,:)
    end do  
  
    

end program example
    