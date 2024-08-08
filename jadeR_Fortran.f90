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
    X_size2 = 3
    allocate(X(X_size1, X_size2))
    X = 3
    call jadeR(X, m, X_size1, X_size2)
    
    
    print *, X
    end program jadeR_Fortran
    
    subroutine jadeR(X, m, X_size1, X_size2)
!  X is tha matrix     m is the number of the modes    
    integer m, n, T, X_size1, X_size2
    real*8 X(X_size1, X_size2)
        
    n = X_size1
    T = X_size2

    
    if (X_size1 < m) then
        write (*,*) "jade -> Do not ask more sources than sensors here!!!"
        return
    else
        write (*,100) "jade -> Looking for", m , "sources"
    end if
    
    write (*,*) "jade -> Removing the mean value"
    
    
100 format(' ',A,' ',I2.2,' ',A)
    end subroutine

