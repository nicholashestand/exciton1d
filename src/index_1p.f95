!*********************************************************************!
!                    index the 1 p k states                           !
!A 1-particle k-state is defined in terms of 1-particle local states  !
!as:                                                                  !
!                                                                     !
!   |k,v> = SUM_n exp(i*2*pi*k*n/N)|n,v> / sqrt(N)                    !
!                                                                     !
!   where the vibrations v are in the shifted potential well,         !  
!   molecule n is excited and all other molecules are in their        !
!   ground electronic and vibrational state                           !
!   In this program, we only work with one k-submatrix at a time      !
!   so only the number of vibratons needs to be indexed               !
!*********************************************************************!
subroutine index_1p()
    use commonvar
    implicit none

    integer vib
    
    !allocate the indexing array and initialize as empty
    allocate( nx_1p ( 0:vibmax ) )
    nx_1p = empty
    
    !index the 1-particle basis states
    !|k,vib> -> nx_1p( vib )
    do vib = 0, vibmax
        kount = kount + 1
        nx_1p( vib ) = kount
    end do
        
end subroutine
