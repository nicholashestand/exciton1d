!*********************************************************************!
!               write the exciton dispersion to a file                !
!    These are the lowest eigenvalues for each k                      !
!*********************************************************************!
subroutine dispersion(k)
    use commonvar
    implicit none

    integer, intent(in) :: k
    integer, parameter :: fno = 999
    
    if ( k == nlbnd ) then
        open ( unit = fno, file = trim(task_title)//'_disp.csv' )
        write( fno, * ) 'k,energy'
        write( fno, '(i4,",",f14.7)') k, eval(1)
        close( fno )
    else
        open ( unit = fno, file = trim(task_title)//'_disp.csv', position='append' )
        write( fno, '(i4,",",f14.7)') k, eval(1)
        close( fno )
    end if
    
end subroutine
