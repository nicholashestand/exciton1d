!*********************************************************************!
!      the kronecker delta function                                   !
!    kd( n, m ) = 1 if n = m, otherwise it is zero                    !
!*********************************************************************!
integer function kd(n,m)
    integer, intent(in) :: n, m

    kd = 0
    if ( n == m ) kd = 1

    return
end function

