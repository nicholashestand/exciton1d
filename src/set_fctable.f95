!********************************************************************!
!     initialize the vibrational overlap tables                      !
!********************************************************************!
subroutine set_fctable()
    use commonvar
    implicit none

    integer vibg, vibn, vibc, viba
    real*8, external :: volap

    !allocate space for the vibrational overlap tables
    allocate( fc_gf ( 0:vibmax, 0:vibmax ) )
    allocate( fc_gc ( 0:vibmax, 0:vibmax ) )
    allocate( fc_ga ( 0:vibmax, 0:vibmax ) )
    allocate( fc_cf ( 0:vibmax, 0:vibmax ) )
    allocate( fc_af ( 0:vibmax, 0:vibmax ) )
    
    !Generate the vibrational overlap tables. Ground state potential
    !well minimum is the reference, all others are shifted by
    !lambda

    ! ground to frenkel
    do vibg = 0,vibmax
    do vibn = 0,vibmax
        fc_gf(vibg, vibn) = volap( 0.d0, vibg, lambda_n, vibn )
    enddo
    enddo

    !ground to cation
    do vibg = 0,vibmax
    do vibc = 0,vibmax
        fc_gc(vibg, vibc) = volap( 0.d0, vibg, lambda_c, vibc )
    enddo
    enddo
        
    !ground to anion
    do vibg = 0,vibmax
    do viba = 0,vibmax
        fc_ga(vibg, viba) = volap( 0.d0, vibg, lambda_a, viba )
    enddo
    enddo

    !cation to frenkel
    do vibn = 0,vibmax
    do vibc = 0,vibmax
        fc_cf(vibc, vibn) = volap( lambda_c, vibc, lambda_n, vibn )
    enddo
    enddo
    
    !anion to frenkel
    do vibn = 0,vibmax
    do viba = 0,vibmax
        fc_af(viba, vibn) = volap( lambda_a, viba, lambda_n, vibn )
    enddo
    enddo
    
end subroutine
!********************************************************************!
!    Return the vibrational overlap for displaced harmonic           !
!    oscillators. The lambda are proportional to the well minimum    !
!    and when squared are equivalent to the Huang-Rhys parameter.    !
!    The vibrational overlap factors are given by the formula        !
!                                                                    !
!    <m|n>=sqrt(m!n!)exp(-lambda^2/2) *                              !
!          SUM_l^(min(m,n))                                          !
!                (-1)^(m-l)/[(m-l)!l!(n-l)!] *                       !
!                lambda^(m+n-2l)                                     !
!                                                                    !
!    which can be derived from the recursion relations found         !
!    on page 167 of MODERN OPTICAL SPECTROSCOPY                      !
!    here lambda is proportional to the equilibrium displacement     !
!********************************************************************!
real*8 function volap( lambda1, vib1, lambda2, vib2 )
    implicit none

    integer, intent (in) :: vib1, vib2
    real*8, intent(in) :: lambda1, lambda2
    integer k
    integer, external :: factorial
    real*8 lambda
    
    !calculate the displacement between the two potential wells
    lambda = lambda2 - lambda1

    volap = 0.d0
    !calculate the vibrational overlap
    !first calculate the summation
    do k = 0, min( vib1, vib2 )
        volap = volap+(-1.d0)**(vib2-k)/        &
                (factorial(vib1-k)*factorial(k)*&
                 factorial(vib2-k))*            &
                 lambda**(vib1+vib2-2*k)
    end do
    volap = volap*dsqrt(1.d0*factorial(vib1)*   &
                             factorial(vib2))*  &
                             dexp(-1.d0*        &
                             lambda**2/2.d0)

end function
!********************************************************************!
!    calcuate factorial                                              !
!********************************************************************!
integer function factorial( n )
    implicit none

    integer, intent(in) :: n
    integer i

    if ( n < 0 ) then
        print*, 'Factorial not calculatable for ', n, ' aborting.'
        stop
    end if

    factorial = 1
    do i = 1, n
        factorial = factorial*i
    end do

end function
