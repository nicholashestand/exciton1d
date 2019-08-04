!*********************************************************************!
!        Calculate the absorption spectrum and write to file          !
!                                                                     !
!   The transition dipole moment is given by:                         !
!                                                                     !
!   <G|u|W_i> = sum_(k,v) c_(k,v) <0|v> if k=0, 0 otherwise           !
!   where <0|v> is the vibrational overlap factor and W_i is the ith  !
!   eigenstate, expanded in terms of local basis states as            !
!                                                                     !
!   W_i = sum_(n,v) c_(n,v)|n,v> + two particle + ct                  !
!                                                                     !
!   ONLY ONE PARTICLE STATES ABSORB!!!                                !
!                                                                     !
!   Absorption to state i is proportional to the transition dipole    !
!   moment squared: |<G|u|W_i>|^2                                     !
!                                                                     !
!   The absorption spectrum is given by:                              !
!       A(E) = sum_(i) |<G|u|W_i>|^2 Gamma(E-E_i, abs_lw)             !
!   Where Gamma is a normal distribution with mean E-E_i and          !
!   standard deviation abs_lw. E_i is the energy of the ith eigenstate!
!                                                                     !
!   The absorption moments are also calculated:                       !
!   The first moment is:                                              !
!       <E> = sum_(i) |<G|u|W_i>|^2*E_(i) / sum_(i) |<G|u|W_i>|^2     !
!   The second moment is:                                             !
!       <E^2> = sum_(i) |<G|u|W_i>|^2*E_(i)^2 / sum_(i) |<G|u|W_i>|^2 !
!   The third moment is:                                              !
!       <E^3> = sum_(i) |<G|u|W_i>|^2*E_(i)^3 / sum_(i) |<G|u|W_i>|^2 !
!                                                                     !
!   The central absorption moments are also calculated:               !
!   The first central moment is:                                      !
!       <E> = sum_(i) |<G|u|W_i>|^2*(E_(i)-<E>) /                     !
!                                           sum_(i) |<G|u|W_i>|^2     !
!   The second central moment is:                                     !
!       <E> = sum_(i) |<G|u|W_i>|^2*(E_(i)-<E>)^2 /                   !
!                                           sum_(i) |<G|u|W_i>|^2     !
!   The third central moment is:                                      !
!       <E> = sum_(i) |<G|u|W_i>|^2*(E_(i)-<E>)^3 /                   !
!                                           sum_(i) |<G|u|W_i>|^2     !
!*********************************************************************!
subroutine absorption()
    use commonvar
    implicit none

    integer vib, h1, state, nsteps, fno, point, m
    complex*16 osc(kount)
    real*8 lineshape, ab, photon_energy, transition_energy, step, &
           moment(3), cmoment(3)
    
    if (.not. one_state ) then
        print*, '>>One particle states are off. '//&
                'Will not calculate the absorption spectrum'
        return
    end if
    
    ! calculate the oscillator strength
    osc = complex_zero
    ! go over all states
    do state = 1, kount
        ! go over all 1p basis states
        do vib = 0, vibmax
            h1 = nx_1p( vib )
            if ( h1 == empty ) cycle
            
            ! assume parallel transition dipole moments
            osc(state) = osc(state) + h(h1,state)*fc_gf(0,vib)
        end do
        osc(state) = osc(state)*dconjg(osc(state))
    end do

    ! calculate the absorption spectrum and write to file
    fno = 999
    open( unit = fno, file = trim(task_title)//'_ab.csv')
    !write( fno, * ) 'energy,absorption'

    ! stet number of spectral points to be evaluated
    nsteps = floor((dabs(maxval(eval) - minval(eval)) &
                + 64.d0*abs_lw)/(10/hw))   !10cm-1 resolution

    ! restrict to a 10000 cm window, in case the ectinf is set very high
    nsteps = min( nsteps, 10000 )
    step = (min(maxval(eval),minval(eval) + 10000/hw) - minval(eval) &
                + 32.d0*abs_lw)/(1.d0*nsteps)
    photon_energy = minval(eval)-16.d0*abs_lw
    do point = 1, nsteps
        photon_energy = photon_energy + step
        ab = 0.d0
    do state = 1, kount
        transition_energy = eval(state)
        ! gaussian lineshape function
        !lineshape = dexp(-(photon_energy - transition_energy)**2/     &
        !                  (2.d0*abs_lw**2))/                          &
        !                  dsqrt(2.d0*abs_lw**2*pi)
        ! lorentzian lineshape function
        lineshape = abs_lw/(1)/                                         &
                    ((photon_energy-transition_energy)**2 + abs_lw**2 )
        ab = ab + osc(state)*lineshape
        ! lorentzian lineshape function
    end do
        write( fno, '(f14.7" ",f14.7)' ) photon_energy*hw, ab
    end do
    close( fno )

    ! calculate the moments of the absorption spectrum
    moment = 0.d0
    do m = 1, 3
        do state = 1, kount
            moment(m) = moment(m) + osc(state)*eval(state)**m
        end do
    end do
    ! divide by the sum of oscillator strengths
    moment = moment / sum(osc)

    ! calculate the central moments of the absorption spectrum
    cmoment = 0.d0
    do m = 1, 3
        do state = 1, kount
            cmoment(m) = cmoment(m) + osc(state)*   &
                         (eval(state) - moment(1))**m
        end do
    end do
    cmoment = cmoment / sum(osc)

    ! write the moments to a file
    open( unit = fno, file = trim(task_title)//'_mom.csv' )
    write(fno, '(a)') 'Moments of the absorption spectrum'
    write(fno, '(a)') 'Moment Number, Moment Value'
    do m = 1, 3
        write( fno, '(i2,",",f14.6)' ) m, moment(m)
    end do
    write(fno, '(a)') 'Central moments of the absorption spectrum'
    write(fno, '(a)') 'Moment Number, Moment Value'
    do m = 1, 3
        write( fno, '(i2,",",f14.6)' ) m, cmoment(m)
    end do

end subroutine
