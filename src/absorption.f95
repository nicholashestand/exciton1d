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
!*********************************************************************!
subroutine absorption()
    use commonvar
    implicit none

    integer vib, h1, state, nsteps, fno, point
    complex*16 osc(kount)
    real*8 lineshape, ab, photon_energy, transition_energy, step
    
    if (.not. one_state ) then
        print*, '>>One particle states are off. '//&
                'Will not calculate the absorption spectrum'
        return
    end if
    
    !calculate the oscillator strength
    osc = complex_zero
    !go over all states
    do state = 1, kount
        !go over all 1p basis states
        do vib = 0, vibmax
            h1 = nx_1p( vib )
            if ( h1 == empty ) cycle
            
            !assume parallel transition dipole moments
            osc(state) = osc(state) + h(h1,state)*fc_gf(0,vib)
        end do
        osc(state) = osc(state)*dconjg(osc(state))
    end do

    !calculate the absorption spectrum and write to file
    fno = 999
    open( unit = fno, file = trim(task_title)//'_ab.csv')
    write( fno, * ) 'energy,absorption'

    nsteps = floor((dabs(maxval(eval) - minval(eval)) + 8.d0*abs_lw)/(10/hw))   !10cm-1 resolution
    step = (maxval(eval) - minval(eval) + 8.d0*abs_lw)/(1.d0*nsteps)
    photon_energy = minval(eval)-4.d0*abs_lw
    do point = 1, nsteps
        photon_energy = photon_energy + step
        ab = 0.d0
    do state = 1, kount
        transition_energy = eval(state)
        !gaussian lineshape function
        lineshape = dexp(-(photon_energy - transition_energy)**2/     &
                         (2.d0*abs_lw**2))/                          &
                         dsqrt(2.d0*abs_lw**2*pi)
        ab = ab + osc(state)*lineshape
    end do
        write( fno, '(f14.7",",f14.7)' ) photon_energy, ab
    end do
    close( fno )
end subroutine
