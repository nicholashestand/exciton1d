!*********************************************************************!
!                  index the 2-particke k states                      !
!A 2-particke k-state is defined in terms of the local states as      !
!                                                                     !
!    |k,v;s,v'> = SUM_n exp(i*2*pi*k*n/N)|n,v;n+s,v'>                 !
!                                                                     !
!where v is the number of vibrations in the shifted potential well    !
!of excited molecule n, and v' is the number of vibrations in the     !
!unshifted well of molecule n+s. All other molecules are in their     !
!ground electronic and vibrational states.                            !
!*********************************************************************!
subroutine index_2p()
    use commonvar
    implicit none

    integer vib, s, svib
    
    !allocate the indexing array and set as empty
    allocate( nx_2p (  0:vibmax, nlbnd:nubnd, 1:vibmax ) )
    nx_2p = empty
    
    !index the two-particle basis states
    !|k,vib;s,svib> -> nx_2p( vib, s, svib )
    do vib = 0, vibmax      !vibration on electronexcited molecule
    do s   = nlbnd, nubnd   !displacement from electronic excited
        if ( s == 0 ) cycle !displacement cannot be zero
    do svib = 1, vibmax     !vibration on ground state molecule
        if ( vib + svib > vibmax ) cycle  !truncate at vibmax
        kount = kount + 1
        nx_2p( vib, s, svib ) = kount
    end do
    end do
    end do
    
end subroutine
