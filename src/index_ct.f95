!*********************************************************************!
!          index the 2-particle charge-transfer  k states             !
!A charge-transfer 2-particle k-state is defined in terms of local    !
!charge transfer states as                                            !
!                                                                     !
!    |k,v;s,v'> = SUM_n exp(i*2*pi*k*n/N)|n,v;n+s,v'> / sqrt(N)       !
!    where n is the cationic molecule with v vibrational quanta       !
!    in the shifted potential well and n+s is the anionic molecule    !
!    with v' vibrational quanta in its shifted well. All other        !
!    molecules are assumed to be in their ground states               !
!    NOTE: Here CT states are restricted to nearest neighbors         !
!    But in general this need not be the case.                        !
!*********************************************************************!
subroutine index_ct()
    use commonvar
    implicit none

    integer cvib, s, avib
    
    !allocate the 2-particle charge-transfer index and 
    !initialize as empty
    allocate( nx_ct (  0:vibmax, nlbnd:nubnd, 0:vibmax ) )
    nx_ct = empty
    
    !index the 2-particle charge-transfer states
    !|k,cvib;s,avib> -> nx_ct( cvib, s, avib )
    do cvib = 0, vibmax    !vibration on cation molecule
    do s    = nlbnd, nubnd !anion displacement from cation
        if ( s > 1 .or. s < -1 ) cycle !restrict to nearest neighbor
        if ( s == 0 ) cycle!displacement cant be zero
    do avib = 0, vibmax    !vibration on anion molecule
        if ( cvib + avib > vibmax ) cycle !truncate at vibmax
        kount = kount + 1
        nx_ct( cvib, s, avib ) = kount
    end do
    end do
    end do
        
end subroutine
