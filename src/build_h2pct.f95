!*********************************************************************!
!                        build the ct2 hamiltonian                    !
!   The coupling element between 2-particle and charge-transfer       !
!   states is given by                                                !
!       <k,v,s,v'|H|k,v'',s',v'''> =                                  !
!       te*fc(v|v'')*fc(v',v''')*kd(s,s')                             !
!      +th*fc(v|v''')*fc(v'|v'')*kd(s,-s')*exp(-2*pi*i*k*s/nmol)      !
!*********************************************************************!
subroutine build_h2pct(k)
    use commonvar
    implicit none
    
    integer, intent(in) :: k
    integer vib, sv, vibv, vibc, sa, viba, h1, h2, kd, &
            bring_inside_nxrange
    
    !choose the 2-particle basis element |k,vib,sv,vibv>
    do vib  = 0, vibmax
    do sv   = nlbnd, nubnd
    do vibv = 1, vibmax

        h1 = nx_2p( vib, sv, vibv ) ! get the basis index and make
        if ( h1 == empty ) cycle    ! sure it is inthe basis

        !choose the charge-transfer basis element |k,vibc,sa,viba>
        do vibc = 0, vibmax
        do sa   = nlbnd, nubnd
        do viba = 0, vibmax
            
            h2 = nx_ct( vibc, sa, viba ) ! get the basis index and make
            if ( h2 == empty ) cycle     ! sure it is in the basis

            !assign the coupling term to the hamiltonian.
            !we only have nearest-neighbor coupling
            !THIS IS EXCHANGE TYPE COUPLING
            if ( abs(sa) == 1 ) then
                h( h1, h2 ) = te*fc_ga(vibv,viba)*fc_cf(vibc,vib) &
                                *kd(sv,sa) +                      &
                              th*fc_gc(vibv,vibc)*fc_af(viba,vib) &
                                *kd(sv,bring_inside_nxrange(-sa)) &
                                *cdexp(2.d0*pi*img*k*sa/(1.d0*nmol))
                h( h2, h1 ) = dconjg(h(h1,h2))
            end if

        end do
        end do
        end do
    end do
    end do
    end do
    
end subroutine
