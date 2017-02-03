!*********************************************************************!
!       build the charge-transfer/2-particle hamiltonian              !
!   The coupling element between 1-particle and charge-transfer       !
!   states is given by                                                !
!       <k,v|H|k,v',s,v''> = te*FC(0|v'')*FC(v'|v) +                  !
!                            th*exp(-2*pi*i*k*s/N)*FC(0|v')*FC(v''|v) !
!*********************************************************************!
subroutine build_h1pct(k)
    use commonvar
    implicit none
    
    integer, intent(in) :: k
    integer vib, vibc, sa, viba, h1, h2
    
    !choose the 1 particle basis element |k,vib>
    do vib = 0, vibmax

        h1 = nx_1p(vib) !get the basis index
        
        !choose the charge-transfer basis element |k,vibc,sa,viba>
        do vibc = 0, vibmax
        do sa   = nlbnd, nubnd
        do viba = 0, vibmax

            h2 = nx_ct( vibc, sa, viba ) !get the basis index and make
            if ( h2 == empty ) cycle     !sure it is in the basis set
            
            !assign the coupling term to the hamiltonian
            !we only have nearest-neighbor coupling
            if ( sa == -1 .or. sa == 1 ) then
                h( h1, h2 ) = te*fc_ga(0,viba)*fc_cf(vibc,vib) +    &
                          th*cdexp(-2.d0*pi*img*k*sa/(1.d0*nmol)) * &
                          fc_gc(0,vibc)*fc_af(viba,vib)
                h( h2, h1 ) = dconjg( h(h1,h2) )
            end if
                    
        end do
        end do
        end do
    end do
        
end subroutine
