!*********************************************************************!
!             build the 2-particle hamiltonian                        !
!There are two types of coupling possible:                            !
!                                                                     !
!   Linker Coupling:                                                  !
!   Here the ground state vibrations must be the same and the exciton !
!   moves s-s' molecules. This is just like coupling 1P 1P states     !
!   except there is now a vibrational excitation as well              !
!   <k,v,s,v'|H|k,v'',s',v'> =                                        !
!               exp(-i*2*pi*k*(s'-s)/N)*JCoul*FC(0|v)*FC(0|v'')       !
!   where it must be the case that abs(s'-s) = 1 (the exciton moves   !
!   one unit)                                                         !
!                                                                     !
!   Exchange Coupling:                                                !
!   Here, the exciton moves to the vibrationally excited molecule     !
!   This is different than the 1-particle, 1-particle coupling        !
!   <k,v,s,v'|H|k,v'',-s,v'''> =                                      ! 
!               exp(-i*2*pi*k*s/N)*JCoul*FC(v'''|v)FC(v'|v'')         !
!*********************************************************************!
subroutine build_h2p(k)
    use commonvar
    implicit none
    
    integer, intent(in) :: k
    integer vib1, s1, vibv1, vib2, s2, vibv2, h1, h2, ds
    complex*16 cpl_sum, modulate
    
    !choose the first basis element |k,vib1,s1,vibv1>
    do vib1  = 0, vibmax
    do s1    = nlbnd, nubnd
    do vibv1 = 1, vibmax

        h1 = nx_2p( vib1, s1, vibv1 )  !get the basis index
        if ( h1 == empty ) cycle       !and check to make sure it is 
                                       !in the basis set

        !Add the vibrational and monomer energy to the diagonal
        h( h1, h1 ) = ( vib1 + vibv1 ) * 1.d0 + ES1

        !choose the second basis element |k,vib2,s2,vibv2>
        do vib2  = 0, vibmax
        do s2    = nlbnd, nubnd
        do vibv2 = 1, vibmax

            h2 = nx_2p( vib2, s2, vibv2 ) !get the basis index
            if ( h2 == empty ) cycle          !and check to make sure
                                              !it is in the basis set

            !calculate the coupling term
            cpl_sum = complex_zero !initialize to zero

            !LINKER TYPE COUPLING
            if ( vibv1 == vibv2 ) then
                !calculate the distance of exciton transfer
                ds = s1 - s2

                !bring inside the aggregate range, if outside
                if ( ds < nlbnd ) ds = ds + nmol
                if ( ds > nubnd ) ds = ds - nmol

                !add the coupling to the coupling sum... There
                !is only nearest neighbor couping and the phase
                !depends on whether the exciton is going "left"
                !or "right"
                if ( ds == 1 .or. ds == -1 ) then
                    modulate = cdexp( -2.d0*pi*img*k*ds/(1.d0*nmol))
                    cpl_sum = cpl_sum + modulate * JCoul *          &
                                    fc_gf(0,vib1)*fc_gf(0,vib2)
                end if      
            end if

            !EXCHANGE TYPE COUPLING
            if ( s1 == -s2 ) then

                !the distance of exciton transfer is s1
                ds = s1

                !add the coupling to the coupling sum... There
                !is only nearest neighbor coupling and the phase
                !depends on whether the exciton is going "left" or
                !"right"
                if ( s1 == -1 .or. s1 == 1 ) then
                    modulate = cdexp( -2.d0*pi*img*k*ds/(1.d0*nmol))
                    cpl_sum = cpl_sum + modulate * JCoul *          &
                                    fc_gf(vibv2,vib1)*fc_gf(vibv1,vib2)
                end if
            end if   
            
            !EXCHANGE TYPE II - SELF EXCHANGE FOR EVEN LATTICES
            !(this will really only affect the dimer case)
            if ( s1 == s2 ) then
                
                !the distance of exciton transfer is s1
                ds = s1

                !This will only happen for even lattices when s1
                !is equal to nmol/2
                if ( mod( nmol, 2 ) == 0 .and. abs(s1) == nmol/2 ) then

                !add the coupling to the coupling sum... There
                !is only nearest neighbor coupling and the phase
                !depends on whether the exciton is going "left" or
                !"right"
                if ( s1 == -1 .or. s1 == 1 ) then
                    modulate = cdexp( -2.d0*pi*img*k*ds/(1.d0*nmol))
                    cpl_sum = cpl_sum + modulate * JCoul *          &
                                    fc_gf(vibv2,vib1)*fc_gf(vibv1,vib2)
                end if

                end if
            end if

            !add the coupling to the Hamiltonian
            h( h1, h2 ) = h( h1, h2 ) + cpl_sum
        end do
        end do
        end do
    end do
    end do
    end do
        
end subroutine
