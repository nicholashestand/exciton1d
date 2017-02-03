!*********************************************************************!
!             build the 1-particle 2-particle hamiltonian             !
!   The 1-particle 2-particle matrix elements are given by            !
!   <k,v|H|k,v',s,v''>=                                               !
!                   exp(i*2*pi*k*s/N)*JCoul*FC(v''|v)*FC(0|v')        !
!
!   These are like the exchange type couplings for the 2-particle     !
!   Hamiltonian.                                                      !
!*********************************************************************!
subroutine build_h1p2p(k)
    use commonvar
    implicit none
    
    integer, intent(in) :: k
    integer vib1, vib2, s2, vibv2, h1, h2, ds
    complex*16 cpl_sum, modulate
    
    !choose the 1-particle basis set |k,vib1>
    do vib1 = 0, vibmax
        
        h1 = nx_1p( vib1 )   ! get the basis index and make sure it 
        if (h1==empty) cycle !is in the basis set

        !choose the 2-particle basis set |k,vib2;s2,vibv2>
        do vib2  = 0, vibmax
        do s2    = nlbnd, nubnd
        do vibv2 = 1, vibmax

            h2 = nx_2p( vib2, s2, vibv2 ) !get the basis index and make
            if ( h2 == empty ) cycle      !sure it is in the basis set

            !the exciton moves -s2
            ds = -s2

            !EXCHANGE TYPE COUPLING
            !only include nearest neighbor coupling
            if ( ds == -1 .or. ds == 1 ) then
                modulate = cdexp( -2.d0*pi*img*k*ds/(1.d0*nmol))
                h( h1, h2 ) = modulate * JCoul *                  &
                              fc_gf(vibv2,vib1)*fc_gf(0,vib2)
            end if
                          
            !Also set the Hermitian Conjugate here
            h( h2, h1 ) = dconjg( h( h1, h2) )
        end do
        end do
        end do
    end do
end subroutine
