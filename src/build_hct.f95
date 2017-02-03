!*********************************************************************!
!                      build the ct hamiltonian                       !
!    Since we only have nearest neighbor charge-transfer states,      !
!    there is no charge-transfer/charge-transfer coupling. All        !
!    we do here is set the charge-transfer energy                     !
!*********************************************************************!
subroutine build_hct(k)
    use commonvar
    implicit none
    
    integer, intent(in) :: k
    integer vibc, sa, viba, h1

    !choose the charge-transfer state |k,vibc,sa,viba>
    do vibc = 0, vibmax
    do sa = nlbnd, nubnd
    do viba = 0, vibmax
        
        h1 = nx_ct( vibc, sa, viba ) !get the basis index and check
        if ( h1 == empty ) cycle     !that it is not empty

        !assign the energy 
        h( h1, h1 ) = ECT + (vibc + viba)*1.d0

        !For systems where displacements can be greater than one
        !this subroutine needs to be extended to account for 
        !charge transfer between these states.

    end do
    end do
    end do
    
end subroutine
