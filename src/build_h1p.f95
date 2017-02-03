!*********************************************************************!
!             build the 1-particle hamiltonian                        !
!   The diagonal energies are given by                                ! 
!   <k,v|H|k,v> = ES1 + v*hw + J(k)*FC(0,v)*FC(0,v)                   !
!   Where J(k) arises from the Coulomb coupling and depends on the    !
!   number of nearest neighbors:                                      !
!   For zero neighbors (i.e. a monomer) J(k) = 0                      !
!   For one neighbor (i.e. a dimer ) J(k) = JCoul*cos(k)              !
!   For two neighbors (nmol >2) J(k) = 2 JCoul*cos(k)                 !
!   and FC(0,v) is the vibrational overlap factor                     !
!   Below J(k) is calculated as                                       !
!                                                                     !
!   Off diagonal entries are given by                                 !
!   <k,v|H|k,v'> = J(k)*FC(0|v)*FC(0|v')                              !
!*********************************************************************!
subroutine build_h1p(k)
    use commonvar
    implicit none
    
    integer, intent(in) :: k
    integer vib1, vib2, h1, h2
    real*8 Jk 

    !choose the first basis element |k,vib1>
    do vib1 = 0, vibmax

        h1 = nx_1p( vib1 ) !get the basis index

        !Add the vibrational and monomer energy to the diagonal
        h( h1, h1 ) = vib1*1.d0 + ES1

        !choose the second basis element |k,vib2> 
        !and calculate J(k)*FC(0,vib1)*FC(0,vib2)
        do vib2 = 0, vibmax

            h2 = nx_1p( vib2 ) !get the basis index

            !calculate Jk
            if ( nmol == 1 ) then
                Jk = 0.d0
            else if ( nmol == 2 ) then
                Jk = JCoul * cos( 2*pi*k/nmol )
            else
                Jk = 2.d0 * JCoul * cos( 2*pi*k/nmol )
            end if

            !multiply by the volap factors and assign
            h( h1, h2 ) = h(h1,h2) + Jk * fc_gf(0,vib1) * fc_gf(0,vib2)
        end do
    end do

end subroutine
