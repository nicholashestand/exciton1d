!*********************************************************************!
!                      build the ct hamiltonian                       !
!    The diagonal energies are:                                       !
!    <k,v;s,v'|H|k,v;s,v'> = ECT(s) + (v+v')*hw                       !
!    where ECT(s) = [ECTInf*(s-1) - ECT]/s                            !
!    ECTInf is the energy of a CT state separated to infinity         !
!                                                                     !
!    The off diagonal matrix elements are:                            !
!    <k,v;s,v'|H|k,v'';s',v'''> =                                     !
!           te*FC(0|v')*FC(0|v''')*kd(v|v'') +                        !
!           th*exp(-i*2*pi*k*(s'-s))*FC(0|v)*FC(0|v'')*kd(v'|v''')    !
!    for |s-s'|=1 and zero otherwise. Here kd is the kronecker delta  !
!*********************************************************************!
subroutine build_hct(k)
    use commonvar
    implicit none
    
    integer, intent(in) :: k
    integer vibc, sa, viba, h1, vibc2, sa2, viba2, h2, s, &
            bring_inside_nxrange, kd

    ! get the diagonal energies
    ! choose the charge-transfer state |k,vibc,sa,viba>
    do vibc = 0, vibmax
    do sa = nlbnd, nubnd
    do viba = 0, vibmax
        
        h1 = nx_ct( vibc, sa, viba ) !get the basis index and check
        if ( h1 == empty ) cycle     !that it is not empty

        ! assign the energy 
        s = abs(sa)
        h( h1, h1 ) = (ECTInf*(s-1) + ECT)/(s*1.d0) + (vibc + viba)*1.d0

        ! choose the second basis element k, vibc2, sa2, viba2
        do vibc2 = 0, vibmax
        do sa2 = nlbnd, nubnd
        do viba2 = 0, vibmax
        
            h2 = nx_ct( vibc2, sa2, viba2 ) ! get the basis index and
            if ( h2 == empty ) cycle        ! check that not empty
            ! calculate the electron/hole displacement
            s = bring_inside_nxrange(sa2-sa)
            if ( abs(s) .ne. 1) cycle  ! only nearest neighbor
                                       ! charge transfer allowed

            ! calculate the matrix element
            h( h1, h2 ) = te * fc_ga(0,viba)*fc_ga(0,viba2)*kd(vibc, vibc2) + &
                          th * fc_gc(0,vibc)*fc_gc(0,vibc2)*kd(viba, viba2) * &
                          cdexp(-2.d0*pi*img*k*s/(1.d0*nmol))
            h( h2, h1 ) = dconjg( h( h1, h2 ) )

        end do
        end do
        end do

    end do
    end do
    end do

end subroutine
!*********************************************************************!
!       bring index inside the periodic range                         !
!*********************************************************************!
integer function bring_inside_nxrange(s)
    use commonvar
    implicit none
    integer, intent(in) :: s

    bring_inside_nxrange = s
    if ( s > nubnd ) bring_inside_nxrange = s - nmol
    if ( s < nlbnd  ) bring_inside_nxrange = s + nmol
    return

end function
!*********************************************************************!
!      the kronecker delta function... move                           !
!*********************************************************************!
integer function kd(n,m)
    integer, intent(in) :: n, m

    kd = 0
    if ( n == m ) kd = 1

    return
end function

