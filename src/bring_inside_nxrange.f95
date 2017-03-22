!*********************************************************************!
!   Helper function to keep indices inside the range [nlbnd, nubnd]   !
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
