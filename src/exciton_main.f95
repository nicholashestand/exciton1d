!*********************************************************************!
!
!   This program calculates the unpolarized absorption spectrum and 
!   band dispersion of a linear, 1D aggregate using the Frenkel/
!   charge-transfer Holstein model. See the README file for more 
!   information.
!                                                   Written By:
!                                                   Nicholas Hestand
!*********************************************************************!
program exciton_main
    use commonvar
    implicit none

    integer k   !the k-index
    
    !read the user input file and set simulation parameters
    call read_in_para()     
   
    !index the multiparticle basis set
    kount = 0
    if ( one_state ) call index_1p()
    if ( two_state ) call index_2p()
    if ( ct_state ) call index_ct()
    
    !make sure the number of requested eigenstates is less than the
    !total number of eigenstates possible
    esnum = min(esnum,kount) 

    !build the franck-condon table for the vibrational overlap factors
    call set_fctable()   
    
    !allocate space for the Hamiltonian matrix and eigenvalue array
    allocate ( h(kount,kount) ) !the Hamiltonian
    allocate ( eval(kount) )    !eigenvalues


    print'(a)', ' Will now build the Hamiltonian, the'//&
                   ' diminsion of each k-block is:'
    print'(i6)',  kount
    print*, '**********************************'//&
            '**********************************'

    !build each k-block of the Hamiltonian diagonalize,
    !and calculate the observables
    do k = nlbnd, nubnd

        !initialize hamiltonian to zero
        h = complex_zero
        
        !build the hamiltonian
        if ( one_state ) call build_h1p(k)
        if ( two_state ) call build_h2p(k)
        if ( one_state .and. two_state ) call build_h1p2p(k)
        if ( ct_state ) call build_hct(k)
        if ( one_state .and. ct_state ) call build_h1pct(k)
        if ( two_state .and. ct_state ) call build_h2pct(k)

        !diagonalize the hamiltonian
        if ( k == 0 .or. esnum == kount ) then
            call diagonalize(h, kount, eval, 'A', kount )
        else
            call diagonalize(h, kount, eval, 'I', esnum)
        end if   
        
        !calculate absorption spectrum
        !(only k=0 absorbes assuming parallel dipoles)
        if ( k == 0 ) then
            call absorption()
        end if
        call dispersion(k)
        
        print*, 'Done with wavevector k: ', k
    end do

    !write the parameter file
    call para_out()
    print*, '**********************************'//&
            '**********************************'
    print*, 'Program exited successfully.'

end program
