!*********************************************************************!
!   subroutine reads parameters from user input file 
!   most variables reside in the commonvar module and are used
!   throughout the program in various subroutines
!*********************************************************************!
subroutine read_in_para()
    use commonvar
    implicit none

    logical         exists
    character*100   buff, label, fname
    integer         fno, ios, line, pos, errstat
    parameter       (fno = 201)
    
    !get the input file name as the first argument from
    !command line, if given otherwise use default parameters
    call get_command_argument( 1, fname, status = errstat )
    if ( errstat .ne. 0 ) then
        print*, 'No control file given. Using default parameters'
        goto 1010
    end if
    !check that the given input file exists, abort if not
    inquire( file = trim( fname ), exist = exists )
    if ( .not. exists ) then
        print*, 'Input file not found...aborting'
        stop
    end if
    
    !open the user input file and read in the parameters
    open( unit = fno, file = fname, status = 'old', action = 'read' )
    ios = 0  !the in/out status
    line = 0 !the current line number
    print*, 'Reading the input file...'
    print*, '**********************************'//&
            '**********************************'
    do while ( ios == 0 ) !continue the loop until end of file
        read( fno, '(a)', iostat = ios ) buff !read a line
        if ( ios == 0 ) then
            line = line + 1
            !parse the line into a label and parameter
            pos   = scan( buff, ' ' )
            label = buff( 1:pos )
            buff  = buff( pos + 1: )
            if ( label(1:1) == '#' ) cycle !treat as a comment
            !find the label and assign the appropriate value to 
            !the variable
            select case ( label )
            case('task_title')
                read( buff, *, iostat=ios ) task_title
                print*, '   Setting task_title to: '//trim(task_title)
            case('nmol')
                read( buff, *, iostat=ios ) nmol
                print'(a,i4)', '    Setting nmol to: ', nmol
            case('vibmax')
                read( buff, *, iostat=ios ) vibmax
                print'(a,i4)', '    Setting vibmax to: ', vibmax
            case('hw')
                read( buff, *, iostat=ios ) hw
                print'(a,f8.2)', &
                    '    Setting vibrational energy to (cm-1): ', hw
            case('lambda')
                read( buff, *, iostat=ios ) lambda_n
                print'(a,f8.2)', '    Setting neutral lambda to: ', lambda_n
            case('lambda+')
                read( buff, *, iostat=ios ) lambda_c
                print'(a,f8.2)', '    Setting cation lambda to: ', &
                                      lambda_c  
            case('lambda-')
                read( buff, *, iostat=ios ) lambda_a
                print'(a,f8.2)', '    Setting anion lambda to: ', &
                                      lambda_a  
            case('JCoul')
                read( buff, *, iostat=ios) JCoul
                print'(a,f8.2)', '    Setting JCoul to (cm-1): ', JCoul   
            case('ES1')
                read( buff, *, iostat=ios) ES1
                print'(a,f8.2)', '    Setting ES1 to ', ES1
            case('te')
                read( buff, *, iostat=ios) te
                print'(a,f8.2)', '    Setting te to (cm-1): ', te   
            case('th')
                read( buff, *, iostat=ios) th
                print'(a,f8.2)', '    Setting th to (cm-1): ', th   
            case('ECT')
                read( buff, *, iostat=ios) ECT
                print'(a,f8.2)', '    Setting ECT to (cm-1): ', ECT
            case('one_state')
                read( buff, *, iostat=ios) one_state
                if ( one_state ) &
                    print*, '   One particle states are turned on.'
                if ( .not.one_state ) &
                    print*, '   One particle states are turned off.'
            case('two_state')
                read( buff, *, iostat=ios) two_state
                if ( two_state ) &
                    print*, '   Two particle states are turned on.'
                if ( .not.two_state ) &
                    print*, '   Two particle states are turned off.'
            case('ct_state')
                read( buff, *, iostat=ios) ct_state
                if ( ct_state ) &
                    print*, '   CT states are turned on.'
                if ( .not.ct_state ) &
                    print*, '   CT states are turned off.'
            case('abs_lw')
                read( buff, *, iostat=ios) abs_lw
                print'(a,f8.2)', &
                    '    Setting the linewidth to (cm-1): ', abs_lw
            case('esnum')
                read( buff, *, iostat=ios) esnum
                print'(a,i8)', '    Setting the number of eigenstates to',&
                                    esnum
            case default
                print*, '    invalid label at line, ', line
                print*, '    press enter to continue or ctrl+c to abort'
                read*
            end select
        end if
    end do
    close( fno )     !close the file
    print*, '**********************************'//&
            '**********************************'
    
1010 continue
    print*, 'Calculating derived parameters in units of hw.'

    !normalize parameters to units of hw    
    JCoul  = JCoul / hw
    ES1    = ES1 / hw
    te     = te / hw
    th     = th / hw
    ECT    = ECT / hw
    abs_lw = abs_lw /hw
    
    !Set all Huang-Rhys factors to zero if vibmax is zero
    !This assumes that the user just wants to calculate the
    !free exciton properties
    if ( vibmax == 0 ) then
        lambda_n = 0.d0
        lambda_c = 0.d0
        lambda_a = 0.d0
    end if

    !set the maximum left and right displacement from a 
    !given molecule given periodic boundary conditions
    nlbnd  =  -nmol/2+(1-1*mod(nmol,2))
    nubnd =   nmol/2
end subroutine
