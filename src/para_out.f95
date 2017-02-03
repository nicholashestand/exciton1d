!*********************************************************************!
!        write parameters to file
!*********************************************************************!
subroutine para_out()
    use commonvar
    implicit none

    integer f_no
    character*100 f_name

    f_no = 26
    f_name =trim(task_title)//'_para.csv'
    open( unit=f_no, file = f_name )
    
    write( f_no, * ) 'parameter, values, all energies in cm-1'
    write( f_no, * ) 'task title, ', trim(task_title)
    write( f_no, * ) 'nmol, ', nmol
    write( f_no, * ) 'vibmax, ', vibmax
    write( f_no, * ) '@@@@@@,@@@@@@'
    write( f_no, * ) 'basis state, on'
    write( f_no, * ) '1p, ', one_state
    write( f_no, * ) '2p, ', two_state
    write( f_no, * ) 'ct, ', ct_state
    write( f_no, * ) 'total,',kount
    write( f_no, * ) '@@@@@@,@@@@@@'
    write( f_no, * ) 'vib energy,', hw
    write( f_no, * ) 'lambda,', lambda_n
    write( f_no, * ) 'lambda+,',lambda_c
    write( f_no, * ) 'lambda-',lambda_a
    write( f_no, * ) 'JCoul,',JCoul*hw
    write( f_no, * ) 'ES1,',ES1*hw
    write( f_no, * ) 'te,',te*hw
    write( f_no, * ) 'th,',th*hw
    write( f_no, * ) 'ECT,',ECT*hw    
    write( f_no, * ) '@@@@@@@@,@@@@@@@@'
    write( f_no, *) 'abs linewidth,', abs_lw*hw
    close( f_no )

end subroutine
