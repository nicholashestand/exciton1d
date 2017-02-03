!*********************************************************************!
!    Module containing variables and parameters used in exciton1D     !
!    subroutines.
!*********************************************************************!
module commonvar
    implicit none

    !=================================================================!
    !   Simulation parameters, all in units of hw                     !
    !=================================================================!

    !Simulation Title
    character*100     task_title

    !Aggregate Size
    integer :: nmol = 1
    
    !Vibrational Parameters
    integer:: vibmax    = 0           !Max vibrations in basis set
    real*8 :: hw        = 1400.d0     !vibration energy
    real*8 :: lambda_n  = 1.d0        !Neutral lambda (harmonic shift)
    real*8 :: lambda_c  = 0.d0        !Cation lambda (harmonic shift)
    real*8 :: lambda_a  = 0.d0        !Anion lambda (harmonic shift)

    !Coupling and Energies
    real*8 :: JCoul     = 0.d0        !Nearest neighbor Coulomb coupling
    real*8 :: ES1       = 0.d0        !Monomer Transition Energy
    real*8 :: te        = 0.d0        !Nearest neighbor electron transfer integral
    real*8 :: th        = 0.d0        !Nearest neighbor hole transfer integral
    real*8 :: ECT       = 0.d0        !Charge Transfer Energy
    
    !multiparticle basis states
    logical :: one_state    =.true.
    logical :: two_state    =.false.
    logical :: ct_state     =.false.
    
    !absorption linewidth
    real*8  :: abs_lw   = 0.1d0
    
    !franck condon tables
    real*8, allocatable :: fc_gf(:,:)
    real*8, allocatable :: fc_gc(:,:)
    real*8, allocatable :: fc_ga(:,:)
    real*8, allocatable :: fc_af(:,:)
    real*8, allocatable :: fc_cf(:,:)

    !constants
    real*8, parameter :: pi = 4.d0*datan(1.d0)
    !cm-1 per electronvolt
    real*8, parameter :: ev = 8065.d0
    !plancks constant times the speed of light in nm*hw 
    real*8, parameter :: hc = 1.23984193d3 * ev 
    !boltzman constant units of cm-1 k
    real*8, parameter :: kbman = 0.6956925d0 
    !reduced planks constant in wavenumber * s
    real*8, parameter :: hbar = 6.58211951440d-16 * ev 

    !basis set counter
    integer :: kount = 0

    !basis set indexes
    integer, allocatable :: nx_1p(:)
    integer, allocatable :: nx_2p(:,:,:)
    integer, allocatable :: nx_ct(:,:,:)

    !the hamiltonian and eigenvalues
    complex*16, allocatable :: h(:,:)
    real*8, allocatable :: eval(:)       

    !empty parameter
    integer, parameter :: empty = -1  
    !parameters for complex numbers
    complex*16, parameter :: complex_zero = ( 0.d0, 0.d0 )
    complex*16, parameter :: img = ( 0.d0, 1.d0 )
    
    !bounds
    integer nlbnd, nubnd
            
    !number of eigenstates to find for each k
    integer :: esnum = 1    
            
end module
