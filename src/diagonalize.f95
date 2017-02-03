!*********************************************************************!
!                        Complex Diagonalization                      !
!    Find only a set of eigenvectors and eigenvalues                  !
!    Makes a call to the lapack routine                               !
!*********************************************************************!
subroutine diagonalize(a, n, w, rrange, iu)
    implicit none

    character*1                :: jobz = 'V'
    character*1, intent(in)    :: rrange
    character*1                :: uplo = 'U'
    integer, intent(in)        :: n 
    complex*16, intent(inout)  :: a(n,n)
    integer	                   :: lda
    real*8                     :: vl = -3.d0
    real*8                     :: vu = 20.d0
    integer                    :: il=1
    integer, intent(in)        :: iu
    real*8                     :: abstol
    integer	                   :: m
    real*8, intent(out)        :: w(n)
    complex*16, allocatable    :: z(:,:)
    integer	                   :: ldz
    integer	                   :: isuppz( 2*n )
    complex*16, allocatable    :: work(:)
    integer	                   :: lwork
    real*8, allocatable        :: rwork(:)
    integer	                   :: lrwork   
    integer, allocatable       :: iwork(:)
    integer	                   :: liwork
    integer	                   :: info
    complex*16                 :: workdim
    real*8                     :: rworkdim
    integer                    :: iworkdim
	integer    i,j

    real*8, external :: dlamch

    abstol = dlamch( 'safe minimum' )
    lda = n
    ldz = n
    allocate( z(n,n) )
    !query for workspace dimensions
    !==============================================!
    !    query for work dimensions
    !==============================================!
    lwork  = -1
    lrwork = -1
    liwork = -1
    call zheevr( jobz, rrange, uplo, n, a,             &
                 lda, vl, vu, il, iu,                  &
                 abstol, m, w, z, ldz,                 &
                 isuppz, workdim, lwork, rworkdim, lrwork,   &
                 iworkdim, liwork, info )
    lwork = workdim
    lrwork = int(rworkdim)
    liwork = int(iworkdim)
    allocate( work(lwork), rwork(lrwork), iwork(liwork) )
    !calculate the eigenspectrum
    call zheevr( jobz, rrange, uplo, n, a,             &
                 lda, vl, vu, il, iu,                  &
                 abstol, m, w, z, ldz,                 &
                 isuppz, work, lwork, rwork, lrwork,   &
                 iwork, liwork, info )
    !send back the eigenvectors
    a(:,1:m) = z(:,1:m)     
    deallocate( work, rwork, iwork, z )
end subroutine
