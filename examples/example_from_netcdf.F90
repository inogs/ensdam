program example_from_netcdf
  use ensdam_score_crps
  use ensdam_score_rcrv
  use ensdam_score_entropy
  use ensdam_score_optimality
  use ensdam_storng
  use ensdam_obserror
  use ensdam_meanstd
  use netcdf

  implicit none

  character(len=1024):: fileNetCDF
  integer ncid, stat, id_prior, id_post, id_obs, counter

  integer m        ! Size of the ensemble
  integer n        ! Size of the state vector
  real(kind=8), parameter :: sigma=0.3   ! Observation error standard deviation

  real(kind=8), allocatable, dimension(:,:) :: prior_ensemble
  real(kind=8), allocatable, dimension(:,:) :: posterior_ensemble

  real(kind=4), allocatable, dimension(:,:) :: prior_ensemble_read
  real(kind=4), allocatable, dimension(:,:) :: posterior_ensemble_read
  real(kind=4), allocatable, dimension(:)  :: observations_read


  real(kind=8), allocatable, dimension(:) :: reference_truth
  real(kind=8), allocatable, dimension(:) :: observations

  real(kind=8) :: crps,crps_reliability,crps_resolution  ! CRPS score
  real(kind=8) :: rcrv_bias, rcrv_spread                 ! RCRV score
  real(kind=8), dimension(2) :: entropy_score            ! entropy score
  real(kind=8) :: optimality                             ! optimality score

  real(kind=8), allocatable, dimension(:) :: ens_mean      ! mean of prior ensemble
  real(kind=8), allocatable, dimension(:) :: ens_var       ! variance of prior ensemble
  real(kind=8), allocatable, dimension(:) :: std_obs       ! observation error standard deviation
  real(kind=8), dimension(2,2) :: binary_pref ! Reference probability distribution for two binary events

  integer :: i, j    ! array indices

  integer :: nproc=1  ! Number of processors
  integer :: iproc=0  ! Current processor index
  real(kind=8) a
  integer, allocatable, dimension(:) :: outcome

#if defined MPI
  ! Definition to make the example compatible with MPI version of the library
  include "mpif.h"
  integer, save :: mpi_code

  ! Initialize parallel computation
  call mpi_init(mpi_code)
  call mpi_comm_size(mpi_comm_world,nproc,mpi_code)
  call mpi_comm_rank(mpi_comm_world,iproc,mpi_code)
#endif

! Read command line
!First, make sure the right number of inputs have been provided
IF(COMMAND_ARGUMENT_COUNT().NE.1)THEN
  WRITE(*,*)'ERROR, ONE COMMAND-LINE ARGUMENT (filename) REQUIRED, STOPPING'
  STOP
ENDIF
CALL GET_COMMAND_ARGUMENT(1,fileNetCDF)

   call getDimension(fileNetCDF,'m',m)
   call getDimension(fileNetCDF,'n',n)
   allocate(prior_ensemble(n,m))
   allocate(posterior_ensemble(n,m))

   allocate(prior_ensemble_read(n,m))
   allocate(posterior_ensemble_read(n,m))
   allocate(observations_read(n))

   allocate(reference_truth(n))
   allocate(observations(n))
   allocate(ens_mean(n))
   allocate(ens_var(n))
   allocate(std_obs(n))
   allocate(outcome(2))

   counter=0
   stat = nf90_open(fileNetCDF, nf90_nowrite, ncid)
   call handle_err1(stat, counter,FileNetCDF)
   stat = nf90_inq_varid (ncid, 'prior_ensemble', id_prior)
   call handle_err2(stat, fileNetCDF,'prior_ensemble')
   call handle_err1(stat, counter,FileNetCDF)
   stat = nf90_get_var (ncid,id_prior,prior_ensemble_read)
   call handle_err2(stat, fileNetCDF,'prior_ensemble')
   call handle_err1(stat, counter,FileNetCDF)

   stat = nf90_inq_varid (ncid, 'posterior_ensemble', id_post)
   call handle_err2(stat, fileNetCDF,'posterior_ensemble')
   call handle_err1(stat, counter,FileNetCDF)
   stat = nf90_get_var (ncid,id_post,posterior_ensemble_read)

   stat = nf90_inq_varid (ncid, 'observations', id_obs)
   call handle_err2(stat, fileNetCDF,'observations')
   call handle_err1(stat, counter,FileNetCDF)
   stat = nf90_get_var (ncid,id_post,observations_read)

   stat = nf90_close(ncid)
   call handle_err1(stat, counter,FileNetCDF)



  do j=1,m
  do i=1,n
    prior_ensemble(i,j) = real(prior_ensemble_read(i,j),8)
  enddo
  enddo
  do j=1,m
  do i=1,n
    posterior_ensemble(i,j) = real(posterior_ensemble_read(i,j),8)
  enddo
  enddo

  reference_truth=0.0
  do i=1, n
     call kiss_gaussian(a)
     reference_truth(i) = a*1.e-2  !1.e-7
  enddo
  do i=1,n
    reference_truth(i) = reference_truth(i) + real(observations_read(i),8)
    !reference_truth(i) = real(observations_read(i),8)
    observations(i) = real(observations_read(i),8)
  enddo



!  do j=1,m
!  do i=1,n
!    call kiss_gaussian(prior_ensemble(i,j))
!  enddo
!  enddo

  ! Sample reference truth from the same distribution
!  do i=1,n
!    call kiss_gaussian(reference_truth(i))
!  enddo

  ! Generate observations by adding N(0,sigma) perturbations to the reference truth
!  do i=1,n
!    call kiss_gaussian(observations(i))
!  enddo
!
!  observations(:) = reference_truth(:) + sigma * observations(:)

  ! Compute the posterior ensemble by conditioning the prior ensemble on observations
  ! i) Generate perttubations to observations for each ensemble member
!  do j=1,m
!    do i=1,n
!      call kiss_gaussian(posterior_ensemble(i,j))
!    enddo
!    posterior_ensemble(:,j) = sigma * posterior_ensemble(:,j) + observations(:)
!  enddo



  ! Compute CRPS score, using reference truth as verification data
  call crps_score(crps,crps_reliability,crps_resolution,prior_ensemble,reference_truth)
!  print '(a,2f8.5)', 'Prior CRPS reliability and resolution:    ',crps_reliability,crps_resolution
   print '(2f12.5)', crps_reliability,crps_resolution


  call crps_score(crps,crps_reliability,crps_resolution,posterior_ensemble,reference_truth)
!  print '(a,2f8.5)', 'Posterior CRPS reliability and resolution:',crps_reliability,crps_resolution
    print '(2f12.5)', crps_reliability,crps_resolution

  ! Compute RCRV score, using reference truth as verification data
  call rcrv_score(rcrv_bias,rcrv_spread,prior_ensemble,reference_truth)
!  print '(a,2e15.5)', 'Prior RCRV bias and spread:    ',rcrv_bias,rcrv_spread
   print '(2f12.5)', rcrv_bias,rcrv_spread


  call rcrv_score(rcrv_bias,rcrv_spread,posterior_ensemble,reference_truth)
!  print '(a,2e15.5)', 'Posterior RCRV bias and spread:',rcrv_bias,rcrv_spread
   print '(2f12.5)', rcrv_bias,rcrv_spread


  do j=1,m
  do i=1,n
    prior_ensemble(i,j) = log10(prior_ensemble(i,j))
  enddo
  enddo
  do j=1,m
  do i=1,n
    posterior_ensemble(i,j) = log10(posterior_ensemble(i,j))
  enddo
  enddo
  do i=1,n
    reference_truth(i) = log10(reference_truth(i))
    observations(i)    = log10(observations(i))
  enddo

  ! Compute CRPS score, using reference truth as verification data
  call crps_score(crps,crps_reliability,crps_resolution,prior_ensemble,reference_truth)
   print '(2f12.5)', crps_reliability,crps_resolution


  call crps_score(crps,crps_reliability,crps_resolution,posterior_ensemble,reference_truth)
    print '(2f12.5)', crps_reliability,crps_resolution

  ! Compute RCRV score, using reference truth as verification data
  call rcrv_score(rcrv_bias,rcrv_spread,prior_ensemble,reference_truth)
   print '(2f12.5)', rcrv_bias,rcrv_spread


  call rcrv_score(rcrv_bias,rcrv_spread,posterior_ensemble,reference_truth)
   print '(2f12.5)', rcrv_bias,rcrv_spread



  STOP













  ! Compute entropy score
  call events_probability(binary_pref,prior_ensemble,binary_event_outcomes)
  print '(a,2f6.3)', 'Prior probability distribution (event 1):    ',binary_pref(1,:)
  print '(a,2f6.3)', 'Prior probability distribution (event 2):    ',binary_pref(2,:)

  call events_probability(binary_pref,posterior_ensemble,binary_event_outcomes)
  print '(a,2f6.3)', 'Posterior probability distribution (event 1):',binary_pref(1,:)
  print '(a,2f6.3)', 'Posterior probability distribution (event 2):',binary_pref(2,:)

  call events_probability(binary_pref,prior_ensemble,binary_event_outcomes)
  call events_score(entropy_score,posterior_ensemble,binary_pref,binary_event_outcomes)
  print '(a,f6.3)', 'Entropy score (posterior vs prior, event 1):',entropy_score(1)
  print '(a,f6.3)', 'Entropy score (posterior vs prior, event 2):',entropy_score(2)

  ! Compute optimality score
  std_obs(:) = sigma

  call optimality_score(optimality,prior_ensemble,observations,std_obs)
  print '(a,f8.5)', 'Prior optimality score:    ',optimality

  call optimality_score(optimality,posterior_ensemble,observations,std_obs)
  print '(a,f8.5)', 'Posterior optimality score:',optimality

contains


  ! Callback routine to the outcome of the user defined events
  ! for a given ensemble member
  subroutine binary_event_outcomes(x,outcome)
  implicit none
  real(kind=8), dimension(:), intent(in) :: x    ! ensemble member
  integer, dimension(:), intent(out) :: outcome  ! events outcome

  if (sum(x*x)/(size(x,1)-1)<1.) then
    outcome(1) = 1
  else
    outcome(1) = 2
  endif

  if (maxval(abs(x))<3.4) then
    outcome(2) = 1
  else
    outcome(2) = 2
  endif

  end subroutine binary_event_outcomes


      !****************************************************************************


        SUBROUTINE getDIMENSION(fileNetCDF,dimname,n)
        use netcdf
        implicit none

        character,intent(in) :: fileNetCDF*(*) ,dimname*(*)
        integer,intent(inout) :: n



        ! local

        integer DIMid,ncid,stat
        character(LEN=100) junk
        integer counter

        counter = 0
        stat = nf90_open(fileNetCDF, nf90_nowrite, ncid)
       call handle_err1(stat, counter,FileNetCDF)
        stat = nf90_inq_dimid (ncid, dimname, DIMid)
       call handle_err1(stat, counter,FileNetCDF)
        stat = nf90_Inquire_Dimension (ncid, DIMid, junk, n)
       call handle_err1(stat, counter,FileNetCDF)
        stat = nf90_close(ncid)
       call handle_err1(stat, counter,FileNetCDF)
        END SUBROUTINE getDIMENSION


        subroutine handle_err1(status,mycount, fileNetCDF)
        USE netcdf
        integer status,mycount
        character fileNetCDF*(*)
        mycount =mycount+1
        if(status .ne. nf90_NoErr)  then
           write(*,*) 'netcdf call',mycount,'with status = ',status
           write(*,*)  'file :', fileNetCDF
           write(*,*) nf90_strerror(status)
           write(*,*) 'Stopped'
           STOP 1
        endif

        end subroutine handle_err1
        subroutine handle_err2(status,fileNetCDF,varname)
        USE netcdf
        integer status
        character fileNetCDF*(*) ,varname*(*)

        if(status .ne. nf90_NoErr)  then
           write(*,*) 'ERROR in Var = ', varname, ' file :', fileNetCDF
        endif

        end subroutine handle_err2


end program example_from_netcdf
