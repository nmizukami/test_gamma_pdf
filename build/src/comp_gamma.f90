MODULE comp_gamma_module

!numeric type
USE nrtype

! privary
implicit none
private

public::check_gamma
public::basinUH

integer(i4b),parameter   :: integerMissing=-9999   ! missing value for integers
real(dp),    parameter   :: realMissing=-9999._dp  ! missing value for real numbers

contains

  SUBROUTINE check_gamma(fshape, tscale, frac_future, ierr, message)
  USE gamma_func_module, ONLY : gammp                   ! interface for the incomplete gamma function
  implicit none
  ! input
  real(dp), intent(in)                   :: fshape         ! shapef parameter in gamma distribution
  real(dp), intent(in)                   :: tscale         ! time scale parameter
  ! output
  real(dp), allocatable, intent(out)     :: frac_future(:) ! fraction in future time steps
  INTEGER(I4B), intent(out)              :: ierr           ! error code
  CHARACTER(*), intent(out)              :: message        ! error message
  ! locals
  real(dp)                               :: tfuture        ! future time (units of dt)
  real(dp)                               :: alamb          ! scale parameter
  integer(i4b)                           :: jTim           ! (loop through future time steps)
  integer(i4b)                           :: ntdh           ! number of values on the time delay histogram
  real(dp)                               :: cumprob        ! cumulative probability at JTIM
  real(dp)                               :: psave          ! cumulative probability at JTIM-1

  ! initialize error control
  ierr=0; message='check_gamma/'

  ntdh=100
  alamb = 1.0_dp/tscale                   ! scale parameter
  !alamb = fshape/tscale                  ! scale parameter

  ! allocate space for the time-delay histogram
  if (.not.allocated(frac_future)) then
    allocate(frac_future(ntdh), stat=ierr)
    if(ierr/=0)then; message=trim(message)//'unable to allocate space for the time delay histogram'; return; endif
  endif
  ! loop through time steps and compute the fraction of runoff in future time steps
  psave = 0.                                                 ! cumulative probability at JTIM-1
  do jTim=1,ntdh
   tfuture            = real(jTim, kind(dp))          ! future time
   cumprob            = gammp(fshape, alamb*tfuture)  ! cumulative probability at JTIM
   frac_future(jTim)  = max(0._DP, cumprob-psave)     ! probability between JTIM-1 and JTIM
   psave              = cumprob                       ! cumulative probability at JTIM-1
   !WRITE(*,'(I5,1X,F20.5,1X,2(F11.5))') JTIM, TFUTURE, frac_future(JTIM), CUMPROB
  end do
  ! ensure that the fractions sum to 1.0 (account for rounding errors, and not enough bins)
  frac_future(:) = frac_future(:) / sum(frac_future(:))

  END SUBROUTINE check_gamma


  SUBROUTINE basinUH(dt, fshape, tscale, frac_future, IERR, MESSAGE)
  USE gamma_func_module, ONLY : gammp                   ! interface for the incomplete gamma function
  IMPLICIT NONE
  ! input
  real(dp), intent(in)                   :: dt             ! model time step
  real(dp), intent(in)                   :: fshape         ! shapef parameter in gamma distribution
  real(dp), intent(in)                   :: tscale         ! time scale parameter
  ! output
  real(dp), allocatable, intent(out)     :: frac_future(:) ! fraction in future time steps
  INTEGER(I4B), INTENT(OUT)              :: IERR           ! error code
  CHARACTER(*), INTENT(OUT)              :: MESSAGE        ! error message
  ! locals
  REAL(DP)                               :: alamb          ! scale parameter
  REAL(DP)                               :: ntdh_min       ! minimum number of time delay points
  REAL(DP)                               :: ntdh_max       ! maximum number of time delay points
  REAL(DP)                               :: ntdh_try       ! trial number of time delay points
  INTEGER(I4B)                           :: itry           ! index of trial value
  INTEGER(I4B), PARAMETER                :: MAXTRY=100     ! maximum number of trial values
  INTEGER(I4B)                           :: ntdh           ! number of values on the time delay histogram
  INTEGER(I4B)                           :: JTIM           ! (loop through future time steps)
  REAL(DP)                               :: TFUTURE        ! future time (units of dt)
  REAL(DP)                               :: X_VALUE        ! xvalue to evaluate using gammp
  REAL(DP)                               :: CUMPROB        ! cumulative probability at JTIM
  REAL(DP)                               :: PSAVE          ! cumulative probability at JTIM-1
  ! ---------------------------------------------------------------------------------------
  ! initialize error control
  ierr=0; message='basinUH/'
  ! use a Gamma distribution with shape parameter, fshape, and time parameter, tscale, input
  !alamb = fshape/tscale                  ! scale parameter
  alamb = 1.0_dp/tscale                  ! scale parameter
  ! find the desired number of future time steps
  ! check if the cummulative Gamma distribution is close to 1.00 for given model time step, tscale and fsahpe.
  x_value = alamb*dt
  cumprob = gammp(fshape, x_value)
  if(cumprob > 0.999_dp) then ! in case if the cumprob is close to 1 in one model time step
   ntdh_try = 1.999_dp
  else
   ntdh_min = 1._dp
   ntdh_max = 1000._dp
   ntdh_try = 0.5_dp*(ntdh_min + ntdh_max)
   do itry=1,maxtry
    x_value = alamb*dt*ntdh_try
    cumprob = gammp(fshape, x_value)
    !print*, tscale, ntdh_try, cumprob, x_value, itry
    if(cumprob < 0.99_dp)  ntdh_min = ntdh_try
    if(cumprob > 0.999_dp) ntdh_max = ntdh_try
    if(cumprob > 0.99_dp .and. cumprob < 0.999_dp) exit
    ntdh_try = 0.5_dp*(ntdh_min + ntdh_max)
    if(itry==maxtry)then; ierr=20; message=trim(message)//'cannot identify the maximum number of bins for the tdh'; return; endif
   end do
  endif
  ntdh = ceiling(ntdh_try)

  ! allocate space for the time-delay histogram
  if (.not.allocated(frac_future)) then
    allocate(frac_future(ntdh), stat=ierr)
    if(ierr/=0)then; message=trim(message)//'unable to allocate space for the time delay histogram'; return; endif
  endif
  ! loop through time steps and compute the fraction of runoff in future time steps
  psave = 0.                                                 ! cumulative probability at JTIM-1
  DO jTim=1,ntdh
   tfuture            = REAL(jTim, kind(dp))*dt       ! future time
   cumprob            = gammp(fshape, alamb*tfuture)  ! cumulative probability at JTIM
   frac_future(jTim)  = MAX(0._DP, cumprob-psave)     ! probability between JTIM-1 and JTIM
   psave              = cumprob                       ! cumulative probability at JTIM-1
   !WRITE(*,'(I5,1X,F20.5,1X,2(F11.5))') JTIM, TFUTURE, frac_future(JTIM), CUMPROB
  END DO
  ! ensure that the fractions sum to 1.0 (account for rounding errors, and not enough bins)
  frac_future(:) = frac_future(:) / SUM(frac_future(:))
  ! ---------------------------------------------------------------------------------------
  END SUBROUTINE basinUH

END MODULE comp_gamma_module
