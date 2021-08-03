PROGRAM main

! ******
! provide access to external data, subroutines
! ****************************************************
! variable types
USE nrtype                                       ! variable types, etc.
USE comp_gamma_module,  only : basinUH                ! construct basin unit hydrograph
USE comp_gamma_module,  only : check_gamma                ! construct basin unit hydrograph
USE read_param_module,  only : read_param

implicit none

integer(i4b), parameter :: iulog=6

! ******
! define variables
! ************************
real(dp), allocatable            :: frac_future(:)
real(dp)                         :: fshape                          ! shape parameter in time delay histogram (=gamma distribution) [-]
real(dp)                         :: tscale                          ! scaling factor for the time delay histogram [sec]
real(dp)                         :: dt                              ! time step [sec]
integer(i4b)                     :: iTim                            !
character(len=strLen),parameter  :: ancil_dir = './'                ! name of the namelist file
character(len=strLen),parameter  :: param_nml = 'param.nml.default' ! name of the namelist file
integer(i4b)                     :: ierr                            ! error code
character(len=strLen)            :: cmessage                        ! error message of downwind routine


! read the routing parameter namelist
call read_param(trim(ancil_dir)//trim(param_nml), fshape, tscale, dt, ierr,cmessage)
if(ierr/=0) call handle_err(ierr, cmessage)

print*, 'fshape = ', fshape
print*, 'tscale = ', tscale
print*, 'dt     = ', dt

! get lag times in the basin unit hydrograph (not sure this is right place...)
!call basinUH(dt, fshape, tscale, frac_future, ierr, cmessage)   ! generate gamma distribution based UH (x axis is second)
call check_gamma(fshape, tscale, frac_future, ierr, cmessage)    ! generate gamma distribution (PDF) without dimension
if(ierr/=0) call handle_err(ierr, cmessage)

open (1, file = 'output.txt', status='replace', action='write')
do iTim = 1, size(frac_future)
  write(1,'(F10.1,x,E15.7)') real(iTim, kind=dp), frac_future(iTim)
end do

print*, 'sum of frac_future = ', sum(frac_future)

CONTAINS

 SUBROUTINE handle_err(err,message)
 implicit none
 integer(i4b),intent(in)::err             ! error code
 character(*),intent(in)::message         ! error message
 if(err/=0)then
   write(iulog,*) 'FATAL ERROR: '//trim(message)
   call flush(6)
   stop
 endif
 END SUBROUTINE handle_err

END PROGRAM main
