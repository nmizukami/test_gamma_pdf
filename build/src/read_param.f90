MODULE read_param_module

 USE nrtype, ONLY : i4b, dp, strLen

 private
 public :: read_param

 CONTAINS

  SUBROUTINE read_param(fname, fshape, tscale, dt, ierr, message)

   USE ascii_util_module, ONLY : file_open          ! open file (performs a few checks as well)

   implicit none
   ! input variables
   character(*),intent(in)    :: fname              ! parameter namelist name
   ! output variables
   real(dp),    intent(out)   :: fshape         ! shapef parameter in gamma distribution
   real(dp),    intent(out)   :: tscale         ! time scale parameter
   real(dp),    intent(out)   :: dt             ! model time step
   integer(i4b),intent(out)   :: ierr               ! error code
   character(*),intent(out)   :: message            ! error message
   ! local variables
   character(len=strLen)      :: cmessage           ! error message from subroutine
   integer(i4b)               :: iunit              ! file unit

   namelist /HSLOPE/fshape, tscale, dt                ! route simulated runoff through the local basin

   ! error control initialization
   ierr=0; message='read_param'

   ! read the name list
   call file_open(trim(fname),iunit, ierr, cmessage)
   if(ierr/=0)then; ierr=30; message=trim(cmessage)//': '//trim(fname); return;endif
   read(iunit, nml=HSLOPE)              ! basin IRF routing parameters
   close(iunit)

  END SUBROUTINE read_param

END MODULE read_param_module
