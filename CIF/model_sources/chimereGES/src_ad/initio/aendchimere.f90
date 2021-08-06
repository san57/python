!-*-f90-*-
subroutine aendchimere 

  use netcdf
  use chimere_consts
  use chimere_common
 
  implicit none

#define NCERR(lnum) if(ncstat/=NF90_NOERR) call nc_err(ncstat,lnum,'aendchimere.f90')

  !*****************************************************************************************
  integer :: ncstat

  !*****************************************************************************************
  !  All final operations                                                 

  ncstat=nf90_close(aoutea_ncid)
  NCERR(__LINE__)
  if (optemisb.ne.0) then
  ncstat=nf90_close(aouteb_ncid)
  NCERR(__LINE__)
  endif
  ncstat=nf90_close(aoutbc_ncid)
  NCERR(__LINE__)
  ncstat=nf90_close(aoutini_ncid)
  NCERR(__LINE__)
  ncstat=nf90_close(aend_ncid)
  NCERR(__LINE__)
!  close(79)

END subroutine aendchimere
