!-*-f90-*-
subroutine endchimere 

  use netcdf
  use chimere_consts
  use chimere_common
  implicit none

#define NCERR(lnum) if(ncstat/=NF90_NOERR) call nc_err(ncstat,lnum,'endchimere.f90')


  !*****************************************************************************************
  integer :: ncstat

  !*****************************************************************************************
  !  All final operations                                                 


  print*,rank,'  Closing opened files',meteo_ncid,emisa_ncid,bounconc_ncid

  ncstat=nf90_close(meteo_ncid)
  NCERR(__LINE__)
  ncstat=nf90_close(emisa_ncid)
  NCERR(__LINE__)
  ncstat=nf90_close(bounconc_ncid)
  NCERR(__LINE__)
  if (ad.eq.0) then
  
  if(NetCDF_output.eq.1)then
    ncstat=nf90_close(out_ncid)
NCERR(__LINE__)
    ncstat=nf90_close(par_ncid)
NCERR(__LINE__)
  endif

  ncstat=nf90_close(end_ncid)
  NCERR(__LINE__)
  ncstat=nf90_close(depo_ncid)
  NCERR(__LINE__)
  
  endif 
  
END subroutine endchimere
