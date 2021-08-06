subroutine inicdf_tl

  ! initializes the netcdf input files that remain open during all the chimere run

  use netcdf
  use chimere_consts
  use chimere_common
  use wholedomain_common
  implicit none

#define NCERR(lnum) if(ncstat/=NF90_NOERR) call nc_err(ncstat,lnum,'inicdf_tl.f90')

  !*****************************************************************************************

  ! local scalars
  integer :: ncstat
  integer :: meteo_time_dimid
  integer :: meteo_date_dimid
  integer :: meteo_zonal_dimid
  integer :: meteo_merid_dimid
  integer :: meteo_layers_dimid

  integer :: meteo_times
  integer :: dl,we,sn,bt
  real(kind=8),dimension(nzonal_domain,nmerid_domain) :: buf2d


!*****************************************************************************************

!  Opening anthropic emission increment file
  ncstat=nf90_open(fnemisaincr,NF90_NOWRITE,emisaincr_ncid)
  NCERR(__LINE__)  
  
  if(optemisb.ne.0) then
  !  Opening biogenic emission increment file
  ncstat=nf90_open(fnemisbincr,NF90_NOWRITE,emisbincr_ncid)
  NCERR(__LINE__)
  endif
  
  ! Opening lateral boundary concentration increments file
  ncstat=nf90_open(fnbounconcincr,NF90_NOWRITE,bounconcincr_ncid)
  NCERR(__LINE__)

end subroutine inicdf_tl
