subroutine inicdf

  ! initializes the netcdf input files that remain open during all the chimere run

  use netcdf
  use chimere_consts
  use chimere_common
  use wholedomain_common
  implicit none

#define NCERR(lnum) if(ncstat/=NF90_NOERR) call nc_err(ncstat,lnum,'inicdf.f90')

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

! lmbb add flags for files opening: emis anthro
  !print*,'  Opening anthropic emission file'
  ncstat=nf90_open(fnemisa,NF90_NOWRITE,emisa_ncid)
  NCERR(__LINE__)
  
! lmbb biogenic file must always be opened   
  if(optemisb.ne.0) then
  !print*,'  Opening biogenic emission file'
  ncstat=nf90_open(fnemisb,NF90_NOWRITE,emisb_ncid)
  NCERR(__LINE__)
  endif
  
  print*,' Opening METEO file',fnmeteo
  ncstat=nf90_open(fnmeteo,NF90_NOWRITE,meteo_ncid)
  !print*,rank,'METEO file number',meteo_ncid
  NCERR(__LINE__)

  ! Check correctness of dimensions
  ! IDs
  ncstat=nf90_inq_dimid(meteo_ncid, 'Time',        meteo_time_dimid)
  NCERR(__LINE__)
  ncstat=nf90_inq_dimid(meteo_ncid, 'DateStrLen',  meteo_date_dimid)
  NCERR(__LINE__)
  ncstat=nf90_inq_dimid(meteo_ncid, 'west_east',   meteo_zonal_dimid)
  NCERR(__LINE__)
  ncstat=nf90_inq_dimid(meteo_ncid, 'south_north', meteo_merid_dimid)
  NCERR(__LINE__)
  ncstat=nf90_inq_dimid(meteo_ncid, 'bottom_top',  meteo_layers_dimid)
  NCERR(__LINE__)

  ! lengths
  ncstat=nf90_inquire_dimension(meteo_ncid, &
       meteo_time_dimid,len=meteo_times)
  NCERR(__LINE__)
  if(meteo_times<nhourrun) call exit1('inicdf : Not enough time slots in METEO file')

  ncstat=nf90_inquire_dimension(meteo_ncid,meteo_date_dimid,len=dl)
  NCERR(__LINE__)
  if (dl/=dlen) call exit1('inicdf : date format error in METEO file')

  ncstat=nf90_inquire_dimension(meteo_ncid,meteo_zonal_dimid,len=we)
  NCERR(__LINE__)
  if (we/=nzonal_domain) call exit1('inicdf : WE dimension error in METEO file')

  ncstat=nf90_inquire_dimension(meteo_ncid,meteo_merid_dimid,len=sn)
  NCERR(__LINE__)
  if (sn/=nmerid_domain) call exit1('inicdf : SN dimension error in METEO file')

  ncstat=nf90_inquire_dimension(meteo_ncid,meteo_layers_dimid,len=bt)
  NCERR(__LINE__)
  if (bt/=nvert_raw) then
    print *, 'bottom_top: ', bt, ' /=  nvert_raw: ', nvert_raw
    call exit1('inicdf : BT dimension error in METEO file')
  endif


  !print*,'  Read longitude/latitude'

  !print*,' Coordinates are read from METEO file (already open)'
  ncstat=nf90_inq_varid(meteo_ncid,'lon',meteo_lon_varid)
  NCERR(__LINE__)

  ncstat=nf90_inq_varid(meteo_ncid,'lat',meteo_lat_varid)
  NCERR(__LINE__)
  ncstat=nf90_inq_varid(meteo_ncid,'a_vcoord',meteo_avcoord_varid)
  NCERR(__LINE__)
  ncstat=nf90_inq_varid(meteo_ncid,'b_vcoord',meteo_bvcoord_varid)
  NCERR(__LINE__)

  ncstat=nf90_get_var(meteo_ncid,meteo_lon_varid,buf2d)
  NCERR(__LINE__)
  xlong=buf2d

  ncstat=nf90_get_var(meteo_ncid,meteo_lat_varid,buf2d)
  NCERR(__LINE__)
  xlati=buf2d

  ncstat=nf90_get_var(meteo_ncid,meteo_avcoord_varid,avcoord)
  NCERR(__LINE__)
  ncstat=nf90_get_var(meteo_ncid,meteo_bvcoord_varid,bvcoord)
  NCERR(__LINE__)
 
! lmbb add flags for files opening: emis anthro

  !print*,' Opening lateral boundary concentrations file'
  ncstat=nf90_open(fnbounconc,NF90_NOWRITE,bounconc_ncid)
  NCERR(__LINE__)

end subroutine inicdf
