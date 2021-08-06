subroutine ainiphys

  !  Initialization of meteo variables at the beginning of the adjoint
  !  There are some changes in units

  use chimere_consts
  use chimere_common
  use netcdf  
  use wholedomain_common
#ifdef IFORT
  !use ifport
  !use ifposix
#endif
  implicit none

#define NCERR(lnum) if(ncstat/=NF90_NOERR) call nc_err(ncstat,lnum,'ainiphys.f90')

#define NOTCONT(vname) if (ncstat/=NF90_NOERR) call not_cont(vname)

#define NFCHAR NF90_CHAR
#define NFLOAT NF90_FLOAT

  !*****************************************************************************************
  integer :: idm
  integer :: idex
  
  
  character(len=dlen) :: datebuf     ! to hold current date in MM5 format
 
 ! netcdf stuff
  integer :: ncstat                  ! return code for netCDF functions
  integer,dimension(1) :: stvec1   ! start vectors for R/W functions
  integer,dimension(1) :: cntvec1  ! count vectors for R/W functions
  integer,dimension(3) :: stvec3   ! start vectors for R/W functions
  integer,dimension(3) :: cntvec3  ! count vectors for R/W functions
  integer,dimension(4) :: stvec4   ! start vectors for R/W functions
  integer,dimension(4) :: cntvec4  ! count vectors for R/W functions
  
! External functions
  integer :: mm5date2numeric
  !************************************************************************************

 ! METEO file (netCDF) is already open

  ! Reading hourly meteo variables.
  ! Times (before-last hour)

  ncstat=nf90_get_var( &
       meteo_ncid, &
       meteo_times_varid, &
       datebuf, &
       (/1,nhourrun/),(/dlen,1/))
  NCERR(__LINE__)
  idex = idate(nhourrun-1)
  idm=mm5date2numeric(datebuf)
  if(idm.ne.idex) then
    print *,'*** ERROR: WRONG EXPECTED DATE IN METEO FILE'
    print *,'IHOURRUN=nhourrun-1,  EXPECTED= ',idex,' METEO= ',idm
    call exit1('Exiting')
  endif


  ! 2D fields
  stvec3=(/1,1,nhourrun/)
  cntvec3=(/nzonal_domain,nmerid_domain,1/)
  ncstat=nf90_get_var(meteo_ncid,meteo_tem2_varid,tem2(:,:,1),stvec3,cntvec3)
  NCERR(__LINE__)
  ncstat=nf90_get_var(meteo_ncid,meteo_sreh_varid,sreh(:,:,1),stvec3,cntvec3)
  NCERR(__LINE__)
  ncstat=nf90_get_var(meteo_ncid,meteo_atte_varid,atte(:,:,1),stvec3,cntvec3)
  NCERR(__LINE__)
  ncstat=nf90_get_var(meteo_ncid,meteo_hght_varid,hght(:,:,1),stvec3,cntvec3)
  NCERR(__LINE__)
  ncstat=nf90_get_var(meteo_ncid,meteo_usta_varid,usta(:,:,1),stvec3,cntvec3)
  NCERR(__LINE__)
  ncstat=nf90_get_var(meteo_ncid,meteo_aerr_varid,aerr(:,:,1),stvec3,cntvec3)
  NCERR(__LINE__)
  ncstat=nf90_get_var(meteo_ncid,meteo_obuk_varid,obuk(:,:,1),stvec3,cntvec3)
  NCERR(__LINE__)
  ncstat=nf90_get_var(meteo_ncid,meteo_wsta_varid,wsta(:,:,1),stvec3,cntvec3)
  NCERR(__LINE__)
  ncstat=nf90_get_var(meteo_ncid,meteo_topc_varid,topc(:,:,1),stvec3,cntvec3)
  NCERR(__LINE__)

  ! 3D fields
  stvec4=(/1,1,1,nhourrun/)
  cntvec4=(/nzonal_domain,nmerid_domain,nvert_raw,1/)

     ncstat=nf90_get_var(meteo_ncid,meteo_winz_varid,winz(:,:,:,1),stvec4,cntvec4)
     NCERR(__LINE__)
     ncstat=nf90_get_var(meteo_ncid,meteo_winm_varid,winm(:,:,:,1),stvec4,cntvec4)
     NCERR(__LINE__)
     ncstat=nf90_get_var(meteo_ncid,meteo_temp_varid,temp(:,:,:,1),stvec4,cntvec4)
     NCERR(__LINE__)
     ncstat=nf90_get_var(meteo_ncid,meteo_sphu_varid,sphu(:,:,:,1),stvec4,cntvec4)
     NCERR(__LINE__)
     ncstat=nf90_get_var(meteo_ncid,meteo_airm_varid,airm(:,:,:,1),stvec4,cntvec4)
     NCERR(__LINE__)
     ncstat=nf90_get_var(meteo_ncid,meteo_clwc_varid,clwc(:,:,:,1),stvec4,cntvec4)
     NCERR(__LINE__)
     ! lmbb deepconv
     ncstat=nf90_get_var(meteo_ncid,meteo_dpeu_varid,dpeu(:,:,:,1),stvec4,cntvec4)
     NCERR(__LINE__)
     ncstat=nf90_get_var(meteo_ncid,meteo_dped_varid,dped(:,:,:,1),stvec4,cntvec4)
     NCERR(__LINE__)
     ncstat=nf90_get_var(meteo_ncid,meteo_dpdu_varid,dpdu(:,:,:,1),stvec4,cntvec4)
     NCERR(__LINE__)
     ncstat=nf90_get_var(meteo_ncid,meteo_dpdd_varid,dpdd(:,:,:,1),stvec4,cntvec4)
     NCERR(__LINE__)
     ncstat=nf90_get_var(meteo_ncid,meteo_alti_varid,hlay(:,:,:,1),stvec4,cntvec4)
     NCERR(__LINE__)
     ncstat=nf90_get_var(meteo_ncid,meteo_kzzz_varid,kzzz(:,:,:,1),stvec4,cntvec4)
     NCERR(__LINE__)


  ! Change of units for 3D variables to obtain molec, cm, s
  winz(:,:,:,1) = winz(:,:,:,1) * 1d2
  winm(:,:,:,1) = winm(:,:,:,1) * 1d2
  sphu(:,:,:,1) = sphu(:,:,:,1) * airm(:,:,:,1) * 1.6
  hlay(:,:,:,1) = hlay(:,:,:,1) * 1d2
  clwc(:,:,:,1) = clwc(:,:,:,1) * airm(:,:,:,1) *29./an
  kzzz(:,:,:,1) = kzzz(:,:,:,1) * 1d4
  
! lmbb * 1e-1 factor for conversion from kg/m2/s into g/cm2/s 
  dpeu(:,:,:,1) = dpeu(:,:,:,1) * an/29.*1.d-1
  dped(:,:,:,1) = dped(:,:,:,1) * an/29.*1.d-1
  dpdu(:,:,:,1) = dpdu(:,:,:,1) * an/29.*1.d-1
  dpdd(:,:,:,1) = dpdd(:,:,:,1) * an/29.*1.d-1
  
  !  Change of units for 2D variables
  hght(:,:,1) = hght(:,:,1)*1d2
  usta(:,:,1) = usta(:,:,1)*1d2
  aerr(:,:,1) = aerr(:,:,1)*1d-2
  obuk(:,:,1) = obuk(:,:,1)*1d2
  wsta(:,:,1) = wsta(:,:,1)*1d2
  where(topc(:,:,1).lt.dzero) topc(:,:,1)=dzero
  
  
  ! Initialization of next hours' values
  ! Times (last hour)
  ncstat=nf90_get_var( &
       meteo_ncid, &
       meteo_times_varid, &
       datebuf, &
       (/1,nhourrun+1/),(/dlen,1/))
  NCERR(__LINE__)
  idex = idate(nhourrun)
  idm=mm5date2numeric(datebuf)
  if(idm.ne.idex) then
     print *,'*** ERROR: WRONG EXPECTED DATE IN METEO FILE'
     print *,'IHOURRUN=nhourrun+1,  EXPECTED= ',idex,' METEO= ',idm
     stop
  endif

  ! 2D fields
  stvec3=(/1,1,nhourrun+1/)
  cntvec3=(/nzonal_domain,nmerid_domain,1/)
  ncstat=nf90_get_var(meteo_ncid,meteo_tem2_varid,tem2(:,:,2),stvec3,cntvec3)
  NCERR(__LINE__)
  ncstat=nf90_get_var(meteo_ncid,meteo_sreh_varid,sreh(:,:,2),stvec3,cntvec3)
  NCERR(__LINE__)
  ncstat=nf90_get_var(meteo_ncid,meteo_atte_varid,atte(:,:,2),stvec3,cntvec3)
  NCERR(__LINE__)
  ncstat=nf90_get_var(meteo_ncid,meteo_hght_varid,hght(:,:,2),stvec3,cntvec3)
  NCERR(__LINE__)
  ncstat=nf90_get_var(meteo_ncid,meteo_usta_varid,usta(:,:,2),stvec3,cntvec3)
  NCERR(__LINE__)
  ncstat=nf90_get_var(meteo_ncid,meteo_aerr_varid,aerr(:,:,2),stvec3,cntvec3)
  NCERR(__LINE__)
  ncstat=nf90_get_var(meteo_ncid,meteo_obuk_varid,obuk(:,:,2),stvec3,cntvec3)
  NCERR(__LINE__)
  ncstat=nf90_get_var(meteo_ncid,meteo_wsta_varid,wsta(:,:,2),stvec3,cntvec3)
  NCERR(__LINE__)
  ncstat=nf90_get_var(meteo_ncid,meteo_topc_varid,topc(:,:,2),stvec3,cntvec3)
  NCERR(__LINE__)

  ! 3D fields
  stvec4=(/1,1,1,nhourrun+1/)
  cntvec4=(/nzonal_domain,nmerid_domain,nvert_raw,1/)

     ncstat=nf90_get_var(meteo_ncid,meteo_winz_varid,winz(:,:,:,2),stvec4,cntvec4)
     NCERR(__LINE__)
     ncstat=nf90_get_var(meteo_ncid,meteo_winm_varid,winm(:,:,:,2),stvec4,cntvec4)
     NCERR(__LINE__)
     ncstat=nf90_get_var(meteo_ncid,meteo_temp_varid,temp(:,:,:,2),stvec4,cntvec4)
     NCERR(__LINE__)
     ncstat=nf90_get_var(meteo_ncid,meteo_sphu_varid,sphu(:,:,:,2),stvec4,cntvec4)
     NCERR(__LINE__)
     ncstat=nf90_get_var(meteo_ncid,meteo_airm_varid,airm(:,:,:,2),stvec4,cntvec4)
     NCERR(__LINE__)
     ncstat=nf90_get_var(meteo_ncid,meteo_clwc_varid,clwc(:,:,:,2),stvec4,cntvec4)
     NCERR(__LINE__)
     ! lmbb deepconv
     ncstat=nf90_get_var(meteo_ncid,meteo_dpeu_varid,dpeu(:,:,:,2),stvec4,cntvec4)
     NCERR(__LINE__)
     ncstat=nf90_get_var(meteo_ncid,meteo_dped_varid,dped(:,:,:,2),stvec4,cntvec4)
     NCERR(__LINE__)
     ncstat=nf90_get_var(meteo_ncid,meteo_dpdu_varid,dpdu(:,:,:,2),stvec4,cntvec4)
     NCERR(__LINE__)
     ncstat=nf90_get_var(meteo_ncid,meteo_dpdd_varid,dpdd(:,:,:,2),stvec4,cntvec4)
     NCERR(__LINE__)
     ncstat=nf90_get_var(meteo_ncid,meteo_alti_varid,hlay(:,:,:,2),stvec4,cntvec4)
     NCERR(__LINE__)
     ncstat=nf90_get_var(meteo_ncid,meteo_kzzz_varid,kzzz(:,:,:,2),stvec4,cntvec4)
     NCERR(__LINE__)


  ! Change of units for 3D variables to obtain molec, cm, s
  
  winz(:,:,:,2) = winz(:,:,:,2)*1d2
  winm(:,:,:,2) = winm(:,:,:,2)*1d2
  sphu(:,:,:,2) = sphu(:,:,:,2)*airm(:,:,:,2)*1.6
  hlay(:,:,:,2) = hlay(:,:,:,2)*1d2
  clwc(:,:,:,2) = clwc(:,:,:,2)*airm(:,:,:,1)*29./an
  kzzz(:,:,:,2) = kzzz(:,:,:,2)*1d4

! lmbb * 1e-1 factor for conversion from kg/m2/s into g/cm2/s 
  dpeu(:,:,:,2) = dpeu(:,:,:,2) * an/29.*1.d-1
  dped(:,:,:,2) = dped(:,:,:,2) * an/29.*1.d-1
  dpdu(:,:,:,2) = dpdu(:,:,:,2) * an/29.*1.d-1
  dpdd(:,:,:,2) = dpdd(:,:,:,2) * an/29.*1.d-1

  ! Change of units for 2D variables
        hght(:,:,2) = max(hght(:,:,2)*1d2, hlay(:,:,1,2)+1d2)
        usta(:,:,2) = usta(:,:,2)*1d2
        aerr(:,:,2) = aerr(:,:,2)*1d-2
        obuk(:,:,2) = obuk(:,:,2)*1d2
        wsta(:,:,2) = wsta(:,:,2)*1d2
        where (topc(:,:,2).lt.dzero) topc(:,:,2)=dzero

end subroutine ainiphys

