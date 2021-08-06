subroutine iniphys

  !  Initialization of time parameters and meteo variables
  !  There are some changes in units

  use chimere_consts
  use chimere_common
  use wholedomain_common
  use netcdf
#ifdef IFORT
  !use ifport
  !use ifposix
#endif
  implicit none

#define NCERR(lnum) if(ncstat/=NF90_NOERR) call nc_err(ncstat,lnum,'iniphys.f90')

#define NOTCONT(vname) if (ncstat/=NF90_NOERR) call not_cont(vname)



#define NFCHAR NF90_CHAR
#define NFLOAT NF90_FLOAT

  !*****************************************************************************************
  integer :: iiy,iim,iid,iih,nh
  integer :: idm,i
  integer :: idex
  integer :: iystart,imstart,idstart,ihstart
  integer :: idyst

  character(len=dlen) :: datebuf       ! to hold current date in MM5 format

  ! netcdf stuff
  integer :: ncstat      ! return code for netCDF functions
  integer :: timedimid   ! Time dimension IDs
  integer :: datedimid   ! Date string dimension ID
  integer :: wedimid     ! WE dimension IDs, non-staggered
  integer :: sndimid     ! SN dimension IDs, non-staggered
  integer :: btdimid     ! BT dimension IDs, non-staggered
  integer :: we,sn,bt,ntimes                  ! dimensions lengths
  integer,dimension(1) :: stvec1   ! start vectors for R/W functions
  integer,dimension(1) :: cntvec1  ! count vectors for R/W functions
  integer,dimension(3) :: stvec3   ! start vectors for R/W functions
  integer,dimension(3) :: cntvec3  ! count vectors for R/W functions
  integer,dimension(4) :: stvec4   ! start vectors for R/W functions
  integer,dimension(4) :: cntvec4  ! count vectors for R/W functions

  character(len=1024) :: history_in
  integer :: hin_len
  integer :: dl
  integer :: idrefe

  ! functions
  integer :: idaytype
  integer :: interdat
  integer :: mm5date2numeric

  !************************************************************************************

  ! Initialization of time variables
  call ddate(idatestart,iystart,imstart,idstart,ihstart)
  call rdate(iystart,1,1,0,idyst)
  djul = interdat(idyst,idatestart)/24.
  do nh=0,nhourrun
     ! definition of 'idate' from the starting date 'idatestart'
     call reldat(idatestart,nh,idate(nh)) ! The YYYYMMDDHH date at all times
     call ddate(idate(nh)       &
          ,iyear(nh),imonth(nh) &         ! Gives date elements        
          ,iday( nh),ihour( nh))
     idtyp(nh) = idaytype(idate(nh))      ! The day type at all times
  enddo
  ! Interval (in hour) between years' winter solstice and start of run
  ! This is required in routine ZENITH
  call ddate(idatestart,iiy,iim,iid,iih) 
  call rdate(iiy,12,21,12,idrefe)
  soltim = interdat(idrefe,idatestart)

  ! METEO file (netCDF) is already open

  ! Read "history" attribute from input file
  ncstat=nf90_inquire_attribute(meteo_ncid,NF90_GLOBAL,'history',len=hin_len)
  NCERR(__LINE__)
  ncstat=nf90_get_att(meteo_ncid, NF90_GLOBAL, 'history', history_in)
  NCERR(__LINE__)
  history_in=history_in(1:hin_len)

  ! Checking dimensions consistency between netCDF file
  ! and chimere parameters
  !   dimensions IDs
  ncstat=nf90_inq_dimid(meteo_ncid,'Time',timedimid)
  NCERR(__LINE__)
  ncstat=nf90_inq_dimid(meteo_ncid,'DateStrLen',datedimid)
  NCERR(__LINE__)
  ncstat=nf90_inq_dimid(meteo_ncid,'west_east',wedimid)
  NCERR(__LINE__)
  ncstat=nf90_inq_dimid(meteo_ncid,'south_north',sndimid)
  NCERR(__LINE__)
  ncstat=nf90_inq_dimid(meteo_ncid,'bottom_top',btdimid)
  NCERR(__LINE__)
  !   lengths
  ncstat=nf90_inquire_dimension(meteo_ncid,timedimid,len=ntimes)
  NCERR(__LINE__)
  ncstat=nf90_inquire_dimension(meteo_ncid,datedimid,len=dl)
  NCERR(__LINE__)
  ncstat=nf90_inquire_dimension(meteo_ncid,wedimid,len=we)
  NCERR(__LINE__)
  ncstat=nf90_inquire_dimension(meteo_ncid,sndimid,len=sn)
  NCERR(__LINE__)
  ncstat=nf90_inquire_dimension(meteo_ncid,btdimid,len=bt)
  NCERR(__LINE__)
  if (dl/=dlen)        call exit1('iniphys : date format error in input file')
  if (we/=nzonal_domain)      call exit1( 'iniphys : west_east dimension inconsistency')
  if (sn/=nmerid_domain)      call exit1( 'iniphys : south_north dimension inconsistency')
  if (bt/=nverti)   call exit1( 'iniphys : bottom-top dimension inconsistency')

!  print*,' Reading variable IDs'
  ncstat=nf90_inq_varid(meteo_ncid,'Times', meteo_times_varid)
  NOTCONT('Times')
  ncstat=nf90_inq_varid(meteo_ncid, 'winz',  meteo_winz_varid)
  NOTCONT('winz')
  ncstat=nf90_inq_varid(meteo_ncid, 'winm',  meteo_winm_varid)
  NOTCONT('winm')
  ncstat=nf90_inq_varid(meteo_ncid, 'temp',  meteo_temp_varid)
  NOTCONT('temp')
  ncstat=nf90_inq_varid(meteo_ncid, 'sphu',  meteo_sphu_varid)
  NOTCONT('sphu')
  ncstat=nf90_inq_varid(meteo_ncid, 'airm',  meteo_airm_varid)
  NOTCONT('airm')
  ncstat=nf90_inq_varid(meteo_ncid, 'clwc',  meteo_clwc_varid)
  NOTCONT('clwc')
  ! lmbb nphourm
  ncstat=nf90_inq_varid(meteo_ncid,'nphourm',meteo_nphourm_varid)
  NOTCONT('nphourm')
  ! lmbb deep conv
  ncstat=nf90_inq_varid(meteo_ncid, 'dpeu',  meteo_dpeu_varid)
  NOTCONT('dpeu')
  ncstat=nf90_inq_varid(meteo_ncid, 'dped',  meteo_dped_varid)
  NOTCONT('dped')
  ncstat=nf90_inq_varid(meteo_ncid, 'dpdu',  meteo_dpdu_varid)
  NOTCONT('dpdu')
  ncstat=nf90_inq_varid(meteo_ncid, 'dpdd',  meteo_dpdd_varid)
  NOTCONT('dpdd')
  ncstat=nf90_inq_varid(meteo_ncid, 'alti',  meteo_alti_varid)
  NOTCONT('alti')
  ncstat=nf90_inq_varid(meteo_ncid, 'kzzz',  meteo_kzzz_varid)
  NOTCONT('kzzz')
  ncstat=nf90_inq_varid(meteo_ncid, 'tem2',  meteo_tem2_varid)
  NOTCONT('tem2')
  ncstat=nf90_inq_varid(meteo_ncid, 'sreh',  meteo_sreh_varid)
  NOTCONT('sreh')
  ncstat=nf90_inq_varid(meteo_ncid, 'atte',  meteo_atte_varid)
  NOTCONT('atte')
  ncstat=nf90_inq_varid(meteo_ncid, 'hght',  meteo_hght_varid)
  NOTCONT('hght')
  ncstat=nf90_inq_varid(meteo_ncid, 'usta',  meteo_usta_varid)
  NOTCONT('usta')
  ncstat=nf90_inq_varid(meteo_ncid, 'aerr',  meteo_aerr_varid)
  NOTCONT('aerr')
  ncstat=nf90_inq_varid(meteo_ncid, 'obuk',  meteo_obuk_varid)
  NOTCONT('obuk')
  ncstat=nf90_inq_varid(meteo_ncid, 'wsta',  meteo_wsta_varid)
  NOTCONT('wsta')
  ncstat=nf90_inq_varid(meteo_ncid, 'topc',  meteo_topc_varid)
  NOTCONT('topc')

!  print*,' Reading hourly meteo variables.'
!  print*,'Times (first hour)'

  ncstat=nf90_get_var( &
       meteo_ncid, &
       meteo_times_varid, &
       datebuf, &
       (/1,1/),(/dlen,1/))
  NCERR(__LINE__)
  idex = idate(0)
  idm=mm5date2numeric(datebuf)
  if(idm.ne.idex) then
     print *,'*** ERROR: WRONG EXPECTED DATE IN METEO FILE'
     print *,'IHOURRUN=0,  EXPECTED= ',idex,' METEO= ',idm
     call exit1('Exiting')
  endif

!  print*,' 1D fields'
  ! lm nphourm
  stvec1=(/1/)
  cntvec1=(/nhourrun+1/)
  ncstat=nf90_get_var(meteo_ncid,meteo_nphourm_varid,nphourm(0:nhourrun),stvec1,cntvec1)
  NCERR(__LINE__)
!  do nh=0,nhourrun
!    print*,nh,idate(nh),nphourm(nh)
!  enddo
  
!  print*,' 2D fields'
  stvec3=(/1,1,1/)
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

!  print*,' 3D fields'
  stvec4=(/1,1,1,1/)
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


!  print*,' Change of units for 3D variables to obtain molec, cm, s'
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
  
  
!  print*," Initialization of next hours' values"
  ! Times (second hour)
  ncstat=nf90_get_var( &
       meteo_ncid, &
       meteo_times_varid, &
       datebuf, &
       (/1,2/),(/dlen,1/))
  NCERR(__LINE__)
  idex = idate(1)
  idm=mm5date2numeric(datebuf)
  if(idm.ne.idex) then
    print *,'*** ERROR: WRONG EXPECTED DATE IN METEO FILE'
    print *,'IHOURRUN=1,  EXPECTED= ',idex,' METEO= ',idm
    call exit1('Exiting')
  endif

  ! 2D fields
  stvec3=(/1,1,2/)
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
  stvec4=(/1,1,1,2/)
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

end subroutine iniphys


subroutine not_cont(vname)
  character*(*) :: vname
    print *
    print *, 'iniphys.f90 : meteo file does not contain'
    print *, 'the mandatory variable ',vname
    call exit1('Exiting')
   end subroutine not_cont
