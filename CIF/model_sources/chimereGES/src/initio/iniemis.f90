subroutine iniemis

  !  Initialization of emissions

  use netcdf
  use chimere_consts
  use chimere_common
  use wholedomain_common
  implicit none

#define NCERR(lnum) if(ncstat/=NF90_NOERR) call nc_err(ncstat,lnum,'iniemis.f90')

  !********************************************************************************
  character(len=15) :: charspec
  character(len=splen) :: name
  integer :: ns,ne,nk,np,i,j,l
  integer :: ide,idex
  integer :: ifnanthro,ifnbiogen
  integer :: dl,spl,mbio,ph

  character(len=dlen) :: datebuf     ! to hold current date in MM5 format

  ! netcdf stuff
  integer :: ncstat                  ! return code for netCDF functions
  integer :: atimedimid,btimedimid   ! Time dimension IDs
  integer :: aspstrdimid,bspstrdimid    ! Species string dimension ID
  integer :: adatestrdimid,bdatestrdimid   ! Date string dimension ID
  integer :: awedimid,bwedimid       ! WE dimension IDs, non-staggered
  integer :: asndimid,bsndimid       ! SN dimension IDs, non-staggered
  integer :: abtdimid,bbtdimid                ! BT dimension IDs, non-staggered
  integer :: nadimid,nbiodimid
  integer :: we,sn,na,ntimes,blev                  ! dimensions lengths
  integer,dimension(4) :: stvec   ! start vectors for R/W functions
  integer,dimension(4) :: cntvec  ! count vectors for R/W functions

  ! External functions
  integer :: mm5date2numeric

   character*15 cspec
   character*100000 usedval,notusedval
  !*****************************************************************************
  !  Initialisation of flag/address arrays

  inemisa = 0
  inemisb = 0

!======================================================================
! lmbb Reading anthropic species names
  notusedval=''
  call opfi(ifnanthro,fnanthro,'f','o',0)
  do ne=1,nemisa
    read(ifnanthro,*)we,charspec
    call findspec(charspec,ns)
    if(ns.gt.0.and.ns.le.nspec) then
        inemisa(ns) = ne
    else
        notusedval=notusedval(1:len_trim(notusedval))  &
	&      //charspec(1:len_trim(charspec))//'; '
5422    continue
    endif
  enddo
  close(ifnanthro)

  print*,'-------------------------------------------'
  print *,'* INIEMIS: Emitted but ignored species '
  if(len_trim(notusedval).gt.1) print*,'   Anthropogenic: ',notusedval(1:len_trim(notusedval)-1)
  

  ! Anthropic emission file (netCDF) is already open

  ! Checking dimensions consistency between file and chimere parameters
  ! dimensions IDs
  ncstat=nf90_inq_dimid(emisa_ncid,'Time',atimedimid)
  NCERR(__LINE__)
  ncstat=nf90_inq_dimid(emisa_ncid,'west_east',awedimid)
  NCERR(__LINE__)
  ncstat=nf90_inq_dimid(emisa_ncid,'south_north',asndimid)
  NCERR(__LINE__)
  ncstat=nf90_inq_dimid(emisa_ncid,'bottom_top',abtdimid)
  NCERR(__LINE__)
  ncstat=nf90_inq_dimid(emisa_ncid,'SpStrLen',aspstrdimid)
  NCERR(__LINE__)
  ncstat=nf90_inq_dimid(emisa_ncid,'DateStrLen',adatestrdimid)
  NCERR(__LINE__)
  ncstat=nf90_inq_dimid(emisa_ncid,'Species',nadimid)
  NCERR(__LINE__)
  ! lengths
  ncstat=nf90_inquire_dimension(emisa_ncid,atimedimid,len=ntimes)
  NCERR(__LINE__)
  ncstat=nf90_inquire_dimension(emisa_ncid,awedimid,len=we)
  NCERR(__LINE__)
  ncstat=nf90_inquire_dimension(emisa_ncid,asndimid,len=sn)
  NCERR(__LINE__)
  ncstat=nf90_inquire_dimension(emisa_ncid,aspstrdimid,len=spl)
  NCERR(__LINE__)
  ncstat=nf90_inquire_dimension(emisa_ncid,adatestrdimid,len=dl)
  NCERR(__LINE__)
  ncstat=nf90_inquire_dimension(emisa_ncid,nadimid,len=na)
  NCERR(__LINE__)
  if (we/=nzonal_domain)      call exit1( 'iniemis : west_east dimension inconsistency')
  if (sn/=nmerid_domain)      call exit1( 'iniemis : south_north dimension inconsistency')
  if (dl/=dlen)        call exit1( 'iniemis : date format error in A-emissions file')
  if (spl/=splen)      call exit1( 'iniemis : species name format error in A-emissions file')
  ! checking time
  ncstat=nf90_inq_varid(emisa_ncid,'Times', emisa_times_varid)
  NCERR(__LINE__)
  !! first hour
  ncstat=nf90_get_var( &
       emisa_ncid, &
       emisa_times_varid, &
       datebuf, &
       (/1,1/),(/dlen,1/))
  NCERR(__LINE__)
  idex=idate(0)
  ide=mm5date2numeric(datebuf)
  if(ide.ne.idex) then
     print *,'*** ERROR: WRONG EXPECTED DATE IN AEMISSIONS FILE'
     print *,'IHOURRUN=0,  EXPECTED= ',idex,' AEMISSIONS= ',ide
     call exit1('Exiting')

  endif
  !! next hour
  ncstat=nf90_get_var( &
       emisa_ncid, &
       emisa_times_varid, &
       datebuf, &
       (/1,2/),(/dlen,1/))
  NCERR(__LINE__)
  idex = idate(1)
  ide=mm5date2numeric(datebuf)
  if(ide.ne.idex) then
     print *,'*** ERROR: WRONG EXPECTED DATE IN AEMISSIONS FILE'
     print *,'IHOURRUN=1,  EXPECTED= ',idex,' AEMISSIONS= ',ide
     call exit1('Exiting')
  endif

  ! Getting the varids of each species in AEMISSIONS file
  call opfi(ifnanthro,fnanthro,'f','o',0)
  do ne=1,nemisa
     read(ifnanthro,*)we,charspec
     ncstat=nf90_inq_varid(emisa_ncid,charspec,emisavarid(ne))
     NCERR(__LINE__)
  end do
  close(ifnanthro)

!============================================================================
  !  Reading biogenic species names

  notusedval=''
  call opfi(ifnbiogen,fnbiogen,'f','o',0)
  do ne=1,nemisb
     read(ifnbiogen,*)we,charspec
     call findspec(charspec,ns)
     if(ns.gt.0.and.ns.le.nspec) then
        inemisb(ns) = ne
     else
	notusedval=notusedval(1:len_trim(notusedval))  &
	&      //charspec(1:len_trim(charspec))//'; '
5423    continue
     endif
  enddo
  close(ifnbiogen)
  print*,'  Biogenic: ',notusedval(1:len_trim(notusedval)-1)
  
  if(optemisb.ne.0) then
  ! Biogenic emission file (netCDF) is already open
  ! Checking dimensions consistency between file and chimere parameters
  ! dimensions IDs
  ncstat=nf90_inq_dimid(emisb_ncid,'Time',btimedimid)
  NCERR(__LINE__)
  ncstat=nf90_inq_dimid(emisb_ncid,'DateStrLen',bdatestrdimid)
  NCERR(__LINE__)
  ncstat=nf90_inq_dimid(emisb_ncid,'west_east',bwedimid)
  NCERR(__LINE__)
  ncstat=nf90_inq_dimid(emisb_ncid,'south_north',bsndimid)
  NCERR(__LINE__)
  ncstat=nf90_inq_dimid(emisb_ncid,'bottom_top',bbtdimid)
  NCERR(__LINE__)
!  ncstat=nf90_inq_dimid(emisb_ncid,'biospecies',nbiodimid)
  ncstat=nf90_inq_dimid(emisb_ncid,'SpStrLen',bspstrdimid)
  NCERR(__LINE__)
  ncstat=nf90_inq_dimid(emisb_ncid,'DateStrLen',bdatestrdimid)
  NCERR(__LINE__)
  ncstat=nf90_inq_dimid(emisb_ncid,'Species',nbiodimid)
  NCERR(__LINE__)
  ! lengths
  ncstat=nf90_inquire_dimension(emisb_ncid,btimedimid,len=ntimes)
  NCERR(__LINE__)
  ncstat=nf90_inquire_dimension(emisb_ncid,bdatestrdimid,len=dl)
  NCERR(__LINE__)
  ncstat=nf90_inquire_dimension(emisb_ncid,bwedimid,len=we)
  NCERR(__LINE__)
  ncstat=nf90_inquire_dimension(emisb_ncid,bsndimid,len=sn)
  NCERR(__LINE__)  
  ncstat=nf90_inquire_dimension(emisb_ncid,bbtdimid,len=blev)
  NCERR(__LINE__)
  ncstat=nf90_inquire_dimension(emisb_ncid,nbiodimid,len=mbio)
  NCERR(__LINE__)
  ncstat=nf90_inquire_dimension(emisb_ncid,bspstrdimid,len=spl)
  NCERR(__LINE__)
  
  if (dl/=dlen)        call exit1( 'iniemis : date format error in input file')
  if (we/=nzonal_domain)      call exit1( 'iniemis : west_east dimension inconsistency')
  if (sn/=nmerid_domain)      call exit1( 'iniemis : south_north dimension inconsistency')
  if(blev/=1)          call exit1('iniemis : number of levels for B-emissions not equal to 1')
  if (spl/=splen)      call exit1( 'iniemis : species name format error in B-emissions file')
  ! checking time
  ncstat=nf90_inq_varid(emisb_ncid,'Times', emisb_times_varid)
  NCERR(__LINE__)
!  ncstat=nf90_inq_varid(emisb_ncid, 'emisb',  emisb_emisb_varid)
!  NCERR(__LINE__)

  !  Reading biogenic emissions
!  cntvec=(/nemisb,we,sn,1/)
  !! first hour
  ncstat=nf90_get_var( &
       emisb_ncid, &
       emisb_times_varid, &
       datebuf, &
       (/1,1/),(/dlen,1/))
  NCERR(__LINE__)
  idex=idate(0)
  ide=mm5date2numeric(datebuf)
  if(ide.ne.idex) then
     print *,'*** ERROR: WRONG EXPECTED DATE IN BEMISSIONS FILE'
     print *,'IHOURRUN=0,  EXPECTED= ',idex,' BEMISSIONS= ',ide
     call exit1('Exiting')
  endif
  !! next hour
  ncstat=nf90_get_var( &
       emisb_ncid, &
       emisb_times_varid, &
       datebuf, &
       (/1,2/),(/dlen,1/))
  NCERR(__LINE__)
  idex = idate(1)
  ide=mm5date2numeric(datebuf)
  if(ide.ne.idex) then
     print *,'*** ERROR: WRONG EXPECTED DATE IN BEMISSIONS FILE'
     print *,'IHOURRUN=1,  EXPECTED= ',idex,' BEMISSIONS= ',ide
     call exit1('Exiting')
  endif
  ! Getting the varids of each species in BEMISSIONS file
  call opfi(ifnbiogen,fnbiogen,'f','o',0)
  do ne=1,nemisb
     read(ifnbiogen,*)ph,charspec
     ncstat=nf90_inq_varid(emisb_ncid,charspec,emisbvarid(ne))
     NCERR(__LINE__)
  end do
  close(ifnbiogen)
  
  cntvec=(/we,sn,1,1/)
  stvec=(/1,1,1,1/)
  do ne=1,nemisb
    ncstat=nf90_get_var( &
       emisb_ncid, &
       emisbvarid(ne), &
       emisb(ne,:,:,1), &
       stvec,cntvec)
  NCERR(__LINE__)
  enddo
  stvec=(/1,1,1,2/)
  do ne=1,nemisb
    ncstat=nf90_get_var( &
       emisb_ncid, &
       emisbvarid(ne), &
       emisb(ne,:,:,2), &
       stvec,cntvec)
  NCERR(__LINE__)
  enddo
  endif

end subroutine iniemis
