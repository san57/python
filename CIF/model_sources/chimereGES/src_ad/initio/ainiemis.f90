subroutine ainiemis

  !  Initialization of emissions

  use netcdf
  use chimere_consts
  use chimere_common
  use wholedomain_common
  
  implicit none

#define NCERR(lnum) if(ncstat/=NF90_NOERR) call nc_err(ncstat,lnum,'iniemis.f90')

  !********************************************************************************
  integer :: ide,idex

  character(len=dlen) :: datebuf     ! to hold current date in MM5 format

  ! netcdf stuff
  integer :: ncstat                  ! return code for netCDF functions
  integer,dimension(4) :: stvec   ! start vectors for R/W functions
  integer,dimension(4) :: cntvec  ! count vectors for R/W functions
  integer :: ifnbiogen,ne,we
  character(len=15) :: charspec

  ! External functions
  integer :: mm5date2numeric

   character*15 cspec
   character*100000 usedval,notusedval

!======================================================================
  ! Anthropic emission file (netCDF) is already open

  !! before-last hour
  ncstat=nf90_get_var( &
       emisa_ncid, &
       emisa_times_varid, &
       datebuf, &
       (/1,nhourrun/),(/dlen,1/))
  NCERR(__LINE__)
  idex=idate(nhourrun-1)
  ide=mm5date2numeric(datebuf)
  if(ide.ne.idex) then
     print *,'*** ERROR: WRONG EXPECTED DATE IN AEMISSIONS FILE'
     print *,'IHOURRUN=nhourrun-1,  EXPECTED= ',idex,' AEMISSIONS= ',ide
     stop
  endif
  !! next hour
  ncstat=nf90_get_var( &
       emisa_ncid, &
       emisa_times_varid, &
       datebuf, &
       (/1,nhourrun+1/),(/dlen,1/))
  NCERR(__LINE__)
  idex = idate(nhourrun)
  ide=mm5date2numeric(datebuf)
  if(ide.ne.idex) then
     print *,'*** ERROR: WRONG EXPECTED DATE IN AEMISSIONS FILE'
     print *,'IHOURRUN=nhourrun,  EXPECTED= ',idex,' AEMISSIONS= ',ide
     stop
  endif


!============================================================================
  if(optemisb.ne.0) then
    !  Reading biogenic emissions
    cntvec=(/nzonal_domain,nmerid_domain,1,1/)
    ! Getting the varids of each species in BEMISSIONS file
    call opfi(ifnbiogen,fnbiogen,'f','o',0)
    do ne=1,nemisb
      read(ifnbiogen,*)we,charspec
      ncstat=nf90_inq_varid(emisb_ncid,charspec,emisbvarid(ne))
      NCERR(__LINE__)
    end do
    close(ifnbiogen)
    !! before-last hour
    ncstat=nf90_get_var( &
       emisb_ncid, &
       emisb_times_varid, &
       datebuf, &
       (/1,nhourrun/),(/dlen,1/))
NCERR(__LINE__)
    idex=idate(nhourrun-1)
    ide=mm5date2numeric(datebuf)
    if(ide.ne.idex) then
      print *,'*** ERROR: WRONG EXPECTED DATE IN BEMISSIONS FILE'
      print *,'IHOURRUN=nhourrun-1,  EXPECTED= ',idex,' BEMISSIONS= ',ide
      stop
    endif
    stvec=(/1,1,1,nhourrun/)
    do ne=1,nemisb
      ncstat=nf90_get_var( &
       emisb_ncid, &
       emisbvarid(ne), &
       emisb(ne,:,:,1), &
       stvec,cntvec)
      NCERR(__LINE__)
    enddo
    !! next hour
    ncstat=nf90_get_var( &
       emisb_ncid, &
       emisb_times_varid, &
       datebuf, &
       (/1,nhourrun+1/),(/dlen,1/))
NCERR(__LINE__)
    idex = idate(nhourrun)
    ide=mm5date2numeric(datebuf)
    if(ide.ne.idex) then
      print *,'*** ERROR: WRONG EXPECTED DATE IN BEMISSIONS FILE'
      print *,'IHOURRUN=nhourrun,  EXPECTED= ',idex,' BEMISSIONS= ',ide
      stop
    endif
  
    stvec=(/1,1,1,nhourrun+1/)
    do ne=1,nemisb
      ncstat=nf90_get_var( &
       emisb_ncid, &
       emisbvarid(ne), &
       emisb(ne,:,:,2), &
       stvec,cntvec)
      NCERR(__LINE__)
    enddo
   
  endif

end subroutine ainiemis
