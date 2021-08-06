subroutine readhour_tl(ksens)
  !----------------------------------------------------------------

  use netcdf
  use chimere_consts
  use chimere_common
  use wholedomain_common
  implicit none

#define NCERR(lnum) if(ncstat/=NF90_NOERR) call nc_err(ncstat,lnum,'readhour_tl.f90')


  !****************************************************************************
  ! subroutine arguments
  integer :: ksens

  ! local variables
  integer :: ide,idr,idm,idex
  integer :: i,j,ne,nl,ns,ispec,ih,ik,ip,ivert
  integer :: izo,ime,ilev
  character(len=dlen) :: datebuf
  real(kind=8),allocatable,dimension(:,:,:)   :: buf3
  real(kind=8),allocatable,dimension(:,:,:)   :: buf4
  ! netCDF stuff
  integer :: ncstat
  integer,dimension(3) :: stvec3  ! start vectors for R/W functions
  integer,dimension(3) :: cntvec3 ! count vectors for R/W functions
  integer,dimension(4) :: stvec4  ! start vectors for R/W functions
  integer,dimension(4) :: cntvec4 ! count vectors for R/W functions

  ! functions
  integer :: mm5date2numeric

  !******************************************************************************
  !  Reading anthropic emissions
  ! Reminder : 
  !     records in the netCDF file are numbered starting with 1
  !       -> record 1 is the initial state
  !     run hours are numbered starting with O
  !       -> ihourrun==0 is the initial state

  !-------------------------------------------
  ncstat=nf90_get_var( &
          emisa_ncid, &
          emisa_times_varid, &
          datebuf, &
          (/1,ihourrun+1/),(/dlen,1/))
  NCERR(__LINE__)
  ide=mm5date2numeric(datebuf)
  idex = idate(ihourrun)
  if(ide.ne.idex) then
        print *,'*** ERROR: WRONG EXPEC. DATE IN AEMISSIONS FILE'
        print *,'IHOURRUN=',ihourrun,' EXPECTED: ' &
             ,idex,' AEMISSIONS: ',ide
        call exit1('Exiting')
  endif
  
  
  ! Gas emissions for each active anthropic species
  allocate(buf3(nzonal_domain,nmerid_domain,nlevemis))
  do ne=1,nemisa
        ncstat=nf90_get_var(                        &
             emisa_ncid,                            &
             emisavarid(ne),               &
             buf3,                                  &
             (/1,     1,     1,       ihourrun+1/), &
             (/nzonal_domain,nmerid_domain,nlevemis,1         /)  &
             )
        NCERR(__LINE__)
        do ilev=1,nlevemis
           do ime=1,nmerid_domain
              do izo=1,nzonal_domain
                 emisaloc(ne,izo,ime,ilev)=buf3(izo,ime,ilev)
              enddo
           enddo
        enddo
  enddo
  deallocate(buf3)
  
  
  ! READING INCREMENTS
  ncstat=nf90_get_var( &
          emisaincr_ncid, &
          emisaincr_times_varid, &
          datebuf, &
          (/1,ihourrun+1/),(/dlen,1/))
  NCERR(__LINE__)
  ide=mm5date2numeric(datebuf)
  idex = idate(ihourrun)
  if(ide.ne.idex) then
        print *,'*** ERROR: WRONG EXPEC. DATE IN AEMISSIONS INCREMENT FILE'
        print *,'IHOURRUN=',ihourrun,' EXPECTED: ' &
             ,idex,' AEMISSIONS INCR: ',ide
        call exit1('Exiting')
  endif
  ! Gas emissions for each active anthropic species
  allocate(buf4(nzonal_domain,nmerid_domain,nlevemis))
  do ne=1,nemisa
        ncstat=nf90_get_var(                        &
             emisaincr_ncid,                            &
             emisaincrvarid(ne),               &
             buf4,                                  &
             (/1,     1,     1,       ihourrun+1/), &
             (/nzonal_domain,nmerid_domain,nlevemis,1         /)  &
             )
        NCERR(__LINE__)
        do ilev=1,nlevemis
           do ime=1,nmerid_domain
              do izo=1,nzonal_domain
                 emisaloc_tl(ne,izo,ime,ilev)=buf4(izo,ime,ilev)
              enddo
           enddo
        enddo
  enddo
  deallocate(buf4)
  
   
  !-------------------------------------------
  !  Reading next-hour biogenic emissions
  if(optemisb.ne.0) then
  ncstat=nf90_get_var( &
       emisb_ncid, &
       emisb_times_varid, &
       datebuf, &
       (/1,ihourrun+2/),(/dlen,1/))
  NCERR(__LINE__)
  ide=mm5date2numeric(datebuf)
  idex = idate(ihourrun+1)
  if(ide.ne.idex) then
    print *,'*** ERROR: WRONG EXPEC. DATE IN BEMISSIONS FILE'
    print *,'IHOURRUN+1=',ihourrun+1,' EXPECTED: ' &
          ,idex,' BEMISSIONS: ',ide
    call exit1('Exiting')
  endif
  do ne=1,nemisb
    ncstat=nf90_get_var( &
       emisb_ncid, &
       emisbvarid(ne), &
       emisb(ne,:,:,ksens), &
       (/1,1,1,ihourrun+2/),(/nzonal_domain,nmerid_domain,1,1/))
NCERR(__LINE__)
  enddo
  ! Reading increments
  ncstat=nf90_get_var( &
       emisbincr_ncid, &
       emisbincr_times_varid, &
       datebuf, &
       (/1,ihourrun+2/),(/dlen,1/))
  NCERR(__LINE__)
  ide=mm5date2numeric(datebuf)
  idex = idate(ihourrun+1)
  if(ide.ne.idex) then
    print *,'*** ERROR: WRONG EXPEC. DATE IN BEMISSIONS INCREMENT FILE'
    print *,'IHOURRUN+1=',ihourrun+1,' EXPECTED: ' &
          ,idex,' INCRBEMISSIONS: ',ide
    call exit1('Exiting')
  endif
  do ne=1,nemisb
  ncstat=nf90_get_var( &
       emisbincr_ncid, &
       emisbincrvarid(ne), &
       emisb_tl(ne,:,:,ksens), &
       (/1,1,1,ihourrun+2/),(/nzonal_domain,nmerid_domain,1,1/))
  NCERR(__LINE__)
  enddo
  
  endif

  !  Reading next-hour meteo parameters

  ncstat=nf90_get_var( &
       meteo_ncid, &
       meteo_times_varid, &
       datebuf, &
       (/1,ihourrun+2/),(/dlen,1/))
  NCERR(__LINE__)
  idm=mm5date2numeric(datebuf)
  idex = idate(ihourrun+1)
  if(idm.ne.idex) then
    print *,'*** ERROR: WRONG EXPECTED DATE IN METEO FILE'
    print *,'IHOURRUN+1=',ihourrun+1,' EXPECTED: ' &
          &         ,idex,' METEO: ',idm
    call exit1('Exiting')
  endif

  ! 2D fields
  stvec3=(/1,1,ihourrun+2/)
  cntvec3=(/nzonal_domain,nmerid_domain,1/)

  ncstat=nf90_get_var(meteo_ncid,meteo_tem2_varid,tem2(:,:,ksens),stvec3,cntvec3)
  NCERR(__LINE__)
  ncstat=nf90_get_var(meteo_ncid,meteo_sreh_varid,sreh(:,:,ksens),stvec3,cntvec3)
  NCERR(__LINE__)
  ncstat=nf90_get_var(meteo_ncid,meteo_atte_varid,atte(:,:,ksens),stvec3,cntvec3)
  NCERR(__LINE__)
  ncstat=nf90_get_var(meteo_ncid,meteo_hght_varid,hght(:,:,ksens),stvec3,cntvec3)
  NCERR(__LINE__)
  ncstat=nf90_get_var(meteo_ncid,meteo_usta_varid,usta(:,:,ksens),stvec3,cntvec3)
  NCERR(__LINE__)
  ncstat=nf90_get_var(meteo_ncid,meteo_aerr_varid,aerr(:,:,ksens),stvec3,cntvec3)
  NCERR(__LINE__)
  ncstat=nf90_get_var(meteo_ncid,meteo_obuk_varid,obuk(:,:,ksens),stvec3,cntvec3)
  NCERR(__LINE__)
  ncstat=nf90_get_var(meteo_ncid,meteo_wsta_varid,wsta(:,:,ksens),stvec3,cntvec3)
  NCERR(__LINE__)
  ncstat=nf90_get_var(meteo_ncid,meteo_topc_varid,topc(:,:,ksens),stvec3,cntvec3)
  NCERR(__LINE__)
  ! 3D fields
  stvec4=(/1,1,1,ihourrun+2/)
  cntvec4=(/nzonal_domain,nmerid_domain,nverti,1/)

  ncstat=nf90_get_var(meteo_ncid,meteo_winz_varid,winz(:,:,:,ksens),stvec4,cntvec4)
  NCERR(__LINE__)
  ncstat=nf90_get_var(meteo_ncid,meteo_winm_varid,winm(:,:,:,ksens),stvec4,cntvec4)
  NCERR(__LINE__)
  ncstat=nf90_get_var(meteo_ncid,meteo_temp_varid,temp(:,:,:,ksens),stvec4,cntvec4)
  NCERR(__LINE__)
  ncstat=nf90_get_var(meteo_ncid,meteo_sphu_varid,sphu(:,:,:,ksens),stvec4,cntvec4)
  NCERR(__LINE__)
  ncstat=nf90_get_var(meteo_ncid,meteo_airm_varid,airm(:,:,:,ksens),stvec4,cntvec4)
  NCERR(__LINE__)
  ncstat=nf90_get_var(meteo_ncid,meteo_clwc_varid,clwc(:,:,:,ksens),stvec4,cntvec4)
  NCERR(__LINE__)
  ! lmbb deep conv
  ncstat=nf90_get_var(meteo_ncid,meteo_dpeu_varid,dpeu(:,:,:,ksens),stvec4,cntvec4)
  NCERR(__LINE__)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  dpeu(4,3,1,ksens)=7.5d-6
!  dpeu(4,3,2,ksens)=1.5d-5
!  dpeu(4,3,3,ksens)=1.8d-5
!  dpeu(4,3,4,ksens)=2.8d-6
!  dpeu(4,3,5,ksens)=4.5d-8
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
  ncstat=nf90_get_var(meteo_ncid,meteo_dped_varid,dped(:,:,:,ksens),stvec4,cntvec4)
  NCERR(__LINE__)
  ncstat=nf90_get_var(meteo_ncid,meteo_dpdu_varid,dpdu(:,:,:,ksens),stvec4,cntvec4)
  NCERR(__LINE__)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  dpdu(4,3,3,ksens)=3.7d-8
!  dpdu(4,3,4,ksens)=6.2d-6
!  dpdu(4,3,5,ksens)=4.d-5
!  dpdu(4,3,6,ksens)=1.1d-5
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
  ncstat=nf90_get_var(meteo_ncid,meteo_dpdd_varid,dpdd(:,:,:,ksens),stvec4,cntvec4)
  NCERR(__LINE__)

  ncstat=nf90_get_var(meteo_ncid,meteo_alti_varid,hlay(:,:,:,ksens),stvec4,cntvec4)
  NCERR(__LINE__)
  ncstat=nf90_get_var(meteo_ncid,meteo_kzzz_varid,kzzz(:,:,:,ksens),stvec4,cntvec4)
  NCERR(__LINE__)

  !  Change of units for 3D variables to obtain molec, cm, s

  winz(:,:,:,ksens) = winz(:,:,:,ksens)*1d2
  winm(:,:,:,ksens) = winm(:,:,:,ksens)*1d2
  sphu(:,:,:,ksens) = sphu(:,:,:,ksens)*airm(:,:,:,ksens)*1.6
  hlay(:,:,:,ksens) = hlay(:,:,:,ksens)*1d2
  clwc(:,:,:,ksens) = clwc(:,:,:,ksens)*airm(:,:,:,ksens)*29./an
  kzzz(:,:,:,ksens) = kzzz(:,:,:,ksens)*1d4

  ! lmbb deep convection: * 1e-1 factor for conversion from kg/m2/s into g/cm2/s 
  dpeu(:,:,:,ksens) = dpeu(:,:,:,ksens) * an/29.*1.d-1
  dped(:,:,:,ksens) = dped(:,:,:,ksens) * an/29.*1.d-1
  dpdu(:,:,:,ksens) = dpdu(:,:,:,ksens) * an/29.*1.d-1
  dpdd(:,:,:,ksens) = dpdd(:,:,:,ksens) * an/29.*1.d-1

  !  Change of units for 2D variables

  hght(:,:,ksens) = hght(:,:,ksens)*1d2
  usta(:,:,ksens) = usta(:,:,ksens)*1d2
  aerr(:,:,ksens) = aerr(:,:,ksens)*1d-2
  obuk(:,:,ksens) = obuk(:,:,ksens)*1d2
  wsta(:,:,ksens) = wsta(:,:,ksens)*1d2
  where (topc(:,:,ksens).lt.dzero) topc(:,:,ksens) = dzero

  !  Reading next-hour boundary conditions
  ncstat=nf90_get_var( &
       bounconc_ncid, bounconc_times_varid, datebuf, &
       (/1,ihourrun+2/),(/dlen,1/))
  NCERR(__LINE__)
  idr=mm5date2numeric(datebuf)
  idex = idate(ihourrun+1)
  if(idr.ne.idex) then
    print *,'*** ERROR: WRONG EXPECTED DATE IN BOUN_CONCS FILE'
    print *,'IHOURRUN+1=',ihourrun+1,'  EXPECTED= ' &
	  ,idex,' BOUN_CONCS= ',idr
    call exit1('Exiting')
  endif
  allocate(buf3(nspecboun, nhbound_domain, nverti))
  stvec4 =(/1	    , 1      , 1     , ihourrun+2 /)
  cntvec4=(/nspecboun, nhbound_domain, nverti, 1	   /)
  ncstat=nf90_get_var(bounconc_ncid, latconc_conc_varid, buf3, stvec4, cntvec4)
  NCERR(__LINE__)

  nl=0
  do ivert=1,nverti
    do ih=1,nhbound_domain
	nl=nl+1
	do ispec=1,nspecboun
	   if(isboun(ispec).gt.0) then
	      boundlat(nl,isboun(ispec),ksens) = buf3(ispec, ih, ivert)
	   end if
	end do
    end do
  end do
  deallocate(buf3)

  !  Reading next-hour top boundary conditions

  allocate(buf3(nspecboun,nzonal_domain,nmerid_domain))
  stvec4 =(/1	    , 1     , 1,      ihourrun+2 /)
  cntvec4=(/nspecboun, nzonal_domain, nmerid_domain, 1	  /)
  ncstat=nf90_get_var(bounconc_ncid,topconc_conc_varid, buf3, stvec4, cntvec4)
  NCERR(__LINE__)
  do ime=1,nmerid_domain
    do izo=1,nzonal_domain
	do ispec=1,nspecboun
	   if(isboun(ispec).gt.0) then
	      boundtop(izo,ime,isboun(ispec),ksens) = buf3(ispec,izo,ime)
	   end if
	end do
    end do
  end do
  deallocate(buf3)
 !  Reading next-hour boundary condition increments
      ncstat=nf90_get_var( &
          bounconcincr_ncid, bounconcincr_times_varid, datebuf, &
          (/1,ihourrun+2/),(/dlen,1/))
      NCERR(__LINE__)
      idr=mm5date2numeric(datebuf)
      idex = idate(ihourrun+1)
      if(idr.ne.idex) then
        print *,'*** ERROR: WRONG EXPECTED DATE IN BOUN_CONCS INCREMENT FILE'
        print *,'IHOURRUN+1=',ihourrun+1,'  EXPECTED= ' &
             ,idex,' BOUN_CONCS= ',idr
        call exit1('Exiting')
      endif
      allocate(buf3(nspecboun, nhbound_domain, nverti))
      stvec4 =(/1       , 1      , 1     , ihourrun+2 /)
      cntvec4=(/nspecboun, nhbound_domain, nverti, 1          /)
      ncstat=nf90_get_var(bounconcincr_ncid, latconcincr_conc_varid, buf3, stvec4, cntvec4)
      NCERR(__LINE__)

      nl=0
      do ivert=1,nverti
        do ih=1,nhbound_domain
           nl=nl+1
           do ispec=1,nspecboun
              if(isboun(ispec).gt.0) then
                 boundlat_tl(nl,isboun(ispec),ksens) = buf3(ispec, ih, ivert)
              end if
           end do
        end do
      end do
      deallocate(buf3)

      !  Reading next-hour top boundary condition increments

      allocate(buf3(nspecboun,nzonal_domain,nmerid_domain))
      stvec4 =(/1       , 1     , 1,      ihourrun+2 /)
      cntvec4=(/nspecboun, nzonal_domain, nmerid_domain, 1          /)
      ncstat=nf90_get_var(bounconcincr_ncid,topconcincr_conc_varid, buf3, stvec4, cntvec4)
      NCERR(__LINE__)
      do ime=1,nmerid_domain
        do izo=1,nzonal_domain
           do ispec=1,nspecboun
              if(isboun(ispec).gt.0) then
                 boundtop_tl(izo,ime,isboun(ispec),ksens) = buf3(ispec,izo,ime)
              end if
           end do
        end do
      end do
      deallocate(buf3)

end subroutine readhour_tl
