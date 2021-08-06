!-*-f90-*-
subroutine iniread_tl

  !  Reading simulation parameters and file names.
  !  This is the first initialization routine
  !  INPUT :  ---
  !  OUTPUT:  NDAYS      Number of simulated days (must be less than NDAYSMAX)
  !           IDATESTART A YYYYMMDDHH date for simulation start (HH=00 mandatory)
  !           FILE NAMES All input file names.

  use chimere_consts
  use chimere_common
  use wholedomain_common
  use netcdf

  implicit none

  !*****************************************************************************************
  integer :: ifn_args,ierr
  character(len=*),parameter :: fniARGS = 'chimere.nml'

  namelist /args/  &
       version,    &
       domain,   &
       nphour_ref, &
       ichemstep ,  &
       nzonal_domain, &
       nmerid_domain, &
       nvert_raw , &
       nlevemis , &
       nspec , &
       nemisa, &
       nemisb, &
       ndep,    &
       nreac , &
       nfam , &
       nspresc , &
       nreactamax , &
       ntemps , &
       ntabmax , &
       nlevphotmax, &
       ntabuzenmax , &
       nphotmax , &
       ivsurf , &
       ideepconv ,&
       idatestart, &
       nhourrun ,   &
       nzdoms , &
       nmdoms , &
       fninit,     &
       fnoutspec,  &
       fnspec,     &
       fnchem,     &
       fnstoi,     &
       fnrates,    &
       fnfamilies, &
       fnphot,     &
       fnanthro,   &
       fnbiogen,   &
       fnemisa,    &
!       fnemisb,    &
       fnbounconc, &
       fnmeteo,    &
       fndepoespe, &
       fndepopars, &
       fnwetd,     &
       fnlanduse,  &
       fnout,      &
       fnothers,   &
       fnconcs,    &
       fndepos,    &
       fniCOOcorn , &
       usechemistry , &
       usedepos , &
       useemissions , &
       usetransmix , &
       usewetdepos , &
       useabsclipconc, &
       dryairout, &
       nvegtype , &
       nlduse, &
       nparammax ,&
       clipconc ,&
       ntyperate, &
       ihoursu, &
       nivout, &
       nsaveconcs, &
       nsavedepos, &
       optemisb, &
       nitgs , &
       nitgssu , &
       fnemisaincr,&
       fnemisbincr, &
       fnbounconcincr,&
       fninitincr, &
!       fnipol, &
!       polar,  &
       NetCDF_output
! 
  !*****************************************************************************************

  call opfi(ifn_args,fniARGS,'f','o')
  
  read(ifn_args,nml=args,iostat=ierr)

  if (ierr/=0) then
    call exit1('iniread_tl: problem reading arguments namelist')
  end if

  !  Check point
  nhoriz_domain=nmerid_domain*nzonal_domain
  nsho=nhoriz_domain/2
  if(nsho.gt.nhoriz_domain) then
    print *,'*** ERROR: NSHO (in CHIMERE.H) MUST BE <=NHORIZ'
    print *,'           NHORIZ is the number of cells per layer'
    call exit1('Exiting')
  endif

  nspectot = nspec + nfam + nspresc
  nverti    = nvert_raw+ivsurf-1
  nhbound_domain = 2 * (nmerid_domain + nzonal_domain) ! Number of lateral boundary cells in one layer
  nlatbound_domain = nhbound_domain*nverti     ! Number of lateral boundary cells


end subroutine iniread_tl
