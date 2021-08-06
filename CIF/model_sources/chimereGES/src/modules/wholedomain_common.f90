! This file contains a module defining various workspaces.

module wholedomain_common

  use chimere_common
  implicit none


 ! IP Observations 
  real(kind=8), allocatable, dimension(:,:):: tabobs_glo
  integer,allocatable, dimension(:) :: nobs
  integer :: nobs_glo
  ! fin IP Observations

  ! Geometry
  integer :: nzonal_domain ! Total number of zonal cells
  integer :: nmerid_domain ! Total number of meridian cells
  integer :: nhbound_domain  ! Number of lateral boundary cells in one layer
  integer :: nlatbound_domain ! Toal number of lateral boundary cells
  integer :: nhoriz_domain ! Num. horiz. cells
  integer :: nverti  ! Number of vertical layers

  real(kind=8),dimension(:,:),allocatable         :: xlati   ! Latitudes
  real(kind=8),dimension(:,:),allocatable         :: xlong   ! Longitudes
  real(kind=8),dimension(:,:),allocatable         :: clati   ! Cosine of latitudes
  real(kind=8),dimension(:,:),allocatable         :: slati   ! Sine of latitudes
  real(kind=8),dimension(:,:),allocatable :: xsize   ! Zonal cell length
  real(kind=8),dimension(:,:),allocatable :: ysize   ! Meridional cell length
  real(kind=8),dimension(:,:),allocatable :: xbasx   ! x coord of local i vector
  real(kind=8),dimension(:,:),allocatable :: xbasy   ! y coord of local i vector
  real(kind=8),dimension(:,:),allocatable :: ybasx   ! x coord of local j vector
  real(kind=8),dimension(:,:),allocatable :: ybasy   ! y coord of local j vector

  real(kind=8),dimension(:),allocatable   :: avcoord ! Sigma coefficient
  real(kind=8),dimension(:),allocatable   :: bvcoord ! Sigma coefficient

  ! Outputs as NetCDF
  integer :: nsho    ! Cell for screen display (NSHO<=NHORIZ)

  ! Chemical workspace
  real(kind=8),dimension(:,:,:),allocatable :: conc   ! Current concentration array
  real(kind=8),dimension(:,:,:),allocatable :: conc_tl
  real(kind=8),dimension(:,:,:,:),allocatable :: aconc  
  real(kind=8),dimension(:,:,:),allocatable :: conco  ! Previous one (for TWOSTEP)
  real(kind=8),dimension(:,:,:),allocatable :: conco_tl
  real(kind=8),dimension(:,:,:,:),allocatable :: aconco  
  real(kind=8),dimension(:,:,:,:),allocatable :: concini  ! Initial concentration array
  real(kind=8),dimension(:,:,:,:),allocatable :: aconcini
  real(kind=8),dimension(:,:,:,:),allocatable :: aconcsave

  ! Boundary conditions
  real(kind=8),dimension(:,:,:),allocatable :: boundlat   ! Hourly lateral boundary concentrations
  real(kind=8),dimension(:,:,:),allocatable :: boundlat_tl
  real(kind=8),dimension(:,:,:),allocatable :: aboundlat 
  real(kind=8),dimension(:,:,:,:),allocatable :: boundtop   ! Hourly top boundary concentrations
  real(kind=8),dimension(:,:,:,:),allocatable :: boundtop_tl
  real(kind=8),dimension(:,:,:,:),allocatable :: aboundtop 

 ! Meteorology
  integer :: ideepconv ! type of deep convection
  real(kind=8),dimension(:,:,:,:),allocatable   :: winz    ! Hourly Zonal wind
  real(kind=8),dimension(:,:,:,:),allocatable   :: winm    ! Hourly Meridional wind
  real(kind=8),dimension(:,:,:,:),allocatable   :: temp    ! Hourly Temperature
  real(kind=8),dimension(:,:,:,:),allocatable   :: sphu    ! Hourly Specific humidity
  real(kind=8),dimension(:,:,:,:),allocatable   :: airm    ! Hourly density
  real(kind=8),dimension(:,:,:,:),allocatable   :: hlay    ! Hourly model layer heights
  real(kind=8),dimension(:,:,:,:),allocatable   :: kzzz    ! Hourly Kz
  real(kind=8),dimension(:,:,:,:),allocatable   :: clwc    ! Hourly Liquid (+ice) water content
  real(kind=8),dimension(:,:,:),allocatable          :: hght    ! Hourly mixing layer heights
  real(kind=8),dimension(:,:,:),allocatable          :: atte    !  Hourly radiation attenuation
  real(kind=8),dimension(:,:,:),allocatable          :: tem2    ! Hourly 2m-temperature
  real(kind=8),dimension(:,:,:),allocatable          :: usta    ! Hourly Ustar
  real(kind=8),dimension(:,:,:),allocatable          :: aerr    ! Hourly aerodynamic resistance
  real(kind=8),dimension(:,:,:),allocatable          :: obuk    ! Hourly Obukov length
  real(kind=8),dimension(:,:,:),allocatable          :: wsta    ! Hourly Wstar
  real(kind=8),dimension(:,:,:),allocatable          :: topc    ! Hourly precipitation
  real(kind=8),dimension(:,:,:),allocatable          :: sreh    ! Hourly surface relative humidity
  integer,dimension(:,:,:),allocatable  :: incloud
  ! lmbb deepconv
  real(kind=8),dimension(:,:,:,:),allocatable   :: dpeu    ! Entrainment in updraft (kg/m2/s)
  real(kind=8),dimension(:,:,:,:),allocatable   :: dped    ! Detrainment in updraft  (kg/m2/s)
  real(kind=8),dimension(:,:,:,:),allocatable   :: dpdu    ! Entrainment in downdraft (kg/m2/s)
  real(kind=8),dimension(:,:,:,:),allocatable   :: dpdd    ! Detrainment in downdraft (kg/m2/s)
  integer,dimension(:,:),allocatable         :: ideep   ! deep convection or not in a column

   ! Meteorology
  real(kind=8),dimension(:,:,:),allocatable   :: airmloc  ! Current density
  real(kind=8),dimension(:,:,:),allocatable    :: winzloc  ! Current zonal wind
  real(kind=8),dimension(:,:,:),allocatable    :: winmloc  ! Current meridional wind
  real(kind=8),dimension(:,:,:),allocatable    :: thlayloc ! Current layer thicknesses
  real(kind=8),dimension(:,:,:),allocatable    :: winvloc  ! Current vertical wind
  real(kind=8),dimension(:,:,:),allocatable    :: sphuloc  ! Current specific humidity
  real(kind=8),dimension(:,:,:),allocatable    :: temploc  ! Current temperature
  real(kind=8),dimension(:,:,:),allocatable    :: winxloc  ! Current mixing velocities
  real(kind=8),dimension(:,:,:),allocatable    :: hlayloc  ! Current model layer heights
  real(kind=8),dimension(:,:,:),allocatable    :: kzzzloc  ! Current top-layer Kz
  real(kind=8),dimension(:,:,:),allocatable    :: clwcloc
  real(kind=8),dimension(:,:,:),allocatable    :: presloc
  real(kind=8),dimension(:,:),allocatable    :: hghtloc ! Current mixing height
  real(kind=8),dimension(:,:),allocatable    :: atteloc ! Current radiation attenuation
  real(kind=8),dimension(:,:),allocatable    :: zeniloc ! Current zenith angle
  real(kind=8),dimension(:,:),allocatable    :: tem2loc ! Current 2m-temperature
  real(kind=8),dimension(:,:),allocatable    :: ustaloc ! Current ustar
  real(kind=8),dimension(:,:),allocatable    :: aerrloc ! Current aerodynamic resistance
  real(kind=8),dimension(:,:),allocatable    :: obukloc ! Current Obukov length
  real(kind=8),dimension(:,:),allocatable    :: wstaloc ! Current Wstar
  real(kind=8),dimension(:,:),allocatable    :: topcloc ! Current total precipitation
  real(kind=8),dimension(:,:,:),allocatable :: dtenloc

  ! Emissions workspace
  real(kind=8),dimension(:,:,:,:),allocatable :: emisaloc ! Current anthropic emissions
  real(kind=8),dimension(:,:,:,:),allocatable :: emisaloc_tl
  real(kind=8),dimension(:,:,:,:,:),allocatable  :: aemisaloc 
  real(kind=8),dimension(:,:,:,:),allocatable        :: emisb    ! Input biogenic emis. at exact hours
  real(kind=8),dimension(:,:,:,:),allocatable        :: emisb_tl
  real(kind=8),dimension(:,:,:,:),allocatable :: aemisb    
 
  ! Deposition and soil processes
  real(kind=8),dimension(:,:,:,:),allocatable :: fveg
  real(kind=8),dimension(:,:,:),allocatable :: drydep
  real(kind=8),dimension(:,:,:),allocatable :: drydep_tl
  real(kind=8),dimension(:,:,:),allocatable :: wetdep
  real(kind=8),dimension(:,:,:),allocatable :: wetdep_tl 
  real(kind=8),dimension(:,:,:),allocatable:: dland
  real(kind=8),dimension(:,:,:),allocatable :: depoloc


contains

  subroutine master_allocwhole
! Along the horizontal dimensions
    allocate(xlong(nzonal_domain,nmerid_domain))
    allocate(xlati(nzonal_domain,nmerid_domain))
    allocate(clati(nzonal_domain,nmerid_domain))
    allocate(slati(nzonal_domain,nmerid_domain))
    allocate(ideep(nzonal_domain,nmerid_domain))
    allocate(xsize(-1:nzonal_domain+2,-1:nmerid_domain+2))
    allocate(xbasx(0:nzonal_domain+1,0:nmerid_domain+1))
    allocate(xbasy(0:nzonal_domain+1,0:nmerid_domain+1))
    allocate(ysize(-1:nzonal_domain+2,-1:nmerid_domain+2))
    allocate(ybasx(0:nzonal_domain+1,0:nmerid_domain+1))
    allocate(ybasy(0:nzonal_domain+1,0:nmerid_domain+1))
    allocate(avcoord(nverti))
    allocate(bvcoord(nverti))    
    allocate(fveg(nzonal_domain,nmerid_domain,nvegtype,nlduse))
    allocate(emisaloc(nemisa,nzonal_domain,nmerid_domain,nlevemis))
    allocate(emisb(nemisb,nzonal_domain,nmerid_domain,2))
    allocate(dland(nzonal_domain,nmerid_domain,nlduse))
    allocate(tem2(nzonal_domain,nmerid_domain,2))
    allocate(sreh(nzonal_domain,nmerid_domain,2))
    allocate(atte(nzonal_domain,nmerid_domain,2))
    allocate(hght(nzonal_domain,nmerid_domain,2))
    allocate(usta(nzonal_domain,nmerid_domain,2))
    allocate(aerr(nzonal_domain,nmerid_domain,2))
    allocate(obuk(nzonal_domain,nmerid_domain,2))
    allocate(wsta(nzonal_domain,nmerid_domain,2))
    allocate(topc(nzonal_domain,nmerid_domain,2))
    allocate(winz(nzonal_domain,nmerid_domain,nverti,2))
    allocate(winm(nzonal_domain,nmerid_domain,nverti,2))
    allocate(temp(nzonal_domain,nmerid_domain,nverti,2))
    allocate(sphu(nzonal_domain,nmerid_domain,nverti,2))
    allocate(airm(nzonal_domain,nmerid_domain,nverti,2))
    allocate(clwc(nzonal_domain,nmerid_domain,nverti,2))
    allocate(dpeu(nzonal_domain,nmerid_domain,nverti,2))
    allocate(dped(nzonal_domain,nmerid_domain,nverti,2))
    allocate(dpdu(nzonal_domain,nmerid_domain,nverti,2))
    allocate(dpdd(nzonal_domain,nmerid_domain,nverti,2))
    allocate(hlay(nzonal_domain,nmerid_domain,nverti,2))
    allocate(kzzz(nzonal_domain,nmerid_domain,nverti,2))
    allocate(presloc(nzonal_domain,nmerid_domain,nverti))
    allocate(incloud(nzonal_domain,nmerid_domain,nverti))
    allocate(dtenloc(nzonal_domain,nmerid_domain,nverti))
    allocate(hghtloc(nzonal_domain,nmerid_domain))
    allocate(airmloc(0:nzonal_domain+1,0:nmerid_domain+1,1:nverti+1))
    allocate(atteloc(nzonal_domain,nmerid_domain))
    allocate(winvloc(nzonal_domain,nmerid_domain,nverti))
    allocate(winxloc(nzonal_domain,nmerid_domain,nverti))
    allocate(winzloc(nzonal_domain,nmerid_domain,nverti))
    allocate(winmloc(nzonal_domain,nmerid_domain,nverti))
    allocate(sphuloc(0:nzonal_domain+1,0:nmerid_domain+1,1:nverti+1))
    allocate(temploc(nzonal_domain,nmerid_domain,nverti))
    allocate(hlayloc(nzonal_domain,nmerid_domain,nverti))
    allocate(kzzzloc(nzonal_domain,nmerid_domain,nverti))
    allocate(clwcloc(nzonal_domain,nmerid_domain,nverti))
    allocate(thlayloc(nzonal_domain,nmerid_domain,nverti))
    allocate(tem2loc(nzonal_domain,nmerid_domain))
    allocate(ustaloc (nzonal_domain,nmerid_domain))
    allocate(aerrloc (nzonal_domain,nmerid_domain))
    allocate(obukloc (nzonal_domain,nmerid_domain))
    allocate(wstaloc (nzonal_domain,nmerid_domain))
    allocate(topcloc (nzonal_domain,nmerid_domain))
    
! For chemistry
    allocate(conc(-2:nzonal_domain+3,-2:nmerid_domain+3,1:nverti+1))
    allocate(conco(nzonal_domain,nmerid_domain,nverti))
    allocate(concini(nspectot,nzonal_domain,nmerid_domain,nverti))
    allocate(drydep(nspec,nzonal_domain,nmerid_domain))
    allocate(wetdep(nspec,nzonal_domain,nmerid_domain))
    allocate(zeniloc(nzonal_domain,nmerid_domain))
    allocate(boundlat(nlatbound_domain,nspec,2)) ! Hourly lateral boundary concentrations
    allocate(boundtop(nzonal_domain,nmerid_domain,nspec,2)) ! Hourly top boundary concentrations
    allocate(depoloc(nspec,nzonal_domain,nmerid_domain))

  end subroutine master_allocwhole

  subroutine master_allocwhole_tl
    allocate(emisaloc_tl(nemisa,nzonal_domain,nmerid_domain,nlevemis))
    allocate(emisb_tl(nemisb,nzonal_domain,nmerid_domain,2))
    allocate(conc_tl(-2:nzonal_domain+3,-2:nmerid_domain+3,1:nverti+1))
    allocate(conco_tl(nzonal_domain,nmerid_domain,nverti))
    allocate(drydep_tl(nspec,nzonal_domain,nmerid_domain))
    allocate(wetdep_tl(nspec,nzonal_domain,nmerid_domain))
    allocate(boundlat_tl(nlatbound_domain,nspec,2))
    allocate(boundtop_tl(nzonal_domain,nmerid_domain,nspec,2))
  end subroutine master_allocwhole_tl
  subroutine master_allocwhole_ad
    allocate(aconc(nspectot,-2:nzonal_domain+3,-2:nmerid_domain+3,1:nverti+1))
    allocate(aconco(nspectot,nzonal_domain,nmerid_domain,nverti))
    allocate(aconcini(nspectot,nzonal_domain,nmerid_domain,nverti))
    allocate(aconcsave(nspectot,nzonal_domain+(6*nzdoms),nmerid_domain+(6*nmdoms),1:nverti+1))

  end subroutine master_allocwhole_ad

end module wholedomain_common
