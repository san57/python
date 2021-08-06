! This file contains a module defining various workspaces.

module chimere_common
  use chimere_consts
  use message_defs

  implicit none

      ! Chimere types
! A structure to hold informations about the fields to write to file "out"
  type :: output_species_type
    integer              :: varid          ! netcdf var id
    integer              :: varid_tl          ! netcdf var id  for TL increment
    integer              :: iaddr          ! address in species(:) table
    character(len=splen) :: name           ! name from species(:) table
    character(len=16)    :: units
  end type output_species_type
  type(output_species_type),dimension(:),allocatable :: output_species

   ! A structure to hold information about active species
  ! lmbb add a new var in module species_type
  type :: species_type
    integer           :: varid
    character(len=16) :: name
    integer           :: transp
    integer           :: transpv
    integer           :: bounddry ! Boundary conditions are in dry mole if 1
    integer           :: varido ! the matching conco
    integer :: varid_tl ! IP for matching TL species
    integer  :: varido_tl ! IP for the conco of the TL species
    real(kind=8)         :: fspec          ! conversion factor
  end type species_type
  type(species_type),dimension(:),allocatable :: species

  ! A structure to hold information about the fields to write to file "par"
  type :: output_par_type
    integer           :: varid
    character(len=16) :: name
    character(len=32) :: long_name
    character(len=16) :: units
    logical           :: printit
  end type output_par_type
  type(output_par_type),dimension(:),allocatable :: output_par



  ! A structure to hold information about the fields to write to file "depo"
  type :: output_depo_type
    integer           :: varid
    character(len=20) :: name
  end type output_depo_type
  type(output_depo_type),dimension(:),allocatable :: output_depo

     ! Outputs as NetCDF
  integer :: NetCDF_output ! 1 if out.nc 0 if None
  integer :: NetCDF_parout ! 1 if par.nc 0 if None
  integer :: NetCDF_outtype      ! Type of outputs in netcdf
  integer :: nparammax   ! Max number of output parameters
  integer :: nivout  ! Number of vertical layers in output files
  integer :: nsaveconcs ! Save concentrations for restart every ... hours
  integer :: nsavedepos ! Save cumulated deposition every ... hours

 ! a structure to hold information about the aemisa fields
  type :: emisa_species_type
    integer :: varid
    character(len=splen) :: name
  end type emisa_species_type
  type(emisa_species_type),dimension(:),allocatable :: emisa_species

  ! a structure to hold information about the aemisb fields
  type :: emisb_species_type
    integer :: varid
    character(len=splen) :: name
  end type emisb_species_type
  type(emisb_species_type),dimension(:),allocatable :: emisb_species
  
! a structure to hold information about the aconc fields
  type :: aconc_species_type
    integer :: varid
    character(len=splen) :: name
  end type aconc_species_type
  type(aconc_species_type),dimension(:),allocatable :: aconc_species

  ! a structure to hold information about the aboundconc fields
  type :: aboun_species_type
    integer :: varid
    character(len=splen) :: name
  end type aboun_species_type
  type(aboun_species_type),dimension(:),allocatable :: aboun_species

  ! The following structure is used by MPI
  type :: dom_type
    integer :: izstart    ! index of the first cell of this subdomain along zonal axis
    integer :: izend      ! index of the last cell
    integer :: nzcount    ! number of cells along zonal axis
    integer :: imstart    ! index of the first cell of this subdomain along meridian axis
    integer :: imend      ! index of the last cell
    integer :: nmcount    ! number of cells along meridian axis
    integer :: i          ! subdomain zonal index 
    integer :: j          ! subdomain meridian index
  end type dom_type
  type(dom_type),dimension(:),allocatable :: dom

    ! MPI
  integer :: fwd
  integer :: tl
  integer :: ad
  integer :: nmdoms,nzdoms ! Number of MPI subdomains in meridian and zonal directions
  integer :: ndoms   ! Total number of subdomains nzdoms*nmdoms
  ! MPI
  integer :: wrk_comm                        ! Communicator local to workers
  integer :: rank                            ! process rank
  integer :: nzonalmax,nmeridmax

    ! Model version
  character(len=240) :: version
  integer :: usechemistry,usedepos,usetransmix,usewetdepos,useemissions,useabsclipconc, dryairout
  real(kind=8) :: clipconc = 1d0         ! Clipping value for the TWOSTEP algorithm

   ! Timing
  integer :: nphour_ref ! Number of physical steps per hour
  integer :: ichemstep ! Number of chemical refined iterations
  integer :: nitgs, nitgssu ! Numbers of Gauss-Seidel iterations
  integer      :: ihourrun                   ! Hour counter
  integer      :: nhourrun                   ! Simulation duration (Hours)
  integer      :: idatestart                 ! Starting date in YYYYMMDDHH format
  real(kind=8) :: soltim                     ! Time interval with year's solstice
  real(kind=8) :: djul                       ! Julian day
  real(kind=8) :: thour                      ! Time (hour) since last exact hour
  integer :: ihoursu  ! Number of spin-up hours
  integer,dimension(:),allocatable :: idate  ! Dates along the run
  integer,dimension(:),allocatable :: ihour  ! UT hours
  integer,dimension(:),allocatable :: iday   ! Day in the month
  integer,dimension(:),allocatable :: imonth ! Month
  integer,dimension(:),allocatable :: iyear  ! Year
  integer,dimension(:),allocatable :: idtyp  ! Day types (week, sat, sun)

  ! adaptative time-step
  integer,dimension(:),allocatable :: nphourm
  integer :: nphour
  real(kind=8) :: dtr
  real(kind=8) :: dtr2

   ! Geometry
  character(len=20) :: domain
  integer :: nvert_raw ! Raw nb of vertical levels, before insert. of sublayer
  integer :: ivsurf  ! Emission sublayer (2) or not (1)
  real(kind=8) :: psurf ! surface pressure in Pa, from VCOORD (b coeff)/METEO


  ! Chemical workspace
  integer :: nreac ! Number of reactions
  integer :: ntyperate ! Max number of reaction types
  integer :: nreactamax ! Max number of reactants/reaction
  !XX faire pareil pour les autres: les enlever d'ici et les mettre dans le yml, section utilisateurs confirmes
  ! avec valeur par defaut = celle qui etait la avant
  integer :: nlevphotmax ! Max number of tabulated photolysis levels
  integer :: ntabuzenmax ! Max number of tabulated zenith angles
  integer :: nphotmax  ! Max number of photolysis reactions
  integer :: nspectot ! Total num. species
  integer :: ntemps ! Number of tabul. temperatures for stoichio.
  integer :: ntabmax ! Max number of rate constants
  integer :: nspresc     ! Number of prescribed species (O2, H2O...)
  integer :: ndep ! Number of deposited species
  integer :: nfam          ! Number of "family" species (for output)

  real(kind=8),dimension(:,:,:),allocatable :: stoi     ! Stoichiometric coeffs.
  real(kind=8),dimension(:),allocatable     :: zetaref  ! Tabulated cos of zenith angles
  real(kind=8),dimension(:,:),allocatable   :: tabrate  ! Reaction rate constants
  real(kind=8),dimension(:),allocatable     :: tabtemp  ! Tabulated temperatures for stoic.
  real(kind=8),dimension(:,:),allocatable   :: dtabrate ! Their slope
  real(kind=8),dimension(:),allocatable     :: altiphot ! Altitudes of photolysis levels
  real(kind=8),dimension(:,:,:),allocatable :: photoj   ! Clear sky photolysis rates
  real(kind=8)                              :: actflx

  integer,dimension(:),allocatable  :: nreactants  ! Number of reactants
  integer,dimension(:),allocatable  :: ityperate   ! Reaction rate types
  integer,dimension(:),allocatable  :: iphoto	   ! Photolysis reaction address
  integer,dimension(:,:),allocatable :: irctt       ! Address of reactants
  integer,dimension(:),allocatable    :: kreacp      ! Number of reactions producing
  integer,dimension(:),allocatable    :: kreacl      ! Number of reactions destructing
  integer,dimension(:,:),allocatable  :: ireacp      ! Addresses of products
  integer,dimension(:,:),allocatable  :: ireacl      ! Addresses of reactants
  integer,dimension(:,:),allocatable  :: ifam        ! Addresses of elements in family
  integer,dimension(:),allocatable    :: ltabrate    ! Their number of constants
  integer,dimension(:),allocatable    :: nelem       ! Num. of elements in families

  integer :: ntabuzen    ! Number of tabulated zenith angles
  integer :: nphot       ! Number of photolysis reactions
  integer :: nlevphot    ! Number of tabulated photolysis levels
  integer :: inNH3       ! NH3  species number
  integer :: inSO2       ! SO2  species number
  integer :: inHNO3      ! HNO3 species number
  integer :: noutspec    ! Number of output species! chemical workspace
  integer :: nspec ! Number of active species


  ! Emissions workspace
  integer :: nlevemis ! Number of emission levels
  integer :: hpulse ! Hour pulse Elise Potier
  integer :: nemisa ! Max number of anthropic emitted species
  integer :: nemisb ! Max number of biogenic emitted species
  integer :: optemisb

  integer,dimension(1000)    :: ispecemip
  integer,dimension(:),allocatable   :: inemisa  ! Addresses of anthropic species
  integer,dimension(:),allocatable   :: inemisb  ! Addresses of biogenic species


  ! Boundary conditions

  integer,allocatable,dimension(:)  :: isboun     ! Addresses of boundary species
  integer                           :: iopboun    ! Time-dependent (1), fixed (0)  bounds
  integer                           :: nspecboun  ! Number of read species at bounds
  integer                           :: nspecbounout  ! Number of out species at bounds

  ! Deposition and soil processes
  integer :: nvegtype
  integer :: nlduse
  integer,dimension(:,:),allocatable           :: inwetd
  real(kind=8),dimension(:),allocatable     :: gmax
  real(kind=8),dimension(:),allocatable     :: fmin
  real(kind=8),dimension(:),allocatable     :: deptmin
  real(kind=8),dimension(:),allocatable     :: deptopt
  real(kind=8),dimension(:),allocatable     :: deptmax
  real(kind=8),dimension(:),allocatable     :: depalph
  real(kind=8),dimension(:),allocatable     :: depvpd1
  real(kind=8),dimension(:),allocatable     :: depvpd2
  real(kind=8),dimension(:),allocatable     :: depsgs
  real(kind=8),dimension(:),allocatable     :: depegs
  real(kind=8),dimension(:),allocatable     :: depsgl
  real(kind=8),dimension(:),allocatable     :: depegl
  real(kind=8),dimension(:),allocatable     :: deplai1
  real(kind=8),dimension(:),allocatable     :: deplai2
  real(kind=8),dimension(:),allocatable     :: depphe0
  real(kind=8),dimension(:),allocatable     :: depphe1
  real(kind=8),dimension(:),allocatable     :: depphe2
  real(kind=8),dimension(:),allocatable     :: zcanopy
  real(kind=8),dimension(:),allocatable     :: RGSO3
  real(kind=8),dimension(:),allocatable     :: RGSSO2
  real(kind=8),dimension(:),allocatable     :: so2rh
  real(kind=8),dimension(:),allocatable        :: dHx
  real(kind=8),dimension(:),allocatable        :: df0
  real(kind=8),dimension(:),allocatable        :: factRb
  real(kind=8),dimension(:),allocatable        :: factD
  real(kind=8),dimension(:),allocatable        :: Rm

  integer,dimension(:),allocatable             :: i0L
  integer,dimension(2)                 :: nwetd
  integer,dimension(:),allocatable           :: indepo
  integer,dimension(:),allocatable           :: inexcha
  integer,dimension(:),allocatable           :: iexcha

  integer                              :: ndepo


  ! File names

  character(len=256) :: fninit      ! Initial conditions
  character(len=256) :: fninitincr  ! Increments on Initial conditions 
  character(len=256) :: fnout       ! Hourly output concentrations
  character(len=256) :: fnexc       ! Hourly output concentrations in soil, vegetation & water
  character(len=256) :: fnothers    ! Hourly output other variables
  character(len=256) :: fnconcs     ! Final backup of all concentrations
  character(len=256) :: fnendex     ! Final backup of concentrations in soil, vegetation and sea
  character(len=256) :: fndepos     ! Backup of cumulated deposition
  character(len=256) :: fnoutspec   ! Names of output species
  character(len=256) :: fnspec      ! Names of active species
  character(len=256) :: fnchem      ! Chemical reactions
  character(len=256) :: fnstoi      ! Stoichiometric coefficients
  character(len=256) :: fnrates     ! Reaction rates
  character(len=256) :: fnfamilies  ! Families list and contents
  character(len=256) :: fnanthro    ! Anthorpic emitted species
  character(len=256) :: fnbiogen    ! Biogenic emitted species
  character(len=256) :: fnemisa     ! Anthropic emissions file
  character(len=256) :: fnemisaincr    ! Anthorpic emission increments file
  character(len=256) :: fnemisb     ! Biogenic emissions file
  character(len=256) :: fnemisbincr ! Biogenic emission increments file
  character(len=256) :: fnmeteo     ! Hourly meteo parameters
  character(len=256) :: fndepoespe  ! Names of deposited species and parameterse
  character(len=256) :: fndepopars  ! Vegetation dependent dry deposition params
  character(len=256) :: fnwetd      ! Wet deposited species and parameters
  character(len=256) :: fnlanduse   ! Land use map
  character(len=256) :: fnbounspec  ! Names of species at the lateral boundarie
  character(len=256) :: fnbounconc  ! Concentrations at the lateral boundaries
  character(len=256) :: fnbounconcincr ! Increments of Concentrations at the lateral boundaries
  character(len=256) :: fnphot      ! Photolysis parameters file
  character(len=256) :: fniCOOcorn  ! Cell corner coordinate file
  character(len=256) :: fnipol  ! Cell corner coordinate file
  character(len=132) :: fnainit
  integer  :: polar ! pole in the domain
  character(len=256) :: afnout ! template nmae for aout.nc files
   
  ! files units and ids
  ! binary files
  integer :: ifnconcs
  integer :: ifndepos
  integer :: iopainit
  integer :: ifnainit ! initial aconcentrations
 


  ! netcdf files
  integer :: emisa_ncid
  integer :: emisaincr_ncid
  integer :: emisa_times_varid
  integer :: emisaincr_times_varid
  integer,dimension(:),allocatable :: emisavarid ! varids of species in AEMISSIONS file
  integer,dimension(:),allocatable :: emisaincrvarid ! varids of species in AEMISSIONS INCREMENT  file
  integer,dimension(:),allocatable :: emisbvarid ! varids of species in BEMISSIONS file
  integer,dimension(:),allocatable :: emisbincrvarid ! varids of species in BEMISSIONS INCREMENT  file

  integer :: emisb_ncid
  integer :: emisbincr_ncid
  integer :: emisb_times_varid
  integer :: emisbincr_times_varid
  integer :: emisb_emisb_varid

  integer :: meteo_ncid
  integer :: meteo_times_varid
  integer :: meteo_lon_varid
  integer :: meteo_lat_varid
  integer :: meteo_avcoord_varid
  integer :: meteo_bvcoord_varid
  integer :: meteo_tem2_varid
  integer :: meteo_sreh_varid
  integer :: meteo_atte_varid
  integer :: meteo_hght_varid
  integer :: meteo_usta_varid
  integer :: meteo_aerr_varid
  integer :: meteo_obuk_varid
  integer :: meteo_wsta_varid
  integer :: meteo_topc_varid
  integer :: meteo_winz_varid
  integer :: meteo_winm_varid
  integer :: meteo_temp_varid
  integer :: meteo_sphu_varid
  integer :: meteo_airm_varid
  integer :: meteo_clwc_varid
  ! lmbb deepconv
  integer :: meteo_dpeu_varid
  integer :: meteo_dped_varid
  integer :: meteo_dpdu_varid
  integer :: meteo_dpdd_varid
  integer :: meteo_nphourm_varid
  
  integer :: meteo_alti_varid
  integer :: meteo_kzzz_varid

  integer :: bounconc_ncid
  integer :: bounconcincr_ncid
  integer :: bounconc_times_varid
  integer :: bounconcincr_times_varid
  integer :: latconc_conc_varid
  integer :: latconcincr_conc_varid
  integer :: topconc_conc_varid
  integer :: topconcincr_conc_varid
  integer :: atopconc_conc_varid
  integer :: alatconc_conc_varid
  
  integer :: out_ncid
  integer :: out_times_varid
  integer :: out_lon_varid
  integer :: out_lat_varid
  integer,dimension(4) :: out_times_varid_ad
  integer,dimension(4) :: out_lon_varid_ad
  integer,dimension(4) :: out_lat_varid_ad
  
  integer :: out_avcoord_varid
  integer :: out_bvcoord_varid
  integer :: out_cutoff_varid
  integer :: out_hlay_varid
  integer :: out_airm_varid
  integer :: out_relh_varid
  integer :: out_temp_varid
  integer,dimension(4) :: out_species_varid
  integer :: aoutea_ncid,aouteb_ncid,aoutbc_ncid,aoutini_ncid ! Final backup of all adjoint outputs
  integer :: exc_ncid
  integer :: exc_times_varid
  integer :: exc_lon_varid
  integer :: exc_lat_varid

  integer :: par_ncid
  integer :: par_times_varid
  integer :: par_lon_varid
  integer :: par_lat_varid

  integer :: end_ncid
  integer :: end_times_varid
  integer :: end_lon_varid
  integer :: end_lat_varid
  integer :: end_species_varid
  integer :: aend_ncid ! final state of the adjoint

  integer :: endex_ncid
  integer :: endex_times_varid
  integer :: endex_lon_varid
  integer :: endex_lat_varid
  integer :: endex_species_varid

  integer :: depo_ncid
  integer :: depo_times_varid
  integer :: depo_lon_varid
  integer :: depo_lat_varid
  integer :: depo_species_varid


    integer :: chckexch

contains

  subroutine master_allocall
 ! print*,' Along the time'
    allocate(idate (0:nhourrun))
    allocate(ihour (0:nhourrun))
    allocate(iday  (0:nhourrun))
    allocate(imonth(0:nhourrun))
    allocate(iyear (0:nhourrun))
    allocate(idtyp (0:nhourrun))
    allocate(nphourm (0:nhourrun))
  !print*,' For chemistry',ntyperate
    allocate(emisa_species(nemisa))
    allocate(emisb_species(nemisb))
    allocate(aconc_species(nspec))
!    print*,'UUUUUUUUUUUUUUU', nspecboun
!    allocate(aboun_species(nspecboun))
    allocate(ltabrate(ntyperate))
    allocate(species(nspectot))
    allocate(emisavarid(nemisa))
    allocate(emisaincrvarid(nemisa))
    allocate(emisbvarid(nemisb))
    allocate(emisbincrvarid(nemisb))
    allocate(output_species(nspectot))
    allocate(output_par(nparammax))
    allocate(output_depo(2*nspec))
    allocate(kreacp(nspectot))
    allocate(kreacl(nspectot))
    allocate(nreactants(nreac))
    allocate(irctt(nreac,nreactamax))
    allocate(ireacl(nspectot,nreac))
    allocate(ireacp(nspectot,nreac))
    allocate(stoi(nspectot,nreac,ntemps))
    allocate(tabrate(ntabmax,nreac))
    allocate(iphoto(nreac))
    allocate(ityperate(nreac))
    allocate(indepo(0:nspec))
    allocate(inemisa(nspec))
    allocate(inemisb(nspec))
    allocate(altiphot(nlevphotmax))
    allocate(photoj(ntabuzenmax,nlevphotmax,nphotmax))
    allocate(zetaref(ntabuzenmax))
    allocate(tabtemp(ntemps))
    allocate(nelem(nfam))
    allocate(ifam(nfam,nspectot*10))

  !print*,'  Deposition and soil processes'
    allocate(factRb(nspec))
    allocate(factD(nspec))
    allocate(Rm(nspec))
    allocate(dHx(nspec))
    allocate(df0(nspec))
    allocate(deptmin(nvegtype))
    allocate(deptmax(nvegtype))
    allocate(deptopt(nvegtype))
    allocate(depvpd1(nvegtype))
    allocate(depvpd2(nvegtype))
    allocate(depegs(nvegtype))
    allocate(depsgs(nvegtype))
    allocate(depegl(nvegtype))
    allocate(depsgl(nvegtype))
    allocate(deplai1(nvegtype))
    allocate(deplai2(nvegtype))
    allocate(depphe0(nvegtype))
    allocate(depphe1(nvegtype))
    allocate(depphe2(nvegtype))
    allocate(depalph(nvegtype))
    allocate(gmax  (nvegtype))
    allocate(fmin  (nvegtype))
    allocate(zcanopy(nvegtype))
    allocate(inwetd(nspec,2))
    allocate(rgso3(nvegtype))
    allocate(so2rh(nvegtype))
    allocate(rgsso2(nvegtype))

  end subroutine master_allocall

  subroutine master_deallocall
    !print*,' Along the time'
    deallocate(idate)
    deallocate(ihour)
    deallocate(iday)
    deallocate(imonth)
    deallocate(iyear)
    deallocate(idtyp)
    deallocate(nphourm)
    !print*,'params'
    deallocate(int_params)
    deallocate(dbl_params)
    ! print*,'Others'
!    deallocate(dom)
!    deallocate(isboun)
  end subroutine master_deallocall

end module chimere_common
