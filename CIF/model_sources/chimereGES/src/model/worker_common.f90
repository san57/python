module worker_common

    use message_defs
    use chimere_consts

    implicit none

    integer :: nitgssu  ! Number of G-S iters. for TWOSTEP during spin-up
    integer :: nitgs  ! Number of G-S iters. for TWOSTEP

    integer :: usechemistry, usedepos, usetransmix, usewetdepos, useemissions, useabsclipconc, dryairout

    ! A structure to hold information about active species
    ! lmbb add a new var in module species_type
    type :: species_type
        integer :: varid
        character(len = 16) :: name
        integer :: transp
        integer :: transpv
        integer :: bounddry ! Boundary conditions are in dry mole if 1
        integer :: varido ! the matching conco
        real(kind = 8) :: fspec          ! conversion factor
    end type species_type
    type(species_type), dimension(:), allocatable :: species

    ! A structure to hold informations about the fields to write to file "out"
    type :: output_species_type
        integer :: varid          ! netcdf var id
        integer :: iaddr          ! address in species(:) table
        character(len = splen) :: name           ! name from species(:) table
        character(len = 16) :: units
    end type output_species_type
    type(output_species_type), dimension(:), allocatable :: output_species

    ! MPI
    integer :: fwd
    integer :: tl
    integer :: ad
    integer :: wrk_comm                        ! Communicator local to workers
    integer :: rank                            ! process rank
    integer :: nmdoms, nzdoms ! Number of MPI subdomains in meridian and zonal directions

    ! Timing
    integer :: nphour_ref ! Number of physical steps per hour
    integer :: ihoursu
    integer :: ichemstep ! Number of chemical refined iterations
    integer :: nhourrun                   ! Simulation duration (Hours)
    real(kind = 8) :: soltim                     ! Time interval with year's solstice
    real(kind = 8) :: djul                       ! Julian day
    integer, dimension(:), allocatable :: ihour  ! UT hours
    integer, dimension(:), allocatable :: imonth ! Month
    integer, dimension(:), allocatable :: idtyp  ! Day types (week, sat, sun)
    integer :: nphour
    real(kind = 8) :: dtr
    real(kind = 8) :: dtr2
    integer :: ihourrun                   ! Hour counter
    real(kind = 8) :: thour                      ! Time (hour) since last exact hour

    ! IP Observations
    real(kind = 8), allocatable, dimension(:, :) :: tabobs
    integer :: nobs
    ! fin IP  Observations
    ! IP adj
    integer :: ifnconcp
    integer :: ifnconch
    ! fin IP adj

    !  Geometry
    integer :: nverti
    integer :: nzonalmax  ! Max Number of zonal cells for all subdomains
    integer :: nmeridmax ! Max Number of meridional cells for all subdomains
    integer :: nzonal                          ! Number of zonal cells
    integer :: nmerid                          ! Number of meridional cells
    integer :: dom_i                           ! subdomain zonal index
    integer :: dom_j                           ! subdomain meridian index
    ! IP
    integer :: imstart, imend, izstart, izend ! limites en im et iz par rapport a la grille complete
    ! fin IP
    real(kind = 8), dimension(:, :), allocatable :: xlong    ! Longitudes
    real(kind = 8), dimension(:, :), allocatable :: xlati    ! Latitudes
    real(kind = 8), dimension(:, :), allocatable :: clati    ! Cosine of latitudes
    real(kind = 8), dimension(:, :), allocatable :: slati    ! Sine of latitudes
    real(kind = 8), dimension(:, :), allocatable :: xsize    ! Zonal cell length
    real(kind = 8), dimension(:, :), allocatable :: ysize    ! Meridional cell length
    real(kind = 8), dimension(:, :), allocatable :: xbasx    ! x coord of local i vector
    real(kind = 8), dimension(:, :), allocatable :: xbasy    ! y coord of local i vector
    real(kind = 8), dimension(:, :), allocatable :: ybasx    ! x coord of local j vector
    real(kind = 8), dimension(:, :), allocatable :: ybasy    ! y coord of local j vector
    real(kind = 8 ) :: psurf ! surface pressure in Pa

    !  Chemical workspace
    real(kind = 8) :: clipconc ! Clipping value for the TWOSTEP algorithm
!    integer :: ntabuzenmax ! Max number of tabulated zenith angles
!    integer :: nphotmax  ! Max number of photolysis reactions
    integer :: nspectot ! Total num. species
    integer :: ntemps ! Number of tabul. temperatures for stoichio.
!    integer :: nlevphotmax ! Max number of tabulated photolysis levels
    integer :: ntabmax ! Max number of rate constants
    integer :: nspec ! Number of active species
    integer :: nreac ! Number of reactions
    integer :: nreactamax ! Max number of reactants/reaction
    integer :: nfam          ! Number of "family" species (for output)
    integer :: ndep ! Number of deposited species
    integer :: noutspec    ! Number of output species! chemical workspace
    integer :: nsaveconcs ! Save concentrations for restart every ... hours

    real(kind = 8), dimension(:, :, :, :), allocatable :: conc     ! Current concentration array
    real(kind = 8), dimension(:, :, :, :), allocatable :: conc_tl
    real(kind = 8), dimension(:, :, :, :), allocatable :: aconc     ! Current concentration array
    real(kind = 8), dimension(:, :, :, :), allocatable :: conco    ! Previous one (for TWOSTEP)
    real(kind = 8), dimension(:, :, :, :), allocatable :: conco_tl
    real(kind = 8), dimension(:, :, :, :), allocatable :: aconco
    real(kind = 8), dimension(:, :, :, :), allocatable :: rate     ! Current rate coefficients (K&J)
    real(kind = 8), dimension(:, :, :, :), allocatable :: phrate   ! Current photolysis rates
    real(kind = 8), dimension(:, :, :), allocatable :: wgstl    ! Low temperature wght for stoich.
    real(kind = 8), dimension(:, :, :), allocatable :: wgsth    ! High temperature wght for stoich.
    integer, dimension(:, :, :), allocatable :: istoit    ! Tabulated T address for stoich.
    integer, dimension(:), allocatable :: nreactants  ! Number of reactants
    integer, dimension(:), allocatable :: ityperate   ! Reaction rate types
    integer, dimension(:), allocatable :: iphoto       ! Photolysis reaction address
    integer, dimension(:, :), allocatable :: irctt       ! Address of reactants
    integer, dimension(:), allocatable :: kreacp      ! Number of reactions producing
    integer, dimension(:), allocatable :: kreacl      ! Number of reactions destructing
    integer, dimension(:, :), allocatable :: ireacp      ! Addresses of products
    integer, dimension(:, :), allocatable :: ireacl      ! Addresses of reactants
    integer, dimension(:, :), allocatable :: ifam        ! Addresses of elements in family
    integer, dimension(:), allocatable :: ltabrate    ! Their number of constants
    integer, dimension(:), allocatable :: nelem       ! Num. of elements in families
    real(kind = 8), dimension(:, :, :), allocatable :: stoi     ! Stoichiometric coeffs.
    real(kind = 8), dimension(:), allocatable :: zetaref  ! Tabulated cos of zenith angles
    real(kind = 8), dimension(:, :), allocatable :: tabrate  ! Reaction rate constants
    real(kind = 8), dimension(:), allocatable :: tabtemp  ! Tabulated temperatures for stoic.
    real(kind = 8), dimension(:), allocatable :: altiphot ! Altitudes of photolysis levels
    real(kind = 8), dimension(:, :, :), allocatable :: photoj   ! Clear sky photolysis rates

    integer :: ntabuzen    ! Number of tabulated zenith angles
    integer :: nphot       ! Number of photolysis reactions
    integer :: nlevphot    ! Number of tabulated photolysis levels
    integer :: inNH3       ! NH3  species number
    integer :: inSO2       ! SO2  species number
    integer :: inHNO3      ! HNO3 species number

    ! Meteorology
    integer :: ideepconv
    real(kind = 8), dimension(:, :, :), allocatable :: airmloc    ! Current density
    real(kind = 8), dimension(:, :, :), allocatable :: winzloc    ! Current zonal wind
    real(kind = 8), dimension(:, :, :), allocatable :: winmloc    ! Current meridional wind
    real(kind = 8), dimension(:, :, :), allocatable :: thlayloc   ! Current layer thicknesses
    real(kind = 8), dimension(:, :, :), allocatable :: winvloc    ! Current vertical wind
    real(kind = 8), dimension(:, :, :), allocatable :: sphuloc    ! Current specific humidity
    real(kind = 8), dimension(:, :, :), allocatable :: temploc    ! Current temperature
    real(kind = 8), dimension(:, :, :), allocatable :: winxloc    ! Current mixing velocities
    real(kind = 8), dimension(:, :, :), allocatable :: hlayloc    ! Current model layer heights
    real(kind = 8), dimension(:, :, :), allocatable :: kzzzloc    ! Current top-layer Kz
    real(kind = 8), dimension(:, :, :), allocatable :: dtenloc
    real(kind = 8), dimension(:, :, :), allocatable :: clwcloc
    real(kind = 8), dimension(:, :, :), allocatable :: presloc
    real(kind = 8), dimension(:, :, :), allocatable :: scinloc
    real(kind = 8), dimension(:, :, :), allocatable :: phloc
    real(kind = 8), dimension(:, :), allocatable :: hghtloc   ! Current mixing height
    real(kind = 8), dimension(:, :), allocatable :: atteloc   ! Current radiation attenuation
    real(kind = 8), dimension(:, :), allocatable :: zeniloc   ! Current zenith angle
    real(kind = 8), dimension(:, :), allocatable :: tem2loc   ! Current 2m-temperature
    real(kind = 8), dimension(:, :), allocatable :: ustaloc   ! Current ustar
    real(kind = 8), dimension(:, :), allocatable :: aerrloc   ! Current aerodynamic resistance
    real(kind = 8), dimension(:, :), allocatable :: obukloc   ! Current Obukov length
    real(kind = 8), dimension(:, :), allocatable :: wstaloc   ! Current Wstar
    real(kind = 8), dimension(:, :), allocatable :: topcloc   ! Current total precipitation
    real(kind = 8), dimension(:, :), allocatable :: srehloc   ! Current surface relative humidity
    real(kind = 8), dimension(:, :), allocatable :: acraloc
    real(kind = 8), dimension(:, :, :), allocatable :: uwest     ! Zonal western wind
    real(kind = 8), dimension(:, :, :), allocatable :: ueast     ! Zonal eastern wind
    real(kind = 8), dimension(:, :, :), allocatable :: unorth    ! Meridional northern wind
    real(kind = 8), dimension(:, :, :), allocatable :: usouth    ! Meridional southern wind
    real(kind = 8), dimension(:, :, :), allocatable :: hwest     ! Western thickness
    real(kind = 8), dimension(:, :, :), allocatable :: heast     ! Eastern thickness
    real(kind = 8), dimension(:, :, :), allocatable :: hnorth    ! Northern thickness
    real(kind = 8), dimension(:, :, :), allocatable :: hsouth    ! Southern thickness
    real(kind = 8), dimension(:, :, :), allocatable :: swest     ! Western cell side length
    real(kind = 8), dimension(:, :, :), allocatable :: seast     ! Eastern cell side length
    real(kind = 8), dimension(:, :, :), allocatable :: snorth    ! Northern cell side length
    real(kind = 8), dimension(:, :, :), allocatable :: ssouth    ! Southern cell side length
    ! lmbb deepconv
    real(kind = 8), dimension(:, :, :), allocatable :: dpeuloc  ! Entrainment in updraft (kg/m2/s)
    real(kind = 8), dimension(:, :, :), allocatable :: dpedloc  ! Detrainment in updraft  (kg/m2/s)
    real(kind = 8), dimension(:, :, :), allocatable :: dpduloc  ! Entrainment in downdraft (kg/m2/s)
    real(kind = 8), dimension(:, :, :), allocatable :: dpddloc  ! Detrainment in downdraft (kg/m2/s)
    real(kind = 8), dimension(:, :, :), allocatable :: flxuloc  ! Current updraught mass flux
    real(kind = 8), dimension(:, :, :), allocatable :: flxdloc  ! Current downdraught mass flux
    ! flux of species in the...
    real(kind = 8), dimension(:, :, :, :), allocatable :: flxuconc  ! ..updraught (molecules/cm2/s)
    real(kind = 8), dimension(:, :, :, :), allocatable :: flxdconc  ! ..downdraught (molecules/cm2/s)
    real(kind = 8), dimension(:, :, :, :), allocatable :: flxeconc  ! ..environment (molecules/cm2/s)
    real(kind = 8), dimension(:, :, :, :), allocatable :: flxuconc_tl  ! ..updraught (molecules/cm2/s)
    real(kind = 8), dimension(:, :, :, :), allocatable :: flxdconc_tl  ! ..downdraught (molecules/cm2/s)
    real(kind = 8), dimension(:, :, :, :), allocatable :: flxeconc_tl  ! ..environment (molecules/cm2/s)
    real(kind = 8), dimension(:, :, :, :), allocatable :: aflxuconc  ! ..updraught (molecules/cm2/s) CS transmix
    real(kind = 8), dimension(:, :, :, :), allocatable :: aflxdconc  ! ..downdraught (molecules/cm2/s) CS transmix
    real(kind = 8), dimension(:, :, :, :), allocatable :: aflxeconc  ! ..environment (molecules/cm2/s) CS transmix
    integer, dimension(:, :), allocatable :: ideep

    real(kind = 8), dimension(:, :, :, :), allocatable :: winz      ! Hourly Zonal wind
    real(kind = 8), dimension(:, :, :, :), allocatable :: winm      ! Hourly Meridional wind
    real(kind = 8), dimension(:, :, :, :), allocatable :: temp      ! Hourly Temperature
    real(kind = 8), dimension(:, :, :, :), allocatable :: sphu      ! Hourly Specific humidity
    real(kind = 8), dimension(:, :, :, :), allocatable :: airm      ! Hourly density
    real(kind = 8), dimension(:, :, :, :), allocatable :: hlay      ! Hourly model layer heights
    real(kind = 8), dimension(:, :, :, :), allocatable :: kzzz      ! Hourly Kz
    real(kind = 8), dimension(:, :, :, :), allocatable :: clwc      ! Hourly Liquid (+ice) water content
    ! lmbb deepconv
    real(kind = 8), dimension(:, :, :, :), allocatable :: dpeu    ! Entrainment in updraft (kg/m2/s)
    real(kind = 8), dimension(:, :, :, :), allocatable :: dped    ! Detrainment in updraft  (kg/m2/s)
    real(kind = 8), dimension(:, :, :, :), allocatable :: dpdu    ! Entrainment in downdraft (kg/m2/s)
    real(kind = 8), dimension(:, :, :, :), allocatable :: dpdd    ! Detrainment in downdraft (kg/m2/s)

    real(kind = 8), dimension(:, :, :), allocatable :: hght      ! Hourly mixing layer heights
    real(kind = 8), dimension(:, :, :), allocatable :: atte      ! Hourly radiation attenuation
    real(kind = 8), dimension(:, :, :), allocatable :: tem2      ! Hourly 2m-temperature
    real(kind = 8), dimension(:, :, :), allocatable :: usta      ! Hourly Ustar
    real(kind = 8), dimension(:, :, :), allocatable :: aerr      ! Hourly aerodynamic resistance
    real(kind = 8), dimension(:, :, :), allocatable :: obuk      ! Hourly Obukov length
    real(kind = 8), dimension(:, :, :), allocatable :: wsta      ! Hourly Wstar
    real(kind = 8), dimension(:, :, :), allocatable :: topc      ! Hourly total precipitation
    real(kind = 8), dimension(:, :, :), allocatable :: sreh      ! Hourly surface relative humidity
    integer, dimension(:, :, :), allocatable :: incloud


    !  Transport

    real(kind = 8), dimension(:, :, :, :), allocatable :: vfluxi   ! Input vertical fluxes
    real(kind = 8), dimension(:, :, :), allocatable :: vfluxo   ! Output vertical fluxes
    real(kind = 8), dimension(:, :, :), allocatable :: fluxw    ! Western fluxes
    real(kind = 8), dimension(:, :, :), allocatable :: fluxe    ! Eastern fluxes
    real(kind = 8), dimension(:, :, :), allocatable :: fluxs    ! Southern fluxes
    real(kind = 8), dimension(:, :, :), allocatable :: fluxn    ! Northern fluxes

    !  Emissions workspace
    integer :: nlevemis ! Number of emission levels
    integer :: hpulse ! Pulse hour Elise Potier
    integer :: nemisa ! Max number of anthropic emitted species
    integer :: nemisb ! Max number of biogenic emitted species
    integer :: optemisb
    integer, dimension(1000) :: ispecemip
    real(kind = 8), dimension(:, :, :), allocatable :: emisbloc ! Current biogenic emissions
    real(kind = 8), dimension(:, :, :), allocatable :: emisbloc_tl ! Current biogenic emissions increments
    real(kind = 8), dimension(:, :, :), allocatable :: aemisbloc
    real(kind = 8), dimension(:, :, :, :), allocatable :: emisaloc ! Current anthropic emissions
    real(kind = 8), dimension(:, :, :, :), allocatable :: emisaloc_tl ! Current anthropic emissions incrments
    real(kind = 8), allocatable, dimension(:, :, :, :, :) :: aemisaloc
    real(kind = 8), dimension(:, :, :, :), allocatable :: emisb    ! Input biogenic emis. at exact hours
    real(kind = 8), allocatable, dimension(:, :, :, :) :: aemisb
    real(kind = 8), dimension(:, :, :, :), allocatable :: emisb_tl ! Input biogenic emis. increments at exact hours
    integer, dimension(:), allocatable :: inemisa  ! Addresses of anthropic species
    integer, dimension(:), allocatable :: inemisb  ! Addresses of biogenic species

    !  Deposition and soil processes
    integer :: nvegtype
    integer :: nlduse
    integer :: ndepo
    integer, dimension(2) :: nwetd
    integer, dimension(:, :), allocatable :: inwetd
    real(kind = 8), dimension(:), allocatable :: gmax
    real(kind = 8), dimension(:), allocatable :: fmin
    real(kind = 8), dimension(:, :, :, :), allocatable :: kparloc    ! Partition coefficient
    real(kind = 8), dimension(:, :, :, :), allocatable :: kvaloc
    real(kind = 8), dimension(:, :, :), allocatable :: khloc
    real(kind = 8), dimension(:, :, :, :), allocatable :: ktrveloc
    real(kind = 8), dimension(:, :, :, :), allocatable :: fveg
    real(kind = 8), dimension(:, :, :, :), allocatable :: drydepi
    real(kind = 8), dimension(:, :, :, :), allocatable :: drydepi_tl
    real(kind = 8), dimension(:, :, :, :), allocatable :: wetdepi
    real(kind = 8), dimension(:, :, :, :), allocatable :: wetdepi_tl
    real(kind = 8), dimension(:, :, :, :), allocatable :: wetdr1
    real(kind = 8), dimension(:, :, :, :), allocatable :: wetdr2
    real(kind = 8), dimension(:, :, :), allocatable :: drydep
    real(kind = 8), dimension(:, :, :), allocatable :: drydep_tl
    real(kind = 8), dimension(:, :, :), allocatable :: wetdep
    real(kind = 8), dimension(:, :, :), allocatable :: wetdep_tl
    real(kind = 8), dimension(:, :, :), allocatable :: dland
    real(kind = 8), dimension(:, :, :), allocatable :: depoloc
    integer, dimension(:), allocatable :: indepo
    real(kind = 8), dimension(:), allocatable :: deptmin
    real(kind = 8), dimension(:), allocatable :: deptopt
    real(kind = 8), dimension(:), allocatable :: deptmax
    real(kind = 8), dimension(:), allocatable :: depalph
    real(kind = 8), dimension(:), allocatable :: depvpd1
    real(kind = 8), dimension(:), allocatable :: depvpd2
    real(kind = 8), dimension(:), allocatable :: depsgs
    real(kind = 8), dimension(:), allocatable :: depegs
    real(kind = 8), dimension(:), allocatable :: depsgl
    real(kind = 8), dimension(:), allocatable :: depegl
    real(kind = 8), dimension(:), allocatable :: deplai1
    real(kind = 8), dimension(:), allocatable :: deplai2
    real(kind = 8), dimension(:), allocatable :: depphe0
    real(kind = 8), dimension(:), allocatable :: depphe1
    real(kind = 8), dimension(:), allocatable :: depphe2
    real(kind = 8), dimension(:), allocatable :: zcanopy
    real(kind = 8), dimension(:), allocatable :: RGSO3
    real(kind = 8), dimension(:), allocatable :: RGSSO2
    real(kind = 8), dimension(:), allocatable :: so2rh
    real(kind = 8), dimension(:), allocatable :: dHx
    real(kind = 8), dimension(:), allocatable :: df0
    real(kind = 8), dimension(:), allocatable :: factRb
    real(kind = 8), dimension(:), allocatable :: factD
    real(kind = 8), dimension(:), allocatable :: Rm

    integer :: chckexch

contains
    subroutine worker_allocall
        ! along time
        allocate(imonth(0:nhourrun))
        allocate(ihour (0:nhourrun))
        allocate(idtyp (0:nhourrun))
        ! On the tile
        allocate(xsize(-1:nzonal + 2, -1:nmerid + 2))
        allocate(ysize(-1:nzonal + 2, -1:nmerid + 2))
        allocate(xbasx(0:nzonal + 1, 0:nmerid + 1))
        allocate(xbasy(0:nzonal + 1, 0:nmerid + 1))
        allocate(ybasx(0:nzonal + 1, 0:nmerid + 1))
        allocate(ybasy(0:nzonal + 1, 0:nmerid + 1))
        allocate(slati(nzonal, nmerid))
        allocate(clati(nzonal, nmerid))
        allocate(xlong(nzonal, nmerid))
        allocate(xlati(nzonal, nmerid))
        allocate(ideep(nzonal, nmerid))
        allocate(fveg(nzonal, nmerid, nvegtype, nlduse))
        allocate(emisaloc(nemisa, nzonal, nmerid, nlevemis))
        allocate(emisb(nemisb, nzonal, nmerid, 2))
        allocate(emisbloc(nemisb, nzonal, nmerid))
        allocate(dland(nzonal, nmerid, nlduse))
        allocate(temp(nzonal, nmerid, nverti, 2))
        allocate(sphu(nzonal, nmerid, nverti, 2))
        allocate(airm(nzonal, nmerid, nverti, 2))
        allocate(clwc(nzonal, nmerid, nverti, 2))
        allocate(dpeu(nzonal, nmerid, nverti, 2))
        allocate(dped(nzonal, nmerid, nverti, 2))
        allocate(dpdu(nzonal, nmerid, nverti, 2))
        allocate(dpdd(nzonal, nmerid, nverti, 2))
        allocate(hlay(nzonal, nmerid, nverti, 2))
        allocate(kzzz(nzonal, nmerid, nverti, 2))
        allocate(tem2(nzonal, nmerid, 2))
        allocate(sreh(nzonal, nmerid, 2))
        allocate(atte(nzonal, nmerid, 2))
        allocate(hght(nzonal, nmerid, 2))
        allocate(usta(nzonal, nmerid, 2))
        allocate(aerr(nzonal, nmerid, 2))
        allocate(obuk(nzonal, nmerid, 2))
        allocate(wsta(nzonal, nmerid, 2))
        allocate(topc(nzonal, nmerid, 2))
        allocate(winz(nzonal, nmerid, nverti, 2))
        allocate(winm(nzonal, nmerid, nverti, 2))
        allocate(hghtloc(nzonal, nmerid))
        allocate(atteloc(nzonal, nmerid))
        allocate(tem2loc(nzonal, nmerid))
        allocate(airmloc(-2:nzonalmax + 3, -2:nmeridmax + 3, 1:nverti + 1))
        allocate(winzloc(0:nzonalmax + 1, 0:nmeridmax + 1, nverti))
        allocate(winmloc(0:nzonalmax + 1, 0:nmeridmax + 1, nverti))
        allocate(presloc(nzonal, nmerid, nverti))
        allocate(incloud(nzonal, nmerid, nverti))
        allocate(dtenloc(nzonal, nmerid, nverti))
        allocate(temploc(nzonal, nmerid, nverti))
        allocate(sphuloc(nzonal, nmerid, nverti))
        allocate(clwcloc(nzonal, nmerid, nverti))
        allocate(hlayloc(nzonal, nmerid, nverti))
        allocate(kzzzloc(nzonal, nmerid, nverti))
        allocate(thlayloc(0:nzonalmax + 1, 0:nmeridmax + 1, nverti))
        allocate(winvloc(nzonal, nmerid, nverti))
        allocate(winxloc(nzonal, nmerid, nverti))
        allocate(uwest(nzonal, nmerid, nverti))
        allocate(hwest(nzonal, nmerid, nverti))
        allocate(swest(nzonal, nmerid, nverti))
        allocate(ueast  (nzonal, nmerid, nverti))
        allocate(unorth (nzonal, nmerid, nverti))
        allocate(usouth (nzonal, nmerid, nverti))
        allocate(heast  (nzonal, nmerid, nverti))
        allocate(hnorth (nzonal, nmerid, nverti))
        allocate(hsouth (nzonal, nmerid, nverti))
        allocate(seast  (nzonal, nmerid, nverti))
        allocate(snorth (nzonal, nmerid, nverti))
        allocate(ssouth (nzonal, nmerid, nverti))
        allocate(fluxw(nzonal, nmerid, nverti))
        allocate(fluxe(nzonal, nmerid, nverti))
        allocate(fluxs(nzonal, nmerid, nverti))
        allocate(fluxn(nzonal, nmerid, nverti))
        allocate(vfluxo(nzonal, nmerid, nverti))
        allocate(vfluxi(nzonal, nmerid, nverti, nverti + 1))
        allocate(dpeuloc(nzonal, nmerid, nverti))
        allocate(dpedloc(nzonal, nmerid, nverti))
        allocate(dpddloc(nzonal, nmerid, nverti))
        allocate(dpduloc(nzonal, nmerid, nverti))
        allocate(flxuloc(nzonal, nmerid, nverti))
        allocate(flxdloc(nzonal, nmerid, nverti))

        ! For chemistry
        allocate(conc(nspectot, -2:nzonal + 3, -2:nmerid + 3, 1:nverti + 1))
        allocate(flxuconc(nspectot, nzonal, nmerid, nverti))
        allocate(flxdconc(nspectot, nzonal, nmerid, nverti))
        allocate(flxeconc(nspectot, nzonal, nmerid, nverti))
        allocate(conco(nspectot, nzonal, nmerid, nverti))
        allocate(rate(nreac, nzonal, nmerid, nverti))
        allocate(phrate(nphot, nzonal, nmerid, nverti))
        allocate(wgstl(nzonal, nmerid, nverti))
        allocate(wgsth(nzonal, nmerid, nverti))
        allocate(drydep(nspec, nzonal, nmerid))
        allocate(wetdep(nspec, nzonal, nmerid))
        allocate(output_species(nspectot))
        allocate(species(nspectot))
        allocate(inemisa(nspec))
        allocate(inemisb(nspec))
        allocate(kreacl(nspectot))
        allocate(kreacp(nspectot))
        allocate(nreactants(nreac))
        allocate(irctt(nreac, nreactamax))
        allocate(ireacl(nspectot, nreac))
        allocate(ireacp(nspectot, nreac))
        allocate(stoi(nspectot, nreac, ntemps))
        allocate(tabrate(ntabmax, nreac))
        allocate(iphoto(nreac))
        allocate(ityperate(nreac))
        allocate(indepo(0:nspec))
        allocate(altiphot(nlevphot))
        allocate(photoj(ntabuzen, nlevphot, nphot))
        allocate(nelem(nfam))
        allocate(ifam(nfam, nspectot * 10))
        allocate(zetaref(ntabuzen))
        allocate(tabtemp(ntemps))
        allocate(zeniloc(nzonal, nmerid))
        allocate(ustaloc(nzonal, nmerid))
        allocate(aerrloc(nzonal, nmerid))
        allocate(obukloc(nzonal, nmerid))
        allocate(wstaloc(nzonal, nmerid))
        allocate(topcloc(nzonal, nmerid))
        allocate(srehloc(nzonal, nmerid))
        allocate(acraloc(nzonal, nmerid))
        allocate(phloc(nzonal, nmerid, nverti))
        allocate(istoit(nzonal, nmerid, nverti))

        !  Deposition and soil processes
        allocate(factRb(nspec))
        allocate(factD(nspec))
        allocate(Rm(nspec))
        allocate(dHx(nspec))
        allocate(df0(nspec))
        allocate(deptmin(nvegtype))
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
        allocate(inwetd(nspec, 2))
        allocate(rgso3(nvegtype))
        allocate(so2rh(nvegtype))
        allocate(rgsso2(nvegtype))
        allocate(depoloc(nspec, nzonal, nmerid))
        allocate(drydepi(nspec, nzonal, nmerid, nverti))
        allocate(wetdepi(nspec, nzonal, nmerid, nverti))
        allocate(wetdr1(nspec,nzonal,nmerid,nverti))
        allocate(wetdr2(nspec,nzonal,nmerid,nverti))

    end subroutine worker_allocall

    subroutine worker_allocall_tl
        allocate(emisaloc_tl(nemisa, nzonal, nmerid, nlevemis))
        allocate(emisb_tl(nemisb, nzonal, nmerid, 2))
        allocate(emisbloc_tl(nemisb, nzonal, nmerid))
        allocate(conc_tl(nspectot, -2:nzonal + 3, -2:nmerid + 3, 1:nverti + 1))
        allocate(flxuconc_tl(nspectot, nzonal, nmerid, nverti))
        allocate(flxdconc_tl(nspectot, nzonal, nmerid, nverti))
        allocate(flxeconc_tl(nspectot, nzonal, nmerid, nverti))
        allocate(conco_tl(nspectot, nzonal, nmerid, nverti))
        allocate(drydep_tl(nspec, nzonal, nmerid))
        allocate(wetdep_tl(nspec, nzonal, nmerid))
        allocate(drydepi_tl(nspec, nzonal, nmerid, nverti))
        allocate(wetdepi_tl(nspec, nzonal, nmerid, nverti))
    end subroutine worker_allocall_tl

    subroutine worker_allocall_ad
        allocate(aconc(nspectot, -2:nzonal + 3, -2:nmerid + 3, 1:nverti + 1)) !Current concentration array
        allocate(aconco(nspectot, nzonal, nmerid, nverti)) !Previous one (for TWOSTEP)
        ! flux of species in the...
        allocate(aflxuconc(nspectot, nzonal, nmerid, nverti))  ! ..updraught (molecules/cm2/s)
        allocate(aflxdconc(nspectot, nzonal, nmerid, nverti))  ! ..downdraught (molecules/cm2/s)
        allocate(aflxeconc(nspectot, nzonal, nmerid, nverti))  ! ..environment (molecules/cm2/s) CS transmix
        allocate(aemisbloc(nemisb, nzonal, nmerid)) !  ADJ Bemissions
    end subroutine  worker_allocall_ad

    subroutine worker_deallocall
        deallocate(int_params)
        deallocate(dbl_params)
    end subroutine worker_deallocall


end module worker_common
