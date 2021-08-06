MODULE SPECIES_NAME

    implicit none

    INTEGER, ALLOCATABLE, DIMENSION(:), save :: restart_ids ! IDs of active species in the restart file

    INTEGER :: ifn_args, ioerr
    CHARACTER(len=256) :: fnacspec    ! Name of active species file
    CHARACTER(len=256) :: fnallspec   ! Name of all species file
    CHARACTER(len=256) :: fnprescr    ! Name of prescribed species file
    CHARACTER(len=256) :: fnprodl     ! Name of prodloss species file
    CHARACTER(len=256) :: fndep       ! Name of deposited species file

    CHARACTER(len=256) :: fnchem      ! Chemical reactions file
    CHARACTER(len=256) :: fnstoi      ! Stoichiometric coefficients file
    CHARACTER(len=256) :: fnrates     ! Reaction rates file
    CHARACTER(len=256) :: fnfamilies  ! Families list and contents file
    CHARACTER(len=256) :: fnjrates    ! Photolysis parameters file

    INTEGER :: iallqmax                 ! Number of species (active + prescribed + prodloss)
    INTEGER :: iqmax                    ! Number of output species (active)
    INTEGER :: iprescrmax               ! Number of prescribed species
    INTEGER :: iprodmax                 ! Number of prod/loss species
    INTEGER :: idepmax                  ! Number of deposited species
    INTEGER :: ijratesmax               ! Number of photolysis reactions with prescribed constant J
    INTEGER :: nreac                    ! Number of reactions
    INTEGER,PARAMETER :: ntyperate   = 50 ! Max number of reaction types
    INTEGER,PARAMETER :: ntabmax     = 10 ! Max number of rate constants


    INTEGER, ALLOCATABLE, DIMENSION(:)            :: kreacp      ! Number of reactions producing
    INTEGER, ALLOCATABLE, DIMENSION(:)            :: kreacl      ! Number of reactions destructing
    INTEGER, ALLOCATABLE, DIMENSION(:,:)          :: idreacp     ! Addresses of products
    INTEGER, ALLOCATABLE, DIMENSION(:,:)          :: idreacl     ! Addresses of reactants
    INTEGER, ALLOCATABLE, DIMENSION(:,:)          :: idspecl     ! Address of reactants
    INTEGER, ALLOCATABLE, DIMENSION(:)            :: idtyperate  ! Reaction rate types
    INTEGER, ALLOCATABLE, DIMENSION(:)            :: nreactants  ! Number of reactants in the reaction
    REAL(kind=8), ALLOCATABLE, DIMENSION(:,:)     :: tabrate     ! Reaction rate constants
    REAL(kind=8), ALLOCATABLE, DIMENSION(:,:)     :: stoi        ! Stoichiometric coeffs.
    INTEGER,DIMENSION(ntyperate)                  :: ltabrate    ! Their number of constants
    REAL, ALLOCATABLE, DIMENSION(:)               :: adv_mass    ! Molar mass of the active species


    LOGICAL :: do_chemistry = .TRUE.


  ! A structure to hold extended information about species
    TYPE :: all_species_type
        INTEGER           :: id
        CHARACTER(len=10) :: name
        CHARACTER(len=2 ) :: type
        LOGICAL           :: dep
        INTEGER           :: iddep
        LOGICAL           :: prod
        INTEGER           :: idprod
    ENDTYPE all_species_type

    ! A structure to hold information about species
    TYPE :: species_type
        INTEGER           :: id
        CHARACTER(len=10) :: name
    ENDTYPE species_type

    ! A structure to hold information about jrates
    TYPE :: idjrates_type
        INTEGER           :: idreac
        CHARACTER(len=10) :: idj
    END TYPE idjrates_type


    TYPE(all_species_type), ALLOCATABLE, DIMENSION(:) :: all_species
    TYPE(species_type), ALLOCATABLE, DIMENSION(:) :: species, prescr_species, prodloss_species, dep_species
    TYPE(idjrates_type), ALLOCATABLE, DIMENSION(:) :: jrates


END MODULE SPECIES_NAME

MODULE CONSTANTS

  REAL, PARAMETER     :: boltz = 1.38044e-16      ! erg/K
  REAL, PARAMETER     :: avogadro = 6.02214199e+23 ! mol-1
  REAL, PARAMETER     :: dry_mass= 28.966    !masse mol air sec
  INTEGER, PARAMETER  :: p_ref=101325.      !pression a 1 atm

END MODULE CONSTANTS

MODULE variables

  IMPLICIT NONE
  SAVE
  INTEGER :: lat_len, lon_len, time_len, lev_len, vect

END MODULE variables