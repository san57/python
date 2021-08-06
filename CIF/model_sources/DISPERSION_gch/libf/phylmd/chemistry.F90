subroutine comp_chemistry(mmr, refprod, refprescr, refjrates, &
        temp, pmid, delt, d_chem, sump, suml, nbprodmax, nblossmax)

  !  Calculation of chemical production and loss terms given the species
  !  (NS) for the current cell of coordinates izo,ime,ivert
  !  INPUT : NS         Species number
  !          izo,ime,ivert coordinates of the current point
  !          KREACL     Number of reactions consuming species NS
  !          IREACL     The Adresses of reactions consuming species NS
  !          RATES      The K or J values for all reactions
  !          NREACTANTS Number of reactants
  !          CONC       Array of instantaneous concensrations
  !          KREACP     Number of reactions producing species NS
  !          IREACP     The Adresses of reactions producing species NS
  !          STOI       Array of Stoichiometric coefficients
  !          WGTSL      Lower weight for interpol. between tabulated tempe
  !          WGTSL      Upper weight for interpol. between tabulated tempe
  !  OUTPUT: CHPR       Chemical production of species NS in box NB
  !          CHLO       Chemical loss       of species NS in box NB

  ! For the moment, all prescribed species are assumed to be prescribed in vmr
  !                 all chemical computations are achieved using molecules/cm3 and then converted to mmr

  USE SPECIES_NAME
  USE CONSTANTS
  USE variables
  USE parallel
  USE dimphy
  USE Write_Field_phy
  IMPLICIT NONE
  INCLUDE "dimensions.h"
  INCLUDE "paramet.h"

  !-----------------------------------------------------------------------
  !      General variables
  !-----------------------------------------------------------------------
  REAL, INTENT(in)    ::  mmr(klon,klev,iqmax)      ! mmr array
  REAL, INTENT(in)    ::  delt                     ! timestep
  REAL, INTENT(in)    ::  refprod(klon,klev,iprodmax)
  REAL, INTENT(in)    ::  refprescr(klon,klev,iprescrmax)
  REAL, INTENT(in)    ::  refjrates(klon,klev,ijratesmax)
  REAL, INTENT(in)    ::  temp(klon,klev)
  REAL, INTENT(in)    ::  pmid(klon,klev)
  INTEGER,INTENT(in)   :: nbprodmax,nblossmax
  REAL, INTENT(out)   ::  d_chem(klon,klev,iqmax)   ! chemical increment
  REAL,INTENT(out)    ::  sump(klev-1,iqmax,nbprodmax,klon)
  REAL,INTENT(out)    ::  suml(klev-1,iqmax,nblossmax,klon)
  !-----------------------------------------------------------------------
  !      Local variables
  !-----------------------------------------------------------------------
  INTEGER :: ns, np
  REAL, ALLOCATABLE, DIMENSION(:,:)   :: invtemp    ! inverse temperature
  REAL, ALLOCATABLE, DIMENSION(:,:)   :: convert  !  vmr to molecule/cm3 conversion factor
  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: prod, loss ! prod/loss terms in VMR
  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: vmr      ! vmr array
  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: molec    ! molecules/cm3 array
  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: rates     ! reaction rates

  INTEGER      :: kr,ire,ir
  INTEGER      :: nr
  INTEGER      :: nk, nl
  REAL, ALLOCATABLE, DIMENSION(:,:) :: trat


  ALLOCATE(invtemp(klon,klev))
  ALLOCATE(prod(klon,klev,iqmax))
  ALLOCATE(loss(klon,klev,iqmax))
  ALLOCATE(vmr(klon,klev,iqmax))
  ALLOCATE(molec(klon,klev,iqmax))
  ALLOCATE(convert(klon,klev))
  ALLOCATE(rates(klon,klev,nreac))
  ALLOCATE(trat(klon,klev))

  loss = 0.
  prod = 0.
  trat = 0.

  ! Compute rates for each reaction
  CALL comp_rates(temp, pmid, refjrates, rates)

  !--- Conversion factor from vmr ---> molecules.cm-3
  convert(:, :) = 10. * pmid(:, :) / ( boltz * temp(:, :) )

  ! =====================================
  ! Main loop over all the active species
  ! =====================================
  DO ns=1,iqmax
    ! Conversion mmr --> vmr
    vmr(:, :, ns) = mmr(:, :, ns) * dry_mass / adv_mass(ns)
    ! Conversion vmr --> molecules/cm3
    molec(:, :, ns) = vmr(:, :, ns) * convert(:, :)
  ENDDO
    
  DO ns=1,iqmax
    ! ========================
    !  Loss Terms by reactions
    ! ========================
    ! Loop over all the loss reactions
    DO kr=1,kreacl(ns)
      ir = idreacl(ns,kr)
      trat(:, :) = rates(:, :, ir)
      ! Loop over all the spec reacting with the active spec
      DO ire=1,nreactants(ir)
          IF (all_species(idspecl(ir,ire))%type == 'ac') THEN
            trat(:, :) = trat(:, :) * molec(:, :, idspecl(ir,ire))
          ELSE IF ((all_species(idspecl(ir,ire))%type == 'pr')) THEN
            trat(:, :) = trat(:, :) * refprescr(:, :, idspecl(ir,ire) - iqmax)
          END IF
      ENDDO

      loss(:, :, ns) = loss(:, :, ns) + trat(:, :)

      DO nk=1,klev-1
        DO nl=1,klon
          suml(nk,kr,kr,nl) = trat(nl, nk)
        ENDDO
      ENDDO

    ENDDO

    ! =============================
    ! Production terms by reactions
    ! =============================
    ! Loop over all the product reactions
    DO kr=1,kreacp(ns)
      ir = idreacp(ns,kr)
      trat(:, :) = rates(:, :, ir)
      ! Loop over all the spec reacting with the active spec
      DO ire=1,nreactants(ir)
        IF (all_species(idspecl(ir,ire))%type == 'ac') THEN
          trat(:, :) = trat(:, :) * molec(:, :, idspecl(ir,ire))
        ELSE IF ((all_species(idspecl(ir,ire))%type == 'pr')) THEN
          trat(:, :) = trat(:, :) * refprescr(:, :, idspecl(ir,ire) - iqmax)
        END IF
      ENDDO

      prod(:, :, ns) = prod(:, :, ns) + trat * stoi(ns,ir)

      DO nk=1,klev-1
        DO nl=1,klon
          sump(nk,kr,kr,nl) = trat(nl, nk) * stoi(ns,ir)
        ENDDO
      ENDDO

    ENDDO

    ! Prod and loss are in molecules.cm-3
    d_chem(:, :, ns) =  prod(:, :, ns) - loss(:, :, ns)

    ! =================================
    ! Production and loss by prodloss3d
    ! =================================
    IF (iprodmax > 0) THEN
      IF (all_species(ns)%prod == .TRUE.) THEN
        ! Assume that refprod is in molecules.cm-3
        d_chem(:, :, ns) = d_chem(:, :, ns) + refprod(:, :, all_species(ns)%idprod)
      END IF
    END IF
    ! Conversion molecules.cm-3 --> mmr
    d_chem(:, :, ns) =  d_chem(:, :, ns) * delt * adv_mass(ns) / (convert(:, :) * dry_mass)

  ENDDO


END subroutine comp_chemistry


! ===============================================================================================================
! ===============================================================================================================


subroutine comp_chemistry_tl(mmr, refprod, refprescr, refjrates, temp, pmid, mmr_tl, &
                        refprod_tl, refprescr_tl, delt,&
                        d_chem, d_chem_tl)

  USE SPECIES_NAME
  USE CONSTANTS
  USE variables
  USE parallel
  USE dimphy
  USE Write_Field_phy
  IMPLICIT NONE
  INCLUDE "dimensions.h"
  INCLUDE "paramet.h"

  !-----------------------------------------------------------------------
  !      General variables
  !-----------------------------------------------------------------------
  REAL, INTENT(in)    ::  mmr(klon,klev,iqmax)      ! mmr array
  REAL, INTENT(in)    ::  delt                     ! timestep
  REAL, INTENT(in)    ::  refprod(klon,klev,iprodmax)
  REAL, INTENT(in)    ::  refprescr(klon,klev,iprescrmax)
  REAL, INTENT(in)    ::  refjrates(klon,klev,ijratesmax)
  REAL, INTENT(in)    ::  temp(klon,klev)
  REAL, INTENT(in)    ::  pmid(klon,klev)
  REAL, INTENT(in)    ::  mmr_tl(klon,klev,iqmax) ! mmr array
  REAL, INTENT(in)    ::  refprod_tl(klon,klev,iprodmax)
  REAL, INTENT(in)    ::  refprescr_tl(klon,klev,iprescrmax)
  REAL, INTENT(out)   ::  d_chem(klon,klev,iqmax)   ! chemical increment
  REAL, INTENT(out)   ::  d_chem_tl(klon,klev,iqmax) ! chemical increment
  !-----------------------------------------------------------------------
  !      Local variables
  !-----------------------------------------------------------------------
  INTEGER :: ns
  REAL, ALLOCATABLE, DIMENSION(:,:)   :: invtemp    ! inverse temperature
  REAL, ALLOCATABLE, DIMENSION(:,:)   :: convert  !  vmr to molecule/cm3 conversion factor
  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: prod, loss ! prod/loss terms in VMR
  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: prod_tl, loss_tl ! prod/loss terms in VMR
  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: vmr, vmr_tl      ! vmr array
  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: molec, molec_tl    ! molecules/cm3 array
  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: rates     ! reaction rates

  INTEGER      :: kr, ire, ir, irea
  INTEGER      :: nr
  INTEGER      :: nk, nl
  REAL, ALLOCATABLE, DIMENSION(:,:) :: mrat, srat, trat

!  !  Initializations

  ALLOCATE(invtemp(klon,klev))
  ALLOCATE(prod(klon,klev,iqmax))
  ALLOCATE(loss(klon,klev,iqmax))
  ALLOCATE(prod_tl(klon,klev,iqmax))
  ALLOCATE(loss_tl(klon,klev,iqmax))
  ALLOCATE(vmr(klon,klev,iqmax))
  ALLOCATE(vmr_tl(klon,klev,iqmax))
  ALLOCATE(molec(klon,klev,iqmax))
  ALLOCATE(molec_tl(klon,klev,iqmax))
  ALLOCATE(convert(klon,klev))
  ALLOCATE(rates(klon,klev,nreac))
  ALLOCATE(srat(klon,klev))
  ALLOCATE(mrat(klon,klev))
  ALLOCATE(trat(klon,klev))

  loss = 0. ; prod = 0.
  loss_tl = 0. ;  prod_tl = 0.

  ! Compute rates for each reaction
  CALL comp_rates(temp, pmid, refjrates, rates)


  !--- Conversion factor from vmr ---> molecules.cm-3
  convert(:, :) = 10. * pmid(:, :) / ( boltz * temp(:, :) )

  ! =====================================
  ! Main loop over all the active species
  ! =====================================

  DO ns=1,iqmax
    ! Conversion mmr --> vmr
    vmr(:, :, ns) = mmr(:, :, ns) * dry_mass / adv_mass(ns)
    vmr_tl(:, :, ns) = mmr_tl(:, :, ns) * dry_mass / adv_mass(ns)
    ! Conversion vmr --> molecules/cm3
    molec(:, :, ns) = vmr(:, :, ns) * convert(:, :)
    molec_tl(:, :, ns) = vmr_tl(:, :, ns) * convert(:, :)
  ENDDO
  
  DO ns=1,iqmax
    ! =====================================
    ! Loss Terms by reactions
    ! =====================================
    ! Loss Terms TL
    DO kr=1,kreacl(ns)
      ir = idreacl(ns,kr)
      srat(:,:) = 0.

      DO ire=1,nreactants(ir)
        mrat(:, :) = rates(:, :, ir)
        DO irea=1,nreactants(ir)
          IF (irea == ire) THEN
            IF (all_species(idspecl(ir,irea))%type == 'ac') THEN
              mrat(:, :) = mrat(:, :) * molec_tl(:, :, idspecl(ir,irea))
            ELSE IF ((all_species(idspecl(ir,irea))%type == 'pr')) THEN
              mrat(:, :) = mrat(:, :) * refprescr_tl(:, :, idspecl(ir,irea) - iqmax)
            ENDIF
          ELSE
            IF (all_species(idspecl(ir,irea))%type == 'ac') THEN
              mrat(:, :) = mrat(:, :) * molec(:, :, idspecl(ir,irea))
            ELSE IF ((all_species(idspecl(ir,irea))%type == 'pr')) THEN
              mrat(:, :) = mrat(:, :) * refprescr(:, :, idspecl(ir,irea) - iqmax)
            ENDIF
          ENDIF
        ENDDO
        srat(:, :) = srat(:, :) + mrat(:,:)
      ENDDO
      loss_tl(:, :, ns) = loss_tl(:, :, ns) + srat(:, :)
    ENDDO

    !  Loss Terms FW
    DO kr=1,kreacl(ns)
      ir = idreacl(ns,kr)
      trat(:, :) = rates(:, :, ir)
      DO ire=1,nreactants(ir)
          IF (all_species(idspecl(ir,ire))%type == 'ac') THEN
            trat(:, :) = trat(:, :) * molec(:, :, idspecl(ir,ire))
          ELSE IF ((all_species(idspecl(ir,ire))%type == 'pr')) THEN
            trat(:, :) = trat(:, :) * refprescr(:, :, idspecl(ir,ire) - iqmax)
          END IF
      ENDDO
      loss(:, :, ns) = loss(:, :, ns) + trat(:, :)
    ENDDO

    ! =====================================
    ! Production Terms by reactions
    ! =====================================
    ! Production Terms TL
    DO kr=1,kreacp(ns)
      ir = idreacp(ns,kr)
      srat(:,:) = 0
      DO ire=1,nreactants(ir)
        mrat(:, :) = rates(:, :, ir)
        DO irea=1,nreactants(ir)
          IF (irea == ire) THEN
            IF (all_species(idspecl(ir,irea))%type == 'ac') THEN
              mrat(:, :) = mrat(:, :) * molec_tl(:, :, idspecl(ir,irea))
            ELSE IF ((all_species(idspecl(ir,irea))%type == 'pr')) THEN
              mrat(:, :) = mrat(:, :) * refprescr_tl(:, :, idspecl(ir,irea) - iqmax)
            ENDIF
          ELSE
            IF (all_species(idspecl(ir,irea))%type == 'ac') THEN
              mrat(:, :) = mrat(:, :) * molec(:, :, idspecl(ir,irea))
            ELSE IF ((all_species(idspecl(ir,irea))%type == 'pr')) THEN
              mrat(:, :) = mrat(:, :) * refprescr(:, :, idspecl(ir,irea) - iqmax)
            ENDIF
          ENDIF
        ENDDO
        srat(:, :) = srat(:, :) + mrat(:, :)
      ENDDO
      prod_tl(:, :, ns) = prod_tl(:, :, ns) + srat(:, :) * stoi(ns,ir)
    ENDDO

    !  Production Terms FW
    DO kr=1,kreacp(ns)
      ir = idreacp(ns,kr)
      trat(:, :) = rates(:, :, ir)
      DO ire=1,nreactants(ir)
        IF (all_species(idspecl(ir,ire))%type == 'ac') THEN
          trat(:, :) = trat(:, :) * molec(:, :, idspecl(ir,ire))
        ELSE IF ((all_species(idspecl(ir,ire))%type == 'pr')) THEN
          trat(:, :) = trat(:, :) * refprescr(:, :, idspecl(ir,ire) - iqmax)
        END IF
      ENDDO
      prod(:, :, ns) = prod(:, :, ns) + trat * stoi(ns,ir)
    ENDDO

    ! Prod and loss are in molecules.cm-3
    d_chem(:, :, ns) =  prod(:, :, ns) - loss(:, :, ns)
    d_chem_tl(:, :, ns) =  prod_tl(:, :, ns) - loss_tl(:, :, ns)

    ! =================================
    ! Production and loss by prodloss3d
    ! =================================
    IF (iprodmax > 0) THEN
      IF (all_species(ns)%prod == .TRUE.) THEN
       ! Assume that refprod is in molecules.cm-3
       d_chem_tl(:, :, ns) = d_chem_tl(:, :, ns) + refprod_tl(:, :, all_species(ns)%idprod)
       d_chem(:, :, ns) = d_chem(:, :, ns) + refprod(:, :, all_species(ns)%idprod)
     END IF
    END IF

    ! Conversion molecules.cm-3 --> mmr
    d_chem(:, :, ns) =  d_chem(:, :, ns) * delt * adv_mass(ns) / (convert(:, :) * dry_mass)
    d_chem_tl(:, :, ns) =  d_chem_tl(:, :, ns) * delt * adv_mass(ns) / (convert(:, :) * dry_mass)


  ENDDO


END subroutine comp_chemistry_tl





subroutine comp_chemistry_ad(mmr, refprescr, refprod, refjrates, temp, pmid, delt, &
                        d_chem_ad, mmr_ad, refprod_ad, refprescr_ad)

  USE SPECIES_NAME
  USE CONSTANTS
  USE variables
  USE parallel
  USE dimphy
  USE Write_Field_phy
  IMPLICIT NONE
  INCLUDE "dimensions.h"
  INCLUDE "paramet.h"

  !-----------------------------------------------------------------------
  !      General variables
  !-----------------------------------------------------------------------
  REAL, INTENT(in)     ::  mmr(klon,klev,iqmax)                ! mmr array
  REAL, INTENT(in)     ::  delt                                ! timestep
  REAL, INTENT(in)     ::  refprod(klon,klev,iprodmax)
  REAL, INTENT(in)     ::  refprescr(klon,klev,iprescrmax)
  REAL, INTENT(in)     ::  refjrates(klon,klev,ijratesmax)
  REAL, INTENT(in)     ::  temp(klon,klev)
  REAL, INTENT(in)     ::  pmid(klon,klev)
  REAL, INTENT(inout)  ::  refprod_ad(klon,klev,iprodmax)      ! prod3D species
  REAL, INTENT(inout)  ::  refprescr_ad(klon,klev,iprescrmax)  ! prescribed species
  REAL, INTENT(inout)  ::  d_chem_ad(klon,klev,iqmax)          ! adjoint increments
  REAL, INTENT(inout)  ::  mmr_ad(klon,klev,iqmax)             ! mmr adjoint array
  !-----------------------------------------------------------------------
  !      Local variables
  !-----------------------------------------------------------------------
  REAL, ALLOCATABLE, DIMENSION(:,:)   :: invtemp            ! inverse temperature
  REAL, ALLOCATABLE, DIMENSION(:,:)   :: convert            ! vmr to molecule/cm3 conversion factor
  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: prod_ad, loss_ad   ! prod/loss terms in VMR
  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: vmr, vmr_ad        ! vmr array
  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: molec, molec_ad    ! vmr array
  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: rates              ! reaction rates
  REAL, ALLOCATABLE, DIMENSION(:,:)   :: nrat               ! sum and multiplication of concentrations and rates


  INTEGER      :: kr,ire,ir, irea
  INTEGER      :: ns, nr

!  !  Initializations

  ALLOCATE(invtemp(klon,klev))
  ALLOCATE(vmr(klon,klev,iqmax))
  ALLOCATE(molec(klon,klev,iqmax))

  ALLOCATE(prod_ad(klon,klev,iqmax))
  ALLOCATE(loss_ad(klon,klev,iqmax))
  ALLOCATE(vmr_ad(klon,klev,iqmax))
  ALLOCATE(molec_ad(klon,klev,iqmax))

  ALLOCATE(convert(klon,klev))
  ALLOCATE(rates(klon,klev,nreac))
  ALLOCATE(nrat(klon,klev))

  ! Compute rates for each reaction
  CALL comp_rates(temp, pmid, refjrates, rates)

  !--- Conversion factor from vmr ---> molecules.cm-3
  convert(:, :) = 10. * pmid(:, :) / ( boltz * temp(:, :) )

  ! ===========
  ! FWD actions
  ! ===========
  DO ns=1,iqmax
    ! Conversion mmr --> vmr
    vmr(:, :, ns) = mmr(:, :, ns) * dry_mass / adv_mass(ns)
    ! Conversion vmr --> molecules/cm3
    molec(:, :, ns) = vmr(:, :, ns) * convert(:, :)
  ENDDO

  ! ====================
  ! Initialization of AD
  ! ====================
  prod_ad(:,:,:) = 0.
  loss_ad(:,:,:) = 0.
  vmr_ad(:,:,:) = 0.
  molec_ad(:,:,:) = 0.
  nrat(:, :) = 0.

  DO ns = iqmax, 1, -1

    d_chem_ad(:, :, ns) =  d_chem_ad(:, :, ns) * delt * adv_mass(ns) / (convert(:, :) * dry_mass)

    ! =================================
    ! Production and loss by prodloss3d
    ! =================================
    IF (iprodmax > 0) THEN
     IF (all_species(ns)%prod == .TRUE.) THEN
       refprod_ad(:, :, all_species(ns)%idprod) = refprod_ad(:, :, all_species(ns)%idprod) + d_chem_ad(:, :, ns)
     END IF
    END IF

    prod_ad(:, :, ns) = prod_ad(:, :, ns) + d_chem_ad(:, :, ns)
    loss_ad(:, :, ns) = loss_ad(:, :, ns) - d_chem_ad(:, :, ns)
    d_chem_ad(:, :, ns) = 0.


    ! =============
    ! Prod Terms AD
    ! =============
    DO kr=kreacp(ns), 1, -1
      ir = idreacp(ns,kr)
      DO ire = nreactants(ir), 1, -1
        nrat(:, :) = rates(:, :, ir)
        DO irea = nreactants(ir), 1, -1
          IF (irea /= ire) THEN
              IF (all_species(idspecl(ir,irea))%type == 'ac') THEN
                nrat(:, :) = nrat(:, :) * molec(:, :, idspecl(ir,irea))
              ELSE IF ((all_species(idspecl(ir,irea))%type == 'pr')) THEN
                nrat(:, :) = nrat(:, :) * refprescr(:, :, idspecl(ir,irea) - iqmax)
              ENDIF
          ENDIF
        ENDDO
        DO irea = nreactants(ir), 1, -1
          IF (irea == ire) THEN
            IF (all_species(idspecl(ir,irea))%type == 'ac') THEN
              molec_ad(:, :, idspecl(ir,irea)) = molec_ad(:, :, idspecl(ir,irea)) &
                      + nrat(:,:) * stoi(ns,ir) * prod_ad(:, :, ns)
            ELSE IF ((all_species(idspecl(ir,irea))%type == 'pr')) THEN
              refprescr_ad(:, :, idspecl(ir,irea) - iqmax) = refprescr_ad(:, :, idspecl(ir,irea) - iqmax) &
                      + nrat(:,:) * stoi(ns,ir) * prod_ad(:, :, ns)
            ENDIF
          ENDIF
        ENDDO
      ENDDO
    ENDDO
    prod_ad(:, :, ns) = 0.
    
    ! =============
    ! Loss Terms AD
    ! =============
    DO kr=kreacl(ns), 1, -1
      ir = idreacl(ns,kr)
      DO ire = nreactants(ir), 1, -1
        nrat(:, :) = rates(:, :, ir)
        DO irea = nreactants(ir), 1, -1
          IF (irea /= ire) THEN
              IF (all_species(idspecl(ir,irea))%type == 'ac') THEN
                nrat(:, :) = nrat(:, :) * molec(:, :, idspecl(ir,irea))
              ELSE IF ((all_species(idspecl(ir,irea))%type == 'pr')) THEN
                nrat(:, :) = nrat(:, :) * refprescr(:, :, idspecl(ir,irea) - iqmax)
              ENDIF
          ENDIF
        ENDDO
        DO irea = nreactants(ir), 1, -1
          IF (irea == ire) THEN
            IF (all_species(idspecl(ir,irea))%type == 'ac') THEN
              molec_ad(:, :, idspecl(ir,irea)) = molec_ad(:, :, idspecl(ir,irea)) &
                      + nrat(:,:) * loss_ad(:, :, ns)
            ELSE IF ((all_species(idspecl(ir,irea))%type == 'pr')) THEN
              refprescr_ad(:, :, idspecl(ir,irea) - iqmax) = refprescr_ad(:, :, idspecl(ir,irea) - iqmax) &
                      + nrat(:,:) * loss_ad(:, :, ns)
            ENDIF
          ENDIF
        ENDDO
      ENDDO
    ENDDO
    
    loss_ad(:, :, ns) = 0.
  
  ENDDO

  DO ns = iqmax, 1, -1
    vmr_ad(:, :, ns) = vmr_ad(:, :, ns) + molec_ad(:, :, ns) * convert(:, :)
    molec_ad(:, :, ns) = 0.
    mmr_ad(:, :, ns) = mmr_ad(:, :, ns) + vmr_ad(:, :, ns) * dry_mass / adv_mass(ns)
    vmr_ad(:, :, ns) = 0.
    
  ENDDO


END subroutine comp_chemistry_ad


