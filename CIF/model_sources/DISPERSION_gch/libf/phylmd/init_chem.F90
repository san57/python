SUBROUTINE findspec(charspec,nospec)

  !  Returns the index of the species, if found among active ones,
  !  and 0 if the species is not found.
  !  INPUT  : CHARSPEC         The species name
  !           SPECIES          Active species names array
  !  OUTPUT : NOSPEC           The species index

    USE SPECIES_NAME
    IMPLICIT NONE

  !*****************************************************************************************
  ! subroutine arguments
  CHARACTER(len=*),INTENT(in) :: charspec
  INTEGER,INTENT(out)             :: nospec
  ! local variables
  INTEGER :: is
  !*****************************************************************************************

  nospec = 0
  DO is=1, iallqmax
     IF(charspec == all_species(is)%name) THEN
        nospec = is
        EXIT
     ENDIF
  ENDDO

END SUBROUTINE findspec





SUBROUTINE init_chem()

    USE SPECIES_NAME
    IMPLICIT NONE

    INTEGER :: ifnallspec, ifnacspec, ifnprescr, ifnprodl, ifndep, ifnchem, ifnstoi, ifnrates, ifnjrates
    INTEGER :: iq, iall, ipr, ipl, idp, ireac, ij, ire, ip, i, j, tr
    INTEGER :: nprods
    INTEGER :: idloss, idprod
    CHARACTER(len = 15), DIMENSION(1000) :: reactants, prods
    REAL, DIMENSION(4) :: coeff


    NAMELIST /args/  fnacspec, fnallspec, fnprescr, fnprodl, fndep , fnchem, &
                     fnstoi, fnrates, fnjrates, fnfamilies,  &
                     iqmax, iallqmax, iprescrmax, iprodmax, idepmax, nreac, ijratesmax


    ! Reading namelist and allocating space
    OPEN(unit=1, file='chemical_scheme.nml', status='OLD')
    READ(1, nml=args, iostat=ioerr)

    if (ioerr/=0) then
        WRITE(*,*) 'readnml: problem reading the chemistry files'
    end if

    ALLOCATE( all_species(iallqmax))
    ALLOCATE( restart_ids(iqmax))
    ALLOCATE( adv_mass(iqmax))
    ALLOCATE( species(iqmax) )
    ALLOCATE( prescr_species(iprescrmax) )
    ALLOCATE( prodloss_species(iprodmax) )
    ALLOCATE( dep_species(idepmax) )
    ALLOCATE( jrates(ijratesmax))
    ALLOCATE( kreacl(iallqmax) )
    ALLOCATE( kreacp(iallqmax) )
    ALLOCATE( nreactants(nreac) )
    ALLOCATE( idreacl(iallqmax, nreac) )
    ALLOCATE( idreacp(iallqmax, nreac) )
    ALLOCATE( stoi(iallqmax, nreac) )
    ALLOCATE( idspecl(nreac, 10) )
    ALLOCATE( idtyperate(nreac) )
    ALLOCATE( tabrate(ntabmax, nreac))

    kreacl = 0.
    kreacp = 0.
    idreacl = 0.
    idreacp = 0.
    idspecl = 0.
    stoi = 1d0
    tabrate = 0.
    
    ! Change ltabrate to properly retrieve the right number of coefficients used for the reaction
    ltabrate(1:25) = (/1, 2, 3, 3, 1, 7, 4, 8, 8, 4, 2, 8, 1, 7, 1, 6, 3, 4, 2, 5, 4, 5, 5, 5, 1/)
    PRINT*, 'If bug around here, think about changing ltabrate coefficients...'
    PRINT*, '... A new chemistry reaction type can be the problem '

    !  Reading active/output species names
    OPEN(unit=330, file=fnacspec, status='OLD',form='formatted', action='read')
    DO iq = 1, iqmax
        READ(330, *, end = 1000)species(iq)%name, restart_ids(iq), adv_mass(iq)
        species(iq)%id = iq
    ENDDO
    1000 CONTINUE
    CLOSE(330)


    !  Reading all species names
    OPEN(unit=330, file=fnallspec, status='OLD',form='formatted', action='read')
    DO iall = 1, iallqmax
        READ(330, *, end = 1001)all_species(iall)%name, all_species(iall)%type
        all_species(iall)%id = iall
        all_species(iall)%dep = .FALSE.
        all_species(iall)%iddep = 0
        all_species(iall)%prod = .FALSE.
        all_species(iall)%idprod = 0
    ENDDO
    1001 CONTINUE
    CLOSE(330)


    !  Reading prescribed species names
    OPEN(unit=330, file=fnprescr, status='OLD',form='formatted', action='read')
    DO ipr = 1, iprescrmax
        READ(330, *, end = 1002)prescr_species(ipr)%name
        prescr_species(ipr)%id = ipr
    ENDDO
    1002 CONTINUE
    CLOSE(330)


    !  Reading prodloss species names
    OPEN(unit=330, file=fnprodl, status='OLD',form='formatted', action='read')
    DO ipl = 1, iprodmax
        READ(330, *, end = 1003)prodloss_species(ipl)%name
        prodloss_species(ipl)%id = ipl
    ENDDO
    1003 CONTINUE
    CLOSE(330)


    !  Reading prodloss species names
    OPEN(unit=330, file=fndep, status='OLD',form='formatted', action='read')
    DO idp = 1, idepmax
        READ(330, *, end = 1004)dep_species(idp)%name
        dep_species(idp)%id = idp
    ENDDO
    1004 CONTINUE
    CLOSE(330)

    ! Make DEP and PROD associations
    DO iall = 1, iallqmax
        DO ipl = 1, iprodmax
            IF (all_species(iall)%name == prodloss_species(ipl)%name) THEN
                all_species(iall)%prod = .TRUE.
                all_species(iall)%idprod = ipl
            END IF
        END DO
        DO idp = 1, idepmax
            IF (all_species(iall)%name == dep_species(idp)%name) THEN
                all_species(iall)%dep = .TRUE.
                all_species(iall)%iddep = idp
            END IF
        END DO
    ENDDO


    !  Reading reaction addressing arrays of chemistry
    OPEN(unit=330, file=fnchem, status='OLD',form='formatted', action='read')

    DO ireac = 1, nreac
        READ(330, *, end = 1005)nreactants(ireac), (reactants(ire), ire = 1, nreactants(ireac)), &
                                    nprods, (prods(ip), ip = 1, nprods)

        !  Addressing the reactant list
        DO ire = 1, nreactants(ireac)
            CALL findspec(reactants(ire), idloss)
            IF (idloss > 0) THEN
                kreacl(idloss) = kreacl(idloss) + 1
                idreacl(idloss, kreacl(idloss)) = ireac
                idspecl(ireac, ire) = idloss
                END IF
        ENDDO

        !  Addressing the product list
        DO ip = 1, nprods
            CALL findspec(prods(ip), idprod)
            IF (idprod > 0) THEN
                kreacp(idprod) = kreacp(idprod) + 1
                idreacp(idprod, kreacp(idprod)) = ireac
            ENDIF
        ENDDO
    ENDDO
    1005 CONTINUE
    CLOSE(330)


    !  Reading and addressing stoichiometric coefficients
    OPEN(unit=330, file=fnstoi, status='OLD',form='formatted', action='read')
    DO i = 1, 1000
        READ(330, *, end = 1006)reactants(1), coeff, ireac
        CALL findspec(reactants(1), idprod)
        IF(idprod > 0) THEN
            stoi(idprod, ireac) = coeff(1)
        ENDIF
    ENDDO
    1006 CONTINUE
    CLOSE(330)


    !  Reading the reaction rates constants
    OPEN(unit=330, file=fnrates, status='OLD',form='formatted', action='read')
    DO ireac = 1, nreac
        READ(330, *, end = 1007)tr, &
                & (tabrate(i, ireac), i = 1, ltabrate(tr))
        idtyperate(ireac) = tr
    ENDDO
    1007 CONTINUE
    CLOSE(330)

    !  Reading the photolysis reaction rates constants
    OPEN(unit=330, file=fnjrates, status='OLD',form='formatted', action='read')
    DO ij = 1, ijratesmax
        READ(330, *, end = 1007) jrates(ij)%idreac, jrates(ij)%idj
    END DO
    
END SUBROUTINE init_chem
