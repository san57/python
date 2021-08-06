SUBROUTINE readchem &
   (irec,refprescr,refprod,refjrates,&
        temp,pmid,convoh)

USE SPECIES_NAME
USE CONSTANTS
USE parallel
IMPLICIT NONE

INCLUDE "netcdf.inc"
INCLUDE "dimensions.h"
INCLUDE "paramet.h"
INCLUDE "ajout.h"

INTEGER :: irec, iq
REAL, DIMENSION(ngridmx,llm,iprescrmax) :: refprescr
REAL, DIMENSION(ngridmx,llm,iprodmax)   :: refprod
REAL, DIMENSION(ngridmx,llm,ijratesmax) :: refjrates
REAL, DIMENSION(ngridmx,llm) :: temp, pmid
REAL, DIMENSION(iim,jjm+1,llm) :: tempo
INTEGER, PARAMETER :: klon=iim*(jjm-1)+2

INTEGER, SAVE :: ncidi ! file ID for kinetics
INTEGER, ALLOCATABLE, DIMENSION(:), SAVE :: ncidiprescr ! file ID for prescribed fluxes
INTEGER, ALLOCATABLE, DIMENSION(:), SAVE :: ncidiprod ! file ID for 3d production
INTEGER, ALLOCATABLE, DIMENSION(:), SAVE :: idprescr ! file ID for prescribed fluxes
INTEGER, ALLOCATABLE, DIMENSION(:), SAVE :: idprod ! file ID for 3d production
INTEGER, ALLOCATABLE, DIMENSION(:), SAVE :: idjrates ! file ID for photolysis constant reactions

INTEGER, SAVE :: idtemp, idpmid

INTEGER, DIMENSION(4) :: start, count

REAL :: rcode
LOGICAL :: convoh
INTEGER :: n,i,j,ig, iret

!
IF (irec == 0) THEN

  ! Allocating variables variables
  ALLOCATE(ncidiprescr(iprescrmax))
  ALLOCATE(ncidiprod(iprodmax))
  ALLOCATE(idprescr(iprescrmax))
  ALLOCATE(idprod(iprodmax))
  ALLOCATE(idjrates(ijratesmax))

  ! open chemical kinetics
  ncidi = NCOPN('kinetic.nc',NCNOWRIT,rcode)
  idtemp = NCVID(ncidi,'temp',rcode)
  PRINT*,'ncidi, idtemp', ncidi, idtemp
  idpmid = NCVID(ncidi,'pmid',rcode)
  PRINT*,'ncidi, idpmid', ncidi, idpmid

  ! open photolysis reactions constants
  DO iq=1,ijratesmax
    idjrates(iq) = NCVID(ncidi, jrates(iq)%idj, rcode)
  END DO
  
  ! open prescribed fields
  DO iq=1,iprescrmax
    ncidiprescr(iq) = NCOPN('prescr_'//trim(prescr_species(iq)%name)//'.nc',NCNOWRIT,rcode)
    idprescr(iq) = NCVID(ncidiprescr(iq),trim(prescr_species(iq)%name), rcode)
  END DO
  
  ! open prod/loss fields
  DO iq=1,iprodmax
    ncidiprod(iq) = NCOPN('prodloss_'//trim(prodloss_species(iq)%name)//'.nc',NCNOWRIT,rcode)
    idprod(iq) = NCVID(ncidiprod(iq),trim(prodloss_species(iq)%name)//'_prod', rcode)
  END DO
  
ELSE
  
  start(1)=1
  start(2)=1
  start(3)=1
  start(4)=irec
  count(1)=iim
  count(2)=jjm+1
  count(3)=llm
  count(4)=1
  
  ! Temperature and pressure
  iret = nf_get_vara_double(ncidi, idtemp, start, count, tempo)
  CALL check_err(iret)
  CALL gr_ecrit_fi(llm,klon,iim,jjm+1,tempo,temp)
  
  iret = nf_get_vara_double(ncidi, idpmid, start, count, tempo)
  CALL check_err(iret)
  CALL gr_ecrit_fi(llm,klon,iim,jjm+1,tempo,pmid)
  
  ! Prescribed species
  DO iq=1,iprescrmax
    iret = nf_get_vara_double(ncidiprescr(iq), idprescr(iq), start, count, tempo)
    CALL check_err(iret)
    CALL gr_ecrit_fi(llm,klon,iim,jjm+1,tempo,refprescr(:,:,iq))
    !  TODO : adapter si VMR = True
    ! Convert to molec.cm-3
    refprescr(:,:,iq) = refprescr(:,:,iq) * 10. * pmid(:,:)  / ( boltz * temp(:,:) )
  END DO

  ! Production and loss
  ! NB: production of HCHO by NMHC already in molec.cm-3.s-1
  DO iq = 1, iprodmax
    iret = nf_get_vara_double(ncidiprod(iq), idprod(iq), start, count, tempo)
    CALL check_err(iret)
    CALL gr_ecrit_fi(llm,klon,iim,jjm+1,tempo,refprod(:,:,iq))
  END DO

  DO iq=1,ijratesmax
    iret = nf_get_vara_double(ncidi, idjrates(iq), start, count, tempo)
    CALL check_err(iret)
    CALL gr_ecrit_fi(llm,klon,iim,jjm+1,tempo,refjrates(:,:,iq))
  END DO

  
ENDIF
RETURN
END

SUBROUTINE check_err(iret)
INTEGER :: iret
INCLUDE "netcdf.inc"
IF (iret /= NF_NOERR) THEN
 PRINT *, iret, nf_strerror(iret)
 STOP
ENDIF
END SUBROUTINE check_err
