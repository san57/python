! ==========================================================================
! LMDZT - Atmospheric tracer transport, forward, tangent-linear and adjoint
!             for use within the PYVAR inversion system
!
! Copyright Laboratoire des Sciences du Climat et de l'Environnement (LSCE)
!           Unite mixte CEA-CNRS-UVSQ
! Initial version from the LMDZ.3.3 model, developed by IPSL
!
! Code manager:
! Frederic Chevallier, LSCE, frederic.chevallier@lsce.ipsl.fr
!
! This software is governed by the CeCILL  license under French law and
! abiding by the rules of distribution of free software.  You can  use,
! modify and/ or redistribute the software under the terms of the CeCILL
! license as circulated by CEA, CNRS and INRIA at the following URL
! "http://www.cecill.info".
!
! As a counterpart to the access to the source code and  rights to copy,
! modify and redistribute granted by the license, users are provided only
! with a limited warranty  and the software's author,  the holder of the
! economic rights,  and the successive licensors  have only  limited
! liability.
!
! In this respect, the user's attention is drawn to the risks associated
! with loading,  using,  modifying and/or developing or reproducing the
! software by the user in light of its specific status of free software,
! that may mean  that it is complicated to manipulate,  and  that  also
! therefore means  that it is reserved for developers  and  experienced
! professionals having in-depth computer knowledge. Users are therefore
! encouraged to load and test the software's suitability as regards their
! requirements in conditions enabling the security of their systems and/or
! data to be ensured and,  more generally, to use and operate it in the
! same conditions as regards security.
!
! The fact that you are presently reading this means that you have had
! knowledge of the CeCILL license and that you accept its terms.
! ==================================================================

SUBROUTINE dynstate0(fichnom,q,time)
  
  USE species_name
  
  IMPLICIT NONE
  
  !=======================================================================
  !
  !   Auteur:  P. Le Van / L.Fairhead
  !   -------
  !
  !   objet:
  !   ------
  !
  !   Lecture de l'etat initial
  !
  !=======================================================================
  !-----------------------------------------------------------------------
  !   Declarations:
  !   -------------
  
 include "dimensions.h"
 include "paramet.h"
 include "temps.h"
 include "comconst.h"
 include "comvert.h"
 include "comgeom.h"
 include "ener.h"
 include "netcdf.inc"
 include "description.h"
 include "serre.h"
 include "logic.h"
  
  !   Arguments:
  !   ----------
  
  CHARACTER(len=*) :: fichnom
  REAL :: vcov(ip1jm,llm),ucov(ip1jmp1,llm),teta(ip1jmp1,llm)
  REAL :: q(ip1jmp1,llm,iqmax),masse(ip1jmp1,llm)
  REAL :: ps(ip1jmp1),phis(ip1jmp1)
  
  REAL :: time
  
  !   Variables 
  !
  INTEGER,parameter :: length=100
  integer :: iq
  REAL :: tab_cntrl(length) ! tableau des parametres du run
  INTEGER :: ierr, nid, nvarid
  CHARACTER(len=3) :: str3
  
  !-----------------------------------------------------------------------
  
  !  Ouverture NetCDF du fichier etat initial
  
  ierr = NF_OPEN (fichnom, NF_NOWRITE,nid)
  IF (ierr/=NF_NOERR) THEN
      WRITE(6,*)' Failure in opening start.nc'
      WRITE(6,*)' ierr = ', ierr
      CALL ABORT
  ENDIF
  
  !
  ierr = NF_INQ_VARID (nid, "controle", nvarid)
  IF (ierr /= NF_NOERR) THEN
      PRINT*, "dynstate0: Field <controle> is absent"
      CALL abort
  ENDIF
#ifdef NC_DOUBLE
  ierr = NF_GET_VAR_DOUBLE(nid, nvarid, tab_cntrl)
#else
  ierr = NF_GET_VAR_REAL(nid, nvarid, tab_cntrl)
#endif
  IF (ierr /= NF_NOERR) THEN
      PRINT*, "dynstate0: Reading failure for <controle>"
      CALL abort
  ENDIF
  
  im         = tab_cntrl(1)
  jm         = tab_cntrl(2)
  lllm       = tab_cntrl(3)
  day_ini    = tab_cntrl(4)
  year_ini   = tab_cntrl(5)
  rad        = tab_cntrl(6)
  omeg       = tab_cntrl(7)
  g          = tab_cntrl(8)
  cpp        = tab_cntrl(9)
  kappa      = tab_cntrl(10)
  daysec     = tab_cntrl(11)
  dtvr       = tab_cntrl(12)
  etot0      = tab_cntrl(13)
  ptot0      = tab_cntrl(14)
  ztot0      = tab_cntrl(15)
  stot0      = tab_cntrl(16)
  ang0       = tab_cntrl(17)
  pa         = tab_cntrl(18)
  preff      = tab_cntrl(19)
  !
  clon       = tab_cntrl(20)
  clat       = tab_cntrl(21)
  grossismx  = tab_cntrl(22)
  grossismy  = tab_cntrl(23)
  !
  IF ( tab_cntrl(24)==1. )  THEN
      fxyhypb  = .TRUE.
      dzoomx   = tab_cntrl(25)
      dzoomy   = tab_cntrl(26)
      taux     = tab_cntrl(28)
      tauy     = tab_cntrl(29)
  ELSE
      fxyhypb = .FALSE.
      ysinus  = .FALSE.
      IF( tab_cntrl(27)==1. ) ysinus = .TRUE. 
  ENDIF
  !   .................................................................
  !
  !
  PRINT*,'rad,omeg,g,cpp,kappa',rad,omeg,g,cpp,kappa
  
  IF(   im/=iim           )  THEN
      PRINT 1,im,iim
      STOP
  ELSE  IF( jm/=jjm       )  THEN
      PRINT 2,jm,jjm
      STOP
  ELSE  IF( lllm/=llm     )  THEN
      PRINT 3,lllm,llm
      STOP
  ENDIF

  ierr = NF_INQ_VARID (nid, "temps", nvarid)
  IF (ierr /= NF_NOERR) THEN
      PRINT*, "dynstate0: Field<temps> is absent"
      CALL abort
  ENDIF
#ifdef NC_DOUBLE
  ierr = NF_GET_VAR_DOUBLE(nid, nvarid, time)
#else
  ierr = NF_GET_VAR_REAL(nid, nvarid, time)
#endif
  IF (ierr /= NF_NOERR) THEN
      PRINT*, "dynstate0: Lecture echouee <temps>"
      CALL abort
  ENDIF


  IF(iqmax>=1) THEN
      DO iq=1,iqmax
        str3(1:1)='q'
        IF(iq>99) THEN
            PRINT*, "Too many tracers"
            CALL abort
        ENDIF
        WRITE(str3(2:3),'(i2.2)') restart_ids(iq)

        WRITE(*,'(i2.2)') restart_ids(iq)
        ierr =  NF_INQ_VARID (nid, str3, nvarid)
        IF (ierr /= NF_NOERR) THEN
            PRINT*, "dynstate0: Field <"//str3//"> is absent"
            PRINT*, "          Initialized at zero"
            CALL initial0(ijp1llm,q(1,1,iq))
        ELSE
#ifdef NC_DOUBLE
            ierr = NF_GET_VAR_DOUBLE(nid, nvarid, q(1,1,iq))
#else
            ierr = NF_GET_VAR_REAL(nid, nvarid, q(1,1,iq))
#endif
            IF (ierr /= NF_NOERR) THEN
                PRINT*, "dynstate0: Reading failed for "//str3
                CALL abort
            ENDIF
        ENDIF
      ENDDO
  ENDIF
  
  ierr = NF_CLOSE(nid)
  
  day_ini=day_ini+INT(time)
  time=time-INT(time)
  
1 FORMAT(//10x,'la valeur de im =',i4,2x,'lue sur le fichier de dem &
     &arrage est differente de la valeur parametree iim =',i4//)
2 FORMAT(//10x,'la valeur de jm =',i4,2x,'lue sur le fichier de dem &
     &arrage est differente de la valeur parametree jjm =',i4//)
3 FORMAT(//10x,'la valeur de lmax =',i4,2x,'lue sur le fichier dema &
     &rrage est differente de la valeur parametree llm =',i4//)
4 FORMAT(//10x,'la valeur de dtrv =',i4,2x,'lue sur le fichier dema &
     &rrage est differente de la valeur  dtinteg =',i4//)
  
  RETURN
END SUBROUTINE dynstate0


SUBROUTINE dynstate1(fichnom)
  
  USE species_name
  
  IMPLICIT NONE
  
  !=======================================================================
  !
  !   Auteur:  P. Le Van / L.Fairhead
  !   -------
  !
  !   objet:
  !   ------
  !
  !   Lecture de l'etat initial
  !
  !=======================================================================
  !-----------------------------------------------------------------------
  !   Declarations:
  !   -------------
  
 include "dimensions.h"
 include "paramet.h"
 include "temps.h"
 include "comconst.h"
 include "comvert.h"
 include "comgeom.h"
 include "ener.h"
 include "netcdf.inc"
 include "description.h"
 include "serre.h"
 include "logic.h"
  
  !   Arguments:
  !   ----------
  
  CHARACTER(len=*) :: fichnom
  REAL :: vcov(ip1jm,llm),ucov(ip1jmp1,llm),teta(ip1jmp1,llm)
  REAL :: q(ip1jmp1,llm,iqmax),masse(ip1jmp1,llm)
  REAL :: ps(ip1jmp1),phis(ip1jmp1)
  
  REAL :: time
  
  !   Variables
  !
  INTEGER,parameter :: length=100
  integer :: iq
  REAL :: tab_cntrl(length) ! tableau des parametres du run
  INTEGER :: ierr, nid, nvarid
  CHARACTER(len=3) :: str3
  
  !-----------------------------------------------------------------------
  
  !  Ouverture NetCDF du fichier etat initial
  
  ierr = NF_OPEN (fichnom, NF_NOWRITE,nid)
  IF (ierr/=NF_NOERR) THEN
      WRITE(6,*)' Failure in opening start.nc'
      WRITE(6,*)' ierr = ', ierr
      CALL ABORT
  ENDIF
  
  !
  ierr = NF_INQ_VARID (nid, "controle", nvarid)
  IF (ierr /= NF_NOERR) THEN
      PRINT*, "dynstate0: Field <controle> is absent"
      CALL abort
  ENDIF
#ifdef NC_DOUBLE
  ierr = NF_GET_VAR_DOUBLE(nid, nvarid, tab_cntrl)
#else
  ierr = NF_GET_VAR_REAL(nid, nvarid, tab_cntrl)
#endif
  IF (ierr /= NF_NOERR) THEN
      PRINT*, "dynstate0: Reading failure for <controle>"
      CALL abort
  ENDIF
  
  im         = tab_cntrl(1)
  jm         = tab_cntrl(2)
  lllm       = tab_cntrl(3)
  day_ini    = tab_cntrl(4)
  year_ini   = tab_cntrl(5)
  rad        = tab_cntrl(6)
  omeg       = tab_cntrl(7)
  g          = tab_cntrl(8)
  cpp        = tab_cntrl(9)
  kappa      = tab_cntrl(10)
  daysec     = tab_cntrl(11)
  dtvr       = tab_cntrl(12)
  etot0      = tab_cntrl(13)
  ptot0      = tab_cntrl(14)
  ztot0      = tab_cntrl(15)
  stot0      = tab_cntrl(16)
  ang0       = tab_cntrl(17)
  pa         = tab_cntrl(18)
  preff      = tab_cntrl(19)
  !
  clon       = tab_cntrl(20)
  clat       = tab_cntrl(21)
  grossismx  = tab_cntrl(22)
  grossismy  = tab_cntrl(23)
  !
  IF ( tab_cntrl(24)==1. )  THEN
      fxyhypb  = .TRUE.
      dzoomx   = tab_cntrl(25)
      dzoomy   = tab_cntrl(26)
      taux     = tab_cntrl(28)
      tauy     = tab_cntrl(29)
  ELSE
      fxyhypb = .FALSE.
      ysinus  = .FALSE.
      IF( tab_cntrl(27)==1. ) ysinus = .TRUE.
  ENDIF
  !   .................................................................
  !
  !
  PRINT*,'rad,omeg,g,cpp,kappa',rad,omeg,g,cpp,kappa
  
  IF(   im/=iim           )  THEN
      PRINT 1,im,iim
      STOP
  ELSE  IF( jm/=jjm       )  THEN
      PRINT 2,jm,jjm
      STOP
  ELSE  IF( lllm/=llm     )  THEN
      PRINT 3,lllm,llm
      STOP
  ENDIF
  
  ierr = NF_INQ_VARID (nid, "rlonu", nvarid)
  IF (ierr /= NF_NOERR) THEN
      PRINT*, "dynstate0: Field<rlonu> is absent"
      CALL abort
  ENDIF
#ifdef NC_DOUBLE
  ierr = NF_GET_VAR_DOUBLE(nid, nvarid, rlonu)
#else
  ierr = NF_GET_VAR_REAL(nid, nvarid, rlonu)
#endif
  IF (ierr /= NF_NOERR) THEN
      PRINT*, "dynstate0: Reading failed for <rlonu>"
      CALL abort
  ENDIF
  
  ierr = NF_INQ_VARID (nid, "rlatu", nvarid)
  IF (ierr /= NF_NOERR) THEN
      PRINT*, "dynstate0: Field<rlatu> is absent"
      CALL abort
  ENDIF
#ifdef NC_DOUBLE
  ierr = NF_GET_VAR_DOUBLE(nid, nvarid, rlatu)
#else
  ierr = NF_GET_VAR_REAL(nid, nvarid, rlatu)
#endif
  IF (ierr /= NF_NOERR) THEN
      PRINT*, "dynstate0: Reading failed for <rlatu>"
      CALL abort
  ENDIF
  
  ierr = NF_INQ_VARID (nid, "rlonv", nvarid)
  IF (ierr /= NF_NOERR) THEN
      PRINT*, "dynstate0: Field<rlonv> is absent"
      CALL abort
  ENDIF
#ifdef NC_DOUBLE
  ierr = NF_GET_VAR_DOUBLE(nid, nvarid, rlonv)
#else
  ierr = NF_GET_VAR_REAL(nid, nvarid, rlonv)
#endif
  IF (ierr /= NF_NOERR) THEN
      PRINT*, "dynstate0: Reading failed for <rlonv>"
      CALL abort
  ENDIF
  
  ierr = NF_INQ_VARID (nid, "rlatv", nvarid)
  IF (ierr /= NF_NOERR) THEN
      PRINT*, "dynstate0: Field<rlatv> is absent"
      CALL abort
  ENDIF
#ifdef NC_DOUBLE
  ierr = NF_GET_VAR_DOUBLE(nid, nvarid, rlatv)
#else
  ierr = NF_GET_VAR_REAL(nid, nvarid, rlatv)
#endif
  IF (ierr /= NF_NOERR) THEN
      PRINT*, "dynstate0: Reading failed for rlatv"
      CALL abort
  ENDIF
  
  ierr = NF_INQ_VARID (nid, "ap", nvarid)
  IF (ierr /= NF_NOERR) THEN
      PRINT*, "dynstate0: Field<ap> is absent"
      CALL abort
  ENDIF
#ifdef NC_DOUBLE
  ierr = NF_GET_VAR_DOUBLE(nid, nvarid, ap)
#else
  ierr = NF_GET_VAR_REAL(nid, nvarid, ap)
#endif
  IF (ierr /= NF_NOERR) THEN
      PRINT*, "dynstate0: Reading failed for <ap>"
      CALL abort
  ENDIF

  ierr = NF_INQ_VARID (nid, "bp", nvarid)
  IF (ierr /= NF_NOERR) THEN
      PRINT*, "dynstate0: Field<bp> is absent"
      CALL abort
  ENDIF
#ifdef NC_DOUBLE
  ierr = NF_GET_VAR_DOUBLE(nid, nvarid, bp)
#else
  ierr = NF_GET_VAR_REAL(nid, nvarid, bp)
#endif
  IF (ierr /= NF_NOERR) THEN
      PRINT*, "dynstate0: Reading failed for <bp>"
      CALL abort
  ENDIF

  
  ierr = NF_INQ_VARID (nid, "cu", nvarid)
  IF (ierr /= NF_NOERR) THEN
      PRINT*, "dynstate0: Field<cu> is absent"
      CALL abort
  ENDIF
#ifdef NC_DOUBLE
  ierr = NF_GET_VAR_DOUBLE(nid, nvarid, cu)
#else
  ierr = NF_GET_VAR_REAL(nid, nvarid, cu)
#endif
  IF (ierr /= NF_NOERR) THEN
      PRINT*, "dynstate0: Reading failed for <cu>"
      CALL abort
  ENDIF
  
  ierr = NF_INQ_VARID (nid, "cv", nvarid)
  IF (ierr /= NF_NOERR) THEN
      PRINT*, "dynstate0: Field<cv> is absent"
      CALL abort
  ENDIF
#ifdef NC_DOUBLE
  ierr = NF_GET_VAR_DOUBLE(nid, nvarid, cv)
#else
  ierr = NF_GET_VAR_REAL(nid, nvarid, cv)
#endif
  IF (ierr /= NF_NOERR) THEN
      PRINT*, "dynstate0: Reading failed for <cv>"
      CALL abort
  ENDIF
  
  ierr = NF_INQ_VARID (nid, "aire", nvarid)
  IF (ierr /= NF_NOERR) THEN
      PRINT*, "dynstate0: Field<aire> is absent"
      CALL abort
  ENDIF
#ifdef NC_DOUBLE
  ierr = NF_GET_VAR_DOUBLE(nid, nvarid, area)
#else
  ierr = NF_GET_VAR_REAL(nid, nvarid, area)
#endif
  IF (ierr /= NF_NOERR) THEN
      PRINT*, "dynstate0: Reading failed for <aire>"
      CALL abort
  ENDIF
  
  ierr = NF_CLOSE(nid)
  
1 FORMAT(//10x,'la valeur de im =',i4,2x,'lue sur le fichier de dem &
     &arrage est differente de la valeur parametree iim =',i4//)
2 FORMAT(//10x,'la valeur de jm =',i4,2x,'lue sur le fichier de dem &
     &arrage est differente de la valeur parametree jjm =',i4//)
3 FORMAT(//10x,'la valeur de lmax =',i4,2x,'lue sur le fichier dema &
     &rrage est differente de la valeur parametree llm =',i4//)
4 FORMAT(//10x,'la valeur de dtrv =',i4,2x,'lue sur le fichier dema &
     &rrage est differente de la valeur  dtinteg =',i4//)
  
  RETURN
END SUBROUTINE dynstate1







