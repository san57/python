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

SUBROUTINE dynrestart0(fichnom,idayref,yearref,phis)
  
  USE SPECIES_NAME
  
  IMPLICIT NONE
  !=======================================================================
  ! Ecriture du fichier de redemarrage sous format NetCDF (initialisation)
  !=======================================================================
  !   Declarations:
  !   -------------
  
  include "dimensions.h"
  include "paramet.h"
  include "comconst.h"
  include "comvert.h"
  include "comgeom.h"
  include "temps.h"
  include "ener.h"
  include "logic.h"
  include "netcdf.inc"
  include "description.h"
  include "serre.h"

  !   Arguments:
  !   ----------
  INTEGER*4 :: idayref,yearref
  REAL :: phis(ip1jmp1)
  CHARACTER(len=*) :: fichnom

  !   Local:
  !   ------
  INTEGER :: iq
  CHARACTER(len=3) :: str3
  INTEGER,parameter :: length=100
  REAL :: tab_cntrl(length) ! tableau des parametres du run
  INTEGER :: ierr
  CHARACTER(len=20) :: modname
  
  !   Variables locales pour NetCDF:
  !
  INTEGER :: dims2(2), dims3(3), dims4(4)
  INTEGER :: idim_index
  INTEGER :: idim_rlonu, idim_rlonv, idim_rlatu, idim_rlatv
  INTEGER :: idim_s, idim_sig
  INTEGER :: idim_tim
  INTEGER :: nid,nvarid
  
  REAL :: hours
  INTEGER :: yyears0,dday0, mmois0
  CHARACTER(len=30) :: unites


  !-----------------------------------------------------------------------
  modname='dynredem'
  
  yyears0 = yearref
  mmois0 = 1
  dday0 = 1
  hours = 0
  
  tab_cntrl(:) = 0.

  tab_cntrl(1)  = FLOAT(iim)
  tab_cntrl(2)  = FLOAT(jjm)
  tab_cntrl(3)  = FLOAT(llm)
  tab_cntrl(4)  = FLOAT(idayref)
  tab_cntrl(5)  = FLOAT(yearref)
  tab_cntrl(6)  = rad
  tab_cntrl(7)  = omeg
  tab_cntrl(8)  = g
  tab_cntrl(9)  = cpp
  tab_cntrl(10) = kappa
  tab_cntrl(11) = daysec
  tab_cntrl(12) = dtvr
  tab_cntrl(13) = etot0
  tab_cntrl(14) = ptot0
  tab_cntrl(15) = ztot0
  tab_cntrl(16) = stot0
  tab_cntrl(17) = ang0
  tab_cntrl(18) = pa
  tab_cntrl(19) = preff
  !
  !    .....    parameters for zooming      ......
  
  tab_cntrl(20)  = clon
  tab_cntrl(21)  = clat
  tab_cntrl(22)  = grossismx
  tab_cntrl(23)  = grossismy
  !
  IF ( fxyhypb )   THEN
      tab_cntrl(24) = 1.
      tab_cntrl(25) = dzoomx
      tab_cntrl(26) = dzoomy
      tab_cntrl(27) = 0.
      tab_cntrl(28) = taux
      tab_cntrl(29) = tauy
  ELSE
      tab_cntrl(24) = 0.
      tab_cntrl(25) = dzoomx
      tab_cntrl(26) = dzoomy
      tab_cntrl(27) = 0.
      tab_cntrl(28) = 0.
      tab_cntrl(29) = 0.
      IF( ysinus )  tab_cntrl(27) = 1.
  ENDIF
  !
  !    .........................................................
  !
  ! Creation du fichier:
  !
  ierr = NF_CREATE(fichnom, NF_CLOBBER, nid)
  IF (ierr/=NF_NOERR) THEN
      WRITE(6,*)" Pb d ouverture du fichier "//fichnom
      WRITE(6,*)' ierr = ', ierr
      CALL ABORT
  ENDIF
  !
  ! Preciser quelques attributs globaux:
  !
  ierr = NF_PUT_ATT_TEXT (nid, NF_GLOBAL, "title", 27, "Fichier demarrage dynamique")
  !
  ! Definir les dimensions du fichiers:
  !
  ierr = NF_DEF_DIM (nid, "index", length, idim_index)
  ierr = NF_DEF_DIM (nid, "rlonu", iip1, idim_rlonu)
  ierr = NF_DEF_DIM (nid, "rlatu", jjp1, idim_rlatu)
  ierr = NF_DEF_DIM (nid, "rlonv", iip1, idim_rlonv)
  ierr = NF_DEF_DIM (nid, "rlatv", jjm, idim_rlatv)
  ierr = NF_DEF_DIM (nid, "sigs", llm, idim_s)
  ierr = NF_DEF_DIM (nid, "sig", llmp1, idim_sig)
  ierr = NF_DEF_DIM (nid, "temps", NF_UNLIMITED, idim_tim)
  !
  ierr = NF_ENDDEF(nid) ! sortir du mode de definition
  !
  ! Definir et enregistrer certains champs invariants:
  !
  ierr = NF_REDEF (nid)
#ifdef NC_DOUBLE
  ierr = NF_DEF_VAR (nid,"controle",NF_DOUBLE,1,idim_index,nvarid)
#else
  ierr = NF_DEF_VAR (nid,"controle",NF_FLOAT,1,idim_index,nvarid)
#endif
  ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 22, "Parametres de controle")
  ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
  ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,tab_cntrl)
#else
  ierr = NF_PUT_VAR_REAL (nid,nvarid,tab_cntrl)
#endif
  !
  ierr = NF_REDEF (nid)
#ifdef NC_DOUBLE
  ierr = NF_DEF_VAR (nid,"rlonu",NF_DOUBLE,1,idim_rlonu,nvarid)
#else
  ierr = NF_DEF_VAR (nid,"rlonu",NF_FLOAT,1,idim_rlonu,nvarid)
#endif
  ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 23, "Longitudes des points U")
  ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
  ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,rlonu)
#else
  ierr = NF_PUT_VAR_REAL (nid,nvarid,rlonu)
#endif
  !
  ierr = NF_REDEF (nid)
#ifdef NC_DOUBLE
  ierr = NF_DEF_VAR (nid,"rlatu",NF_DOUBLE,1,idim_rlatu,nvarid)
#else
  ierr = NF_DEF_VAR (nid,"rlatu",NF_FLOAT,1,idim_rlatu,nvarid)
#endif
  ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 22, "Latitudes des points U")
  ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
  ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,rlatu)
#else
  ierr = NF_PUT_VAR_REAL (nid,nvarid,rlatu)
#endif
  !
  ierr = NF_REDEF (nid)
#ifdef NC_DOUBLE
  ierr = NF_DEF_VAR (nid,"rlonv",NF_DOUBLE,1,idim_rlonv,nvarid)
#else
  ierr = NF_DEF_VAR (nid,"rlonv",NF_FLOAT,1,idim_rlonv,nvarid)
#endif
  ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 23, "Longitudes des points V")
  ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
  ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,rlonv)
#else
  ierr = NF_PUT_VAR_REAL (nid,nvarid,rlonv)
#endif
  !
  ierr = NF_REDEF (nid)
#ifdef NC_DOUBLE
  ierr = NF_DEF_VAR (nid,"rlatv",NF_DOUBLE,1,idim_rlatv,nvarid)
#else
  ierr = NF_DEF_VAR (nid,"rlatv",NF_FLOAT,1,idim_rlatv,nvarid)
#endif
  ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 22, "Latitudes des points V")
  ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
  ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,rlatv)
#else
  ierr = NF_PUT_VAR_REAL (nid,nvarid,rlatv)
#endif
  !
  ierr = NF_REDEF (nid)
#ifdef NC_DOUBLE
  ierr = NF_DEF_VAR (nid,"nivsigs",NF_DOUBLE,1,idim_s,nvarid)
#else
  ierr = NF_DEF_VAR (nid,"nivsigs",NF_FLOAT,1,idim_s,nvarid)
#endif
  ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 28, "Numero naturel des couches s")
  ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
  ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,nivsigs)
#else
  ierr = NF_PUT_VAR_REAL (nid,nvarid,nivsigs)
#endif
  !
  ierr = NF_REDEF (nid)
#ifdef NC_DOUBLE
  ierr = NF_DEF_VAR (nid,"nivsig",NF_DOUBLE,1,idim_sig,nvarid)
#else
  ierr = NF_DEF_VAR (nid,"nivsig",NF_FLOAT,1,idim_sig,nvarid)
#endif
  ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 32, "Numero naturel des couches sigma")
  ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
  ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,nivsig)
#else
  ierr = NF_PUT_VAR_REAL (nid,nvarid,nivsig)
#endif
  !
  ierr = NF_REDEF (nid)
#ifdef NC_DOUBLE
  ierr = NF_DEF_VAR (nid,"ap",NF_DOUBLE,1,idim_sig,nvarid)
#else
  ierr = NF_DEF_VAR (nid,"ap",NF_FLOAT,1,idim_sig,nvarid)
#endif
  ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 26, "Coefficient A pour hybride")
  ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
  ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,ap)
#else
  ierr = NF_PUT_VAR_REAL (nid,nvarid,ap)
#endif
  !
  ierr = NF_REDEF (nid)
#ifdef NC_DOUBLE
  ierr = NF_DEF_VAR (nid,"bp",NF_DOUBLE,1,idim_sig,nvarid)
#else
  ierr = NF_DEF_VAR (nid,"bp",NF_FLOAT,1,idim_sig,nvarid)
#endif
  ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 26, "Coefficient B pour hybride")
  ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
  ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,bp)
#else
  ierr = NF_PUT_VAR_REAL (nid,nvarid,bp)
#endif
  !
  ierr = NF_REDEF (nid)
#ifdef NC_DOUBLE
  ierr = NF_DEF_VAR (nid,"presnivs",NF_DOUBLE,1,idim_s,nvarid)
#else
  ierr = NF_DEF_VAR (nid,"presnivs",NF_FLOAT,1,idim_s,nvarid)
#endif
  ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
  ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,presnivs)
#else
  ierr = NF_PUT_VAR_REAL (nid,nvarid,presnivs)
#endif
  !
  ! Coefficients de passage cov. <-> contra. <--> naturel
  !
  ierr = NF_REDEF (nid)
  dims2(1) = idim_rlonu
  dims2(2) = idim_rlatu
#ifdef NC_DOUBLE
  ierr = NF_DEF_VAR (nid,"cu",NF_DOUBLE,2,dims2,nvarid)
#else
  ierr = NF_DEF_VAR (nid,"cu",NF_FLOAT,2,dims2,nvarid)
#endif
  ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 29, "Coefficient de passage pour U")
  ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
  ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,cu)
#else
  ierr = NF_PUT_VAR_REAL (nid,nvarid,cu)
#endif
  !
  ierr = NF_REDEF (nid)
  dims2(1) = idim_rlonv
  dims2(2) = idim_rlatv
#ifdef NC_DOUBLE
  ierr = NF_DEF_VAR (nid,"cv",NF_DOUBLE,2,dims2,nvarid)
#else
  ierr = NF_DEF_VAR (nid,"cv",NF_FLOAT,2,dims2,nvarid)
#endif
  ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 29, "Coefficient de passage pour V")
  ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
  ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,cv)
#else
  ierr = NF_PUT_VAR_REAL (nid,nvarid,cv)
#endif
  !
  ! Aire de chaque maille:
  !
  ierr = NF_REDEF (nid)
  dims2(1) = idim_rlonv
  dims2(2) = idim_rlatu
#ifdef NC_DOUBLE
  ierr = NF_DEF_VAR (nid,"aire",NF_DOUBLE,2,dims2,nvarid)
#else
  ierr = NF_DEF_VAR (nid,"aire",NF_FLOAT,2,dims2,nvarid)
#endif
  ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 22, "Aires de chaque maille")
  ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
  ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,area)
#else
  ierr = NF_PUT_VAR_REAL (nid,nvarid,area)
#endif
  !
  ! Geopentiel au sol:
  !
  ierr = NF_REDEF (nid)
  dims2(1) = idim_rlonv
  dims2(2) = idim_rlatu
#ifdef NC_DOUBLE
  ierr = NF_DEF_VAR (nid,"phisinit",NF_DOUBLE,2,dims2,nvarid)
#else
  ierr = NF_DEF_VAR (nid,"phisinit",NF_FLOAT,2,dims2,nvarid)
#endif
  ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 19, "Geopotentiel au sol")
  ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
  ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,phis)
#else
  ierr = NF_PUT_VAR_REAL (nid,nvarid,phis)
#endif
  !
  ! Definir les variables pour pouvoir les enregistrer plus tard:
  !
  ierr = NF_REDEF (nid) ! entrer dans le mode de definition
  !
#ifdef NC_DOUBLE
  ierr = NF_DEF_VAR (nid,"temps",NF_DOUBLE,1,idim_tim,nvarid)
#else
  ierr = NF_DEF_VAR (nid,"temps",NF_FLOAT,1,idim_tim,nvarid)
#endif
  ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 19, "Temps de simulation")
  WRITE(unites,200)yyears0,mmois0,dday0
200 FORMAT('days since ',i4,'-',i2.2,'-',i2.2,' 00:00:00')
  ierr = NF_PUT_ATT_TEXT (nid, nvarid, "units", 30, &
     unites)
  
  !
  dims4(1) = idim_rlonu
  dims4(2) = idim_rlatu
  dims4(3) = idim_s
  dims4(4) = idim_tim
#ifdef NC_DOUBLE
  ierr = NF_DEF_VAR (nid,"ucov",NF_DOUBLE,4,dims4,nvarid)
#else
  ierr = NF_DEF_VAR (nid,"ucov",NF_FLOAT,4,dims4,nvarid)
#endif
  ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 9, "Vitesse U")
  !
  dims4(1) = idim_rlonv
  dims4(2) = idim_rlatv
  dims4(3) = idim_s
  dims4(4) = idim_tim
#ifdef NC_DOUBLE
  ierr = NF_DEF_VAR (nid,"vcov",NF_DOUBLE,4,dims4,nvarid)
#else
  ierr = NF_DEF_VAR (nid,"vcov",NF_FLOAT,4,dims4,nvarid)
#endif
  ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 9, "Vitesse V")
  !
  dims4(1) = idim_rlonv
  dims4(2) = idim_rlatu
  dims4(3) = idim_s
  dims4(4) = idim_tim
#ifdef NC_DOUBLE
  ierr = NF_DEF_VAR (nid,"teta",NF_DOUBLE,4,dims4,nvarid)
#else
  ierr = NF_DEF_VAR (nid,"teta",NF_FLOAT,4,dims4,nvarid)
#endif
  ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 11, "Temperature")
  !
  dims4(1) = idim_rlonv
  dims4(2) = idim_rlatu
  dims4(3) = idim_s
  dims4(4) = idim_tim
  
  IF(iqmax>=1) THEN
      DO iq=1,iqmax
        IF (iq>99) THEN
            PRINT*, "Too many tracers"
            CALL abort
        ELSE
            str3(1:1)='q'
            WRITE(str3(2:3),'(i2.2)') restart_ids(iq)
#ifdef NC_DOUBLE
            ierr = NF_DEF_VAR (nid,str3,NF_DOUBLE,4,dims4,nvarid)
#else
            ierr = NF_DEF_VAR (nid,str3,NF_FLOAT,4,dims4,nvarid)
#endif
            ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 12, "Traceurs "//str3)
        ENDIF
      ENDDO
  ENDIF
  !
  dims4(1) = idim_rlonv
  dims4(2) = idim_rlatu
  dims4(3) = idim_s
  dims4(4) = idim_tim
#ifdef NC_DOUBLE
  ierr = NF_DEF_VAR (nid,"masse",NF_DOUBLE,4,dims4,nvarid)
#else
  ierr = NF_DEF_VAR (nid,"masse",NF_FLOAT,4,dims4,nvarid)
#endif
  ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 12, "C est quoi ?")
  !
  dims3(1) = idim_rlonv
  dims3(2) = idim_rlatu
  dims3(3) = idim_tim
#ifdef NC_DOUBLE
  ierr = NF_DEF_VAR (nid,"ps",NF_DOUBLE,3,dims3,nvarid)
#else
  ierr = NF_DEF_VAR (nid,"ps",NF_FLOAT,3,dims3,nvarid)
#endif
  ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 15, "Pression au sol")
  !
  ierr = NF_ENDDEF(nid) ! sortir du mode de definition
  ierr = NF_CLOSE(nid) ! fermer le fichier
  
  PRINT*,'iim,jjm,llm,idayref',iim,jjm,llm,idayref
  PRINT*,'rad,omeg,g,cpp,kappa', rad,omeg,g,cpp,kappa
  
  RETURN
END SUBROUTINE dynrestart0


SUBROUTINE dynrestart1(fichnom,time, vcov,ucov,teta,q,masse,ps)
  
  USE SPECIES_NAME
  
  IMPLICIT NONE
  !=================================================================
  !  Ecriture du fichier de redemarrage sous FORMAT NetCDF
  !=================================================================
  include "dimensions.h"
  include "paramet.h"
  include "description.h"
  include "netcdf.inc"
  include "comvert.h"
  include "comgeom.h"
  
  REAL :: vcov(ip1jm,llm),ucov(ip1jmp1,llm)
  REAL :: teta(ip1jmp1,llm)                   
  REAL :: ps(ip1jmp1),masse(ip1jmp1,llm)                   
  REAL :: q(ip1jmp1,llm,iqmax)
  CHARACTER(len=*) :: fichnom
  
  REAL :: time
  INTEGER :: nid, nvarid
  INTEGER :: ierr
  INTEGER :: iq
  CHARACTER(len=3) :: str3
  CHARACTER(len=20) :: modname
  !
  INTEGER,SAVE :: nb
  DATA nb / 0 /
  
  nb = 0
  modname = 'dynredem1'
  ierr = NF_OPEN(fichnom, NF_WRITE, nid)
  IF (ierr /= NF_NOERR) THEN
      PRINT*, "Pb in opening: "//fichnom
      CALL abort
  ENDIF
  
  !  Writing/extension of the time variable
  
  nb = nb + 1
  ierr = NF_INQ_VARID(nid, "temps", nvarid)
  IF (ierr /= NF_NOERR) THEN
      PRINT *, NF_STRERROR(ierr)
      WRITE(*,*) 'Variable temps is not defined'
      STOP
  ENDIF
#ifdef NC_DOUBLE
  ierr = NF_PUT_VAR1_DOUBLE (nid,nvarid,nb,time)
#else
  ierr = NF_PUT_VAR1_REAL (nid,nvarid,nb,time)
#endif
  PRINT*, "Saving for ", nb, time
  
  !  Ecriture des champs
  !
  ierr = NF_INQ_VARID(nid, "ucov", nvarid)
  IF (ierr /= NF_NOERR) THEN
      PRINT*, "Variable ucov is not defined"
      CALL abort
  ENDIF
#ifdef NC_DOUBLE
  ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,ucov)
#else
  ierr = NF_PUT_VAR_REAL (nid,nvarid,ucov)
#endif
  
  ierr = NF_INQ_VARID(nid, "vcov", nvarid)
  IF (ierr /= NF_NOERR) THEN
      PRINT*, "Variable vcov is not defined"
      CALL abort
  ENDIF
#ifdef NC_DOUBLE
  ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,vcov)
#else
  ierr = NF_PUT_VAR_REAL (nid,nvarid,vcov)
#endif
  
  ierr = NF_INQ_VARID(nid, "teta", nvarid)
  IF (ierr /= NF_NOERR) THEN
      PRINT*, "Variable teta is not defined"
      CALL abort
  ENDIF
#ifdef NC_DOUBLE
  ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,teta)
#else
  ierr = NF_PUT_VAR_REAL (nid,nvarid,teta)
#endif
  
  IF(iqmax>=1) THEN
      DO iq=1,iqmax
        IF (iq>99) THEN
            PRINT*, "Too many tracers"
            CALL abort
        ELSE
            str3(1:1)='q'
            WRITE(str3(2:3),'(i2.2)') restart_ids(iq)
            ierr = NF_INQ_VARID(nid, str3, nvarid)
            IF (ierr /= NF_NOERR) THEN
                PRINT*, "Variable "//str3//" is not defined"
                CALL abort
            ENDIF
#ifdef NC_DOUBLE
            ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,q(1,1,iq))
#else
            ierr = NF_PUT_VAR_REAL (nid,nvarid,q(1,1,iq))
#endif
        ENDIF
      ENDDO
  ENDIF
  !
  ierr = NF_INQ_VARID(nid, "masse", nvarid)
  IF (ierr /= NF_NOERR) THEN
      PRINT*, "Variable masse is not defined"
      CALL abort
  ENDIF
#ifdef NC_DOUBLE
  ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,masse)
#else
  ierr = NF_PUT_VAR_REAL (nid,nvarid,masse)
#endif
  !
  ierr = NF_INQ_VARID(nid, "ps", nvarid)
  IF (ierr /= NF_NOERR) THEN
      PRINT*, "Variable ps is not defined"
      CALL abort
  ENDIF
#ifdef NC_DOUBLE
  ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,ps)
#else
  ierr = NF_PUT_VAR_REAL (nid,nvarid,ps)
#endif
  
  ierr = NF_CLOSE(nid)
  !
  RETURN
END SUBROUTINE dynrestart1

