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

SUBROUTINE nflxtr(pdtime,pmfu,pmfd,pde_u,pen_d, paprs,x,dx) 
  USE dimphy
  USE Write_Field_phy
  IMPLICIT NONE 
  !=====================================================================
  ! Objet : Melange convectif de traceurs a partir des flux de masse 
  ! Date : 13/12/1996 -- 13/01/97
  ! Auteur: O. Boucher (LOA) sur inspiration de Z. X. Li (LMD),
  !         Brinkop et Sausen (1996) et Boucher et al. (1996).
  ! ATTENTION : meme si cette routine se veut la plus generale possible, 
  !             elle a herite de certaines notations et conventions du 
  !             schema de Tiedtke (1993). 
  ! --En particulier, les couches sont numerotees de haut en bas !!!
  !   Ceci est valable pour les flux
  !   mais pas pour les entrees x, pplay, paprs !!!!
  ! --pmfu est positif, pmfd est negatif 
  ! --Tous les flux d'entrainements et de detrainements sont positifs 
  !   contrairement au schema de Tiedtke d'ou les changements de signe!!!! 
  !=====================================================================
  !
  INCLUDE "YOMCST.h"
  INCLUDE "YOECUMF.h" 
  !
  REAL :: pdtime
  !--les flux sont definis au 1/2 niveaux
  !--pmfu(klev+1) et pmfd(klev+1) sont implicitement nuls
  REAL :: pmfu(klon,klev)  ! flux de masse dans le panache montant 
  REAL :: pmfd(klon,klev)  ! flux de masse dans le panache descendant
  REAL :: pde_u(klon,klev) ! flux detraine dans le panache montant
  REAL :: pen_d(klon,klev) ! flux entraine dans le panache descendant
  
  REAL :: paprs(klon,klev+1)  ! pression aux 1/2 couches (bas en haut)
  REAL :: x(klon,klev)        ! q de traceur (bas en haut) 
  REAL :: dx(klon,klev)     ! tendance de traceur  (bas en haut)
  !
  !--flux convectifs mais en variables locales
  REAL :: zmfu(klon,klev+1) 
  REAL :: zmfd(klon,klev+1) 
  REAL :: zen_u(klon,klev) 
  REAL :: zde_u(klon,klev)
  REAL :: zen_d(klon,klev) 
  REAL :: zde_d(klon,klev)
  REAL :: zmfe
  !
  !--variables locales      
  !--les flux de x sont definis aux 1/2 niveaux 
  !--xu et xd sont definis aux niveaux complets
  REAL :: xu(klon,klev)        ! q de traceurs dans le panache montant
  REAL :: xd(klon,klev)        ! q de traceurs dans le panache descendant
  REAL :: zmfux(klon,klev+1)   ! flux de x dans le panache montant
  REAL :: zmfdx(klon,klev+1)   ! flux de x dans le panache descendant
  REAL :: zmfex(klon,klev+1)   ! flux de x dans l'environnement 
  INTEGER :: i, k 
  REAL,PARAMETER :: zmfmin=1.E-10

  
  ! =========================================
  !
  !
  !   Extension des flux UP et DN sur klev+1 niveaux
  ! =========================================
  DO k=1,klev
    DO i=1,klon
      zmfu(i,k)=pmfu(i,k)
      zmfd(i,k)=pmfd(i,k)
    ENDDO
  ENDDO
  DO i=1,klon
    zmfu(i,klev+1)=0.
    zmfd(i,klev+1)=0.
  ENDDO
  
  !--modif pour diagnostiquer les detrainements
  ! =========================================
  !   on privilegie l'ajustement de l'entrainement dans l'ascendance.
  
  DO k=1, klev
    DO i=1, klon
      zen_d(i,k)=pen_d(i,k)
      zde_u(i,k)=pde_u(i,k)
      zde_d(i,k) =-zmfd(i,k+1)+zmfd(i,k)+zen_d(i,k)
      zen_u(i,k) = zmfu(i,k+1)-zmfu(i,k)+zde_u(i,k)
    ENDDO
  ENDDO
  !
  !--calcul des flux dans le panache montant
  ! =========================================
  !
  ! Dans la premiere couche, on prend q comme valeur de qu
  !
  DO i=1, klon
    zmfux(i,1)=0.0 
  ENDDO
  !
  ! Autres couches
  DO k=1,klev
    DO i=1, klon
      IF ((zmfu(i,k+1)+zde_u(i,k))<zmfmin) THEN
          xu(i,k)=x(i,k)
      ELSE
          xu(i,k)=(zmfux(i,k)+zen_u(i,k)*x(i,k)) &
             /(zmfu(i,k+1)+zde_u(i,k))
      ENDIF
      zmfux(i,k+1)=zmfu(i,k+1)*xu(i,k)
    ENDDO
  ENDDO
  !
  !--calcul des flux dans le panache descendant
  ! =========================================
  !   
  DO i=1, klon
    zmfdx(i,klev+1)=0.0 
  ENDDO
  !
  DO k=klev,1,-1
    DO i=1, klon
      IF ((zde_d(i,k)-zmfd(i,k))<zmfmin) THEN
          xd(i,k)=x(i,k)
      ELSE
          xd(i,k)=(zmfdx(i,k+1)-zen_d(i,k)*x(i,k)) / &
             (zmfd(i,k)-zde_d(i,k))
      ENDIF
      zmfdx(i,k)=zmfd(i,k)*xd(i,k)
    ENDDO
  ENDDO
  !
  !--introduction du flux de retour dans l'environnement
  ! =========================================
  !
  do k=2, klev
    do i=1, klon
      zmfe=-zmfu(i,k)-zmfd(i,k)
      if (zmfe<=0.) then
          zmfex(i,k)= zmfe*x(i,k)
      else
          zmfex(i,k)= zmfe*x(i,k-1)
      endif
    enddo
  enddo
  
  DO i=1, klon
    zmfex(i,1)=0.
    zmfex(i,klev+1)=0.
  ENDDO
  !
  !--calcul final des tendances
  !
  DO k=1, klev
    DO i=1, klon
      dx(i,k)=RG/(paprs(i,k)-paprs(i,k+1))*pdtime* &
         ( zmfux(i,k) - zmfux(i,k+1) + &
         zmfdx(i,k) - zmfdx(i,k+1) + &
         zmfex(i,k) - zmfex(i,k+1) )
    ENDDO
  ENDDO
  !
  RETURN 
END SUBROUTINE nflxtr
!
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!
SUBROUTINE nflxtr_tl(pdtime,pmfu,pmfd,pde_u,pen_d, paprs,x,dx,x_tl,dx_tl) 
  USE dimphy
  IMPLICIT NONE 
  !=====================================================================
  ! Objet : Melange convectif de traceurs a partir des flux de masse 
  ! Date : 13/12/1996 -- 13/01/97
  ! Auteur: O. Boucher (LOA) sur inspiration de Z. X. Li (LMD),
  !         Brinkop et Sausen (1996) et Boucher et al. (1996).
  ! ATTENTION : meme si cette routine se veut la plus generale possible, 
  !             elle a herite de certaines notations et conventions du 
  !             schema de Tiedtke (1993). 
  ! --En particulier, les couches sont numerotees de haut en bas !!!
  !   Ceci est valable pour les flux
  !   mais pas pour les entrees x, pplay, paprs !!!!
  ! --pmfu est positif, pmfd est negatif 
  ! --Tous les flux d'entrainements et de detrainements sont positifs 
  !   contrairement au schema de Tiedtke d'ou les changements de signe!!!! 
  !=====================================================================
  !
  INCLUDE "YOMCST.h"
  INCLUDE "YOECUMF.h" 
  !
  REAL :: pdtime
  !--les flux sont definis au 1/2 niveaux
  !--pmfu(klev+1) et pmfd(klev+1) sont implicitement nuls
  REAL :: pmfu(klon,klev)  ! flux de masse dans le panache montant 
  REAL :: pmfd(klon,klev)  ! flux de masse dans le panache descendant
  REAL :: pde_u(klon,klev) ! flux detraine dans le panache montant
  REAL :: pen_d(klon,klev) ! flux entraine dans le panache descendant
  
  REAL :: paprs(klon,klev+1)  ! pression aux 1/2 couches (bas en haut)
  REAL :: x(klon,klev)        ! q de traceur (bas en haut) 
  REAL :: dx(klon,klev)     ! tendance de traceur  (bas en haut)
  REAL :: x_tl(klon,klev)        ! q de traceur (bas en haut) 
  REAL :: dx_tl(klon,klev)     ! tendance de traceur  (bas en haut)
  !
  !--flux convectifs mais en variables locales
  REAL :: zmfu(klon,klev+1) 
  REAL :: zmfd(klon,klev+1) 
  REAL :: zen_u(klon,klev) 
  REAL :: zde_u(klon,klev)
  REAL :: zen_d(klon,klev) 
  REAL :: zde_d(klon,klev)
  REAL :: zmfe
  !
  !--variables locales      
  !--les flux de x sont definis aux 1/2 niveaux 
  !--xu et xd sont definis aux niveaux complets
  REAL :: xu(klon,klev)        ! q de traceurs dans le panache montant
  REAL :: xu_tl(klon,klev)        ! q de traceurs dans le panache montant
  REAL :: xd(klon,klev)        ! q de traceurs dans le panache descendant
  REAL :: xd_tl(klon,klev)        ! q de traceurs dans le panache descendant
  REAL :: zmfux(klon,klev+1)   ! flux de x dans le panache montant
  REAL :: zmfux_tl(klon,klev+1)   ! flux de x dans le panache montant
  REAL :: zmfdx(klon,klev+1)   ! flux de x dans le panache descendant
  REAL :: zmfdx_tl(klon,klev+1)   ! flux de x dans le panache descendant
  REAL :: zmfex(klon,klev+1)   ! flux de x dans l'environnement 
  REAL :: zmfex_tl(klon,klev+1)   ! flux de x dans l'environnement 
  INTEGER :: i, k 
  REAL,PARAMETER ::  zmfmin=1.E-10
  
  ! =========================================
  !
  !
  !   Extension des flux UP et DN sur klev+1 niveaux
  ! =========================================
  DO k=1,klev
    DO i=1,klon
      zmfu(i,k)=pmfu(i,k)
      zmfd(i,k)=pmfd(i,k)
    ENDDO
  ENDDO
  DO i=1,klon
    zmfu(i,klev+1)=0.
    zmfd(i,klev+1)=0.
  ENDDO
  
  !--modif pour diagnostiquer les detrainements
  ! =========================================
  !   on privilegie l'ajustement de l'entrainement dans l'ascendance.
  
  DO k=1, klev
    DO i=1, klon
      zen_d(i,k)=pen_d(i,k)
      zde_u(i,k)=pde_u(i,k)
      zde_d(i,k) =-zmfd(i,k+1)+zmfd(i,k)+zen_d(i,k)
      zen_u(i,k) = zmfu(i,k+1)-zmfu(i,k)+zde_u(i,k)
    ENDDO
  ENDDO
  !
  !--calcul des flux dans le panache montant
  ! =========================================
  !
  ! Dans la premiere couche, on prend q comme valeur de qu
  !
  DO i=1, klon
    zmfux_tl(i,1)=0.0 
    zmfux(i,1)=0.0 
  ENDDO
  !
  ! Autres couches
  DO k=1,klev
    DO i=1, klon
      IF ((zmfu(i,k+1)+zde_u(i,k))<zmfmin) THEN
          xu_tl(i,k)=x_tl(i,k)
          xu(i,k)=x(i,k)
      ELSE
          xu_tl(i,k)=(zmfux_tl(i,k)+zen_u(i,k)*x_tl(i,k)) &
             /(zmfu(i,k+1)+zde_u(i,k))
          xu(i,k)=(zmfux(i,k)+zen_u(i,k)*x(i,k)) &
             /(zmfu(i,k+1)+zde_u(i,k))
      ENDIF
      zmfux_tl(i,k+1)=zmfu(i,k+1)*xu_tl(i,k)
      zmfux(i,k+1)=zmfu(i,k+1)*xu(i,k)
    ENDDO
  ENDDO
  
  !
  !--calcul des flux dans le panache descendant
  ! =========================================
  !   
  DO i=1, klon
    zmfdx_tl(i,klev+1)=0.0 
    zmfdx(i,klev+1)=0.0 
  ENDDO
  !
  DO k=klev,1,-1
    DO i=1, klon
      IF ((zde_d(i,k)-zmfd(i,k))<zmfmin) THEN
          xd_tl(i,k)=x_tl(i,k)
          xd(i,k)=x(i,k)
      ELSE
          xd_tl(i,k)=(zmfdx_tl(i,k+1)-zen_d(i,k)*x_tl(i,k)) / &
             (zmfd(i,k)-zde_d(i,k))
          xd(i,k)=(zmfdx(i,k+1)-zen_d(i,k)*x(i,k)) / &
             (zmfd(i,k)-zde_d(i,k))
      ENDIF
      zmfdx_tl(i,k)=zmfd(i,k)*xd_tl(i,k)
      zmfdx(i,k)=zmfd(i,k)*xd(i,k)
    ENDDO
  ENDDO
  !
  !--introduction du flux de retour dans l'environnement
  ! =========================================
  !
  DO k=2, klev
    DO i=1, klon
      zmfe=-zmfu(i,k)-zmfd(i,k)
      IF (zmfe<=0.) THEN
          zmfex_tl(i,k)= zmfe*x_tl(i,k)
          zmfex(i,k)= zmfe*x(i,k)
      ELSE
          zmfex_tl(i,k)= zmfe*x_tl(i,k-1)
          zmfex(i,k)= zmfe*x(i,k-1)
      ENDIF
    ENDDO
  ENDDO
  
  DO i=1, klon
    zmfex_tl(i,1)=0.
    zmfex(i,1)=0.
    zmfex_tl(i,klev+1)=0.
    zmfex(i,klev+1)=0.
  ENDDO
  !
  !--calcul final des tendances
  !
  DO k=1, klev
    DO i=1, klon
      dx_tl(i,k)=RG/(paprs(i,k)-paprs(i,k+1))*pdtime* &
         ( zmfux_tl(i,k) - zmfux_tl(i,k+1) + &
         zmfdx_tl(i,k) - zmfdx_tl(i,k+1) + &
         zmfex_tl(i,k) - zmfex_tl(i,k+1) )
      dx(i,k)=RG/(paprs(i,k)-paprs(i,k+1))*pdtime* &
         ( zmfux(i,k) - zmfux(i,k+1) + &
         zmfdx(i,k) - zmfdx(i,k+1) + &
         zmfex(i,k) - zmfex(i,k+1) )
    ENDDO
  ENDDO
  !
  RETURN 
END SUBROUTINE nflxtr_tl
!
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!
SUBROUTINE nflxtr_ad(pdtime,pmfu,pmfd,pde_u,pen_d,paprs,x,dx,x_ad,dx_ad) 
  USE dimphy
  IMPLICIT NONE 
  !=====================================================================
  ! Objet : Melange convectif de traceurs a partir des flux de masse 
  ! Date : 13/12/1996 -- 13/01/97
  ! Auteur: O. Boucher (LOA) sur inspiration de Z. X. Li (LMD),
  !         Brinkop et Sausen (1996) et Boucher et al. (1996).
  ! ATTENTION : meme si cette routine se veut la plus generale possible, 
  !             elle a herite de certaines notations et conventions du 
  !             schema de Tiedtke (1993). 
  ! --En particulier, les couches sont numerotees de haut en bas !!!
  !   Ceci est valable pour les flux
  !   mais pas pour les entrees x, pplay, paprs !!!!
  ! --pmfu est positif, pmfd est negatif 
  ! --Tous les flux d'entrainements et de detrainements sont positifs 
  !   contrairement au schema de Tiedtke d'ou les changements de signe!!!! 
  !=====================================================================
  !
  INCLUDE "YOMCST.h"
  INCLUDE "YOECUMF.h" 
  !
  REAL :: pdtime
  !--les flux sont definis au 1/2 niveaux
  !--pmfu(klev+1) et pmfd(klev+1) sont implicitement nuls
  REAL :: pmfu(klon,klev)  ! flux de masse dans le panache montant 
  REAL :: pmfd(klon,klev)  ! flux de masse dans le panache descendant
  REAL :: pde_u(klon,klev) ! flux detraine dans le panache montant
  REAL :: pen_d(klon,klev) ! flux entraine dans le panache descendant
  
  REAL :: paprs(klon,klev+1)  ! pression aux 1/2 couches (bas en haut)
  REAL :: x(klon,klev)        ! q de traceur (bas en haut) 
  REAL :: dx(klon,klev)     ! tendance de traceur  (bas en haut)
  REAL :: x_ad(klon,klev)        ! q de traceur (bas en haut) 
  REAL :: dx_ad(klon,klev)     ! tendance de traceur  (bas en haut)
  !
  !--flux convectifs mais en variables locales
  REAL :: zmfu(klon,klev+1) 
  REAL :: zmfd(klon,klev+1) 
  REAL :: zen_u(klon,klev) 
  REAL :: zde_u(klon,klev)
  REAL :: zen_d(klon,klev) 
  REAL :: zde_d(klon,klev)
  REAL :: zmfe
  !
  !--variables locales      
  !--les flux de x sont definis aux 1/2 niveaux 
  !--xu et xd sont definis aux niveaux complets
  REAL :: xu(klon,klev)        ! q de traceurs dans le panache montant
  REAL :: xu_ad(klon,klev)        ! q de traceurs dans le panache montant
  REAL :: xd(klon,klev)        ! q de traceurs dans le panache descendant
  REAL :: xd_ad(klon,klev)        ! q de traceurs dans le panache descendant
  REAL :: zmfux(klon,klev+1)   ! flux de x dans le panache montant
  REAL :: zmfux_ad(klon,klev+1)   ! flux de x dans le panache montant
  REAL :: zmfdx(klon,klev+1)   ! flux de x dans le panache descendant
  REAL :: zmfdx_ad(klon,klev+1)   ! flux de x dans le panache descendant
  REAL :: zmfex(klon,klev+1)   ! flux de x dans l'environnement 
  REAL :: zmfex_ad(klon,klev+1)   ! flux de x dans l'environnement 
  INTEGER :: i, k 
  REAL,PARAMETER :: zmfmin=1.E-10
  
  ! =========================================
  !
  !
  zmfux_ad(:,:) = 0.
  zmfdx_ad(:,:) = 0.
  zmfex_ad(:,:) = 0.
  xu_ad(:,:)=0.
  
  !   Extension des flux UP et DN sur klev+1 niveaux
  ! =========================================
  DO k=1,klev
    DO i=1,klon
      zmfu(i,k)=pmfu(i,k)
      zmfd(i,k)=pmfd(i,k)
    ENDDO
  ENDDO
  DO i=1,klon
    zmfu(i,klev+1)=0.
    zmfd(i,klev+1)=0.
  ENDDO
  
  !--modif pour diagnostiquer les detrainements
  ! =========================================
  !   on privilegie l'ajustement de l'entrainement dans l'ascendance.
  
  do k=1, klev
    do i=1, klon
      zen_d(i,k)=pen_d(i,k)
      zde_u(i,k)=pde_u(i,k)
      zde_d(i,k) =-zmfd(i,k+1)+zmfd(i,k)+zen_d(i,k)
      zen_u(i,k) = zmfu(i,k+1)-zmfu(i,k)+zde_u(i,k)
    enddo 
  enddo 
  !
  !--calcul des flux dans le panache montant
  ! =========================================
  !
  ! Dans la premiere couche, on prend q comme valeur de qu
  !
  DO i=1, klon
    zmfux(i,1)=0.0 
  ENDDO
  !
  ! Autres couches
  DO k=1,klev
    DO i=1, klon
      IF ((zmfu(i,k+1)+zde_u(i,k))<zmfmin) THEN
          xu(i,k)=x(i,k)
      ELSE
          xu(i,k)=(zmfux(i,k)+zen_u(i,k)*x(i,k)) &
             /(zmfu(i,k+1)+zde_u(i,k))
      ENDIF
      zmfux(i,k+1)=zmfu(i,k+1)*xu(i,k)
    ENDDO
  ENDDO
  !
  !--calcul des flux dans le panache descendant
  ! =========================================
  !   
  DO i=1, klon
    zmfdx(i,klev+1)=0.0 
  ENDDO
  !
  DO k=klev,1,-1
    DO i=1, klon
      IF ((zde_d(i,k)-zmfd(i,k))<zmfmin) THEN
          xd(i,k)=x(i,k)
      ELSE
          xd(i,k)=(zmfdx(i,k+1)-zen_d(i,k)*x(i,k)) / &
             (zmfd(i,k)-zde_d(i,k))
      ENDIF
      zmfdx(i,k)=zmfd(i,k)*xd(i,k)
    ENDDO
  ENDDO
  !
  !--introduction du flux de retour dans l'environnement
  ! =========================================
  !
  DO k=2, klev
    DO i=1, klon
      zmfe=-zmfu(i,k)-zmfd(i,k)
      IF (zmfe<=0.) THEN
          zmfex(i,k)= zmfe*x(i,k)
      ELSE
          zmfex(i,k)= zmfe*x(i,k-1)
      ENDIF
    ENDDO
  ENDDO
  
  DO i=1, klon
    zmfex(i,1)=0.
    zmfex(i,klev+1)=0.
  ENDDO
  !
  !--calcul final des tendances
  !
  DO k=1, klev
    DO i=1, klon
      dx(i,k)=RG/(paprs(i,k)-paprs(i,k+1))*pdtime* &
         ( zmfux(i,k) - zmfux(i,k+1) + &
         zmfdx(i,k) - zmfdx(i,k+1) + &
         zmfex(i,k) - zmfex(i,k+1) )
      !Start of AD
      zmfux_ad(i,k) = zmfux_ad(i,k) + dx_ad(i,k) &
         *RG/(paprs(i,k)-paprs(i,k+1))*pdtime
      zmfux_ad(i,k+1) = zmfux_ad(i,k+1) - dx_ad(i,k) &
         *RG/(paprs(i,k)-paprs(i,k+1))*pdtime
      zmfdx_ad(i,k) = zmfdx_ad(i,k) + dx_ad(i,k) &
         *RG/(paprs(i,k)-paprs(i,k+1))*pdtime
      zmfdx_ad(i,k+1) = zmfdx_ad(i,k+1) - dx_ad(i,k) &
         *RG/(paprs(i,k)-paprs(i,k+1))*pdtime
      zmfex_ad(i,k) = zmfex_ad(i,k) + dx_ad(i,k) &
         *RG/(paprs(i,k)-paprs(i,k+1))*pdtime
      zmfex_ad(i,k+1) = zmfex_ad(i,k+1) - dx_ad(i,k) &
         *RG/(paprs(i,k)-paprs(i,k+1))*pdtime
      dx_ad(i,k) = 0.
    ENDDO
  ENDDO

  DO i=1, klon
    zmfex_ad(i,1)=0.
    zmfex_ad(i,klev+1)=0.
  ENDDO
  
  DO k=2, klev
    DO i=1, klon
      zmfe=-zmfu(i,k)-zmfd(i,k)
      IF (zmfe<=0.) THEN
          x_ad(i,k) = x_ad(i,k) + zmfe*zmfex_ad(i,k)
      ELSE
          x_ad(i,k-1) = x_ad(i,k-1) + zmfe*zmfex_ad(i,k)
      ENDIF
      zmfex_ad(i,k) = 0.
    ENDDO
  ENDDO
  
  xd_ad(:,:) = 0.
  DO k=1,klev
    DO i=1, klon
      xd_ad(i,k) = xd_ad(i,k) + zmfd(i,k)*zmfdx_ad(i,k)
      zmfdx_ad(i,k) = 0.
      IF ((zde_d(i,k)-zmfd(i,k))<zmfmin) THEN
          x_ad(i,k) = x_ad(i,k) + xd_ad(i,k)
      ELSE
          zmfdx_ad(i,k+1) = zmfdx_ad(i,k+1) + xd_ad(i,k)/ &
             (zmfd(i,k)-zde_d(i,k))
          x_ad(i,k) = x_ad(i,k) - zen_d(i,k) * xd_ad(i,k)/ &
             (zmfd(i,k)-zde_d(i,k))
      ENDIF
      xd_ad(i,k) = 0.
    ENDDO
  ENDDO
  DO i=1, klon
    zmfdx_ad(i,klev+1)=0.0 
  ENDDO
  
  DO k=klev,1,-1
    DO i=1, klon
      xu_ad(i,k) = xu_ad(i,k) + zmfu(i,k+1)*zmfux_ad(i,k+1)
      zmfux_ad(i,k+1) = 0.
      IF ((zmfu(i,k+1)+zde_u(i,k))<zmfmin) THEN
          x_ad(i,k)=x_ad(i,k) + xu_ad(i,k)
      ELSE
          zmfux_ad(i,k) = zmfux_ad(i,k) + xu_ad(i,k) &
             /(zmfu(i,k+1)+zde_u(i,k))
          x_ad(i,k) = x_ad(i,k) + zen_u(i,k) * xu_ad(i,k) &
             /(zmfu(i,k+1)+zde_u(i,k))
      ENDIF
      xu_ad(i,k) = 0.
    ENDDO
  ENDDO
  DO i=1, klon
    zmfux_ad(i,1)=0.0 
  ENDDO
  !
  RETURN 
END SUBROUTINE nflxtr_ad
