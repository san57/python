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

SUBROUTINE cltrac(dtime,coef,t,tr,flux,paprs,pplay,delp, &
   d_tr)
  USE dimphy
  IMPLICIT NONE
  !======================================================================
  ! Auteur(s): O. Boucher (LOA/LMD) date: 19961127
  !            inspire de clvent
  ! Objet: diffusion verticale de traceurs avec flux fixe a la surface
  !        ou/et flux du TYPE c-drag
  !======================================================================
  ! Arguments:
  ! dtime----input-R- intervalle du temps (en second)
  ! coef-----input-R- le coefficient d'echange (m**2/s) l>1
  ! t--------input-R- temperature (K)
  ! tr-------input-R- la q. de traceurs
  ! flux-----input-R- le flux de traceurs a la surface
  ! paprs----input-R- pression a inter-couche (Pa)
  ! pplay----input-R- pression au milieu de couche (Pa)
  ! delp-----input-R- epaisseur de couche (Pa)
  ! cdrag----input-R- cdrag pour le flux de surface (non active)
  ! tr0------input-R- traceurs a la surface ou dans l'ocean (non active)
  ! d_tr-----output-R- le changement de tr
  ! flux_tr--output-R- flux de tr
  !======================================================================
  REAL :: dtime
  REAL :: coef(klon,klev)
  REAL :: t(klon,klev), tr(klon,klev)
  REAL :: paprs(klon,klev+1), pplay(klon,klev), delp(klon,klev)
  REAL :: d_tr(klon,klev)
  REAL :: flux(klon), cdrag(klon), tr0(klon)
  !======================================================================
  INCLUDE "YOMCST.h"
  !======================================================================
  INTEGER :: i, k
  REAL :: zx_ctr(klon,2:klev)
  REAL :: zx_dtr(klon,2:klev)
  REAL :: zx_buf(klon)
  REAL :: zx_coef(klon,klev)
  REAL :: local_tr(klon,klev)
  REAL :: zx_alf1(klon), zx_alf2(klon), zx_flux(klon)
  !======================================================================
  local_tr(:,:) = tr(:,:)
  !======================================================================
  
  DO i = 1, klon
    zx_alf1(i) = (paprs(i,1)-pplay(i,2))/(pplay(i,1)-pplay(i,2))
    zx_alf2(i) = 1.0 - zx_alf1(i)
    zx_flux(i) =  -flux(i)*dtime*RG
    !--pour le moment le flux est prescrit
    !--cdrag et zx_coef(1) vaut 0
    cdrag(i) = 0.0 
    tr0(i) = 0.0
    zx_coef(i,1) = cdrag(i)*dtime*RG 
  ENDDO
  !======================================================================
  DO k = 2, klev
    DO i = 1, klon
      zx_coef(i,k) = coef(i,k)*RG/(pplay(i,k-1)-pplay(i,k)) &
         *(paprs(i,k)*2/(t(i,k)+t(i,k-1))/RD)**2
      zx_coef(i,k) = zx_coef(i,k)*dtime*RG
    ENDDO
  ENDDO
  !======================================================================
  DO i = 1, klon
    zx_buf(i) = delp(i,1) + zx_coef(i,1)*zx_alf1(i) + zx_coef(i,2)
    zx_ctr(i,2) = (local_tr(i,1)*delp(i,1)+ &
       zx_coef(i,1)*tr0(i)-zx_flux(i))/zx_buf(i)
    zx_dtr(i,2) = (zx_coef(i,2)-zx_alf2(i)*zx_coef(i,1)) /  &
       zx_buf(i)
  ENDDO
  !
  DO k = 3, klev
    DO i = 1, klon 
      zx_buf(i) = delp(i,k-1) + zx_coef(i,k) &
         + zx_coef(i,k-1)*(1.-zx_dtr(i,k-1))
      zx_ctr(i,k) = (local_tr(i,k-1)*delp(i,k-1) &
         +zx_coef(i,k-1)*zx_ctr(i,k-1) )/zx_buf(i)
      zx_dtr(i,k) = zx_coef(i,k)/zx_buf(i)
    ENDDO
  ENDDO
  DO i = 1, klon
    local_tr(i,klev) = ( local_tr(i,klev)*delp(i,klev) &
       +zx_coef(i,klev)*zx_ctr(i,klev) ) &
       / ( delp(i,klev) + zx_coef(i,klev) &
       -zx_coef(i,klev)*zx_dtr(i,klev) )
  ENDDO
  DO k = klev-1, 1, -1
    DO i = 1, klon
      local_tr(i,k) = zx_ctr(i,k+1) + zx_dtr(i,k+1)*local_tr(i,k+1)
    ENDDO
  ENDDO
  DO k = 1, klev
    DO i = 1, klon
      d_tr(i,k) = local_tr(i,k) - tr(i,k)
    ENDDO
  ENDDO
  !
  RETURN
END SUBROUTINE cltrac
!
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!
SUBROUTINE cltrac_tl(dtime,coef,t,tr,flux,paprs,pplay,delp, &
   d_tr,flux_tl,tr_tl,d_tr_tl)
  USE dimphy
  IMPLICIT NONE
  !======================================================================
  ! Auteur(s): O. Boucher (LOA/LMD) date: 19961127
  !            inspire de clvent
  ! Objet: diffusion verticale de traceurs avec flux fixe a la surface
  !        ou/et flux du TYPE c-drag
  !======================================================================
  ! Arguments:
  ! dtime----input-R- intervalle du temps (en second)
  ! coef-----input-R- le coefficient d'echange (m**2/s) l>1
  ! t--------input-R- temperature (K)
  ! tr-------input-R- la q. de traceurs
  ! flux-----input-R- le flux de traceurs a la surface
  ! paprs----input-R- pression a inter-couche (Pa)
  ! pplay----input-R- pression au milieu de couche (Pa)
  ! delp-----input-R- epaisseur de couche (Pa)
  ! cdrag----input-R- cdrag pour le flux de surface (non active)
  ! tr0------input-R- traceurs a la surface ou dans l'ocean (non active)
  ! d_tr-----output-R- le changement de tr
  ! flux_tr--output-R- flux de tr
  !======================================================================
  REAL :: dtime
  REAL :: coef(klon,klev)
  REAL :: t(klon,klev), tr(klon,klev)
  REAL :: tr_tl(klon,klev)
  REAL :: paprs(klon,klev+1), pplay(klon,klev), delp(klon,klev)
  REAL :: d_tr(klon,klev)
  REAL :: d_tr_tl(klon,klev)
  REAL :: flux(klon), cdrag(klon), tr0(klon)
  REAL :: flux_tl(klon), tr0_tl(klon)
  !======================================================================
  INCLUDE "YOMCST.h"
  !======================================================================
  INTEGER :: i, k
  REAL :: zx_ctr(klon,2:klev)
  REAL :: zx_ctr_tl(klon,2:klev)
  REAL :: zx_dtr(klon,2:klev)
  REAL :: zx_buf(klon)
  REAL :: zx_coef(klon,klev)
  REAL :: local_tr(klon,klev)
  REAL :: local_tr_tl(klon,klev)
  REAL :: zx_alf1(klon), zx_alf2(klon), zx_flux(klon)
  REAL :: zx_flux_tl(klon)
  !======================================================================
  local_tr_tl(:,:) = tr_tl(:,:)
  local_tr(:,:)    = tr(:,:)
  !======================================================================
  DO i = 1, klon
    zx_alf1(i) = (paprs(i,1)-pplay(i,2))/(pplay(i,1)-pplay(i,2))
    zx_alf2(i) = 1.0 - zx_alf1(i)
    zx_flux_tl(i) =  -flux_tl(i)*dtime*RG
    zx_flux(i) =  -flux(i)*dtime*RG
    !--pour le moment le flux est prescrit
    !--cdrag et zx_coef(1) vaut 0
    cdrag(i) = 0.0 
    tr0_tl(i) = 0.0
    tr0(i) = 0.0
    zx_coef(i,1) = cdrag(i)*dtime*RG 
  ENDDO
  !======================================================================
  DO k = 2, klev
    DO i = 1, klon
      zx_coef(i,k) = coef(i,k)*RG/(pplay(i,k-1)-pplay(i,k)) &
         *(paprs(i,k)*2/(t(i,k)+t(i,k-1))/RD)**2
      zx_coef(i,k) = zx_coef(i,k)*dtime*RG
    ENDDO
  ENDDO
  !======================================================================
  DO i = 1, klon
    zx_buf(i) = delp(i,1) + zx_coef(i,1)*zx_alf1(i) + zx_coef(i,2)
    zx_ctr_tl(i,2) = (local_tr_tl(i,1)*delp(i,1)+ &
       zx_coef(i,1)*tr0_tl(i)-zx_flux_tl(i))/zx_buf(i)
    zx_ctr(i,2) = (local_tr(i,1)*delp(i,1)+ &
       zx_coef(i,1)*tr0(i)-zx_flux(i))/zx_buf(i)
    zx_dtr(i,2) = (zx_coef(i,2)-zx_alf2(i)*zx_coef(i,1)) /  &
       zx_buf(i)
  ENDDO
  !
  DO k = 3, klev
    DO i = 1, klon
      zx_buf(i) = delp(i,k-1) + zx_coef(i,k) &
         + zx_coef(i,k-1)*(1.-zx_dtr(i,k-1))
      zx_ctr_tl(i,k) = (local_tr_tl(i,k-1)*delp(i,k-1) &
         +zx_coef(i,k-1)*zx_ctr_tl(i,k-1) )/zx_buf(i)
      zx_ctr(i,k) = (local_tr(i,k-1)*delp(i,k-1) &
         +zx_coef(i,k-1)*zx_ctr(i,k-1) )/zx_buf(i)
      zx_dtr(i,k) = zx_coef(i,k)/zx_buf(i)
    ENDDO
  ENDDO
  DO i = 1, klon
    local_tr_tl(i,klev) = ( local_tr_tl(i,klev)*delp(i,klev) &
       +zx_coef(i,klev)*zx_ctr_tl(i,klev) ) &
       / ( delp(i,klev) + zx_coef(i,klev) &
       -zx_coef(i,klev)*zx_dtr(i,klev) )
    local_tr(i,klev) = ( local_tr(i,klev)*delp(i,klev) &
       +zx_coef(i,klev)*zx_ctr(i,klev) ) &
       / ( delp(i,klev) + zx_coef(i,klev) &
       -zx_coef(i,klev)*zx_dtr(i,klev) )
  ENDDO
  DO k = klev-1, 1, -1
    DO i = 1, klon
      local_tr_tl(i,k) = zx_ctr_tl(i,k+1) + zx_dtr(i,k+1)*local_tr_tl(i,k+1)
      local_tr(i,k) = zx_ctr(i,k+1) + zx_dtr(i,k+1)*local_tr(i,k+1)
    ENDDO
  ENDDO
  DO k = 1, klev
    DO i = 1, klon
      d_tr_tl(i,k) = local_tr_tl(i,k) - tr_tl(i,k)
      d_tr(i,k) = local_tr(i,k) - tr(i,k)
    ENDDO
  ENDDO
  !
  RETURN
END SUBROUTINE cltrac_tl
!
!------------------------------------------------------------------------
!------------------------------------------------------------------------
! 
SUBROUTINE cltrac_ad(dtime,coef,t,tr,flux,paprs,pplay,delp, &
   d_tr,flux_ad,tr_ad,d_tr_ad)
  USE dimphy
  IMPLICIT NONE
  !======================================================================
  ! Auteur(s): O. Boucher (LOA/LMD) date: 19961127
  !            inspire de clvent
  ! Objet: diffusion verticale de traceurs avec flux fixe a la surface
  !        ou/et flux du TYPE c-drag
  !======================================================================
  ! Arguments:
  ! dtime----input-R- intervalle du temps (en second)
  ! coef-----input-R- le coefficient d'echange (m**2/s) l>1
  ! t--------input-R- temperature (K)
  ! tr-------input-R- la q. de traceurs
  ! flux-----input-R- le flux de traceurs a la surface
  ! paprs----input-R- pression a inter-couche (Pa)
  ! pplay----input-R- pression au milieu de couche (Pa)
  ! delp-----input-R- epaisseur de couche (Pa)
  ! cdrag----input-R- cdrag pour le flux de surface (non active)
  ! tr0------input-R- traceurs a la surface ou dans l'ocean (non active)
  ! d_tr-----output-R- le changement de tr
  ! flux_tr--output-R- flux de tr
  !======================================================================
  REAL :: dtime
  REAL :: coef(klon,klev)
  REAL :: t(klon,klev), tr(klon,klev)
  REAL :: tr_ad(klon,klev)
  REAL :: paprs(klon,klev+1), pplay(klon,klev), delp(klon,klev)
  REAL :: d_tr(klon,klev)
  REAL :: d_tr_ad(klon,klev)
  REAL :: flux(klon), cdrag(klon), tr0(klon)
  REAL :: flux_ad(klon)
  !======================================================================
  INCLUDE "YOMCST.h"
  !======================================================================
  INTEGER :: i, k
  REAL :: zx_ctr(klon,2:klev)
  REAL :: zx_ctr_ad(klon,2:klev)
  REAL :: zx_dtr(klon,2:klev)
  REAL :: zx_buf(klon)
  REAL :: zx_coef(klon,klev)
  REAL :: local_tr(klon,klev)
  REAL :: local_tr_ad(klon,klev)
  REAL :: zx_alf1(klon), zx_alf2(klon), zx_flux(klon)
  REAL :: zx_flux_ad(klon)
  !======================================================================
  
  zx_ctr_ad(:,:) = 0.
  local_tr_ad(:,:) = 0.
  flux_ad(:) = 0.
  zx_flux_ad(:) = 0.
  
  local_tr(:,:) = tr(:,:)
  !
  !======================================================================
  DO i = 1, klon
    zx_alf1(i) = (paprs(i,1)-pplay(i,2))/(pplay(i,1)-pplay(i,2))
    zx_alf2(i) = 1.0 - zx_alf1(i)
    zx_flux(i) =  -flux(i)*dtime*RG
    !--pour le moment le flux est prescrit
    !--cdrag et zx_coef(1) vaut 0
    cdrag(i) = 0.0 
    tr0(i) = 0.0
    zx_coef(i,1) = cdrag(i)*dtime*RG 
  ENDDO
  !======================================================================
  DO k = 2, klev
    DO i = 1, klon
      zx_coef(i,k) = coef(i,k)*RG/(pplay(i,k-1)-pplay(i,k)) &
         *(paprs(i,k)*2/(t(i,k)+t(i,k-1))/RD)**2
      zx_coef(i,k) = zx_coef(i,k)*dtime*RG
    ENDDO
  ENDDO
  !======================================================================
  DO i = 1, klon
    zx_buf(i) = delp(i,1) + zx_coef(i,1)*zx_alf1(i) + zx_coef(i,2)
    zx_ctr(i,2) = (local_tr(i,1)*delp(i,1)+ &
       zx_coef(i,1)*tr0(i)-zx_flux(i))/zx_buf(i)
    zx_dtr(i,2) = (zx_coef(i,2)-zx_alf2(i)*zx_coef(i,1)) /  &
       zx_buf(i)
  ENDDO
  !
  DO k = 3, klev
    DO i = 1, klon
      zx_buf(i) = delp(i,k-1) + zx_coef(i,k) &
         + zx_coef(i,k-1)*(1.-zx_dtr(i,k-1))
      zx_ctr(i,k) = (local_tr(i,k-1)*delp(i,k-1) &
         +zx_coef(i,k-1)*zx_ctr(i,k-1) )/zx_buf(i)
      zx_dtr(i,k) = zx_coef(i,k)/zx_buf(i)
    ENDDO
  ENDDO
  DO i = 1, klon
    local_tr(i,klev) = ( local_tr(i,klev)*delp(i,klev) &
       +zx_coef(i,klev)*zx_ctr(i,klev) ) &
       / ( delp(i,klev) + zx_coef(i,klev) &
       -zx_coef(i,klev)*zx_dtr(i,klev) )
  ENDDO
  DO k = klev-1, 1, -1
    DO i = 1, klon
      local_tr(i,k) = zx_ctr(i,k+1) + zx_dtr(i,k+1)*local_tr(i,k+1)
    ENDDO
  ENDDO
  DO k = 1, klev
    DO i = 1, klon
      d_tr(i,k) = local_tr(i,k) - tr(i,k)
      !Start of AD
      local_tr_ad(i,k) = local_tr_ad(i,k) + d_tr_ad(i,k)
      tr_ad(i,k) = tr_ad(i,k) - d_tr_ad(i,k)
      d_tr_ad(i,k) = 0.
    ENDDO
  ENDDO
  DO k = 1,klev-1
    DO i = 1, klon
      zx_ctr_ad(i,k+1) = zx_ctr_ad(i,k+1) + local_tr_ad(i,k)
      local_tr_ad(i,k+1) = local_tr_ad(i,k+1) + zx_dtr(i,k+1)*local_tr_ad(i,k)
      local_tr_ad(i,k) = 0.
    ENDDO
  ENDDO
  DO i = 1, klon
    zx_ctr_ad(i,klev) = zx_ctr_ad(i,klev) + zx_coef(i,klev)*local_tr_ad(i,klev) &
       / ( delp(i,klev) + zx_coef(i,klev) -zx_coef(i,klev)*zx_dtr(i,klev) )
    local_tr_ad(i,klev) = local_tr_ad(i,klev)*delp(i,klev) &
       / ( delp(i,klev) + zx_coef(i,klev) -zx_coef(i,klev)*zx_dtr(i,klev) )
  ENDDO
  DO k = klev,3,-1
    DO i = 1, klon
      zx_buf(i) = delp(i,k-1) + zx_coef(i,k) &
         + zx_coef(i,k-1)*(1.-zx_dtr(i,k-1))
      zx_ctr_ad(i,k-1) = zx_ctr_ad(i,k-1) + zx_coef(i,k-1)*zx_ctr_ad(i,k) &
         / zx_buf(i)
      local_tr_ad(i,k-1) = local_tr_ad(i,k-1) + zx_ctr_ad(i,k)*delp(i,k-1) &
         / zx_buf(i)
      zx_ctr_ad(i,k) = 0.
    ENDDO
  ENDDO
  DO i = 1, klon
    zx_buf(i) = delp(i,1) + zx_coef(i,1)*zx_alf1(i) + zx_coef(i,2)
    zx_flux_ad(i) = zx_flux_ad(i) - zx_ctr_ad(i,2)/zx_buf(i)
    local_tr_ad(i,1) = local_tr_ad(i,1) + zx_ctr_ad(i,2) * delp(i,1)/zx_buf(i)
    zx_ctr_ad(i,2) = 0.
  ENDDO
  DO i = 1, klon
    flux_ad(i) = flux_ad(i) - zx_flux_ad(i) *dtime*RG
    zx_flux_ad(i) = 0.
  ENDDO
  tr_ad(:,:) = tr_ad(:,:) + local_tr_ad(:,:) 
  local_tr_ad(:,:) = 0.
  !
  RETURN
END SUBROUTINE cltrac_ad
