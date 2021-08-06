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

SUBROUTINE cvltr(pdtime,da, phi,mp,paprs,x,upd,dnd,wght,dx)
  USE dimphy
  USE Write_Field_phy
  IMPLICIT NONE 
!=====================================================================
! Objet : convection des traceurs / KE
! Auteurs: M-A Filiberti and J-Y Grandpeix
!=====================================================================

  include "YOMCST.h"
  include "YOECUMF.h" 


! Entree
  REAL  :: pdtime
  REAL  :: da(klon,klev)
  REAL  :: phi(klon,klev,klev)
  REAL  :: mp(klon,klev)
  REAL  :: wght(klon,klev)
  REAL  :: paprs(klon,klev+1) ! pression aux 1/2 couches (bas en haut)
  REAL  :: x(klon,klev)     ! q de traceur (bas en haut) 
  REAL  :: upd(klon,klev)   ! saturated updraft mass flux
  REAL  :: dnd(klon,klev)   ! saturated downdraft mass flux

! Sortie
  REAL  :: dx(klon,klev) ! tendance de traceur  (bas en haut)

! Variables locales     
  REAL  :: zmd(klon,klev,klev)
  REAL  :: za(klon,klev,klev)
  REAL  :: zmfd(klon,klev)
  REAL  :: zmfa(klon,klev)    
  REAL  :: zmfp(klon,klev)
  REAL  :: zmfu(klon,klev)    
  REAL  :: deltap(klon,klev)
  INTEGER  :: i,k,j 
  REAL     :: pdtimeRG
  REAL     :: qs(klon)
  REAL     :: conserv

!TEST - DA NAN
!  DO j=1,klev
!  DO i=1,klon
!   if (da(i,j) /= da(i,j)) then
!    print*,'NANDA'
!   end if  
!  END DO
!  END DO


! =========================================
! calcul des tendances liees au downdraft
! =========================================
  DO i=1,klon
    qs(i)=0.
  END DO
!cdir collapse
  DO j=1,klev
  DO i=1,klon
    zmfd(i,j)=0.
    zmfa(i,j)=0.
    zmfu(i,j)=0.
    zmfp(i,j)=0.
  END DO
  END DO
!cdir collapse
  DO k=1,klev
  DO j=1,klev
  DO i=1,klon
    zmd(i,j,k)=0.
    za (i,j,k)=0.
  END DO
  END DO
  END DO

! calcul de la matrice d echange
! matrice de distribution de la masse entrainee en k

  DO k=1,klev-1
     DO i=1,klon
        zmd(i,k,k)=max(0.,mp(i,k)-mp(i,k+1))
     END DO
  END DO
  DO k=2,klev
     DO j=k-1,1,-1
        DO i=1,klon
!          if(mp(i,j+1).ne.0) then !Changed on 01/08/2018 to fit the on-line code
           if(mp(i,j+1).gt.1.e-10) then
              zmd(i,j,k)=zmd(i,j+1,k)*min(1.,mp(i,j)/mp(i,j+1))
           ENDif
        END DO
     END DO
  END DO
  DO k=1,klev
     DO j=1,klev-1
        DO i=1,klon
           za(i,j,k)=max(0.,zmd(i,j+1,k)-zmd(i,j,k))
        END DO
     END DO
  END DO
!
! rajout du terme lie a l ascendance induite
!
  DO j=2,klev
     DO i=1,klon
        za(i,j,j-1)=za(i,j,j-1)+mp(i,j)
     END DO
  END DO
!
! tendances
!            
  DO k=1,klev
     DO j=1,klev
        DO i=1,klon
           zmfd(i,j)=zmfd(i,j)+za(i,j,k)*(x(i,k)-x(i,j))
        END DO
     END DO
  END DO
!
! =========================================
! calcul des tendances liees aux flux satures
! =========================================

!RL!
!determination concentration dans la couche dalimentation

  DO j=1,klev
     DO i=1,klon
        qs(i)=qs(i)+wght(i,j)*x(i,j)
     END DO
  END DO
!RL!

  DO j=1,klev
     DO i=1,klon
        zmfa(i,j)=da(i,j)*(qs(i)-x(i,j))     !RL!
     END DO
  END DO
  DO k=1,klev
     DO j=1,klev
        DO i=1,klon
           zmfp(i,j)=zmfp(i,j)+phi(i,j,k)*(x(i,k)-x(i,j))
        END DO
     END DO
  END DO

  DO j=1,klev-1
     DO i=1,klon
        zmfu(i,j)=max(0.,upd(i,j+1)+dnd(i,j+1))*(x(i,j+1)-x(i,j))
     END DO
  END DO
  DO j=2,klev
     DO i=1,klon
        zmfu(i,j)=zmfu(i,j)+min(0.,upd(i,j)+dnd(i,j))*(x(i,j)-x(i,j-1))   
     END DO
  END DO

! =========================================
! calcul final des tendances
! =========================================
  DO k=1, klev
     DO i=1, klon
        deltap(i,k)=paprs(i,k)-paprs(i,k+1)
     ENDDO
  ENDDO
  pdtimeRG=pdtime*RG
!cdir collapse
  DO k=1, klev
     DO i=1, klon
        dx(i,k)=(zmfd(i,k)+zmfu(i,k)       &
                +zmfa(i,k)+zmfp(i,k))*pdtimeRG/deltap(i,k)
     ENDDO
  ENDDO

!test conservation
! conserv=0.
!  DO k=1,klev
!   DO i=1,klon
!    conserv = conserv + dx(i,k)*deltap(i,k)/RG
!   ENDDO
!  ENDDO
!  PRINT*,'CONSERVATION-CVLTR',conserv

     
END SUBROUTINE cvltr
!---------------------------------------------------------------------
!---------------------------------------------------------------------
SUBROUTINE cvltr_tl(pdtime,da,phi,mp,paprs,x,x_tl,upd,dnd,wght,dx,dx_tl)

  USE dimphy
  IMPLICIT NONE 
!=====================================================================
! Objet : convection des traceurs / KE
! Auteurs: M-A Filiberti and J-Y Grandpeix
!=====================================================================

  include "YOMCST.h"
  include "YOECUMF.h" 

! Entree
  REAL  :: pdtime
  REAL  :: da(klon,klev)
  REAL  :: phi(klon,klev,klev)
  REAL  :: mp(klon,klev)
  REAL  :: wght(klon,klev)
  REAL  :: paprs(klon,klev+1) ! pression aux 1/2 couches (bas en haut)
  REAL  :: x(klon,klev)     ! q de traceur (bas en haut) 
  REAL  :: x_tl(klon,klev)
  REAL  :: upd(klon,klev)   ! saturated updraft mass flux
  REAL  :: dnd(klon,klev)   ! saturated downdraft mass flux

! Sortie
  REAL  :: dx(klon,klev) ! tendance de traceur  (bas en haut)
  REAL  :: dx_tl(klon,klev)

! Variables locales     
  REAL  :: zmd(klon,klev,klev)
  REAL  :: za(klon,klev,klev)
  REAL  :: zmfd(klon,klev)
  REAL  :: zmfd_tl(klon,klev)
  REAL  :: zmfa(klon,klev)
  REAL  :: zmfa_tl(klon,klev)
  REAL  :: zmfp(klon,klev)
  REAL  :: zmfp_tl(klon,klev)
  REAL  :: zmfu(klon,klev)
  REAL  :: zmfu_tl(klon,klev)
  REAL  :: deltap(klon,klev)
  INTEGER  :: i,k,j
  REAL     :: pdtimeRG
  REAL     :: qs(klon)
  REAL     :: qs_tl(klon)


! =========================================
! calcul des tendances liees au downdraft
! =========================================
  DO i=1,klon
    qs_tl(i)=0.
    qs(i)=0.
  END DO
!cdir collapse
  DO j=1,klev
  DO i=1,klon
    zmfd_tl(i,j)=0.
    zmfd(i,j)=0.
    zmfa_tl(i,j)=0.
    zmfa(i,j)=0.
    zmfu_tl(i,j)=0.
    zmfu(i,j)=0.
    zmfp_tl(i,j)=0.
    zmfp(i,j)=0.
  END DO
  END DO

!cdir collapse
  DO k=1,klev
  DO j=1,klev
  DO i=1,klon
    zmd(i,j,k)=0.
    za (i,j,k)=0.
  END DO
  END DO
  END DO

! calcul de la matrice d echange
! matrice de distribution de la masse entrainee en k

  DO k=1,klev-1
     DO i=1,klon
        zmd(i,k,k)=max(0.,mp(i,k)-mp(i,k+1))
     END DO
  END DO
  DO k=2,klev
     DO j=k-1,1,-1
        DO i=1,klon
!          if(mp(i,j+1).ne.0) then !Changed on 01/08/2018 to fit the on-line code
           if(mp(i,j+1).gt.1.e-10) then
              zmd(i,j,k)=zmd(i,j+1,k)*min(1.,mp(i,j)/mp(i,j+1))
           ENDif
        END DO
     END DO
  END DO
  DO k=1,klev
     DO j=1,klev-1
        DO i=1,klon
           za(i,j,k)=max(0.,zmd(i,j+1,k)-zmd(i,j,k))
        END DO
     END DO
  END DO
!
! rajout du terme lie a l ascendance induite
!
  DO j=2,klev
     DO i=1,klon
        za(i,j,j-1)=za(i,j,j-1)+mp(i,j)
     END DO
  END DO
!
! tendances
!            
  DO k=1,klev
     DO j=1,klev
        DO i=1,klon
           zmfd_tl(i,j)=zmfd_tl(i,j)+za(i,j,k)*(x_tl(i,k)-x_tl(i,j))
           zmfd(i,j)=zmfd(i,j)+za(i,j,k)*(x(i,k)-x(i,j))
        END DO
     END DO
  END DO
!
! =========================================
! calcul des tendances liees aux flux satures
! =========================================


!determination concentration dans la couche dalimentation
  DO j=1,klev
     DO i=1,klon
        qs_tl(i)=qs_tl(i)+wght(i,j)*x_tl(i,j)
        qs(i)=qs(i)+wght(i,j)*x(i,j)
     END DO
  END DO


  DO j=1,klev
     DO i=1,klon
        zmfa_tl(i,j)=da(i,j)*(qs_tl(i)-x_tl(i,j)) 
        zmfa(i,j)=da(i,j)*(qs(i)-x(i,j))    
     END DO
  END DO
  DO k=1,klev
     DO j=1,klev
        DO i=1,klon
           zmfp_tl(i,j)=zmfp_tl(i,j)+phi(i,j,k)*(x_tl(i,k)-x_tl(i,j))
           zmfp(i,j)=zmfp(i,j)+phi(i,j,k)*(x(i,k)-x(i,j))
        END DO
     END DO
  END DO

  DO j=1,klev-1
     DO i=1,klon
        zmfu_tl(i,j)=max(0.,upd(i,j+1)+dnd(i,j+1))*(x_tl(i,j+1)-x_tl(i,j))
        zmfu(i,j)=max(0.,upd(i,j+1)+dnd(i,j+1))*(x(i,j+1)-x(i,j))
     END DO
  END DO
  DO j=2,klev
     DO i=1,klon
        zmfu_tl(i,j)=zmfu_tl(i,j)+min(0.,upd(i,j)+dnd(i,j))*(x_tl(i,j)-x_tl(i,j-1))
        zmfu(i,j)=zmfu(i,j)+min(0.,upd(i,j)+dnd(i,j))*(x(i,j)-x(i,j-1))
     END DO
  END DO


! =========================================
! calcul final des tendances
! =========================================
  DO k=1, klev
     DO i=1, klon
        deltap(i,k)=paprs(i,k)-paprs(i,k+1)
     ENDDO
  ENDDO
  pdtimeRG=pdtime*RG
!cdir collapse
  DO k=1, klev
     DO i=1, klon
        dx_tl(i,k)=(zmfd_tl(i,k)+zmfu_tl(i,k)       &
                +zmfa_tl(i,k)+zmfp_tl(i,k))*pdtimeRG/deltap(i,k)
        dx(i,k)=(zmfd(i,k)+zmfu(i,k)       &
                +zmfa(i,k)+zmfp(i,k))*pdtimeRG/deltap(i,k)
     ENDDO
  ENDDO

END SUBROUTINE cvltr_tl
!
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!
SUBROUTINE cvltr_ad(pdtime,da, phi, mp,paprs,x,x_ad,upd,dnd,wght,dx,dx_ad)

  USE dimphy
  IMPLICIT NONE 
!=====================================================================
! Objet : convection des traceurs / KE
! Auteurs: M-A Filiberti and J-Y Grandpeix
!=====================================================================

  include "YOMCST.h"
  include "YOECUMF.h" 

! Entree
  REAL  :: pdtime
  REAL  :: da(klon,klev)
  REAL  :: phi(klon,klev,klev)
  REAL  :: mp(klon,klev)
  REAL  :: wght(klon,klev)
  REAL  :: paprs(klon,klev+1) ! pression aux 1/2 couches (bas en haut)
  REAL  :: x(klon,klev)     ! q de traceur (bas en haut) 
  REAL  :: x_ad(klon,klev)
  REAL  :: upd(klon,klev)   ! saturated updraft mass flux
  REAL  :: dnd(klon,klev)   ! saturated downdraft mass flux
    
! Sortie
  REAL  :: dx(klon,klev) ! tendance de traceur  (bas en haut)
  REAL  :: dx_ad(klon,klev)
  
! Variables locales     
  REAL  :: zmd(klon,klev,klev)
  REAL  :: za(klon,klev,klev)
  REAL  :: zmfd(klon,klev)
  REAL  :: zmfd_ad(klon,klev)
  REAL  :: zmfa(klon,klev)
  REAL  :: zmfa_ad(klon,klev)
  REAL  :: zmfp(klon,klev)
  REAL  :: zmfp_ad(klon,klev)
  REAL  :: zmfu(klon,klev)
  REAL  :: zmfu_ad(klon,klev)
  REAL  :: deltap(klon,klev)
  INTEGER  :: i,k,j
  REAL     :: pdtimeRG
  REAL     :: qs(klon)
  REAL     :: qs_ad(klon)


!=====================================================================
!Initialisation des variables adjointes locales
 zmfd_ad(:,:) = 0.
 zmfa_ad(:,:) = 0.
 zmfp_ad(:,:) = 0.
 zmfu_ad(:,:) = 0.
 qs_ad(:) = 0. 

! =========================================
! calcul des tendances liees au downdraft
! =========================================
  DO i=1,klon
    qs(i)=0.
  END DO
!cdir collapse
  DO j=1,klev
  DO i=1,klon
    zmfd(i,j)=0.
    zmfa(i,j)=0.
    zmfu(i,j)=0.
    zmfp(i,j)=0.
  END DO
  END DO
!cdir collapse
  DO k=1,klev
  DO j=1,klev
  DO i=1,klon
    zmd(i,j,k)=0.
    za (i,j,k)=0.
  END DO
  END DO
  END DO

! calcul de la matrice d echange
! matrice de distribution de la masse entrainee en k

  DO k=1,klev-1
     DO i=1,klon
        zmd(i,k,k)=max(0.,mp(i,k)-mp(i,k+1))
     END DO
  END DO
  DO k=2,klev
     DO j=k-1,1,-1
        DO i=1,klon
!          if(mp(i,j+1).ne.0) then !Changed on 01/08/2018 to fit the on-line code
           if(mp(i,j+1).gt.1.e-10) then
              zmd(i,j,k)=zmd(i,j+1,k)*min(1.,mp(i,j)/mp(i,j+1))
           ENDif
        END DO
     END DO
  END DO
  DO k=1,klev
     DO j=1,klev-1
        DO i=1,klon
           za(i,j,k)=max(0.,zmd(i,j+1,k)-zmd(i,j,k))
        END DO
     END DO
  END DO
!
! rajout du terme lie a l ascendance induite
!
  DO j=2,klev
     DO i=1,klon
        za(i,j,j-1)=za(i,j,j-1)+mp(i,j)
     END DO
  END DO
!
! tendances
!            
  DO k=1,klev
     DO j=1,klev
        DO i=1,klon
           zmfd(i,j)=zmfd(i,j)+za(i,j,k)*(x(i,k)-x(i,j))
        END DO
     END DO
  END DO
!
! =========================================
! calcul des tendances liees aux flux satures
! =========================================

!RL!
!determination concentration dans la couche dalimentation

  DO j=1,klev
     DO i=1,klon
        qs(i)=qs(i)+wght(i,j)*x(i,j)
     END DO
  END DO
!RL!

  DO j=1,klev
     DO i=1,klon
        zmfa(i,j)=da(i,j)*(qs(i)-x(i,j))     !RL!
     END DO
  END DO
  DO k=1,klev
     DO j=1,klev
        DO i=1,klon
           zmfp(i,j)=zmfp(i,j)+phi(i,j,k)*(x(i,k)-x(i,j))
        END DO
     END DO
  END DO

  DO j=1,klev-1
     DO i=1,klon
        zmfu(i,j)=max(0.,upd(i,j+1)+dnd(i,j+1))*(x(i,j+1)-x(i,j))
     END DO
  END DO
  DO j=2,klev
     DO i=1,klon
        zmfu(i,j)=zmfu(i,j)+min(0.,upd(i,j)+dnd(i,j))*(x(i,j)-x(i,j-1))   
     END DO
  END DO

! =========================================
! calcul final des tendances
! =========================================
  DO k=1, klev
     DO i=1, klon
        deltap(i,k)=paprs(i,k)-paprs(i,k+1)
     END DO
  END DO
  pdtimeRG=pdtime*RG
!cdir collapse
  DO k=1, klev
     DO i=1, klon
        dx(i,k)=(zmfd(i,k)+zmfu(i,k)       &
                +zmfa(i,k)+zmfp(i,k))*pdtimeRG/deltap(i,k)
!debut modele adjoint
        zmfd_ad(i,k)=zmfd_ad(i,k)+pdtimeRG/deltap(i,k)*dx_ad(i,k)
        zmfu_ad(i,k)=zmfu_ad(i,k)+pdtimeRG/deltap(i,k)*dx_ad(i,k)
        zmfa_ad(i,k)=zmfa_ad(i,k)+pdtimeRG/deltap(i,k)*dx_ad(i,k)
        zmfp_ad(i,k)=zmfp_ad(i,k)+pdtimeRG/deltap(i,k)*dx_ad(i,k)
        dx_ad(i,k)=0.
     END DO
  END DO

  DO j=klev,2,-1     !DO j=2,klev
     DO i=1,klon
        x_ad(i,j)=x_ad(i,j)+min(0.,upd(i,j)+dnd(i,j))*zmfu_ad(i,j)
        x_ad(i,j-1)=x_ad(i,j-1)-min(0.,upd(i,j)+dnd(i,j))*zmfu_ad(i,j)
     END DO
  END DO
  DO j=klev-1,1,-1     !DO j=1,klev-1,-1
     DO i=1,klon
        x_ad(i,j+1)=x_ad(i,j+1)+max(0.,upd(i,j+1)+dnd(i,j+1))*zmfu_ad(i,j)
        x_ad(i,j)=x_ad(i,j)-max(0.,upd(i,j+1)+dnd(i,j+1))*zmfu_ad(i,j)
     END DO
  END DO     
  zmfu_ad(:,:)=0.

  DO k=klev,1,-1
     DO j=klev,1,-1
        DO i=1,klon
           x_ad(i,k)=x_ad(i,k)+phi(i,j,k)*zmfp_ad(i,j)
           x_ad(i,j)=x_ad(i,j)-phi(i,j,k)*zmfp_ad(i,j)
        END DO
     END DO
  END DO
  zmfp_ad(:,:)=0. 

  DO j=1,klev
     DO i=1,klon
        qs_ad(i)=qs_ad(i)+da(i,j)*zmfa_ad(i,j)
        x_ad(i,j)=x_ad(i,j)-da(i,j)*zmfa_ad(i,j)
        zmfa_ad(i,j)=0.
     END DO
  END DO

  DO j=1,klev
     DO i=1,klon
        x_ad(i,j)=x_ad(i,j)+wght(i,j)*qs_ad(i)
     END DO
  END DO
  qs_ad(:)=0.

  DO k=klev,1,-1
     DO j=klev,1,-1
        DO i=1,klon
           x_ad(i,k)=x_ad(i,k)+za(i,j,k)*zmfd_ad(i,j)
           x_ad(i,j)=x_ad(i,j)-za(i,j,k)*zmfd_ad(i,j)
        END DO
     END DO
  END DO
  zmfd_ad(:,:)=0.


END SUBROUTINE cvltr_ad
