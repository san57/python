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

SUBROUTINE fxyhyper ( yzoom, grossy, dzoomy,tauy  , &
   xzoom, grossx, dzoomx,taux  , &
   rlatu,yprimu,rlatv,yprimv,rlatu1,  yprimu1,  rlatu2,  yprimu2  , & 
   rlonu,xprimu,rlonv,xprimv,rlonm025,xprimm025,rlonp025,xprimp025)

  IMPLICIT NONE
  !
  !      Auteur :  P. Le Van .
  !
  !      d'apres  formulations de R. Sadourny .
  !
  !
  !     Ce spg calcule les latitudes( routine fyhyp ) et longitudes( fxhyp )
  !            par des  fonctions  a tangente hyperbolique .
  !
  !     Il y a 3 parametres ,en plus des coordonnees du centre du zoom (xzoom
  !                      et  yzoom )   :  
  !
  !     a) le grossissement du zoom  :  grossy  ( en y ) et grossx ( en x )
  !     b) l' extension     du zoom  :  dzoomy  ( en y ) et dzoomx ( en x )
  !     c) la raideur de la transition du zoom  :   taux et tauy   
  !
  !  N.B : Il vaut mieux avoir   :   grossx * dzoomx <  pi    ( radians )
  ! ******
  !                  et              grossy * dzoomy <  pi/2  ( radians )
  !
 include "dimensions.h"
 include "paramet.h"
  
  
  !
  !   .....  Arguments  ...
  REAL :: xzoom,yzoom,grossx,grossy,dzoomx,dzoomy,taux,tauy
  REAL :: rlatu(jjp1), yprimu(jjp1),rlatv(jjm), yprimv(jjm), &
     rlatu1(jjm), yprimu1(jjm), rlatu2(jjm), yprimu2(jjm)
  REAL :: rlonu(iip1),xprimu(iip1),rlonv(iip1),xprimv(iip1), &
     rlonm025(iip1),xprimm025(iip1), rlonp025(iip1),xprimp025(iip1)
  REAL*8 :: dxmin, dxmax , dymin, dymax
  
  !   ....   var. locales   .....
  !
  INTEGER :: i,j
  !
  
  CALL fyhyp ( yzoom, grossy, dzoomy,tauy  , &
     rlatu, yprimu,rlatv,yprimv,rlatu2,yprimu2,rlatu1,yprimu1 , &
     dymin,dymax                                               )
  
  CALL fxhyp(xzoom,grossx,dzoomx,taux,rlonm025,xprimm025,rlonv, &
     xprimv,rlonu,xprimu,rlonp025,xprimp025 , dxmin,dxmax         )
  
  
  DO i = 1, iip1
    IF(rlonp025(i)<rlonv(i))  THEN
        WRITE(6,*) ' Attention !  rlonp025 < rlonv',i
        STOP
    ENDIF
    
    IF(rlonv(i)<rlonm025(i))  THEN 
        WRITE(6,*) ' Attention !  rlonm025 > rlonv',i
        STOP
    ENDIF
    
    IF(rlonp025(i)>rlonu(i))  THEN
        WRITE(6,*) ' Attention !  rlonp025 > rlonu',i
        STOP
    ENDIF
  ENDDO
  
  WRITE(6,*) '  *** TEST DE COHERENCE  OK    POUR   FX **** '
  
  !
  DO j = 1, jjm
    !
    IF(rlatu1(j)<=rlatu2(j))   THEN
        WRITE(6,*)'Attention ! rlatu1 < rlatu2 ',rlatu1(j), rlatu2(j),j
        STOP 13
    ENDIF
    !
    IF(rlatu2(j)<=rlatu(j+1))  THEN
        WRITE(6,*)'Attention ! rlatu2 < rlatup1 ',rlatu2(j),rlatu(j+1),j
        STOP 14
    ENDIF
    !
    IF(rlatu(j)<=rlatu1(j))    THEN
        WRITE(6,*)' Attention ! rlatu < rlatu1 ',rlatu(j),rlatu1(j),j
        STOP 15
    ENDIF
    !
    IF(rlatv(j)<=rlatu2(j))    THEN
        WRITE(6,*)' Attention ! rlatv < rlatu2 ',rlatv(j),rlatu2(j),j
        STOP 16
    ENDIF
    !
    IF(rlatv(j)>=rlatu1(j))    THEN
        WRITE(6,*)' Attention ! rlatv > rlatu1 ',rlatv(j),rlatu1(j),j
        STOP 17
    ENDIF
    !
    IF(rlatv(j)>=rlatu(j))     THEN
        WRITE(6,*) ' Attention ! rlatv > rlatu ',rlatv(j),rlatu(j),j
        STOP 18
    ENDIF
    !
  ENDDO
  !
  WRITE(6,*) '  *** TEST DE COHERENCE  OK    POUR   FY **** '
  !
  WRITE(6,18)
  WRITE(6,*) '  Latitudes  '
  WRITE(6,*) ' *********** '
  WRITE(6,18)
  WRITE(6,3)  dymin, dymax
  WRITE(6,*) ' Si cette derniere est trop lache , modifiez les par &
     &ametres  grossism , tau , dzoom pour Y et repasser ! '
  !
  WRITE(6,18)
  WRITE(6,*) '  Longitudes  '
  WRITE(6,*) ' ************ '
  WRITE(6,18)
  WRITE(6,3)  dxmin, dxmax
  WRITE(6,*) ' Si cette derniere est trop lache , modifiez les par &
     &ametres  grossism , tau , dzoom pour Y et repasser ! '
  WRITE(6,18)
  !
3 FORMAT(1x, ' Au centre du zoom , la longueur de la maille est', &
     ' d environ ',f8.2 ,' degres  ', &
     ' alors que la maille en dehors de la zone du zoom est d environ &
     &', f8.2,' degres ' )
18 FORMAT(/)
  
  RETURN
END SUBROUTINE fxyhyper

