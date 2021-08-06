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

SUBROUTINE  exner_hyb ( ngrid, ps, p,alpha,beta, pks, pk, pkf )
  !
  !     Auteurs :  P.Le Van  , Fr. Hourdin  .
  !    ..........
  !
  !    ....  ngrid, ps,p             sont des argum.d'entree  au sous-prog ...
  !    .... alpha,beta, pks,pk,pkf   sont des argum.de sortie au sous-prog ...
  !
  !   ************************************************************************
  !    Calcule la fonction d'Exner pk = Cp * p ** kappa , aux milieux des 
  !    couches .   Pk(l) sera calcule aux milieux  des couches l ,entre les
  !    pressions p(l) et p(l+1) ,definis aux interfaces des llm couches .
  !   ************************************************************************
  !  .. N.B : Au sommet de l'atmosphere,  p(llm+1) = 0. , et ps et pks sont
  !    la pression et la fonction d'Exner  au  sol  .
  !
  !                                 -------- z                                   
  !    A partir des relations  ( 1 ) p*dz(pk) = kappa *pk*dz(p)      et
  !                            ( 2 ) pk(l) = alpha(l)+ beta(l)*pk(l-1)
  !    ( voir note de Fr.Hourdin )  ,
  !
  !    on determine successivement , du haut vers le bas des couches, les 
  !    coef. alpha(llm),beta(llm) .,.,alpha(l),beta(l),,,alpha(2),beta(2), 
  !    puis pk(ij,1). Ensuite ,on calcule,du bas vers le haut des couches,  
  !     pk(ij,l)  donne  par la relation (2),  pour l = 2 a l = llm .
  !
  !
  USE parallel
  IMPLICIT NONE
  !
 include "dimensions.h"
 include "paramet.h"
 include "comconst.h"
 include "comgeom.h"
 include "comvert.h"
 include "serre.h"
  
  INTEGER :: ngrid
  REAL :: p(ngrid,llmp1),pk(ngrid,llm),pkf(ngrid,llm)
  REAL :: ps(ngrid),pks(ngrid), alpha(ngrid,llm),beta(ngrid,llm)
  
  !    .... variables locales   ...
  
  INTEGER :: l, ij
  REAL :: unpl2k,dellta
  
  REAL :: ppn(iim),pps(iim)
  REAL :: xpn, xps
  REAL :: SSUM
  EXTERNAL :: SSUM
  INTEGER :: ije,ijb,jje,jjb
  !
  !$OMP MASTER           
  unpl2k    = 1.+ 2.* kappa
  !
  ijb=ij_begin
  ije=ij_end
  
  DO   ij  = ijb, ije
    pks(ij) = cpp * ( ps(ij)/preff ) ** kappa
  ENDDO
  
  IF (north_pole) THEN
      DO  ij   = 1, iim
        ppn(ij) = area(   ij   ) * pks(  ij     )
      ENDDO
      xpn      = SSUM(iim,ppn,1) /apoln
      
      DO ij   = 1, iip1
        pks(   ij     )  =  xpn
      ENDDO
  ENDIF
  
  IF (south_pole) THEN
      DO  ij   = 1, iim
        pps(ij) = area(ij+ip1jm) * pks(ij+ip1jm )
      ENDDO
      xps      = SSUM(iim,pps,1) /apols 
      
      DO ij   = 1, iip1
        pks( ij+ip1jm )  =  xps
      ENDDO
  ENDIF

  !
  !
  !    .... Calcul des coeff. alpha et beta  pour la couche l = llm ..
  !
  DO     ij      = ijb,ije
    alpha(ij,llm) = 0.
    beta (ij,llm) = 1./ unpl2k
  ENDDO
  !
  !     ... Calcul des coeff. alpha et beta  pour l = llm-1  a l = 2 ...
  !
  DO l = llm -1 , 2 , -1
    !
    DO ij = ijb, ije
      dellta = p(ij,l)* unpl2k + p(ij,l+1)* ( beta(ij,l+1)-unpl2k )
      alpha(ij,l)  = - p(ij,l+1) / dellta * alpha(ij,l+1)
      beta (ij,l)  =   p(ij,l  ) / dellta
    ENDDO
    !
  ENDDO
  
  !
  !  ***********************************************************************
  !     .....  Calcul de pk pour la couche 1 , pres du sol  ....
  !
  
  DO   ij   = ijb, ije
    pk(ij,1) = ( p(ij,1)*pks(ij) - 0.5*alpha(ij,2)*p(ij,2) )  / &
       (  p(ij,1)* (1.+kappa) + 0.5*( beta(ij,2)-unpl2k )* p(ij,2) )
  ENDDO
  !
  !    ..... Calcul de pk(ij,l) , pour l = 2 a l = llm  ........
  !
  DO l = 2, llm
    DO   ij   = ijb, ije
      pk(ij,l) = alpha(ij,l) + beta(ij,l) * pk(ij,l-1)
    ENDDO
  ENDDO
  !
  !
  !      CALL SCOPY   ( ngrid * llm, pk, 1, pkf, 1 )
  pkf(ijb:ije,1:llm)=pk(ijb:ije,1:llm)
  !$OMP END MASTER
  !$OMP BARRIER
  
  jjb=jj_begin
  jje=jj_end
  CALL filtreg( pkf,jjb,jje, jmp1, llm, 2, 1, .TRUE., 1 )
  
  RETURN
END SUBROUTINE exner_hyb
