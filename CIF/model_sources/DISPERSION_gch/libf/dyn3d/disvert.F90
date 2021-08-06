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

SUBROUTINE disvert(pa,preff,ap,bp,dpres,presnivs,nivsigs,nivsig)
  
  !    Auteur :  P. Le Van .
  !
  IMPLICIT NONE
  
 include "dimensions.h"
 include "paramet.h"
  !
  !=======================================================================
  !
  !
  !    s = sigma ** kappa   :  coordonnee  verticale
  !    dsig(l)            : epaisseur de la couche l ds la coord.  s
  !    sig(l)             : sigma a l'interface des couches l et l-1
  !    ds(l)              : distance entre les couches l et l-1 en coord.s
  !
  !=======================================================================
  !
  REAL :: pa,preff
  REAL :: ap(llmp1),bp(llmp1),dpres(llm),nivsigs(llm),nivsig(llmp1)
  REAL :: presnivs(llm)
  !
  !   declarations:
  !   -------------
  !
  REAL :: sig(llm+1),dsig(llm)
  !
  INTEGER :: l
  REAL :: snorm
  REAL :: alpha,beta,gama,delta,deltaz,h
  INTEGER :: np,ierr
  REAL :: pi,x
  
  REAL :: SSUM
  EXTERNAL :: SSUM
  !
  !-----------------------------------------------------------------------
  !
  pi=2.*ASIN(1.)
  
  DO l=1,llm
    nivsigs(l) = FLOAT(l)
  ENDDO
  
  DO l=1,llmp1
    nivsig(l)= FLOAT(l)
  ENDDO
  

  !
  !   .....  ap, bp, pa et preff sont lus  sur les fichiers start .....
  !

  
  
  PRINT *,' BP '
  PRINT *,  bp
  PRINT *,' AP '
  PRINT *,  ap
  
  DO l = 1, llm
    dpres(l) = bp(l) - bp(l+1)
    presnivs(l) = 0.5 *( ap(l)+bp(l)*preff + ap(l+1)+bp(l+1)*preff )
  ENDDO
  
  PRINT *,' PRESNIVS '
  PRINT *,presnivs
  
  RETURN
END SUBROUTINE disvert
