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

SUBROUTINE diverg(klevel,x,y,div)
  !
  !     P. Le Van
  !
  !  *********************************************************************
  !  ... calcule la divergence a tous les niveaux d'1 vecteur de compos. 
  !     x et y...
  !              x et y  etant des composantes covariantes   ...
  !  *********************************************************************
  IMPLICIT NONE
  !
  !      x  et  y  sont des arguments  d'entree pour le s-prog
  !        div      est  un argument  de sortie pour le s-prog
  !
  !
  !   ---------------------------------------------------------------------
  !
  !    ATTENTION : pendant ce s-pg , ne pas toucher au COMMON/scratch/  .
  !
  !   ---------------------------------------------------------------------
 include "dimensions.h"
 include "paramet.h"
 include "comgeom.h"
  !
  !    ..........          variables en arguments    ...................
  !
  INTEGER :: klevel
  REAL :: x( ip1jmp1,klevel ),y( ip1jm,klevel ),div( ip1jmp1,klevel )
  INTEGER :: l,ij
  !
  !    ...............     variables  locales   .........................

  REAL :: aiy1( iip1 ) , aiy2( iip1 )
  REAL :: sumypn,sumyps
  !    ...................................................................
  !
  EXTERNAL :: SSUM
  REAL :: SSUM
  !
  !
  DO l = 1,klevel
    !
    DO  ij = iip2, ip1jm - 1
      div( ij + 1, l )     = & 
         cvusurcu( ij+1 ) * x( ij+1,l ) - cvusurcu( ij ) * x( ij , l) + &
         cuvsurcv(ij-iim) * y(ij-iim,l) - cuvsurcv(ij+1) * y(ij+1,l) 
    ENDDO
    !
    !     ....  correction pour  div( 1,j,l)  ......
    !     ....   div(1,j,l)= div(iip1,j,l) ....
    !
    !DIR$ IVDEP
    DO  ij = iip2,ip1jm,iip1
      div( ij,l ) = div( ij + iim,l )
    ENDDO
    !
    !     ....  calcul  aux poles  .....
    !
    DO  ij  = 1,iim
      aiy1(ij) =    cuvsurcv(    ij       ) * y(     ij     , l )
      aiy2(ij) =    cuvsurcv( ij+ ip1jmi1 ) * y( ij+ ip1jmi1, l )
    ENDDO
    sumypn = SSUM ( iim,aiy1,1 ) / apoln
    sumyps = SSUM ( iim,aiy2,1 ) / apols
    !
    DO  ij = 1,iip1
      div(     ij    , l ) = - sumypn
      div( ij + ip1jm, l ) =   sumyps
    ENDDO
  ENDDO
  !
  
  DO l = 1, klevel
    DO ij = iip2,ip1jm
      div(ij,l) = div(ij,l) * unsarea(ij)
    ENDDO
  ENDDO
  !
  RETURN
END SUBROUTINE diverg
