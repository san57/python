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

SUBROUTINE gr_dyn_fi(nfield,im,jm,ngrid,pdyn,pfi)
  IMPLICIT NONE
  !=======================================================================
  !   passage d'un champ de la grille scalaire a la grille physique
  !=======================================================================

  !-----------------------------------------------------------------------
  !   declarations:
  !   -------------

  INTEGER :: im,jm,ngrid,nfield
  REAL :: pdyn(im,jm,nfield)
  REAL :: pfi(ngrid,nfield)
  
  INTEGER :: j,ifield,ig
  EXTERNAL :: SCOPY
  
  !-----------------------------------------------------------------------
  !   calcul:
  !   -------
  
  IF(ngrid/=2+(jm-2)*(im-1)) STOP 'probleme de dim'
  !   traitement des poles
  CALL SCOPY(nfield,pdyn,im*jm,pfi,ngrid)
  CALL SCOPY(nfield,pdyn(1,jm,1),im*jm,pfi(ngrid,1),ngrid)
  
  !   traitement des point normaux
  DO ifield=1,nfield
    DO j=2,jm-1
      ig=2+(j-2)*(im-1)
      CALL SCOPY(im-1,pdyn(1,j,ifield),1,pfi(ig,ifield),1)
    ENDDO
  ENDDO

  RETURN
END SUBROUTINE gr_dyn_fi
SUBROUTINE gr_dyn_fi_ad(nfield,im,jm,ngrid,pdyn,pfi)
  ! All real variables are AD variables
  IMPLICIT NONE
  INTEGER :: im,jm,ngrid,nfield
  REAL :: pdyn(im,jm,nfield)
  REAL :: pfi(ngrid,nfield)
  
  INTEGER :: i,j,ifield,ig
  EXTERNAL :: SCOPY
  
  IF(ngrid/=2+(jm-2)*(im-1)) STOP 'probleme de dim'
  
  pdyn = 0.
  DO ifield=nfield,1,-1
    DO j=jm-1,2,-1
      ig=2+(j-2)*(im-1)
      DO i=im-1,1,-1
        pdyn(i,j,ifield) = pdyn(i,j,ifield) + pfi(ig-1+i,ifield)
        pfi(ig-1+i,ifield) = 0.
      ENDDO
    ENDDO
  ENDDO

  pdyn(1,jm,:) = pdyn(1,jm,:) + pfi(ngrid,:)
  pdyn(1,1,:) = pdyn(1,1,:) + pfi(1,:)
  
  pfi(ngrid,:)=0.
  pfi(1,:)=0.
  RETURN
END SUBROUTINE gr_dyn_fi_ad

SUBROUTINE gr_fi_dyn_ad(nfield,ngrid,im,jm,pfi,pdyn)
  ! All real variables are AD variables
  IMPLICIT NONE
  
  INTEGER :: im,jm,ngrid,nfield
  REAL :: pdyn(im,jm,nfield)
  REAL :: pfi(ngrid,nfield)
  
  INTEGER  :: i,j,ifield,ig
  EXTERNAL :: SCOPY
  
  pfi = 0.
  DO ifield=nfield,1,-1
    DO j=jm-1,2,-1
      ig=2+(j-2)*(im-1)
      pdyn(1,j,ifield) = pdyn(1,j,ifield) + pdyn(im,j,ifield)
      pdyn(im,j,ifield) = 0.
      DO i=im-1,1,-1
        pfi(ig-1+i,ifield) = pfi(ig-1+i,ifield) + pdyn(i,j,ifield)
        pdyn(i,j,ifield) = 0.
      ENDDO
    ENDDO
    DO i=1,im
      pfi(1,ifield) = pfi(1,ifield) + pdyn(i,1,ifield)
      pdyn(i,1,ifield) = 0.
      pfi(ngrid,ifield) = pfi(ngrid,ifield) + pdyn(i,jm,ifield)
      pdyn(i,jm,ifield) = 0.
    ENDDO
  ENDDO
  
  RETURN
END SUBROUTINE gr_fi_dyn_ad

