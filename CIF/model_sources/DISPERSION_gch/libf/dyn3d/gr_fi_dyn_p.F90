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

SUBROUTINE gr_fi_dyn_p(nfield,ngrid,im,jm,pfi,pdyn)
  USE mod_interface_dyn_phys
  USE dimphy
  USE parallel
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
  
  INTEGER :: i,j,ifield,ig
  
  !-----------------------------------------------------------------------
  !   calcul:
  !   -------
  !$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
  DO ifield=1,nfield
    
    DO ig=1,klon
      i=index_i(ig)
      j=index_j(ig)
      pdyn(i,j,ifield)=pfi(ig,ifield)
      IF (i==1) pdyn(im,j,ifield)=pdyn(i,j,ifield)
    ENDDO
    
    !   traitement des poles
    IF (north_pole) THEN
        DO i=1,im
          pdyn(i,1,ifield)=pdyn(1,1,ifield)
        ENDDO
    ENDIF
    
    IF (south_pole) THEN
        DO i=1,im
          pdyn(i,jm,ifield)=pdyn(1,jm,ifield)
        ENDDO
    ENDIF
    
  ENDDO
  !$OMP END DO NOWAIT
  RETURN
END SUBROUTINE gr_fi_dyn_p
