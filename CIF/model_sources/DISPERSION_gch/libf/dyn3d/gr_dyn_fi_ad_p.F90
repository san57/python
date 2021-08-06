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

SUBROUTINE gr_dyn_fi_ad_p(nfield,im,jm,ngrid,pdyn,pfi)
  ! All real variables are AD variables
  USE mod_interface_dyn_phys
  USE dimphy
  USE parallel
  IMPLICIT NONE
  INTEGER :: im,jm,ngrid,nfield
  REAL :: pdyn(im,jm,nfield)
  REAL :: pfi(ngrid,nfield)

  INTEGER :: i,j,ifield,ig
  
  pdyn(:,:,:)=0.
  DO ifield=1,nfield
    
    DO ig=1,ngrid
      i=index_i(ig)
      j=index_j(ig)
      pdyn(i,j,ifield)=pfi(ig,ifield)
    ENDDO
    
    IF (north_pole) THEN
        pdyn(:,1,ifield)=0
        pdyn(1,1,ifield)=pdyn(1,1,ifield)+pfi(1,ifield)
    ENDIF
    IF (south_pole) THEN
        pdyn(:,jm,ifield)=0
        pdyn(1,jm,ifield)=pdyn(1,jm,ifield)+pfi(ngrid,ifield)
    ENDIF
  END DO
  
  RETURN
END SUBROUTINE gr_dyn_fi_ad_p
