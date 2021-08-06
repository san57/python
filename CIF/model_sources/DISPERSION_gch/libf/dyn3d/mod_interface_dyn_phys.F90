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

MODULE mod_interface_dyn_phys
  INTEGER,SAVE,dimension(:),allocatable :: index_i
  INTEGER,SAVE,dimension(:),allocatable :: index_j
  
  
CONTAINS
  
  SUBROUTINE Init_interface_dyn_phys
    USE mod_phys_lmdz_mpi_data
    IMPLICIT NONE
    include 'dimensions.h'    
    
    INTEGER :: i,j,k
    
    ALLOCATE(index_i(klon_mpi))
    ALLOCATE(index_j(klon_mpi))
    
    k=1
    IF (is_north_pole) THEN
      index_i(k)=1
      index_j(k)=1
      k=2
    ELSE
      DO i=ii_begin,iim
	index_i(k)=i
	index_j(k)=jj_begin
	k=k+1
       ENDDO
    ENDIF
    
    DO j=jj_begin+1,jj_end-1
      DO i=1,iim
	index_i(k)=i
	index_j(k)=j
	k=k+1
      ENDDO
    ENDDO
    
    IF (is_south_pole) THEN
      index_i(k)=1
      index_j(k)=jj_end
    ELSE
      DO i=1,ii_end
	index_i(k)=i
	index_j(k)=jj_end
	k=k+1
       ENDDO
    ENDIF
  
  END SUBROUTINE Init_interface_dyn_phys 

END MODULE mod_interface_dyn_phys
