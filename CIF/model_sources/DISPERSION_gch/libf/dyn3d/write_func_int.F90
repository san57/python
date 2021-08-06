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

  SUBROUTINE write_func_int(wfunc_int,nsec)
  USE parallel
  IMPLICIT NONE
 include "dimensions.h"
 include "paramet.h"
 include 'mpif.h' 
      INTEGER :: nsec
      REAL, DIMENSION(llm,nsec,iip1,jjp1)  :: wfunc_int
      REAL, DIMENSION(llm,nsec,iip1,jjp1)  :: wfunc_int_glo
      INTEGER, DIMENSION(0:mpi_size-1) :: recv_count
      INTEGER, DIMENSION(0:mpi_size-1) :: displ
      INTEGER :: n,ierr
      
      DO n=0,mpi_size-1
        recv_count(n)=llm*nsec*iip1*jj_nb_para(n)
        displ(n)=llm*nsec*iip1*(jj_begin_para(n)-1)
      ENDDO
      
      
      PRINT*, recv_count
      PRINT*,displ
      PRINT*,llm*nsec*iip1*jj_nb
      CALL MPI_GATHERV(wfunc_int(1,1,1,jj_begin),llm*nsec*iip1*jj_nb,MPI_REAL8,   &
                        wfunc_int_glo,recv_count,displ,MPI_REAL8,0,COMM_LMDZ,ierr)

      IF (mpi_rank==0) THEN
         open(99,file='wfunc.bin',form='unformatted')
         write(99) wfunc_int_glo
         close(99)
      ENDIF
        
  END SUBROUTINE write_func_int
