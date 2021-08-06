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

  SUBROUTINE gather_tabobs(tabobs,nobs,nobs_para,nobs_glo,tabobs_glo)
    USE parallel
    IMPLICIT NONE
    INCLUDE 'mpif.h'
    
    INTEGER :: nobs,nobs_glo
    INTEGER, DIMENSION(0:mpi_size-1) :: nobs_para
    REAL,DIMENSION(nobs,10) :: tabobs
    REAL,DIMENSION(10,nobs) :: tmp
    REAL,DIMENSION(nobs_glo,10) :: tabobs_glo
    REAL,DIMENSION(10,nobs_glo) :: tabobs_tmp
    
    INTEGER, DIMENSION(0:mpi_size-1) :: recv_count
    INTEGER, DIMENSION(0:mpi_size-1) :: displ
    INTEGER :: n,i,j, ierr
    
    ! Definir arguments pour mpi_gatherv
    DO n=1,mpi_size-1
      recv_count(n)=10*nobs_para(n)
      displ(n)=SUM(10*nobs_para(0:n-1))
    ENDDO
    recv_count(0)=10*nobs_para(0)
    displ(0)=0
    
    tabobs_tmp(:,:)=0.
    
    tmp=transpose(tabobs)
    CALL mpi_gatherv(tmp,10*nobs,mpi_real8,tabobs_tmp,recv_count, &
       displ,mpi_real8,0,comm_lmdz,ierr)
    
    ! Remettre les donnees dans le bon ordre depuis tabobs_tmp a tabobs_glo sur noeud de rank 0
    IF (mpi_rank == 0) THEN
        j=0
        DO n=0,mpi_size-1 
          DO i=1,nobs_glo
            IF (tabobs_glo(i,3)>=jj_begin_para(n) .AND. tabobs_glo(i,3)<=jj_end_para(n)) THEN
                j=j+1
                tabobs_glo(i,:)=tabobs_tmp(:,j)
            ENDIF
          ENDDO
        ENDDO
    ENDIF
    
  END SUBROUTINE gather_tabobs

