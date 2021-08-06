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

MODULE mod_phys_lmdz_omp_data

  INTEGER,SAVE :: omp_size
  INTEGER,SAVE :: omp_rank
  LOGICAL,SAVE :: is_omp_root
  LOGICAL,SAVE :: is_ok_omp
  
  INTEGER,SAVE,DIMENSION(:),ALLOCATABLE :: klon_omp_para_nb
  INTEGER,SAVE,DIMENSION(:),ALLOCATABLE :: klon_omp_para_begin
  INTEGER,SAVE,DIMENSION(:),ALLOCATABLE :: klon_omp_para_end    
  
  INTEGER,SAVE :: klon_omp
  INTEGER,SAVE :: klon_omp_begin
  INTEGER,SAVE :: klon_omp_end
!$OMP  THREADPRIVATE(omp_rank,klon_omp,is_omp_root,klon_omp_begin,klon_omp_end)

CONTAINS
  
  SUBROUTINE Init_phys_lmdz_omp_data(klon_mpi)
    USE dimphy
    IMPLICIT NONE
    INTEGER, INTENT(in) :: klon_mpi

    INTEGER :: i

#ifdef _OPENMP    
    INTEGER :: OMP_GET_NUM_THREADS
    EXTERNAL OMP_GET_NUM_THREADS
    INTEGER :: OMP_GET_THREAD_NUM
    EXTERNAL OMP_GET_THREAD_NUM
#endif  

#ifdef _OPENMP
!$OMP MASTER
        is_ok_omp=.TRUE.
        omp_size=OMP_GET_NUM_THREADS()
!$OMP END MASTER
        omp_rank=OMP_GET_THREAD_NUM()    
#else    
    is_ok_omp=.FALSE.
    omp_size=1
    omp_rank=0
#endif

   is_omp_root=.FALSE.
!$OMP MASTER
   IF (omp_rank==0) THEN
     is_omp_root=.TRUE.
   ELSE
     PRINT *,'ANORMAL : OMP_MASTER /= 0'
     STOP
   ENDIF
!$OMP END MASTER


!$OMP MASTER 
    ALLOCATE(klon_omp_para_nb(0:omp_size-1))
    ALLOCATE(klon_omp_para_begin(0:omp_size-1))
    ALLOCATE(klon_omp_para_end(0:omp_size-1))
    
    DO i=0,omp_size-1
      klon_omp_para_nb(i)=klon_mpi/omp_size
      IF (i<MOD(klon_mpi,omp_size)) klon_omp_para_nb(i)=klon_omp_para_nb(i)+1
    ENDDO
    
    klon_omp_para_begin(0) = 1
    klon_omp_para_end(0) = klon_omp_para_nb(0)
    
    DO i=1,omp_size-1
      klon_omp_para_begin(i)=klon_omp_para_end(i-1)+1
      klon_omp_para_end(i)=klon_omp_para_begin(i)+klon_omp_para_nb(i)-1
    ENDDO
!$OMP END MASTER
!$OMP BARRIER
   
    klon_omp=klon_omp_para_nb(omp_rank)
    klon_omp_begin=klon_omp_para_begin(omp_rank)
    klon_omp_end=klon_omp_para_end(omp_rank)
    
    CALL Print_module_data
    
  END SUBROUTINE Init_phys_lmdz_omp_data

  SUBROUTINE Print_module_data
  IMPLICIT NONE

!$OMP CRITICAL  
  PRINT *,'--------> TASK ',omp_rank
  PRINT *,'omp_size =',omp_size
  PRINT *,'omp_rank =',omp_rank
  PRINT *,'is_omp_root =',is_omp_root
  PRINT *,'klon_omp_para_nb =',klon_omp_para_nb
  PRINT *,'klon_omp_para_begin =',klon_omp_para_begin
  PRINT *,'klon_omp_para_end =',klon_omp_para_end    
  PRINT *,'klon_omp =',klon_omp
  PRINT *,'klon_omp_begin =',klon_omp_begin
  PRINT *,'klon_omp_end =',klon_omp_end    
!$OMP END CRITICAL

  END SUBROUTINE Print_module_data
END MODULE mod_phys_lmdz_omp_data
