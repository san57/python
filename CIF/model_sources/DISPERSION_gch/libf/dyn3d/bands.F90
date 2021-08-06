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

MODULE Bands
  
  INTEGER, PARAMETER :: bands_caldyn=1
  INTEGER, PARAMETER :: bands_vanleer=2
  INTEGER, PARAMETER :: bands_dissip=3
  
  INTEGER,DIMENSION(:),ALLOCATABLE :: jj_Nb_Caldyn
  INTEGER,DIMENSION(:),ALLOCATABLE :: jj_Nb_vanleer
  INTEGER,DIMENSION(:),ALLOCATABLE :: jj_Nb_vanleer2
  INTEGER,DIMENSION(:),ALLOCATABLE :: jj_Nb_dissip
  INTEGER,DIMENSION(:),ALLOCATABLE :: jj_Nb_physic
  INTEGER,DIMENSION(:),ALLOCATABLE :: jj_Nb_physic_bis
  INTEGER,DIMENSION(:),ALLOCATABLE :: distrib_phys
  
CONTAINS
  
  SUBROUTINE AllocateBands
    USE parallel
    IMPLICIT NONE
    
    ALLOCATE(jj_Nb_Caldyn(0:MPI_Size-1))
    ALLOCATE(jj_Nb_vanleer(0:MPI_Size-1))
    ALLOCATE(jj_Nb_vanleer2(0:MPI_Size-1))
    ALLOCATE(jj_Nb_dissip(0:MPI_Size-1))
    ALLOCATE(jj_Nb_physic(0:MPI_Size-1))
    ALLOCATE(jj_Nb_physic_bis(0:MPI_Size-1))
    ALLOCATE(distrib_phys(0:MPI_Size-1))
    
  END SUBROUTINE AllocateBands
  
  SUBROUTINE Read_distrib
    USE parallel
    IMPLICIT NONE
    
    INCLUDE "dimensions.h"
    INTEGER :: i,j
    CHARACTER (len=4) :: siim,sjjm,sllm,sproc
    CHARACTER (len=255) :: filename
    INTEGER :: unit_number=10
    INTEGER :: ierr
    
    CALL AllocateBands
    WRITE(siim,'(i3)') iim
    WRITE(sjjm,'(i3)') jjm
    WRITE(sllm,'(i3)') llm
    WRITE(sproc,'(i3)') mpi_size
    filename='Bands_'//TRIM(ADJUSTL(siim))//'x'//TRIM(ADJUSTL(sjjm))//'x'//TRIM(ADJUSTL(sllm))//'_'  &
       //TRIM(ADJUSTL(sproc))//'prc.dat'    
    
    OPEN(UNIT=unit_number,FILE=TRIM(filename),STATUS='old',FORM='formatted',IOSTAT=ierr)
    
    IF (ierr==0) THEN
        
        DO i=0,mpi_size-1
          READ (unit_number,*) j,jj_nb_caldyn(i)
        ENDDO
        
        DO i=0,mpi_size-1
          READ (unit_number,*) j,jj_nb_vanleer(i)
        ENDDO
        
        DO i=0,mpi_size-1
          READ (unit_number,*) j,jj_nb_dissip(i)
        ENDDO
        
        DO i=0,mpi_size-1
          READ (unit_number,*) j,distrib_phys(i)
        ENDDO
	
	CLOSE(unit_number)  
        
    ELSE
        DO i=0,mpi_size-1
          jj_nb_caldyn(i)=(jjm+1)/mpi_size
	  IF (i<MOD(jjm+1,mpi_size)) jj_nb_caldyn(i)=jj_nb_caldyn(i)+1
        ENDDO
        
        jj_nb_vanleer(:)=jj_nb_caldyn(:)
        jj_nb_dissip(:)=jj_nb_caldyn(:)
        
        !	do i=0,mpi_size-1
        !	  distrib_phys(i)=(iim*(jjm-1)+2)/mpi_size
        !	  IF (i<MOD(iim*(jjm-1)+2,mpi_size)) distrib_phys(i)=distrib_phys(i)+1
        !	enddo
	DO i=0,mpi_size-1
	  distrib_phys(i)=jj_nb_caldyn(i)*iim
	ENDDO
        
        distrib_phys(0)=distrib_phys(0)-iim+1
        distrib_phys(mpi_size-1)=distrib_phys(mpi_size-1)-iim+1
        
    ENDIF
    
  END SUBROUTINE Read_distrib
  
  
  SUBROUTINE  Set_Bands 
    USE parallel
    USE mod_phys_lmdz_para, ONLY : jj_para_begin,jj_para_end
    IMPLICIT NONE
    INCLUDE 'dimensions.h'    
    INTEGER :: i
    
    DO i=0,mpi_size-1
      jj_nb_vanleer2(i)=(jjm+1)/mpi_size
      IF (i<MOD(jjm+1,mpi_size)) jj_nb_vanleer2(i)=jj_nb_vanleer2(i)+1
    ENDDO
    
    DO i=0,MPI_Size-1
      jj_Nb_physic(i)=jj_para_end(i)-jj_para_begin(i)+1
      IF (i/=0) THEN
          IF (jj_para_begin(i)==jj_para_end(i-1)) THEN
              jj_Nb_physic(i-1)=jj_Nb_physic(i-1)-1
          ENDIF
      ENDIF
    ENDDO
    
    DO i=0,MPI_Size-1
      jj_Nb_physic_bis(i)=jj_para_end(i)-jj_para_begin(i)+1
      IF (i/=0) THEN
          IF (jj_para_begin(i)==jj_para_end(i-1)) THEN
              jj_Nb_physic_bis(i)=jj_Nb_physic_bis(i)-1
          ELSE
              jj_Nb_physic_bis(i-1)=jj_Nb_physic_bis(i-1)+1
              jj_Nb_physic_bis(i)=jj_Nb_physic_bis(i)-1
	  ENDIF
      ENDIF
    ENDDO
    
  END SUBROUTINE Set_Bands
  
  
  SUBROUTINE AdjustBands_caldyn
    USE times
    USE parallel
    IMPLICIT NONE
    
    REAL :: minvalue,maxvalue
    INTEGER :: min_proc,max_proc
    INTEGER :: i,j
    REAL,ALLOCATABLE,DIMENSION(:) :: value
    INTEGER,ALLOCATABLE,DIMENSION(:) :: index
    REAL :: tmpvalue
    INTEGER :: tmpindex
    
    ALLOCATE(value(0:mpi_size-1))
    ALLOCATE(INDEX(0:mpi_size-1))
    
    
    CALL allgather_timer_average
    
    DO i=0,mpi_size-1
      value(i)=timer_average(jj_nb_caldyn(i),timer_caldyn,i)
      INDEX(i)=i
    ENDDO
    
    DO i=0,mpi_size-2
      DO j=i+1,mpi_size-1
        IF (value(i)>value(j)) THEN
	    tmpvalue=value(i)
	    value(i)=value(j)
	    value(j)=tmpvalue
	    
	    tmpindex=INDEX(i)
	    INDEX(i)=INDEX(j)
	    INDEX(j)=tmpindex
        ENDIF
      ENDDO
    ENDDO
    
    maxvalue=value(mpi_size-1)
    max_proc=INDEX(mpi_size-1)           
    
    DO i=0,mpi_size-2
      minvalue=value(i)
      min_proc=INDEX(i)
      IF (jj_nb_caldyn(max_proc)>3) THEN
          IF (timer_iteration(jj_nb_caldyn(min_proc)+1,timer_caldyn,min_proc)<=1 ) THEN
              jj_nb_caldyn(min_proc)=jj_nb_caldyn(min_proc)+1
              jj_nb_caldyn(max_proc)=jj_nb_caldyn(max_proc)-1
              EXIT
          ELSE
              IF (timer_average(jj_nb_caldyn(min_proc)+1,timer_caldyn,min_proc)                 &
                 -timer_delta(jj_nb_caldyn(min_proc)+1,timer_caldyn,min_proc) < maxvalue) THEN
                  jj_nb_caldyn(min_proc)=jj_nb_caldyn(min_proc)+1
                  jj_nb_caldyn(max_proc)=jj_nb_caldyn(max_proc)-1
                  EXIT
              ENDIF
          ENDIF
      ENDIF
    ENDDO
    
    DEALLOCATE(value)
    DEALLOCATE(index)
    
  END SUBROUTINE AdjustBands_caldyn
  
  SUBROUTINE AdjustBands_vanleer
    USE times
    USE parallel
    IMPLICIT NONE
    
    REAL :: minvalue,maxvalue
    INTEGER :: min_proc,max_proc
    INTEGER :: i,j
    REAL,ALLOCATABLE,DIMENSION(:) :: value
    INTEGER,ALLOCATABLE,DIMENSION(:) :: index
    REAL :: tmpvalue
    INTEGER :: tmpindex
    
    ALLOCATE(value(0:mpi_size-1))
    ALLOCATE(INDEX(0:mpi_size-1))
    
    
    CALL allgather_timer_average
    
    DO i=0,mpi_size-1
      value(i)=timer_average(jj_nb_vanleer(i),timer_vanleer,i)
      INDEX(i)=i
    ENDDO
    
    DO i=0,mpi_size-2
      DO j=i+1,mpi_size-1
        IF (value(i)>value(j)) THEN
	    tmpvalue=value(i)
	    value(i)=value(j)
	    value(j)=tmpvalue
	    
	    tmpindex=INDEX(i)
	    INDEX(i)=INDEX(j)
	    INDEX(j)=tmpindex
        ENDIF
      ENDDO
    ENDDO
    
    maxvalue=value(mpi_size-1)
    max_proc=INDEX(mpi_size-1)           
    
    DO i=0,mpi_size-2
      minvalue=value(i)
      min_proc=INDEX(i)
      
      IF (jj_nb_vanleer(max_proc)>3) THEN
          IF (timer_average(jj_nb_vanleer(min_proc)+1,timer_vanleer,min_proc)==0. .OR. &
             timer_average(jj_nb_vanleer(max_proc)-1,timer_vanleer,max_proc)==0.) THEN
              jj_nb_vanleer(min_proc)=jj_nb_vanleer(min_proc)+1
              jj_nb_vanleer(max_proc)=jj_nb_vanleer(max_proc)-1
              EXIT
          ELSE
              IF (timer_average(jj_nb_vanleer(min_proc)+1,timer_vanleer,min_proc) < maxvalue) THEN
                  jj_nb_vanleer(min_proc)=jj_nb_vanleer(min_proc)+1
                  jj_nb_vanleer(max_proc)=jj_nb_vanleer(max_proc)-1
                  EXIT
              ENDIF
          ENDIF
      ENDIF
    ENDDO
    
    DEALLOCATE(value)
    DEALLOCATE(index)
    
  END SUBROUTINE AdjustBands_vanleer
  
  SUBROUTINE AdjustBands_dissip
    USE times
    USE parallel
    IMPLICIT NONE
    
    REAL :: minvalue,maxvalue
    INTEGER :: min_proc,max_proc
    INTEGER :: i,j
    REAL,ALLOCATABLE,DIMENSION(:) :: value
    INTEGER,ALLOCATABLE,DIMENSION(:) :: index
    REAL :: tmpvalue
    INTEGER :: tmpindex
    
    ALLOCATE(value(0:mpi_size-1))
    ALLOCATE(INDEX(0:mpi_size-1))
    
    
    CALL allgather_timer_average
    
    DO i=0,mpi_size-1
      value(i)=timer_average(jj_nb_dissip(i),timer_dissip,i)
      INDEX(i)=i
    ENDDO
    
    DO i=0,mpi_size-2
      DO j=i+1,mpi_size-1
        IF (value(i)>value(j)) THEN
	    tmpvalue=value(i)
	    value(i)=value(j)
	    value(j)=tmpvalue
	    
	    tmpindex=INDEX(i)
	    INDEX(i)=INDEX(j)
	    INDEX(j)=tmpindex
        ENDIF
      ENDDO
    ENDDO
    
    maxvalue=value(mpi_size-1)
    max_proc=INDEX(mpi_size-1)           
    
    DO i=0,mpi_size-2
      minvalue=value(i)
      min_proc=INDEX(i)
      
      IF (jj_nb_dissip(max_proc)>3) THEN
          IF (timer_iteration(jj_nb_dissip(min_proc)+1,timer_dissip,min_proc)<=1) THEN
              jj_nb_dissip(min_proc)=jj_nb_dissip(min_proc)+1
              jj_nb_dissip(max_proc)=jj_nb_dissip(max_proc)-1
              EXIT
          ELSE
              IF (timer_average(jj_nb_dissip(min_proc)+1,timer_dissip,min_proc)         &
                 - timer_delta(jj_nb_dissip(min_proc)+1,timer_dissip,min_proc) < maxvalue) THEN
                  jj_nb_dissip(min_proc)=jj_nb_dissip(min_proc)+1
                  jj_nb_dissip(max_proc)=jj_nb_dissip(max_proc)-1
                  EXIT
              ENDIF
          ENDIF
      ENDIF
    ENDDO
    
    DEALLOCATE(value)
    DEALLOCATE(index)
    
  END SUBROUTINE AdjustBands_dissip
  
  SUBROUTINE AdjustBands_physic
    USE times
    USE mod_phys_lmdz_para, ONLY : klon_mpi_para_nb
    USE parallel
    IMPLICIT NONE
    
    INTEGER :: i,Index
    REAL,ALLOCATABLE,DIMENSION(:) :: value
    INTEGER,ALLOCATABLE,DIMENSION(:) :: Inc
    REAL :: medium
    INTEGER :: NbTot,sgn
    
    ALLOCATE(value(0:mpi_size-1))
    ALLOCATE(Inc(0:mpi_size-1))
    
  
    CALL allgather_timer_average
    
    medium=0
    DO i=0,mpi_size-1
      value(i)=timer_average(jj_nb_physic(i),timer_physic,i)
      medium=medium+value(i)
    ENDDO
    
    medium=medium/mpi_size      
    NbTot=0
    DO i=0,mpi_size-1
      Inc(i)=NINT(klon_mpi_para_nb(i)*(medium-value(i))/value(i))
      NbTot=NbTot+Inc(i)  
    ENDDO
    
    IF (NbTot>=0) THEN
        Sgn=1
    ELSE
        Sgn=-1
	NbTot=-NbTot
    ENDIF
    
    Index=0
    DO i=1,NbTot
      Inc(Index)=Inc(Index)-Sgn
      Index=Index+1
      IF (Index>mpi_size-1) Index=0
    ENDDO
    
    DO i=0,mpi_size-1
      distrib_phys(i)=klon_mpi_para_nb(i)+inc(i)
    ENDDO
    
  END SUBROUTINE AdjustBands_physic
  
  SUBROUTINE WriteBands
    USE parallel
    IMPLICIT NONE
    INCLUDE "dimensions.h"
    
    INTEGER :: i
    CHARACTER (len=4) :: siim,sjjm,sllm,sproc
    CHARACTER (len=255) :: filename
    INTEGER :: unit_number=10
    INTEGER :: ierr
    
    WRITE(siim,'(i3)') iim
    WRITE(sjjm,'(i3)') jjm
    WRITE(sllm,'(i3)') llm
    WRITE(sproc,'(i3)') mpi_size
    
    filename='Bands_'//TRIM(ADJUSTL(siim))//'x'//TRIM(ADJUSTL(sjjm))//'x'//TRIM(ADJUSTL(sllm))//'_'  &
       //TRIM(ADJUSTL(sproc))//'prc.dat'    
    
    OPEN(UNIT=unit_number,FILE=TRIM(filename),STATUS='replace',FORM='formatted',IOSTAT=ierr)
    
    IF (ierr==0) THEN
        
        !	write (unit_number,*) '*** Bandes caldyn ***'
	DO i=0,mpi_size-1
          WRITE (unit_number,*) i,jj_nb_caldyn(i)
        ENDDO
        
        !	write (unit_number,*) '*** Bandes vanleer ***' 
        DO i=0,mpi_size-1
          WRITE (unit_number,*) i,jj_nb_vanleer(i)
        ENDDO
        
        !        write (unit_number,*) '*** Bandes dissip ***'
        DO i=0,mpi_size-1
          WRITE (unit_number,*) i,jj_nb_dissip(i)
        ENDDO
        
	DO i=0,mpi_size-1
          WRITE (unit_number,*) i,distrib_phys(i)
        ENDDO
	
        CLOSE(unit_number)   
    ELSE 
        PRINT *,'probleme lors de l ecriture des bandes'
    ENDIF
    
  END SUBROUTINE WriteBands
  
END MODULE Bands

  
