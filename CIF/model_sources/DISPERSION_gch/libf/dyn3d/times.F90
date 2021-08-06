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

module times
  integer,private,save :: Last_Count=0
  real, private,save :: Last_cpuCount=0
  logical, private,save :: AllTimer_IsActive=.false.
  
  integer, parameter :: nb_timer = 4
  integer, parameter :: timer_caldyn  = 1
  integer, parameter :: timer_vanleer = 2
  integer, parameter :: timer_dissip = 3
  integer, parameter :: timer_physic = 4
  integer, parameter :: stopped = 1
  integer, parameter :: running = 2
  integer, parameter :: suspended = 3 
  
  integer :: max_size
  real,    allocatable, dimension(:,:,:) :: timer_table
  real,    allocatable, dimension(:,:,:) :: timer_table_sqr 
  integer, allocatable, dimension(:,:,:) :: timer_iteration
  real,    allocatable, dimension(:,:,:) :: timer_average
  real,    allocatable, dimension(:,:,:) :: timer_delta
  real,    allocatable,dimension(:) :: timer_running, last_time
  integer, allocatable,dimension(:) :: timer_state
  
  contains
  
  subroutine init_timer
    use parallel
    implicit none
    include "dimensions.h"
    include "paramet.h"
    
    max_size=jjm+1
    allocate(timer_table(max_size,nb_timer,0:mpi_size-1))
    allocate(timer_table_sqr(max_size,nb_timer,0:mpi_size-1))
    allocate(timer_iteration(max_size,nb_timer,0:mpi_size-1))
    allocate(timer_average(max_size,nb_timer,0:mpi_size-1))
    allocate(timer_delta(max_size,nb_timer,0:mpi_size-1))
    allocate(timer_running(nb_timer))
    allocate(timer_state(nb_timer))
    allocate(last_time(nb_timer))
    
    timer_table(:,:,:)=0
    timer_table_sqr(:,:,:)=0
    timer_iteration(:,:,:)=0
    timer_average(:,:,:)=0
    timer_delta(:,:,:)=0
    timer_state(:)=stopped      
  end subroutine init_timer
  
  subroutine start_timer(no_timer)
    implicit none
    integer :: no_timer
    
    if (AllTimer_IsActive) then
    
      if (timer_state(no_timer)/=stopped) then
        stop 'start_timer :: timer is already running or suspended'
      else
        timer_state(no_timer)=running
      endif
      
      timer_running(no_timer)=0
      call cpu_time(last_time(no_timer))
    
    endif
    
  end subroutine start_timer
  
  subroutine suspend_timer(no_timer)
    implicit none
    integer :: no_timer
     
    if (AllTimer_IsActive) then   
      if (timer_state(no_timer)/=running) then
        stop 'suspend_timer :: timer is not running'
      else
        timer_state(no_timer)=suspended
      endif
    
      timer_running(no_timer)=timer_running(no_timer)-last_time(no_timer)
      call cpu_time(last_time(no_timer))
      timer_running(no_timer)=timer_running(no_timer)+last_time(no_timer)
    endif
  end subroutine suspend_timer
  
  subroutine resume_timer(no_timer)
    implicit none
    integer :: no_timer
     
    if (AllTimer_IsActive) then   
      if (timer_state(no_timer)/=suspended) then
        stop 'resume_timer :: timer is not suspended'
      else
        timer_state(no_timer)=running
      endif
      
      call cpu_time(last_time(no_timer))
    endif
    
  end subroutine resume_timer

  subroutine stop_timer(no_timer)
    use parallel
    implicit none
    integer :: no_timer
    integer :: N
    real :: V,V2
    
    if (AllTimer_IsActive) then
       
      if (timer_state(no_timer)/=running) then
        stop 'stop_timer :: timer is not running'
      else
        timer_state(no_timer)=stopped
      endif
    
      timer_running(no_timer)=timer_running(no_timer)-last_time(no_timer)
      call cpu_time(last_time(no_timer))
      timer_running(no_timer)=timer_running(no_timer)+last_time(no_timer)
    
      timer_table(jj_nb,no_timer,mpi_rank)=timer_table(jj_nb,no_timer,mpi_rank)+timer_running(no_timer)
      timer_table_sqr(jj_nb,no_timer,mpi_rank)=timer_table_sqr(jj_nb,no_timer,mpi_rank)+timer_running(no_timer)**2
      timer_iteration(jj_nb,no_timer,mpi_rank)=timer_iteration(jj_nb,no_timer,mpi_rank)+1
      timer_average(jj_nb,no_timer,mpi_rank)=timer_table(jj_nb,no_timer,mpi_rank)/timer_iteration(jj_nb,no_timer,mpi_rank)
      if (timer_iteration(jj_nb,no_timer,mpi_rank)>=2) then
        N=timer_iteration(jj_nb,no_timer,mpi_rank)
	V2=timer_table_sqr(jj_nb,no_timer,mpi_rank)
	V=timer_table(jj_nb,no_timer,mpi_rank)
	timer_delta(jj_nb,no_timer,mpi_rank)=sqrt(ABS(V2-V*V/N)/(N-1)) 
      else
        timer_delta(jj_nb,no_timer,mpi_rank)=0
      endif
    endif
    
  end subroutine stop_timer
   
  subroutine allgather_timer
    use parallel
    implicit none
    include 'mpif.h'
    integer :: ierr
    integer :: data_size
    real, allocatable,dimension(:,:) :: tmp_table
    
    if (AllTimer_IsActive) then
    
    allocate(tmp_table(max_size,nb_timer))
    
    data_size=max_size*nb_timer
    
    tmp_table(:,:)=timer_table(:,:,mpi_rank)
    call mpi_allgather(tmp_table(1,1),data_size,MPI_REAL8,timer_table(1,1,mpi_rank),data_size,MPI_REAL8,COMM_LMDZ,ierr)

    tmp_table(:,:)=timer_table_sqr(:,:,mpi_rank)
    call mpi_allgather(tmp_table(1,1),data_size,MPI_REAL8,timer_table_sqr(1,1,mpi_rank),data_size,MPI_REAL8,COMM_LMDZ,ierr)
    
    deallocate(tmp_table)
    
    endif
    
  end subroutine allgather_timer
  
  subroutine allgather_timer_average
    use parallel
    implicit none
    include 'mpif.h'
    integer :: ierr
    integer :: data_size
    real, allocatable,dimension(:,:),target :: tmp_table
    integer, allocatable,dimension(:,:),target :: tmp_iter
    
    if (AllTimer_IsActive) then
    
    allocate(tmp_table(max_size,nb_timer))
    allocate(tmp_iter(max_size,nb_timer))
   
    data_size=max_size*nb_timer

    tmp_table(:,:)=timer_average(:,:,mpi_rank)
    call mpi_allgather(tmp_table(1,1),data_size,MPI_REAL8,timer_average(1,1,0),data_size,MPI_REAL8,COMM_LMDZ,ierr)

    tmp_table(:,:)=timer_delta(:,:,mpi_rank)
    call mpi_allgather(tmp_table(1,1),data_size,MPI_REAL8,timer_delta(1,1,0),data_size,MPI_REAL8,COMM_LMDZ,ierr)

    tmp_iter(:,:)=timer_iteration(:,:,mpi_rank)
    call mpi_allgather(tmp_iter(1,1),data_size,MPI_INTEGER,timer_iteration(1,1,0),data_size,MPI_INTEGER,COMM_LMDZ,ierr)
    
    deallocate(tmp_table)
    
    endif
  end subroutine allgather_timer_average
  
  subroutine InitTime
  implicit none
    integer :: count,count_rate,count_max
    
    AllTimer_IsActive=.TRUE.
    if (AllTimer_IsActive) then
      call system_clock(count,count_rate,count_max)
      call cpu_time(Last_cpuCount)
      Last_Count=count
    endif
  end subroutine InitTime
  
  function DiffTime
  implicit none
    double precision :: DiffTime
    integer :: count,count_rate,count_max
  
    call system_clock(count,count_rate,count_max)
    if (Count>=Last_Count) then
      DiffTime=(1.*(Count-last_Count))/count_rate
    else
      DiffTime=(1.*(Count-last_Count+Count_max))/count_rate
    endif
    Last_Count=Count 
  end function DiffTime
  
  function DiffCpuTime
  implicit none
    real :: DiffCpuTime
    real :: Count
    
    call cpu_time(Count)
    DiffCpuTime=Count-Last_cpuCount
    Last_cpuCount=Count 
  end function DiffCpuTime

end module times
