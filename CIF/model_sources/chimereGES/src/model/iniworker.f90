subroutine iniworker

  use chimere_consts
  use worker_common
  use message_defs
  use worker_message_subs

  implicit none

  !****************************************************************************
  integer :: ns, ierr

  call worker_recv_once
  call mpi_comm_rank(wrk_comm,rank,ierr)

  ! IP recpetion obs
  call mpi_recv(nobs,1,mpi_integer,0,ias_nobs,  &
                mpi_comm_world,mpi_status_ignore,ierr)
  print*,'worker numero',rank,' nb obs=',nobs
  allocate(tabobs(nobs,12))
  call mpi_recv(tabobs,nobs*12,mpi_double_precision,0,ias_tabobs,   &
                 mpi_comm_world,mpi_status_ignore,ierr)
!  chckexch=20
!  call mpi_send(chckexch,1,mpi_integer,0,ias_chckexch,mpi_comm_world,ierr)
  ! fin IP
  call mpi_barrier(mpi_comm_world,ierr)

  conc = 1d-17
  ns=0
  do while (ns /= 10000)
    call worker_recv_ns(ns)
    if (ns /= 10000) call worker_recv_conc(ns)
    call mpi_barrier(mpi_comm_world,ierr)
  end do
end subroutine iniworker
