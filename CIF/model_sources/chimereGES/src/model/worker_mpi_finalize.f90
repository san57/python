subroutine worker_mpi_finalize
    use worker_common
    implicit none
    include 'mpif.h'
    integer :: ierr

    call worker_deallocall
    call mpi_barrier(mpi_comm_world, ierr)
    call mpi_finalize(ierr)

end subroutine worker_mpi_finalize


