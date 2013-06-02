module mpitools
  use mpi
  integer :: ierr, rank, nodes_returned
contains
  subroutine get_mpi_rank_nodes
    implicit none
    call mpi_comm_size(mpi_comm_world,nodes_returned,ierr)
    if (ierr /= mpi_success) call mpi_abort(mpi_comm_world,ierr,ierr)

    call mpi_comm_rank(mpi_comm_world,rank,ierr)
    if (ierr /= mpi_success) call mpi_abort(mpi_comm_world,ierr,ierr)
    return
  end subroutine get_mpi_rank_nodes
  
  subroutine mpi_initialize
    implicit none
    call mpi_init(ierr)
    if (ierr /= mpi_success) call mpi_abort(mpi_comm_world,ierr,ierr)
  end subroutine mpi_initialize
 
  subroutine mpi_end_process
    implicit none
    call mpi_finalize(ierr)
  end subroutine mpi_end_process
end module mpitools
