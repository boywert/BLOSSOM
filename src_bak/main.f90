program test
  use common_vars
  use read_parameters
  use runprocs
  implicit none
  character(len=STRINGLEN) :: filename
  integer :: i
  
  call mpi_initialize
  call get_mpi_rank_nodes

  !filename = 'inputs/Config'
  call getarg(1,filename)
  !print*, iargc(), filename
  call read_config(filename)
  call read_zlist

  do i=1,ZLIST_MAXSIZE
     if(zlist(i) > 0) then
        call massfunction(zlist(i))
        !call gen_los_rr(zlist(i))
        call gen_los_rg(zlist(i))
        !call makecorrelation_rr(zlist(i))
        call makecorrelation_rg(zlist(i))
        !call absorpstats(zlist(i))
     end if
  end do

  call mpi_end_process
end program test

