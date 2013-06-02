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
#ifdef SUBCELL
        call call_subcell_corr(zlist(i))
#endif
#ifdef MASSFUNCTION        
        call massfunction(zlist(i))
#endif
#ifdef GENLOS_RR
        call gen_los_rr(zlist(i))
#endif
#ifdef GENLOS_RG
        call gen_los_rg(zlist(i))
#endif
#ifdef CALCORR_RR
        call makecorrelation_rr(zlist(i))
#endif
#ifdef CALCORR_RG
        call makecorrelation_rg(zlist(i))
#endif
#ifdef ABSORPTION
        call makeobserveline(zlist(i))
#endif
     end if
  end do

  call mpi_end_process
end program test

