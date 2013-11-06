program test
  use common_vars
  use read_parameters
  use runprocs
  use omp_lib
  implicit none
  character(len=STRINGLEN) :: filename
  integer :: i
  real(kind=8) :: main_start,main_stop
  call mpi_initialize
  call get_mpi_rank_nodes
  if(rank==0) main_start = omp_get_wtime() 
  !filename = 'inputs/Config'
  call getarg(1,filename)
  !print*, iargc(), filename
  call read_config(filename)
  call read_zlist

  do i=1,ZLIST_MAXSIZE
     if(zlist(i) > 0) then
#ifdef GENRANDOM
	call gen_random(zlist(i))
#endif
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
#ifdef ABSORPTION_RR
        call makeobservedlines_rr(zlist(i))
#endif
#ifdef ABSORPTION_RG
        call makeobservedlines_rg(zlist(i))
#endif

     end if
  end do
  if(rank==0) then
     main_stop = omp_get_wtime() 
     print*,"#####################################################"
     print*,"This process used", main_stop-main_start,"s"
     print*,"This process used", nodes_returned*omp_get_max_threads()*(main_stop-main_start)/3600.,"SUs"
     print*, "Terminate process"
     print*,"#####################################################"
     
  end if
  call mpi_end_process
end program test

