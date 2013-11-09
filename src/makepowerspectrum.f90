#ifdef RR
subroutine makepowerspectrum_rr(z)
#else
subroutine makepowerspectrum_rg(z)
#endif
  use mpi
  use mpitools
  use common_vars
  use omp_lib
  use io_tools
  implicit none
  include 'fftw3.f'
  real(kind=8) :: z
  integer, parameter :: N=1000
  integer(kind=8) :: plan(0:2)
  real(kind=8) :: in(N,0:2)
  complex(kind=8) :: out(N/2+1,0:2)
  integer :: i,j,iret

  do i=1,N
     in(i,0:2) = real(i)**2.
  end do
  print*,""
  do i=0,100
     call dfftw_plan_dft_r2c_1d(plan(omp_get_thread_num()),N,in(:,omp_get_thread_num()),out(:,omp_get_thread_num()),FFTW_ESTIMATE)
     call dfftw_execute(plan(omp_get_thread_num()))
     call dfftw_destroy_plan(plan(omp_get_thread_num()))
  end do
  ! do i=1,N/2+1
  !    print*, real(in(i)), real(out(i))
  ! end do


  
#ifdef RR
end subroutine makepowerspectrum_rr
#else
end subroutine makepowerspectrum_rg
#endif
