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
  integer(kind=8) :: plan
  real(kind=8) :: in(N)
  complex(kind=8) :: out(N/2+1)
  integer :: i,j,iret
  call dfftw_init_threads(iret)
  call dfftw_plan_with_nthreads(omp_get_max_threads())
  do i=1,N
     in(i) = real(i)**2.
  end do
  print*,""
  call dfftw_plan_dft_r2c_1d(plan,N,in,out,FFTW_ESTIMATE)
  call dfftw_execute(plan)
  call dfftw_destroy_plan(plan)
  do i=1,N/2+1
     print*, real(in(i)), real(out(i))
  end do


  
#ifdef RR
end subroutine makepowerspectrum_rr
#else
end subroutine makepowerspectrum_rg
#endif
