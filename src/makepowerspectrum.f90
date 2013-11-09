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
  include "fftw3.f"

  double complex in, out
  dimension in(N), out(N)
  integer*8 plan

  call dfftw_plan_dft_1d(plan,N,in,out,FFTW_FORWARD,FFTW_ESTIMATE)
  call dfftw_execute_dft(plan, in, out)
  call dfftw_destroy_plan(plan)

  
#ifdef RR
end subroutine makepowerspectrum_rr
#else
end subroutine makepowerspectrum_rg
#endif
