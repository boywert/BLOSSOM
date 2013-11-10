#ifdef RR
subroutine makepowerspectrum_rr(z)
#else
subroutine makepowerspectrum_rg(z)
#endif
  use mpi
  use mpitools
  use common_vars
  use conversiontools
  use omp_lib
  use io_tools
  implicit none
  include 'fftw3.f'
  real(kind=8) :: z
  integer, parameter :: N=1000
  integer(kind=8),allocatable :: plan(:)
  real(kind=8),allocatable :: in(:,:)
  complex(kind=8),allocatable :: out(:,:)
  integer :: i,j,iret
  real(kind=8) :: d0
  real(kind=8) :: nu_max,nu_min,max_observe,min_observe
  integer :: freq_nbins
  real(kind=8), allocatable :: frequency_value(:),tmp_distance_value(:)

  !Calculate frequency and distance range
  d0 = z_to_d(z)
  max_observe = d0 + convert_length2physical(real(Boxsize*line_length_factor/2.,8),z) 
  nu_min = d_to_nu(max_observe)
  min_observe = d0 - convert_length2physical(real(Boxsize*line_length_factor/2.,8),z) 
  nu_max = d_to_nu(min_observe)
  freq_nbins = ceiling((nu_max-nu_min)/maxfreqresolution)
  
  !Allocate arrays
  allocate(frequency_value(0:freq_nbins))
  allocate(tmp_distance_value(0:freq_nbins))
  do i=0,freq_nbins
     frequency_value(i) = nu_max - i*obsfreqresolution
     tmp_distance_value(i) = nu_to_d(frequency_value(i))
  end do
  do i=freq_nbins,1,-1
     tmp_distance_value(i) = (tmp_distance_value(i)-tmp_distance_value(i-1))/Mpc
  end do
  tmp_distance_value(0) = 0.
  do i=0,freq_nbins
     print*,tmp_distance_value(i)
  end do
  ! do j= first_l,last_l
  !    call dfftw_plan_dft_r2c_1d(plan(omp_get_thread_num()),N,in(:,omp_get_thread_num()),out(:,omp_get_thread_num()),FFTW_ESTIMATE)
  !    call dfftw_execute(plan(omp_get_thread_num()))
  !    call dfftw_destroy_plan(plan(omp_get_thread_num()))
  ! end do
     
  deallocate(frequency_value,tmp_distance_value)
  
#ifdef RR
end subroutine makepowerspectrum_rr
#else
end subroutine makepowerspectrum_rg
#endif
