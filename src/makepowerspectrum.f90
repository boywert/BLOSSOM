#ifdef RR
subroutine makepowerspectrum_rr(z)
#else
subroutine makepowerspectrum_rg(z)
#endif
  use mpi
  use mpitools
  use common_vars
  use conversiontools
  use arraytools
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
  real(kind=8) :: d0,delta_x
  real(kind=8) :: nu_max,nu_min,max_observe,min_observe
  integer :: freq_nbins,x_nbins
  real(kind=8), allocatable :: frequency_value(:),tmp_distance_value(:),tmp_signal(:)
  real(kind=8), allocatable :: x_array(:), y_array(:)
  !Calculate frequency and distance range
  d0 = z_to_d(z)
  max_observe = d0 + convert_length2physical(real(Boxsize*line_length_factor/2.,8),z) 
  nu_min = d_to_nu(max_observe)
  min_observe = d0 - convert_length2physical(real(Boxsize*line_length_factor/2.,8),z) 
  nu_max = d_to_nu(min_observe)
  freq_nbins = ceiling((nu_max-nu_min)/obsfreqresolution)
  
  !Allocate arrays
  allocate(frequency_value(0:freq_nbins))
  allocate(tmp_distance_value(0:freq_nbins))
  allocate(tmp_signal(0:freq_nbins))
  do i=0,freq_nbins
     frequency_value(i) = nu_max - i*obsfreqresolution
     tmp_distance_value(i) = nu_to_d(frequency_value(i))
     tmp_signal(i) = real(i,8)+1.0
  end do
  do i=freq_nbins,1,-1
     tmp_distance_value(i) = (tmp_distance_value(i)-tmp_distance_value(0))/kpc
  end do
  tmp_distance_value(0) = 0.
  
  !Set up new x arrays
  delta_x = tmp_distance_value(1)              !high frequency has higher delta_x
  x_nbins = ceiling(tmp_distance_value(freq_nbins)/delta_x)-1
  allocate(x_array(0:x_nbins))
  allocate(y_array(0:x_nbins))
  do i=0,x_nbins
     x_array(i) = i*delta_x
  end do
  call array_intrpol(tmp_distance_value,tmp_signal,freq_nbins+1,x_array,y_array,x_nbins+1)
  do i=0,x_nbins
     print*,y_array(i)
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
