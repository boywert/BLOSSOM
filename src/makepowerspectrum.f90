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
  integer(kind=8) :: plan
  integer :: i,j
  real(kind=8) :: d0,delta_x
  real(kind=8) :: nu_max,nu_min,max_observe,min_observe
  integer :: freq_nbins,x_nbins
  real(kind=8), allocatable :: frequency_value(:),tmp_distance_value(:),tmp_signal(:)
  real(kind=8), allocatable :: x_array(:), y_array(:)
  character(len=100) :: str_rank,z_s,str_line
  real(kind=8) :: M0,impact_param,nu_dist,nu_undist,this_absorp,delta_nu,width_real
  complex, allocatable :: fft_result(:)
  real(kind=8), allocatable :: sum_delta_sq(:)
  ! Prepare strings
  write(z_s,'(f10.3)') z
  z_s = adjustl(z_s)
  write(str_rank,'(i10)') rank
  str_rank = adjustl(str_rank)

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
  end do
  do i=freq_nbins,1,-1
     tmp_distance_value(i) = (tmp_distance_value(i)-tmp_distance_value(0))/kpc
  end do
  tmp_distance_value(0) = 0.
  
  !Set up new x arrays
  delta_x = tmp_distance_value(1)              !high frequency has higher delta_x
  x_nbins = ceiling(tmp_distance_value(freq_nbins)/delta_x)-1
  allocate(x_array(1:x_nbins+1))
  allocate(y_array(1:x_nbins+1))
  allocate(fft_result(1:(x_nbins+1)/2+1))
  allocate(sum_delta_sq(1:(x_nbins+1)/2+1))
  do i=1,x_nbins+1
     x_array(i) = (i)*delta_x
  end do
  sum_delta_sq(:) = 0.0
  do j= first_l,last_l
     write(str_line,'(i10)') j
     str_line = adjustl(str_line)

     tmp_signal(:) = 1.0
     open (unit=10, &
          file=trim(result_path)//z_s(1:len_trim(z_s))//'/RG/0.000/'//trim(adjustl(str_rank))//'/sout.'//trim(adjustl(str_line)), &
          form='binary')
     do
        read(10,end=327) M0,impact_param,nu_dist,nu_undist,this_absorp,delta_nu
        width_real = delta_nu/nu0*nu_dist
        tmp_signal = tmp_signal * (1.-this_absorp* exp(-0.5*((nu_dist-frequency_value)/width_real)**2.0))
     end do
327  close(10)

     call array_intrpol(tmp_distance_value,tmp_signal,freq_nbins+1,x_array,y_array,x_nbins+1)

     call dfftw_plan_dft_r2c_1d(plan,x_nbins+1,y_array,fft_result,FFTW_ESTIMATE)
     call dfftw_execute(plan)
     call dfftw_destroy_plan(plan)

     do i=1,(x_nbins+1)/2+1
        print*, fft_result(i)
        !sum_delta_sq(i) = sum_delta_sq(i) + real(fft_result(i),8)**2.
     end do
  end do
  do i=1,(x_nbins+1)/2+1
     !print*,sum_delta_sq(i)/(last_l-first_l+1)
  end do

  call exit
     
  deallocate(frequency_value,tmp_distance_value)
  
#ifdef RR
end subroutine makepowerspectrum_rr
#else
end subroutine makepowerspectrum_rg
#endif
