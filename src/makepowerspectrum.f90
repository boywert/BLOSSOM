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
  integer(kind=8) :: plan,plan_rev
  integer :: i,j
  real(kind=8) :: mean_den
  real(kind=8) :: d0,delta_x
  real(kind=8) :: nu_max,nu_min,max_observe,min_observe,max_box
  integer :: obs_freq_nbins, fine_freq_nbins, x_nbins, obs_to_fine
  real(kind=8), allocatable :: obs_frequency_value(:),fine_frequency_value(:),tmp_distance_value(:),tmp_signal_fine(:),tmp_signal_obs(:)
  real(kind=8), allocatable :: x_array(:), y_array(:),y_tmp(:)
  character(len=100) :: str_rank,z_s,str_line
  real(kind=8) :: M0,impact_param,nu_dist,nu_undist,this_absorp,delta_nu,width_real
  complex(kind=8), allocatable :: fft_result(:)
  real(kind=8), allocatable :: sum_delta_sq(:),ps_3D(:),ps_1D(:),k_1D(:),k_3D(:)
  ! Prepare strings
  write(z_s,'(f10.3)') z
  z_s = adjustl(z_s)
  write(str_rank,'(i10)') rank
  str_rank = adjustl(str_rank)

  !Calculate frequency and distance range
  d0 = z_to_d(z)
  max_observe = d0 !+ convert_length2physical(real(Boxsize*line_length_factor/2.,8),z) 
  nu_min = d_to_nu(max_observe)-obsfreqresolution !add buffer space
  min_observe = d0 - convert_length2physical(real(Boxsize*line_length_factor/2.,8),z) 
  nu_max = d_to_nu(min_observe)+obsfreqresolution !add buffer space
  obs_freq_nbins = ceiling((nu_max-nu_min)/obsfreqresolution)
  fine_freq_nbins = ceiling((nu_max-nu_min)/maxfreqresolution)
  obs_to_fine = int((obsfreqresolution+0.001)/maxfreqresolution)

  !Allocate arrays
  allocate(fine_frequency_value(0:fine_freq_nbins),obs_frequency_value(0:obs_freq_nbins))
  allocate(tmp_distance_value(0:obs_freq_nbins))
  allocate(tmp_signal_obs(0:obs_freq_nbins),tmp_signal_fine(0:fine_freq_nbins))
  do i=0,obs_freq_nbins
     obs_frequency_value(i) = nu_max - i*obsfreqresolution
     tmp_distance_value(i) = nu_to_d(obs_frequency_value(i))
  end do
  do i=0,fine_freq_nbins
     fine_frequency_value(i) = nu_max - i*maxfreqresolution
  end do
  do i=obs_freq_nbins,1,-1
     tmp_distance_value(i) = (tmp_distance_value(i)-tmp_distance_value(0))/Mpc*(1.+z)*h
  end do
  tmp_distance_value(0) = 0.
  
  !Set up new x arrays
  delta_x = tmp_distance_value(1)              !high frequency has higher delta_x
  
  x_nbins = ceiling(tmp_distance_value(obs_freq_nbins)/delta_x)-1
  allocate(x_array(0:x_nbins-1))
  allocate(y_array(0:x_nbins-1))
  allocate(fft_result(0:(x_nbins)/2))
  allocate(sum_delta_sq(0:(x_nbins)/2))
  do i=0,x_nbins-1
     x_array(i) = (i)*delta_x
  end do
  max_box = x_array(x_nbins-1)
  sum_delta_sq(:) = 0.0
  do j= first_l, last_l
     write(str_line,'(i10)') j
     str_line = adjustl(str_line)

     tmp_signal_fine(:) = 1.0
     tmp_signal_obs(:) = 0.0
     print*, trim(result_path)//z_s(1:len_trim(z_s))//'/RG/0.000/'//trim(adjustl(str_rank))//'/sins.'//trim(adjustl(str_line)//".1"
     open (unit=10, &
          file=trim(result_path)//z_s(1:len_trim(z_s))//'/RG/0.000/'//trim(adjustl(str_rank))//'/sins.'//trim(adjustl(str_line)//".1"), &
          form='binary')
     do
        read(10,end=327) M0,impact_param,nu_dist,nu_undist,this_absorp,delta_nu
        width_real = delta_nu/nu0*nu_dist
        tmp_signal_fine = tmp_signal_fine * (1.-this_absorp* exp(-0.5*((nu_dist-fine_frequency_value)/width_real)**2.0))
     end do
327  close(10)
     do i=0,fine_freq_nbins-1
        tmp_signal_obs(i/obs_to_fine) = tmp_signal_obs(i/obs_to_fine) + tmp_signal_fine(i)
     end do
     
     tmp_signal_obs = tmp_signal_obs/obs_to_fine
 
     
     call array_intrpol(tmp_distance_value(0:obs_freq_nbins-1),tmp_signal_obs(0:obs_freq_nbins-1),obs_freq_nbins,x_array(0:x_nbins-1),y_array(0:x_nbins-1),x_nbins)


     y_array = 1.- y_array 
     mean_den = sum(y_array)/x_nbins
     y_array = y_array/mean_den -1.
     

     call dfftw_plan_dft_r2c_1d(plan,x_nbins,y_array,fft_result,FFTW_ESTIMATE)
     call dfftw_execute_dft_r2c(plan,y_array,fft_result)

  
     ! call dfftw_plan_dft_c2r_1d(plan_rev,x_nbins,fft_result,x_array,FFTW_ESTIMATE)
     ! call dfftw_execute(plan_rev)

     sum_delta_sq = sum_delta_sq + (abs(fft_result)/real(x_nbins,8))**2.

     ! call dfftw_destroy_plan(plan_rev)
     call dfftw_destroy_plan(plan)
  end do
  allocate(ps_1D(0:x_nbins/2))
  allocate(ps_3D(1:x_nbins/2))
  allocate(k_1D(0:x_nbins/2))
  allocate(k_3D(1:x_nbins/2))
  ps_1D = sum_delta_sq/real(last_l-first_l+1,8)
  do i=0,x_nbins/2
     k_1D(i) = real(i)/max_box
  end do
  do i=1,x_nbins/2
     ps_3d(i) = -1.* (ps_1D(i)-ps_1D(i-1))/(k_1D(i)-k_1D(i-1))*2.*pi/(k_1D(i))
     print*,k_1D(i),ps_1D(i),ps_3D(i)
  end do
  call exit
     
  !deallocate(frequency_value,tmp_distance_value)
  
#ifdef RR
end subroutine makepowerspectrum_rr
#else
end subroutine makepowerspectrum_rg
#endif
