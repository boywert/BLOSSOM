subroutine absorpstats(z)
  use mpitools
  use omp_lib
  use common_vars
  use io_tools
  use absorptiontools
  use conversiontools

  implicit none

  ! should have used halfbox instead
  real(kind=4) :: Halfbox, BinSize
  integer(kind=4),parameter :: max_haloperlinehisto = 200


  integer(kind=8) :: i,j,k,n_point,totalbin,totalpoint,omp_thread

  integer(kind=4) :: fh_lineid,fh_online,fh_toline,fh_haloid
  integer(kind=8) :: curHalo,innerHalo,block,curhalo_id
  integer(kind=mpi_offset_kind) :: filesize
  integer(kind=4),allocatable :: lineid(:),linelinkedlist(:),haloid(:)
  real(kind=4), allocatable :: online(:),toline(:)
  integer(kind=8) :: haloperline(max_line),headofline(max_line)
  real(kind=8),allocatable :: minihalohisto(:),diskhisto(:)
  real(kind=8),allocatable :: sub_minihalohisto(:,:),sub_diskhisto(:,:)
  real(kind=8),allocatable :: minihalohisto_final(:),diskhisto_final(:)
  real(kind=8),allocatable :: correlation(:)
  integer(kind=4) :: n
  real(kind=8) :: tau, absorp,r(0:max_size),rho(0:max_size)
  integer(kind=8),allocatable :: n_histo(:), n_histo_final(:), sub_n_histo(:,:)
  real(kind=8) :: z,d,d0,nu
  character(len=100) :: str_rank, z_s
  real(kind=8) :: fwhm, sum_fwhm, sum_fwhm_sq, final_fwhm, final_fwhm_sq
  real(kind=8),allocatable :: sub_sum_fwhm(:), sub_sum_fwhm_sq(:)
  logical :: element_flag(1:17)
  real(kind=4),allocatable :: halodata(:,:)
  integer(kind=4) :: halonumber,n_elements
  integer(kind=8),allocatable :: sub_n_point(:)
  integer(kind=4) :: fh_result
  real(kind=8) :: nu_min, nu_max, nu_Range, nu_self, nu_index
  real(kind=8) :: nu_min_fine, nu_max_fine, nu_Range_fine
  real(kind=8) :: M0,sigma_V,impact_param,r_t
  real(kind=8), parameter :: nu_binsize = 1000., nu_binsize_fine = 100.
  real(kind=8), allocatable :: minihaloline_perline(:,:), diskline_perline(:,:)
  character(len=256) :: dummy_string


  !@ define some parameters
  Halfbox = Boxsize/2.
  BinSize = 4.
  omp_thread = omp_get_max_threads()
  if(rank==0) then
     print*, 'Start running absorption line'
     print*,'OMP threads = ',omp_thread
     print*,'MPI threads = ',nodes_returned
     call run_unit_test(z)
  end if



  element_flag(1:17) = .FALSE.
  element_flag(15) = .TRUE. ! 1 => gridmass


  halonumber = get_halo_number(real(z,4))
  n_elements = get_n_elements(element_flag)

  allocate(halodata(1:n_elements,1:halonumber))
  call mpi_read_halo(real(z,4),element_flag,n_elements,halonumber, halodata)
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  call mpi_bcast(halodata,halonumber*n_elements,mpi_real,0,mpi_comm_world,ierr)
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)


  write(z_s,'(f10.3)') z
  z_s = adjustl(z_s)

  write(str_rank,'(i10)') rank
  omp_thread = omp_get_max_threads()


  call mpi_file_open(mpi_comm_self, &
       trim(los_path)//trim(z_s)//'/LINEID/'//trim(adjustl(str_rank)), &
       MPI_MODE_RDONLY, &
       mpi_info_null,fh_lineid,ierr)

  if(rank ==0) print*,'Get N-point ....'

  call mpi_file_get_size (fh_lineid, filesize, ierr)
  n_point = filesize/4

  if(rank ==0) print*, 'total point',n_point
  if(rank ==0) print*,'Allocate lineid,online,linkedlist ....'

  allocate(lineid(1:n_point))
  allocate(haloid(1:n_point))
  allocate(online(1:n_point))
  allocate(toline(1:n_point))
  allocate(linelinkedlist(1:n_point))

  if(rank ==0) print*,'Read data into lineid ....'

  call mpi_file_read(fh_lineid,lineid,n_point,mpi_integer,mpi_status_ignore,ierr)

  if(rank ==0) print*,'Close file lineid ....'

  call mpi_file_close(fh_lineid,ierr)

  if(rank ==0) print*,'Open file for online ....'

  call mpi_file_open(mpi_comm_self, &
       trim(los_path)//trim(z_s)//'/ONLINE/'//trim(adjustl(str_rank)), &
       MPI_MODE_RDONLY, &
       mpi_info_null,fh_online,ierr)

  if(rank ==0) print*,'Read data into online ....'

  call mpi_file_read(fh_online,online,n_point,mpi_real,mpi_status_ignore,ierr)
  if(rank ==0) then
     print*,'Close file online ....'
  endif
  call mpi_file_close(fh_online,ierr)

  if(rank ==0) print*,'Open file for online ....'

  call mpi_file_open(mpi_comm_self, &
       trim(los_path)//trim(z_s)//'/TOLINE/'//trim(adjustl(str_rank)), &
       MPI_MODE_RDONLY, &
       mpi_info_null,fh_toline,ierr)

  if(rank ==0) print*,'Read data into toline ....'

  call mpi_file_read(fh_toline,toline,n_point,mpi_real,mpi_status_ignore,ierr)
  if(rank ==0) then
     print*,'Close file online ....'
  endif
  call mpi_file_close(fh_toline,ierr)


  call mpi_file_open(mpi_comm_self, &
       trim(los_path)//trim(z_s)//'/HALOID/'//trim(adjustl(str_rank)), &
       MPI_MODE_RDONLY, &
       mpi_info_null,fh_haloid,ierr)

  if(rank ==0) print*,'Read data into toline ....'

  call mpi_file_read(fh_haloid,haloid,n_point,mpi_integer,mpi_status_ignore,ierr)
  if(rank ==0) then
     print*,'Close file online ....'
  endif
  call mpi_file_close(fh_haloid,ierr)


  ! make linked list
  if(rank ==0) then
     print*,'Making linkedlist ....'
     print*,'... Allocating.....'
  endif
  haloperline(1:max_line) = 0
  linelinkedlist(1:n_point) = 0
  headofline(1:max_line) = 0


  if(rank ==0) print*,'... Looping for linkedlist....'
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  do i=1,n_point
     !if(rank==0 .and. (lineid(i) > max_line .or. lineid(i) < 1)) &
     !     print*,i,lineid(i)
     !if(rank==1) print*,i,lineid(i)
     haloperline(lineid(i)) = haloperline(lineid(i))+1
     linelinkedlist(i) = headofline(lineid(i))
     headofline(lineid(i)) = i
  end do

  call n_cal(n,r,rho)

  d0 = z_to_d(z)
  nu_max = real(ceiling((d_to_nu(d0-convert_length2physical(Boxsize/2.d0,z)))/nu_binsize),8)*nu_binsize
  nu_min = real(floor((d_to_nu(d0+convert_length2physical(Boxsize/2.d0,z)))/nu_binsize),8)*nu_binsize
  nu_range = nu_max-nu_min

  nu_max_fine = real(ceiling((d_to_nu(d0-convert_length2physical(Boxsize/2.d0,z)))/nu_binsize_fine),8)*nu_binsize_fine
  nu_min_fine = real(floor((d_to_nu(d0+convert_length2physical(Boxsize/2.d0,z)))/nu_binsize_fine),8)*nu_binsize_fine
  nu_range_fine = nu_max_fine-nu_min_fine

  if(rank==0) then
     print*, ""
     print*,"############################################################"
     print*,"Start calculating absorption statistics"
     print*,"############################################################"
     print*,""
     print*, "Max Freq:",nu_max
     print*, "Min Freq:",nu_min
     print*, "Freq range:",nu_range
     print*, "Badwidth",nu_binsize

     print*,""

     print*, "Max Freq (fine):",nu_max_fine
     print*, "Min Freq (fine):",nu_min_fine
     print*, "Freq range (fine):",nu_range_fine
     print*, "Badwidth (fine):",nu_binsize_fine
  end if


  if(rank==0) print*, 'Allocate and set the initials'
  allocate(minihalohisto(1:ceiling(nu_range_fine/nu_binsize_fine)))
  allocate(diskhisto(1:ceiling(nu_range_fine/nu_binsize_fine)))
  allocate(sub_minihalohisto(1:ceiling(nu_range_fine/nu_binsize_fine),0:omp_thread-1))
  allocate(sub_diskhisto(1:ceiling(nu_range_fine/nu_binsize_fine),0:omp_thread-1))
  !allocate(sub_n_histo(0:ceiling(nu_range/nu_binsize),0:omp_thread-1))
  !allocate(n_histo(0:ceiling(nu_range/nu_binsize)))
  !allocate(sub_sum_fwhm(0:omp_thread-1))
  !allocate(sub_sum_fwhm_sq(0:omp_thread-1))
  !allocate(sub_n_point(0:omp_thread-1))




  minihalohisto(1:size(minihalohisto)-1) = 0.
  diskhisto(1:size(diskhisto)-1) = 0.
  sub_minihalohisto(1:ceiling(nu_range_fine/nu_binsize_fine),0:omp_thread-1) = 0.
  sub_diskhisto(1:ceiling(nu_range_fine/nu_binsize_fine),0:omp_thread-1) = 0.
  !n_histo(0:size(n_histo)-1) = 0
  !sub_n_histo(0:ceiling(nu_range/nu_binsize),0:omp_thread-1) = 0
  !sum_fwhm_sq = 0.d0
  !sum_fwhm = 0.d0
  !sub_n_point(0:omp_thread-1) = 0
  !sub_sum_fwhm_sq(0:omp_thread-1) = 0.d0
  !sub_sum_fwhm(0:omp_thread-1) = 0.d0

  allocate(minihaloline_perline(1:ceiling(nu_range_fine/nu_binsize_fine),0:omp_thread-1))
  allocate(diskline_perline(1:ceiling(nu_range_fine/nu_binsize_fine),0:omp_thread-1))

  if(rank==0) print*, 'Start calculation'

  !$omp parallel private(j,nu_self,curHalo_id,curHalo,block,tau,absorp,fwhm,M0,sigma_V,r_t,impact_param) 
  !$omp do
  do i=1,max_line
     curHalo = headofline(i)

#ifdef DEBUG
     if(curHalo < 0 .or. curHalo > n_point) then
        print*, &
             'rank=',rank,'mp=',omp_get_thread_num(), &
             'curHalo=',curhalo
        call abort
     endif
#endif

     minihaloline_perline(1:ceiling(nu_range_fine/nu_binsize_fine),omp_get_thread_num()) = 0.d0
     diskline_perline(1:ceiling(nu_range_fine/nu_binsize_fine),omp_get_thread_num()) = 0.d0
     do while(curHalo /= 0)

        curHalo_id = haloid(curHalo)
        nu_self = d_to_nu(d0+convert_length2physical(real((BoxSize/2.-online(Curhalo)),8),z))

        M0 = convert_mass2physical(real(halodata(1,curHalo_id),8))/M_sol
        impact_param = convert_length2physical(real(toline(curHalo),8),z)
        if(M0 <= 4.d7) then
           call tau_cal(M0,z,impact_param,n,r,rho,sigma_V,tau)
           fwhm = voigt_FWHM(A10,sigma_V,nu_self)    
           minihaloline_perline(1:ceiling(nu_range_fine/nu_binsize_fine),omp_get_thread_num()) = &
                minihaloline_perline(1:ceiling(nu_range_fine/nu_binsize_fine),omp_get_thread_num()) + &
                absorption_line(nu_min_fine,nu_max_fine,nu_binsize_fine,ceiling(nu_range_fine/nu_binsize_fine),nu_self,fwhm,tau)
        else
           call tau_cal(M0,z,impact_param,n,r,rho,sigma_V,tau)
           fwhm = voigt_FWHM(A10,sigma_V,nu_self)  
           diskline_perline(1:ceiling(nu_range_fine/nu_binsize_fine),omp_get_thread_num()) = &
                diskline_perline(1:ceiling(nu_range_fine/nu_binsize_fine),omp_get_thread_num()) + &
                absorption_line(nu_min_fine,nu_max_fine,nu_binsize_fine,ceiling(nu_range_fine/nu_binsize_fine),nu_self,fwhm,tau)
        endif

        !sub_sum_fwhm(omp_get_thread_num()) = sub_sum_fwhm(omp_get_thread_num()) + fwhm
        !sub_sum_fwhm_sq(omp_get_thread_num()) = sub_sum_fwhm_sq(omp_get_thread_num()) + fwhm**2.d0
        !sub_tauhisto(block,omp_get_thread_num()) = sub_tauhisto(block,omp_get_thread_num()) + tau
        !sub_absorphisto(block,omp_get_thread_num()) = sub_absorphisto(block,omp_get_thread_num()) + (1.d0-exp(-1.d0*tau))/nu_binsize
        !sub_n_histo(block,omp_get_thread_num()) = sub_n_histo(block,omp_get_thread_num()) + 1
        !sub_n_point(omp_get_thread_num()) = sub_n_point(omp_get_thread_num()) + 1


        curHalo = linelinkedlist(curHalo)
     end do
     do j=1,ceiling(nu_range_fine/nu_binsize_fine)
        minihaloline_perline(j,omp_get_thread_num()) = &
             dmin1(1.d0,minihaloline_perline(j,omp_get_thread_num()))
        diskline_perline(j,omp_get_thread_num()) = &
             dmin1(1.d0,diskline_perline(j,omp_get_thread_num()))
     end do

     sub_minihalohisto(1:ceiling(nu_range_fine/nu_binsize_fine),omp_get_thread_num()) = &
          sub_minihalohisto(1:ceiling(nu_range_fine/nu_binsize_fine),omp_get_thread_num()) + &
          minihaloline_perline(1:ceiling(nu_range_fine/nu_binsize_fine),omp_get_thread_num())

     sub_diskhisto(1:ceiling(nu_range_fine/nu_binsize_fine),omp_get_thread_num()) = &
          sub_diskhisto(1:ceiling(nu_range_fine/nu_binsize_fine),omp_get_thread_num()) + &
          diskline_perline(1:ceiling(nu_range_fine/nu_binsize_fine),omp_get_thread_num())

  end do
  !$omp end do
  !$omp end parallel


  call MPI_BARRIER(MPI_COMM_WORLD,ierr)


  if(rank==0) then
     print*, "Start to save single line"
     i = line_use_to_save
     curHalo = headofline(i)

#ifdef DEBUG
     if(curHalo < 0 .or. curHalo > n_point) then
        print*, &
             'rank=',rank,'mp=',omp_get_thread_num(), &
             'curHalo=',curhalo
        call abort
     endif
#endif

     minihaloline_perline(1:ceiling(nu_range_fine/nu_binsize_fine),0) = 0.d0
     diskline_perline(1:ceiling(nu_range_fine/nu_binsize_fine),0) = 0.d0
     do while(curHalo /= 0)

        curHalo_id = haloid(curHalo)
        nu_self = d_to_nu(d0+convert_length2physical(real((BoxSize/2.-online(Curhalo)),8),z))

        M0 = convert_mass2physical(real(halodata(curHalo_id,1),8))/M_sol
        impact_param = convert_length2physical(real(toline(curHalo),8),z)
        if(M0 <= 4.d7) then
           call tau_cal(M0,z,impact_param,n,r,rho,sigma_V,tau)
           fwhm = voigt_FWHM(A10,sigma_V,nu_self)    
           minihaloline_perline(1:ceiling(nu_range_fine/nu_binsize_fine),0) = &
                minihaloline_perline(1:ceiling(nu_range_fine/nu_binsize_fine),0) + &
                absorption_line(nu_min_fine,nu_max_fine,nu_binsize_fine,ceiling(nu_range_fine/nu_binsize_fine),nu_self,fwhm,tau)
        else
           call tau_cal(M0,z,impact_param,n,r,rho,sigma_V,tau)
           fwhm = voigt_FWHM(A10,sigma_V,nu_self)  
           diskline_perline(1:ceiling(nu_range_fine/nu_binsize_fine),0) = &
                diskline_perline(1:ceiling(nu_range_fine/nu_binsize_fine),0) + &
                absorption_line(nu_min_fine,nu_max_fine,nu_binsize_fine,ceiling(nu_range_fine/nu_binsize_fine),nu_self,fwhm,tau)
        endif

        !sub_sum_fwhm(omp_get_thread_num()) = sub_sum_fwhm(omp_get_thread_num()) + fwhm
        !sub_sum_fwhm_sq(omp_get_thread_num()) = sub_sum_fwhm_sq(omp_get_thread_num()) + fwhm**2.d0
        !sub_tauhisto(block,omp_get_thread_num()) = sub_tauhisto(block,omp_get_thread_num()) + tau
        !sub_absorphisto(block,omp_get_thread_num()) = sub_absorphisto(block,omp_get_thread_num()) + (1.d0-exp(-1.d0*tau))/nu_binsize
        !sub_n_histo(block,omp_get_thread_num()) = sub_n_histo(block,omp_get_thread_num()) + 1
        !sub_n_point(omp_get_thread_num()) = sub_n_point(omp_get_thread_num()) + 1


        curHalo = linelinkedlist(curHalo)
     end do
     do j=1,ceiling(nu_range_fine/nu_binsize_fine)
        minihaloline_perline(j,0) = &
             dmin1(1.d0,minihaloline_perline(j,0))
        diskline_perline(j,0) = &
             dmin1(1.d0,diskline_perline(j,0))
     end do



     print*, "save file"
     dummy_string = "save_single_line_minihalo.dat"
     call save_absorptionline(minihaloline_perline(1:ceiling(nu_range_fine/nu_binsize_fine),0),ceiling(nu_range_fine/nu_binsize_fine),nu_min_fine,nu_binsize_fine,z_s,dummy_string)

     dummy_string = "save_single_line_galacticdisk.dat"
     call save_absorptionline(diskline_perline(1:ceiling(nu_range_fine/nu_binsize_fine),0),ceiling(nu_range_fine/nu_binsize_fine),nu_min_fine,nu_binsize_fine,z_s,dummy_string)

  endif





  deallocate(minihaloline_perline,diskline_perline)
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  if(rank==0) print*,'Start combining results in nodes'
  n_point = 0
  !if(rank==0) print*, sub_n_histo
  do i=0, omp_thread-1
     !sum_fwhm_sq = sum_fwhm_sq + sub_sum_fwhm_sq(i)
     !sum_fwhm = sum_fwhm + sub_sum_fwhm(i)
     !n_point = n_point + sub_n_point(i)
     !n_histo(0:ceiling(nu_range/nu_binsize)) = n_histo(0:ceiling(nu_range/nu_binsize)) + sub_n_histo(0:ceiling(nu_range/nu_binsize),i)
     minihalohisto(1:ceiling(nu_range_fine/nu_binsize_fine)) = minihalohisto(1:ceiling(nu_range_fine/nu_binsize_fine)) + sub_minihalohisto(1:ceiling(nu_range_fine/nu_binsize_fine),i)
     diskhisto(1:ceiling(nu_range_fine/nu_binsize_fine)) = diskhisto(1:ceiling(nu_range_fine/nu_binsize_fine)) + sub_diskhisto(1:ceiling(nu_range_fine/nu_binsize_fine),i)
  end do

  !deallocate(sub_sum_fwhm, sub_sum_fwhm_sq)
  deallocate(sub_minihalohisto,sub_diskhisto)
  !deallocate(sub_n_point)
  deallocate(halodata)
  deallocate(lineid,online,toline,linelinkedlist,haloid)
  !print*,'rank',rank, n_point, sum_fwhm, sum_fwhm_sq
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  if(rank==0) then
     allocate(minihalohisto_final(1:ceiling(nu_range_fine/nu_binsize_fine)))
     allocate(diskhisto_final(1:ceiling(nu_range_fine/nu_binsize_fine)))
     !allocate(n_histo_final(0:ceiling(nu_range/nu_binsize)))  
  endif
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  call mpi_reduce(minihalohisto, &
       minihalohisto_final,size(minihalohisto), &
       mpi_real8,mpi_sum,0,mpi_comm_world,ierr)
  call mpi_reduce(diskhisto, &
       diskhisto_final,size(diskhisto), &
       mpi_real8,mpi_sum,0,mpi_comm_world,ierr)
  !call mpi_reduce(n_histo, &
  !     n_histo_final,size(n_histo), &
  !     mpi_integer8,mpi_sum,0,mpi_comm_world,ierr)


  !call mpi_reduce(sum_fwhm, &
  !     final_fwhm,1, &
  !     mpi_real8,mpi_sum,0,mpi_comm_world,ierr)
  !call mpi_reduce(sum_fwhm_sq, &
  !     final_fwhm_sq,1, &
  !     mpi_real8,mpi_sum,0,mpi_comm_world,ierr)
  !call mpi_reduce(n_point, &
  !     totalpoint,1, &
  !     mpi_integer8,mpi_sum,0,mpi_comm_world,ierr)

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  deallocate(minihalohisto,diskhisto)
  !deallocate(n_histo)
  if(rank==0) then
     call system("mkdir -p "//trim(result_path)//z_s(1:len_trim(z_s))//'/')
     minihalohisto_final(1:ceiling(nu_range_fine/nu_binsize_fine)) = minihalohisto_final(1:ceiling(nu_range_fine/nu_binsize_fine))/real(max_line,8)
     diskhisto_final(1:ceiling(nu_range_fine/nu_binsize_fine)) = diskhisto_final(1:ceiling(nu_range_fine/nu_binsize_fine))/real(max_line,8)
     dummy_string = "minihalohisto.dat"
     call save_absorptionline(minihalohisto_final,size(minihalohisto_final),nu_min_fine,nu_binsize_fine,z_s,dummy_string)
     dummy_string = "galacticdiskhisto.dat"
     call save_absorptionline(diskhisto_final,size(diskhisto_final),nu_min_fine,nu_binsize_fine,z_s,dummy_string)
     !call system("rm -rf "//trim(result_path)//z_s(1:len_trim(z_s))//'/fwhm.dat')
     !open (unit=50,file=trim(result_path)//z_s(1:len_trim(z_s))//'/'//'fwhm.dat',STATUS = 'NEW')
     !write(50,*) 'mean:', final_fwhm/real(totalpoint,8)
     !print*,  'mean:', final_fwhm/real(totalpoint,8)
     !write(50,*) 'error:', sqrt(final_fwhm_sq/real(totalpoint,8)-(final_fwhm/real(totalpoint,8))**2.d0)/sqrt(real(totalpoint,8))
     !print*, 'error:', sqrt(final_fwhm_sq/real(totalpoint,8)-(final_fwhm/real(totalpoint,8))**2.d0)/sqrt(real(totalpoint,8))
     !write(50,*) 'N:',totalpoint
     !print*, 'N:',totalpoint
     !close(50)
     !call system("rm -rf "//trim(result_path)//z_s(1:len_trim(z_s))//'/transmission.dat')
     !open (unit=50,file=trim(result_path)//z_s(1:len_trim(z_s))//'/'//'transmission.dat',STATUS = 'NEW')
     !do i=0, size(n_histo_final)-1
     !   write(50,*)  nu_min+nu_Binsize*real(i,8)+nu_Binsize/2., n_histo_final(i), 1.d0-absorphisto_final(i)/real(max_line*nodes_returned,8)
     !end do
     !close(50)
     deallocate(minihalohisto_final,diskhisto_final)
     !deallocate(n_histo_final)
  endif

end subroutine absorpstats
