subroutine makeobserveline(z)
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
  real(kind=8) :: nu_self_dist, nu_self_undist
  real(kind=8) :: M0,sigma_V,impact_param,r_t
  real(kind=8), parameter :: nu_binsize = 1000., nu_binsize_fine = 100.
  real(kind=8), allocatable :: minihaloline_perline(:,:), diskline_perline(:,:)
  character(len=256) :: dummy_string
  integer :: status_checkfiles

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
  element_flag(7:9) = .TRUE. ! 1:3 => velocity
  element_flag(15) = .TRUE.
  n_elements = get_n_elements(element_flag)

  !## Read data to node 0
  status_checkfiles = checkfiles(real(z,4))
  if(status_checkfiles == 1) then
     if(rank == 0) then
        call ascii_read_halo(real(z,4),element_flag,n_elements,halonumber)
        allocate(halodata(1:n_elements,1:halonumber))
        halodata(1:n_elements,1:halonumber) = tmpfloat(1:n_elements,1:halonumber)      
        deallocate(tmpfloat)
     end if
     call MPI_BARRIER(MPI_COMM_WORLD,ierr)
     call mpi_bcast(halonumber,1,mpi_integer,0,mpi_comm_world,ierr)
     !print*, "total halos",halonumber
  else if(status_checkfiles ==2) then
     if(rank == 0) then
        halonumber = get_halo_number(real(z,4))
     end if
     call MPI_BARRIER(MPI_COMM_WORLD,ierr)
     call mpi_bcast(halonumber,1,mpi_integer,0,mpi_comm_world,ierr)
     call MPI_BARRIER(MPI_COMM_WORLD,ierr)
     !print*, "total halos",halonumber
     if(rank == 0) allocate(halodata(1:n_elements,1:halonumber))
     call mpi_read_halo(real(z,4),element_flag,n_elements,halonumber, halodata)
  else
     if(rank==0) print*, "No file to read"
     call abort
  end if
  !## End read data to node 0


  if(rank /= 0) allocate(halodata(1:n_elements,1:halonumber))
  call MPI_BARRIER( MPI_COMM_WORLD,ierr)
  !call abort
  if(rank ==0) print*,"Broadcasting from node 0"
  if(rank ==0) print*,halodata(1:n_elements,1)
  call mpi_bcast(halodata,halonumber*n_elements,mpi_real,0,mpi_comm_world,ierr)
  if(rank ==0) print*,"Finish broadcasting from node 0"
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
     haloperline(lineid(i)) = haloperline(lineid(i))+1
     linelinkedlist(i) = headofline(lineid(i))
     headofline(lineid(i)) = i
  end do

  d0 = z_to_d(z)

  do i=1,n_point
 !    curHalo = haloid(i)
 !    comov_dist(i) = online(i)
 !    undistorted_dist(i) = d0+convert_length2physical(real((BoxSize/2.-online(i)),8),z)
 !    undistorted_z = d_to_z(undistorted_dist(i))
 !    distorted_z = undistorted_z+convert_vel2physical(real(dotproduct(-1.*direction(1:3,i),halodata(1:3,curHalo)),8),z)/c
 !    distorted_dist(i) = z_to_d(distorted_z)

 !    haloperline(lineid(i)) = haloperline(lineid(i))+1
 !    if(convert_mass2physical(real(halodata(4,curHalo),8))/M_sol > mass_limit) then
 !       haloperline_large(lineid(i)) =  haloperline_large(lineid(i)) +1
 !    else
 !       haloperline_small(lineid(i)) =  haloperline_small(lineid(i)) + 1
 !    end if
 !    linelinkedlist(i) = headofline(lineid(i))
 !    headofline(lineid(i)) = i
  end do


  
  !M0 = convert_mass2physical(real(halodata(1,curHalo_id),8))/M_sol
  do i=1,max_line
     curHalo = headofline(i)
     do while (curHalo /= 0)
        nu_self_undist = d_to_nu(d0+convert_length2physical(real((BoxSize/2.-online(Curhalo)),8),z))
  !      nu_self_dist = 
        curHalo = linelinkedlist(curHalo)
     end do
  end do
   
end subroutine makeobserveline
