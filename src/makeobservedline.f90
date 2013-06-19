#ifdef RR
subroutine makeobservedlines_rr(z)
#else
subroutine makeobservedlines_rg(z)
#endif
  use mpi
  use mpitools
  use common_vars
  use conversiontools
  use absorptiontools
  use omp_lib
  use io_tools
  use vectortools
  implicit none

  ! should have used halfbox instead
  real(kind=8) :: BinSize
  integer(kind=4),parameter :: max_haloperlinehisto = 200


  integer(kind=8) :: i,j,k,n_point,totalbin,totalpoint,omp_thread
  integer(kind=4) :: fh_hitpoint,fh_direction,fh_haloid
  integer(kind=4) :: fh_lineid,fh_online,fh_toline, line_with_max_halo, max_halo_so_far
  integer(kind=8) :: curHalo,curHaloid,innerHalo,block
  integer(kind=mpi_offset_kind) :: filesize
  integer(kind=4),allocatable :: lineid(:),linelinkedlist(:),haloid(:),haloperline(:),headofline(:)
  real(kind=4), allocatable :: online(:),toline(:),direction(:,:)

  real(kind=8) :: Halfbox,mass_limit,lambda,line_centre
  character(len=100) :: str_rank,z_s
  real(kind=8) :: z,d0,d_self
  integer(kind=4) :: halonumber,n_elements,rlogsteps
  real(kind=4), allocatable :: halodata(:,:)
  real(kind=8), allocatable :: distorted_dist(:),undistorted_dist(:),comov_dist(:)
  logical :: element_flag(1:17)
  real(kind=8) :: undistorted_z,distorted_z,sumtest
  integer(kind=8) :: usedlines, totusedlines
  real(kind=8) :: rmax,rmin,rminlog,rmaxlog,deltalogr
  integer :: status_checkfiles
  real(kind=8) :: nu_dist, nu_undist

  integer(kind=4) :: n
  real(kind=8) :: tau, area_tau,absorp,extend_absorp, r(0:max_size), rho(0:max_size)
  real(kind=8) :: max_observe, min_observe, nu_min, nu_max, d_source,nu_source
  real(kind=8) :: M0,impact_param,sigma_V,gaussian_sd
 


  Halfbox = real(BoxSize,8)/2.d0
  line_centre = real(Boxsize*line_length_factor,8)/2.d0

  write(z_s,'(f10.3)') z
  z_s = adjustl(z_s)
  write(str_rank,'(i10)') rank
  str_rank = adjustl(str_rank)
  omp_thread = omp_get_max_threads()

  mass_limit = 1.d8

  if(rank ==0) then
     print*,'OMP threads = ',omp_thread
     print*,'MPI threads = ',nodes_returned
  endif

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
 
  if(rank ==0) print*,"Start reading LOS data"


#ifdef RR
  call mpi_file_open(mpi_comm_self, &
       trim(los_path)//'/'//z_s(1:len_trim(z_s))//'/RR/LINEID/'//trim(adjustl(str_rank)), &
       MPI_MODE_RDONLY, &
       mpi_info_null,fh_lineid,ierr)
#else
  call mpi_file_open(mpi_comm_self, &
       trim(los_path)//'/'//z_s(1:len_trim(z_s))//'/RG/LINEID/'//trim(adjustl(str_rank)), &
       MPI_MODE_RDONLY, &
       mpi_info_null,fh_lineid,ierr)
#endif
  if(rank ==0) print*,'Get N-point ....'

  call mpi_file_get_size (fh_lineid, filesize, ierr)
  n_point = filesize/4

  if(rank ==0) print*, 'total point',n_point
  if(rank ==0) print*,'Allocate lineid,online,linkedlist ....'

  allocate(lineid(1:n_point))
  allocate(online(1:n_point))
  allocate(toline(1:n_point))
  allocate(haloid(1:n_point))
  allocate(direction(1:3,1:n_point))


  if(rank ==0) print*,'Read data into lineid ....'
  call mpi_file_read(fh_lineid,lineid,n_point,mpi_integer,mpi_status_ignore,ierr)

  if(rank ==0) print*,'Close file lineid ....'

  call mpi_file_close(fh_lineid,ierr)

  if(rank ==0) print*,'Open file for online ....'

#ifdef RR
  call mpi_file_open(mpi_comm_self, &
        trim(los_path)//'/'//z_s(1:len_trim(z_s))//'/RR/ONLINE/'//trim(adjustl(str_rank)), &
       MPI_MODE_RDONLY, &
       mpi_info_null,fh_online,ierr)
#else
  call mpi_file_open(mpi_comm_self, &
        trim(los_path)//'/'//z_s(1:len_trim(z_s))//'/RG/ONLINE/'//trim(adjustl(str_rank)), &
       MPI_MODE_RDONLY, &
       mpi_info_null,fh_online,ierr)
#endif
  if(rank ==0) print*,'Read data into online ....'

  call mpi_file_read(fh_online,online,n_point,mpi_real,mpi_status_ignore,ierr)
  if(rank ==0) then
     print*,'Close file online ....'
  endif
  call mpi_file_close(fh_online,ierr)


  if(rank ==0) print*,'Open file for toline ....'

#ifdef RR
  call mpi_file_open(mpi_comm_self, &
        trim(los_path)//'/'//z_s(1:len_trim(z_s))//'/RR/TOLINE/'//trim(adjustl(str_rank)), &
       MPI_MODE_RDONLY, &
       mpi_info_null,fh_toline,ierr)
#else
  call mpi_file_open(mpi_comm_self, &
        trim(los_path)//'/'//z_s(1:len_trim(z_s))//'/RG/TOLINE/'//trim(adjustl(str_rank)), &
       MPI_MODE_RDONLY, &
       mpi_info_null,fh_toline,ierr)
#endif
  if(rank ==0) print*,'Read data into toline ....'

  call mpi_file_read(fh_toline,toline,n_point,mpi_real,mpi_status_ignore,ierr)
  if(rank ==0) then
     print*,'Close file toline ....'
  endif
  call mpi_file_close(fh_toline,ierr)


  if(rank ==0) print*,'Open file for direction ....'

#ifdef RR
  call mpi_file_open(mpi_comm_self, &
       trim(los_path)//'/'//z_s(1:len_trim(z_s))//'/RR/DIRECTION/'//trim(adjustl(str_rank)), &
       MPI_MODE_RDONLY, &
       mpi_info_null,fh_direction,ierr)
#else
  call mpi_file_open(mpi_comm_self, &
        trim(los_path)//'/'//z_s(1:len_trim(z_s))//'/RG/DIRECTION/'//trim(adjustl(str_rank)), &
       MPI_MODE_RDONLY, &
       mpi_info_null,fh_direction,ierr)
#endif
  if(rank ==0) print*,'Read data into direciton ....'

  call mpi_file_read(fh_direction,direction,3*n_point,mpi_real,mpi_status_ignore,ierr)
  if(rank ==0) then
     print*,'Close file online ....'
  endif
  call mpi_file_close(fh_direction,ierr)


#ifdef RR
  call mpi_file_open(mpi_comm_self, &
       trim(los_path)//trim(z_s)//'/RR/HALOID/'//trim(adjustl(str_rank)), &
       MPI_MODE_RDONLY, &
       mpi_info_null,fh_haloid,ierr)
#else
  call mpi_file_open(mpi_comm_self, &
       trim(los_path)//trim(z_s)//'/RG/HALOID/'//trim(adjustl(str_rank)), &
       MPI_MODE_RDONLY, &
       mpi_info_null,fh_haloid,ierr)
#endif
  if(rank ==0) print*,'Read data into haloid ....'

  call mpi_file_read(fh_haloid,haloid,n_point,mpi_integer,mpi_status_ignore,ierr)
  if(rank ==0) then
     print*,'Close file online ....'
  endif
  call mpi_file_close(fh_haloid,ierr)

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  ! make linked list
  if(rank ==0) then
     print*,'Making linkedlist ....'
     print*,'... Allocating.....'
  endif

  allocate(distorted_dist(1:n_point))
  allocate(comov_dist(1:n_point))
  allocate(undistorted_dist(1:n_point))
  allocate(linelinkedlist(1:n_point))
  allocate(headofline(1:max_line))
  allocate(haloperline(1:max_line))

  if(rank ==0) then
     print*,'Set varaiables to 0 ...'
  endif
  haloperline(1:max_line) = 0
  linelinkedlist(1:n_point) = 0
  headofline(1:max_line) = 0


  if(rank ==0) print*,'... Looping for linkedlist....'
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  d0 = z_to_d(z) !physical distance from observer to the centre of line

  do i=1,n_point
     curHaloid = haloid(i)
     comov_dist(i) = online(i)
     undistorted_dist(i) = d0+convert_length2physical(real(online(i),8),z) - convert_length2physical(real(Boxsize*line_length_factor/2.,8),z)
     undistorted_z = d_to_z(undistorted_dist(i))
     distorted_z = undistorted_z+convert_vel2physical(real(dotproduct(direction(1:3,i),halodata(1:3,curHaloid)),8),z)/c
     distorted_dist(i) = z_to_d(distorted_z)
     !nu_dist = d_to_nu(distorted_dist(i))
     !nu_undist = d_to_nu(undistorted_dist(i))

     !@ Use only minihalos
     M0 = convert_mass2physical(real(halodata(4,curHaloid),8))/M_sol
     if(M0 > 1.e5 .and. M0 < 1.e8) then
        linelinkedlist(i) = headofline(lineid(i))
        headofline(lineid(i)) = i
        haloperline(lineid(i)) = haloperline(lineid(i)) + 1
     end if
  end do

  !call n_cal(n,r,rho)
  
  max_observe = d0 + convert_length2physical(real(Boxsize*line_length_factor/2.,8),z) 
  nu_min = d_to_nu(max_observe)
  min_observe = d0 - convert_length2physical(real(Boxsize*line_length_factor/2.,8),z) 
  nu_max = d_to_nu(min_observe)

#ifndef RR
  d_source = d0 + convert_length2physical(real(Boxsize*(real(line_length_factor)-0.5),8),z) 
  nu_source = d_to_nu(d_source)
#endif

#ifdef DEBUG
  if(rank==0) then
     print*, 'Start simulating observe data'
  end if
#endif
  
  call n_cal(n,r,rho)
  if(rank==0) call system("rm -f "//trim(result_path)//z_s(1:len_trim(z_s))//'/'//'out.dat')
  if(rank==0) open(unit=53,file=trim(result_path)//z_s(1:len_trim(z_s))//'/'//'out.dat',STATUS = 'NEW')
  !@ use openmp in tau_cal instead since it's troublesome 
  !@ to set up multiple arrays to collect information
  
  do i=1,100 !max_line
#ifdef DEBUG
     if(rank==0) then
        print*, 'line =',i
     end if
#endif
     curHalo = headofline(i)
     do while(curHalo /= 0)
        curHaloid = haloid(curHalo)
        nu_dist = d_to_nu(distorted_dist(curhalo))
        nu_undist = d_to_nu(undistorted_dist(curhalo))
        M0 = convert_mass2physical(real(halodata(4,curHaloid),8))/M_sol
        impact_param = convert_length2physical(real(toline(curHalo),8),z)
        call tau_cal(M0,z,impact_param,n,r,rho,sigma_V,area_tau,tau)
        absorp = 1.d0 - exp(-1*tau)
        extend_absorp =  1.d0 - exp(-1*area_tau)
#ifdef DEBUG
        if(rank==0) then
           print*,'absorp',absorp,'impact',impact_param
        end if
#endif
        !@ convert to SI
        gaussian_sd = sqrt(2.)*sigma_V/c*nu_dist
        if(rank==0) write(53,*) int(i), real(nu_dist), real(absorp),real(extend_absorp),real(gaussian_sd)
        !gaussian_sd = sqrt(2.)*sigma_V/c*nu_undist
        !if(rank==0) write(53,*) int(i), real(nu_undist), real(absorp),real(extend_absorp),real(gaussian_sd)
        
        !@ next halo
        curHalo = linelinkedlist(curHalo)
     end do
  end do
  if(rank==0) close(53)
#ifdef RR
end subroutine makeobservedlines_rr
#else
end subroutine makeobservedlines_rg
#endif
