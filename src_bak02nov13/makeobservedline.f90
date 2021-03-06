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
  integer(kind=4) :: fh_hitpoint,fh_direction,fh_haloid,fh_record(21)
  integer(kind=4) :: fh_lineid,fh_online,fh_toline, line_with_max_halo, max_halo_so_far
  integer(kind=8) :: curHalo,curHaloid,innerHalo,block
  integer(kind=mpi_offset_kind) :: filesize
  integer(kind=4),allocatable :: lineid(:),linelinkedlist(:),haloid(:),haloperline(:),headofline(:)
  real(kind=4), allocatable :: online(:),toline(:),direction(:,:)
  integer :: n_cache,mass_index,r_index,tag,status
  real(kind=8) :: Halfbox,mass_limit,lambda,line_centre
  character(len=100) :: str_rank,z_s,str_line
  real(kind=8) :: z,d0,d_self,radius,r0,Delta_c,rho_crit_z
  integer(kind=4) :: n_elements,rlogsteps
  real(kind=8), allocatable :: distorted_dist(:),undistorted_dist(:),rev_undistorted_dist(:),comov_dist(:),rev_distorted_dist(:)
  logical :: element_flag(1:17)
  real(kind=8) :: undistorted_z,distorted_z,sumtest,rev_distorted_z,rev_undistorted_z
  integer(kind=8) :: usedlines, totusedlines
  real(kind=8) :: rmax,rmin,rminlog,rmaxlog,deltalogr
  integer :: status_checkfiles
  real(kind=8) :: nu_dist, nu_undist, nu_undist_rev
  real(kind=8),allocatable :: tau_cache(:,:),areatau_cache(:),delta_nu_cache(:)
  integer(kind=4) :: n,firstline,lastline,overlap_index
  real(kind=8) :: tau, area_tau,absorp,extend_absorp, r(0:max_size), rho(0:max_size),spherepart(0:10)
  real(kind=8) :: max_observe, min_observe, nu_min, nu_max, d_source,nu_source,source_radius,source_diameter,block_area,block_ratio
  real(kind=8) :: M0,impact_param,sigma_V,gaussian_sd,delta_nu,theta,this_absorp,tmp_tau,this_absorp_extend,point_distance
 
#ifdef DEBUG
  if(rank ==0) call system('free')
#endif
  n_cache = 100000
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
  !element_flag(14) = .TRUE.
  n_elements = get_n_elements(element_flag)

  !## Read data to node 0
  status_checkfiles = checkfiles(real(z,4))
  if(status_checkfiles == 1) then
     if(rank == 0) then
        call ascii_read_halo(real(z,4),element_flag,n_elements)
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
        print*, "halonumber=",halonumber
     end if
     call MPI_BARRIER(MPI_COMM_WORLD,ierr)
     call mpi_bcast(halonumber,1,mpi_integer,0,mpi_comm_world,ierr)
     call MPI_BARRIER(MPI_COMM_WORLD,ierr)
#ifdef DEBUG
     print*, '#DEBUG: rank',rank,"halonumber",halonumber
#endif
     if(rank == 0) allocate(halodata(1:n_elements,1:halonumber))
     call mpi_read_halo(real(z,4),element_flag,n_elements)
  else
     if(rank==0) print*, "No file to read"
     call abort
  end if
  !## End read data to node 0

  call MPI_BARRIER( MPI_COMM_WORLD,ierr)
  if(rank /= 0) allocate(halodata(1:n_elements,1:halonumber))
  call MPI_BARRIER( MPI_COMM_WORLD,ierr)
  !call abort
  if(rank ==0) print*,"Broadcasting from node 0"
  if(rank ==0) print*,halodata(1:n_elements,1)
#ifdef DEBUG
  if(rank ==0) call system('free')
  if(rank==0) print*, 'halonumber',halonumber,'n_element',n_elements
#endif
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  call mpi_bcast(halodata,halonumber*n_elements,mpi_real,0,mpi_comm_world,ierr)
  if(rank ==0) print*,"Finish broadcasting from node 0"
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
#ifdef DEBUG
  if(rank ==0) call system('free')
#endif
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
#ifdef DEBUG
  if(rank ==0) call system('free')
#endif
  if(rank ==0) print*, 'total point',n_point
  if(rank ==0) print*,'Allocate lineid,online,linkedlist ....'

  allocate(lineid(1:n_point))
  allocate(online(1:n_point))
  allocate(toline(1:n_point))
  allocate(haloid(1:n_point))
  allocate(direction(1:3,1:n_point))
#ifdef DEBUG
  if(rank ==0) call system('free')
#endif

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
#ifdef DEBUG
  if(rank ==0) call system('free')
#endif

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
#ifdef DEBUG
  if(rank ==0) call system('free')
#endif

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
#ifdef DEBUG
  if(rank ==0) call system('free')
#endif

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
#ifdef DEBUG
  if(rank ==0) call system('free')
#endif
  ! make linked list
  if(rank ==0) then
     print*,'Making linkedlist ....'
     print*,'... Allocating.....'
  endif

  allocate(distorted_dist(1:n_point))
  allocate(rev_undistorted_dist(1:n_point))
  allocate(rev_distorted_dist(1:n_point))
  allocate(comov_dist(1:n_point))
  allocate(undistorted_dist(1:n_point))
  allocate(linelinkedlist(1:n_point))
  allocate(headofline(1:max_line))
  !allocate(haloperline(1:max_line))
#ifdef DEBUG
  if(rank ==0) call system('free')
#endif
  if(rank ==0) then
     print*,'Set varaiables to 0 ...'
  endif
  !haloperline(1:max_line) = 0
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

     rev_undistorted_dist(i) = d0+convert_length2physical(real(Boxsize*line_length_factor-online(i),8),z) - convert_length2physical(real(Boxsize*line_length_factor/2.,8),z)
     rev_undistorted_z = d_to_z(rev_undistorted_dist(i))
     rev_distorted_z = undistorted_z+convert_vel2physical(real(dotproduct(direction(1:3,i)*-1.,halodata(1:3,curHaloid)),8),z)/c
     rev_distorted_dist(i) = z_to_d(rev_distorted_z)
     !nu_dist = d_to_nu(distorted_dist(i))
     !nu_undist = d_to_nu(undistorted_dist(i))

     !@ Use only minihalos
     M0 = convert_mass2physical(real(halodata(4,curHaloid),8))/M_sol
#ifndef INCLUDEPROTOGALACTIC
     if(M0 > 1.e5 .and. M0 < 1.e8) then
#endif
        linelinkedlist(i) = headofline(lineid(i))
        headofline(lineid(i)) = i
        !haloperline(lineid(i)) = haloperline(lineid(i)) + 1
#ifndef INCLUDEPROTOGALACTIC
     end if
#endif
  end do
  deallocate(direction)
  !call n_cal(n,r,rho)
  
  max_observe = d0 + convert_length2physical(real(Boxsize*line_length_factor/2.,8),z) 
  nu_min = d_to_nu(max_observe)
  min_observe = d0 - convert_length2physical(real(Boxsize*line_length_factor/2.,8),z) 
  nu_max = d_to_nu(min_observe)
  print*, 'nu_max=',nu_max,'numin=',nu_min

#ifndef RR
  d_source = d0 + convert_length2physical(real(Boxsize*(real(line_length_factor)-0.5),8),z) 
  nu_source = d_to_nu(d_source)
#endif

#ifdef DEBUG
  if(rank==0) then
     print*, 'Start simulating observe data'
  end if
#endif
  
  call n_cal(z,n,r,rho)
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  call make_spherecut(n,r,rho,spherepart)
  if(rank==0) call system('free')
  allocate(tau_cache(0:max_size,1:n_cache))
  allocate(areatau_cache(1:n_cache))
  allocate(delta_nu_cache(1:n_cache))
  
  if(mod(n_cache,nodes_returned) /= 0) then
     print*, 'nodes',nodes_returned,'ncache',n_cache
     call abort
  end if
  
  do i=rank*n_cache/nodes_returned+1,(rank+1)*n_cache/nodes_returned
     M0 = (1.d8)/n_cache*i 
     print*,'rank',rank,i,M0
     call cache_tau_table(M0,z,n,r,rho,delta_nu_cache(i),tau_cache(0:max_size,i),areatau_cache(i))
  end do
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  do i=1,nodes_returned-1
     tag = i
     
     if (rank == 0) then
        print*,'tranfering tau from',i
        call mpi_recv(tau_cache(0:max_size,i*n_cache/nodes_returned+1:(i+1)*n_cache/nodes_returned),(max_size+1)*(n_cache/nodes_returned),mpi_real8, &
             i,tag,mpi_comm_world,status,ierr)
     elseif (rank == i) then
        call mpi_send(tau_cache(0:max_size,i*n_cache/nodes_returned+1:(i+1)*n_cache/nodes_returned),(max_size+1)*(n_cache/nodes_returned),mpi_real8, &
             0,tag,mpi_comm_world,ierr)
     endif
     call mpi_barrier(mpi_comm_world,ierr)

     if (rank == 0) then
        print*,'tranfering area tau from',i
        call mpi_recv(areatau_cache(i*n_cache/nodes_returned+1:(i+1)*n_cache/nodes_returned),(n_cache/nodes_returned),mpi_real8, &
             i,tag,mpi_comm_world,status,ierr)
     elseif (rank == i) then
        call mpi_send(areatau_cache(i*n_cache/nodes_returned+1:(i+1)*n_cache/nodes_returned),(n_cache/nodes_returned),mpi_real8, &
             0,tag,mpi_comm_world,ierr)
     endif
     call mpi_barrier(mpi_comm_world,ierr)

     if (rank == 0) then
        print*,'tranfering delta nu from',i
        call mpi_recv(delta_nu_cache(i*n_cache/nodes_returned+1:(i+1)*n_cache/nodes_returned),(n_cache/nodes_returned),mpi_real8, &
             i,tag,mpi_comm_world,status,ierr)
     elseif (rank == i) then
        call mpi_send(delta_nu_cache(i*n_cache/nodes_returned+1:(i+1)*n_cache/nodes_returned),(n_cache/nodes_returned),mpi_real8, &
             0,tag,mpi_comm_world,ierr)
     endif
     call mpi_barrier(mpi_comm_world,ierr)
  enddo
  call mpi_barrier(mpi_comm_world,ierr)
  if(rank==0) call system('free')
  call mpi_bcast(delta_nu_cache,n_cache,mpi_real8,0,mpi_comm_world,ierr)
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  call mpi_bcast(areatau_cache,n_cache,mpi_real8,0,mpi_comm_world,ierr)
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  call mpi_bcast(tau_cache,n_cache*(max_size+1),mpi_real8,0,mpi_comm_world,ierr)
  !call abort
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
#ifdef RR
  print*,"mkdir -p "//trim(result_path)//z_s(1:len_trim(z_s))//'/RR/0.400/'//trim(adjustl(str_rank))
  call system("mkdir -p "//trim(result_path)//z_s(1:len_trim(z_s))//'/RR/0.000/'//trim(adjustl(str_rank)))
  call system("mkdir -p "//trim(result_path)//z_s(1:len_trim(z_s))//'/RR/0.400/'//trim(adjustl(str_rank)))
  call system("mkdir -p "//trim(result_path)//z_s(1:len_trim(z_s))//'/RR/0.200/'//trim(adjustl(str_rank)))
  call system("mkdir -p "//trim(result_path)//z_s(1:len_trim(z_s))//'/RR/0.100/'//trim(adjustl(str_rank)))
  call system("mkdir -p "//trim(result_path)//z_s(1:len_trim(z_s))//'/RR/0.050/'//trim(adjustl(str_rank)))
  call system("mkdir -p "//trim(result_path)//z_s(1:len_trim(z_s))//'/RR/0.030/'//trim(adjustl(str_rank)))
  call system("mkdir -p "//trim(result_path)//z_s(1:len_trim(z_s))//'/RR/0.003/'//trim(adjustl(str_rank)))
#else
  print*,"mkdir -p "//trim(result_path)//z_s(1:len_trim(z_s))//'/RG/0.400/'//trim(adjustl(str_rank))
  call system("mkdir -p "//trim(result_path)//z_s(1:len_trim(z_s))//'/RG/0.000/'//trim(adjustl(str_rank)))
  call system("mkdir -p "//trim(result_path)//z_s(1:len_trim(z_s))//'/RG/0.400/'//trim(adjustl(str_rank)))
  call system("mkdir -p "//trim(result_path)//z_s(1:len_trim(z_s))//'/RG/0.200/'//trim(adjustl(str_rank)))
  call system("mkdir -p "//trim(result_path)//z_s(1:len_trim(z_s))//'/RG/0.100/'//trim(adjustl(str_rank)))
  call system("mkdir -p "//trim(result_path)//z_s(1:len_trim(z_s))//'/RG/0.050/'//trim(adjustl(str_rank)))
  call system("mkdir -p "//trim(result_path)//z_s(1:len_trim(z_s))//'/RG/0.030/'//trim(adjustl(str_rank)))
  call system("mkdir -p "//trim(result_path)//z_s(1:len_trim(z_s))//'/RG/0.003/'//trim(adjustl(str_rank)))
#endif


  firstline = first_l
  lastline = last_l 


  
  rho_crit_z = rho_crit_0*(lambda0+Omega_0*(z+1.)**3. &    
       +(1.-lambda0-Omega_0)*(1.+z)**2.)
  Delta_c = 18.*pi**2 


#ifndef USERHO178
  Delta_c = (etaSUS/etaTIS)**3.*Delta_c
#endif

  if(rank==0) then
     print*, 'ref',1.e5,(3.*1.e5*M_sol/(4.*pi*Delta_c*rho_crit_z))**(1./3.)
     print*, 'ref',1.e8,(3.*1.e8*M_sol/(4.*pi*Delta_c*rho_crit_z))**(1./3.)
     print*, 'ref',1.e11,(3.*1.e11*M_sol/(4.*pi*Delta_c*rho_crit_z))**(1./3.)
  end if

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)



  !@!$omp parallel private(std_cputime,str_line,fh_record,ierr,curHalo,curhaloid,nu_dist,nu_undist,M0,impact_param,mass_index,radius,r0,r_index,tau,absorp,area_tau,extend_absorp,delta_nu,this_absorp,source_diameter,source_radius,block_ratio,block_area,overlap_index,theta)
  !@!$omp do
  do i=firstline,lastline
     std_cputime = omp_get_wtime()
     write(str_line,'(i10)') i
     str_line = adjustl(str_line)
! #ifdef RR
!      print*,"rm -f "//trim(result_path)//z_s(1:len_trim(z_s))//'/RR/status.'//str_rank(1:len_trim(str_rank))
!      call system("rm -f "//trim(result_path)//z_s(1:len_trim(z_s))//'/RR/status.'//str_rank(1:len_trim(str_rank)))
!      open(29,file=trim(result_path)//z_s(1:len_trim(z_s))//'/RR/status.'//str_rank(1:len_trim(str_rank)), status='new')
! #else
!      print*,"rm -f "//trim(result_path)//z_s(1:len_trim(z_s))//'/RG/status.'//str_rank(1:len_trim(str_rank))
!      call system("rm -f "//trim(result_path)//z_s(1:len_trim(z_s))//'/RG/status.'//str_rank(1:len_trim(str_rank)))
!      open(29,file=trim(result_path)//z_s(1:len_trim(z_s))//'/RG/status.'//str_rank(1:len_trim(str_rank)), status='new')
! #endif
!      write(29,*) i
!      close(29)   
     
#ifdef DEBUG
     if(rank==0) then
        print*, 'line =',i
     end if
#endif

#ifdef RR
     call system("rm -f "//trim(result_path)//z_s(1:len_trim(z_s))//'/RR/0.000/'//trim(adjustl(str_rank))//'/sout.'//trim(adjustl(str_line)) )
     call mpi_file_open(mpi_comm_self, &
          trim(result_path)//z_s(1:len_trim(z_s))//'/RR/0.000/'//trim(adjustl(str_rank))//'/sout.'//trim(adjustl(str_line)), &
          MPI_MODE_WRONLY + MPI_MODE_CREATE, &
          mpi_info_null,fh_record(1),ierr)
     call system("rm -f "//trim(result_path)//z_s(1:len_trim(z_s))//'/RR/0.400/'//trim(adjustl(str_rank))//'/sout.'//trim(adjustl(str_line)) )
     call mpi_file_open(mpi_comm_self, &
          trim(result_path)//z_s(1:len_trim(z_s))//'/RR/0.400/'//trim(adjustl(str_rank))//'/sout.'//trim(adjustl(str_line)), &
          MPI_MODE_WRONLY + MPI_MODE_CREATE, &
          mpi_info_null,fh_record(4),ierr)
     call system("rm -f "//trim(result_path)//z_s(1:len_trim(z_s))//'/RR/0.200/'//trim(adjustl(str_rank))//'/sout.'//trim(adjustl(str_line)) )
     call mpi_file_open(mpi_comm_self, &
          trim(result_path)//z_s(1:len_trim(z_s))//'/RR/0.200/'//trim(adjustl(str_rank))//'/sout.'//trim(adjustl(str_line)), &
          MPI_MODE_WRONLY + MPI_MODE_CREATE, &
          mpi_info_null,fh_record(7),ierr)
     call system("rm -f "//trim(result_path)//z_s(1:len_trim(z_s))//'/RR/0.100/'//trim(adjustl(str_rank))//'/sout.'//trim(adjustl(str_line)) )
     call mpi_file_open(mpi_comm_self, &
          trim(result_path)//z_s(1:len_trim(z_s))//'/RR/0.100/'//trim(adjustl(str_rank))//'/sout.'//trim(adjustl(str_line)), &
          MPI_MODE_WRONLY + MPI_MODE_CREATE, &
          mpi_info_null,fh_record(10),ierr)
     call system("rm -f "//trim(result_path)//z_s(1:len_trim(z_s))//'/RR/0.050/'//trim(adjustl(str_rank))//'/sout.'//trim(adjustl(str_line)) )
     call mpi_file_open(mpi_comm_self, &
          trim(result_path)//z_s(1:len_trim(z_s))//'/RR/0.050/'//trim(adjustl(str_rank))//'/sout.'//trim(adjustl(str_line)), &
          MPI_MODE_WRONLY + MPI_MODE_CREATE, &
          mpi_info_null,fh_record(13),ierr)

     !minihalo range
     call system("rm -f "//trim(result_path)//z_s(1:len_trim(z_s))//'/RR/0.030/'//trim(adjustl(str_rank))//'/sout.'//trim(adjustl(str_line)) )
     call mpi_file_open(mpi_comm_self, &
          trim(result_path)//z_s(1:len_trim(z_s))//'/RR/0.030/'//trim(adjustl(str_rank))//'/sout.'//trim(adjustl(str_line)), &
          MPI_MODE_WRONLY + MPI_MODE_CREATE, &
          mpi_info_null,fh_record(16),ierr)

     call system("rm -f "//trim(result_path)//z_s(1:len_trim(z_s))//'/RR/0.003/'//trim(adjustl(str_rank))//'/sout.'//trim(adjustl(str_line)) )
     call mpi_file_open(mpi_comm_self, &
          trim(result_path)//z_s(1:len_trim(z_s))//'/RR/0.003/'//trim(adjustl(str_rank))//'/sout.'//trim(adjustl(str_line)), &
          MPI_MODE_WRONLY + MPI_MODE_CREATE, &
          mpi_info_null,fh_record(19),ierr)
#else
     ! point source
     call system("rm -f "//trim(result_path)//z_s(1:len_trim(z_s))//'/RG/0.000/'//trim(adjustl(str_rank))//'/sout.'//trim(adjustl(str_line)) )
     call system("rm -f "//trim(result_path)//z_s(1:len_trim(z_s))//'/RG/0.000/'//trim(adjustl(str_rank))//'/sins.'//trim(adjustl(str_line))//'.*' )
     call mpi_file_open(mpi_comm_self, &
          trim(result_path)//z_s(1:len_trim(z_s))//'/RG/0.000/'//trim(adjustl(str_rank))//'/sout.'//trim(adjustl(str_line)), &
          MPI_MODE_WRONLY + MPI_MODE_CREATE, &
          mpi_info_null,fh_record(1),ierr)
     call mpi_file_open(mpi_comm_self, &
          trim(result_path)//z_s(1:len_trim(z_s))//'/RG/0.000/'//trim(adjustl(str_rank))//'/sins.'//trim(adjustl(str_line))//'.0', &
          MPI_MODE_WRONLY + MPI_MODE_CREATE, &
          mpi_info_null,fh_record(2),ierr)
     call mpi_file_open(mpi_comm_self, &
          trim(result_path)//z_s(1:len_trim(z_s))//'/RG/0.000/'//trim(adjustl(str_rank))//'/sins.'//trim(adjustl(str_line))//'.1', &
          MPI_MODE_WRONLY + MPI_MODE_CREATE, &
          mpi_info_null,fh_record(3),ierr)

     !0.4 diameter
     call system("rm -f "//trim(result_path)//z_s(1:len_trim(z_s))//'/RG/0.400/'//trim(adjustl(str_rank))//'/sout.'//trim(adjustl(str_line)) )
     call system("rm -f "//trim(result_path)//z_s(1:len_trim(z_s))//'/RG/0.400/'//trim(adjustl(str_rank))//'/sins.'//trim(adjustl(str_line))//'.*' )
     call mpi_file_open(mpi_comm_self, &
          trim(result_path)//z_s(1:len_trim(z_s))//'/RG/0.400/'//trim(adjustl(str_rank))//'/sout.'//trim(adjustl(str_line)), &
          MPI_MODE_WRONLY + MPI_MODE_CREATE, &
          mpi_info_null,fh_record(4),ierr)
     call mpi_file_open(mpi_comm_self, &
          trim(result_path)//z_s(1:len_trim(z_s))//'/RG/0.400/'//trim(adjustl(str_rank))//'/sins.'//trim(adjustl(str_line))//'.0', &
          MPI_MODE_WRONLY + MPI_MODE_CREATE, &
          mpi_info_null,fh_record(5),ierr)
     call mpi_file_open(mpi_comm_self, &
          trim(result_path)//z_s(1:len_trim(z_s))//'/RG/0.400/'//trim(adjustl(str_rank))//'/sins.'//trim(adjustl(str_line))//'.1', &
          MPI_MODE_WRONLY + MPI_MODE_CREATE, &
          mpi_info_null,fh_record(6),ierr)

     !0.2 diameter
     call system("rm -f "//trim(result_path)//z_s(1:len_trim(z_s))//'/RG/0.200/'//trim(adjustl(str_rank))//'/sout.'//trim(adjustl(str_line)) )
     call system("rm -f "//trim(result_path)//z_s(1:len_trim(z_s))//'/RG/0.200/'//trim(adjustl(str_rank))//'/sins.'//trim(adjustl(str_line))//'.*' )
     call mpi_file_open(mpi_comm_self, &
          trim(result_path)//z_s(1:len_trim(z_s))//'/RG/0.200/'//trim(adjustl(str_rank))//'/sout.'//trim(adjustl(str_line)), &
          MPI_MODE_WRONLY + MPI_MODE_CREATE, &
          mpi_info_null,fh_record(7),ierr)
     call mpi_file_open(mpi_comm_self, &
          trim(result_path)//z_s(1:len_trim(z_s))//'/RG/0.200/'//trim(adjustl(str_rank))//'/sins.'//trim(adjustl(str_line))//'.0', &
          MPI_MODE_WRONLY + MPI_MODE_CREATE, &
          mpi_info_null,fh_record(8),ierr)
     call mpi_file_open(mpi_comm_self, &
          trim(result_path)//z_s(1:len_trim(z_s))//'/RG/0.200/'//trim(adjustl(str_rank))//'/sins.'//trim(adjustl(str_line))//'.1', &
          MPI_MODE_WRONLY + MPI_MODE_CREATE, &
          mpi_info_null,fh_record(9),ierr)

     !0.1 diameter
     call system("rm -f "//trim(result_path)//z_s(1:len_trim(z_s))//'/RG/0.100/'//trim(adjustl(str_rank))//'/sout.'//trim(adjustl(str_line)) )
     call system("rm -f "//trim(result_path)//z_s(1:len_trim(z_s))//'/RG/0.100/'//trim(adjustl(str_rank))//'/sins.'//trim(adjustl(str_line))//'.*' )
     call mpi_file_open(mpi_comm_self, &
          trim(result_path)//z_s(1:len_trim(z_s))//'/RG/0.100/'//trim(adjustl(str_rank))//'/sout.'//trim(adjustl(str_line)), &
          MPI_MODE_WRONLY + MPI_MODE_CREATE, &
          mpi_info_null,fh_record(10),ierr)
     call mpi_file_open(mpi_comm_self, &
          trim(result_path)//z_s(1:len_trim(z_s))//'/RG/0.100/'//trim(adjustl(str_rank))//'/sins.'//trim(adjustl(str_line))//'.0', &
          MPI_MODE_WRONLY + MPI_MODE_CREATE, &
          mpi_info_null,fh_record(11),ierr)
     call mpi_file_open(mpi_comm_self, &
          trim(result_path)//z_s(1:len_trim(z_s))//'/RG/0.100/'//trim(adjustl(str_rank))//'/sins.'//trim(adjustl(str_line))//'.1', &
          MPI_MODE_WRONLY + MPI_MODE_CREATE, &
          mpi_info_null,fh_record(12),ierr)
     
     
     call system("rm -f "//trim(result_path)//z_s(1:len_trim(z_s))//'/RG/0.050/'//trim(adjustl(str_rank))//'/sout.'//trim(adjustl(str_line)) )
     call system("rm -f "//trim(result_path)//z_s(1:len_trim(z_s))//'/RG/0.050/'//trim(adjustl(str_rank))//'/sins.'//trim(adjustl(str_line))//'.*' )
     call mpi_file_open(mpi_comm_self, &
          trim(result_path)//z_s(1:len_trim(z_s))//'/RG/0.050/'//trim(adjustl(str_rank))//'/sout.'//trim(adjustl(str_line)), &
          MPI_MODE_WRONLY + MPI_MODE_CREATE, &
          mpi_info_null,fh_record(13),ierr)
     call mpi_file_open(mpi_comm_self, &
          trim(result_path)//z_s(1:len_trim(z_s))//'/RG/0.050/'//trim(adjustl(str_rank))//'/sins.'//trim(adjustl(str_line))//'.0', &
          MPI_MODE_WRONLY + MPI_MODE_CREATE, &
          mpi_info_null,fh_record(14),ierr)
     call mpi_file_open(mpi_comm_self, &
          trim(result_path)//z_s(1:len_trim(z_s))//'/RG/0.050/'//trim(adjustl(str_rank))//'/sins.'//trim(adjustl(str_line))//'.1', &
          MPI_MODE_WRONLY + MPI_MODE_CREATE, &
          mpi_info_null,fh_record(15),ierr)


     call system("rm -f "//trim(result_path)//z_s(1:len_trim(z_s))//'/RG/0.030/'//trim(adjustl(str_rank))//'/sout.'//trim(adjustl(str_line)) )
     call system("rm -f "//trim(result_path)//z_s(1:len_trim(z_s))//'/RG/0.030/'//trim(adjustl(str_rank))//'/sins.'//trim(adjustl(str_line))//'.*' )
     call mpi_file_open(mpi_comm_self, &
          trim(result_path)//z_s(1:len_trim(z_s))//'/RG/0.030/'//trim(adjustl(str_rank))//'/sout.'//trim(adjustl(str_line)), &
          MPI_MODE_WRONLY + MPI_MODE_CREATE, &
          mpi_info_null,fh_record(16),ierr)
     call mpi_file_open(mpi_comm_self, &
          trim(result_path)//z_s(1:len_trim(z_s))//'/RG/0.030/'//trim(adjustl(str_rank))//'/sins.'//trim(adjustl(str_line))//'.0', &
          MPI_MODE_WRONLY + MPI_MODE_CREATE, &
          mpi_info_null,fh_record(17),ierr)
     call mpi_file_open(mpi_comm_self, &
          trim(result_path)//z_s(1:len_trim(z_s))//'/RG/0.030/'//trim(adjustl(str_rank))//'/sins.'//trim(adjustl(str_line))//'.1', &
          MPI_MODE_WRONLY + MPI_MODE_CREATE, &
          mpi_info_null,fh_record(18),ierr)


     call system("rm -f "//trim(result_path)//z_s(1:len_trim(z_s))//'/RG/0.003/'//trim(adjustl(str_rank))//'/sout.'//trim(adjustl(str_line)) )
     call system("rm -f "//trim(result_path)//z_s(1:len_trim(z_s))//'/RG/0.003/'//trim(adjustl(str_rank))//'/sins.'//trim(adjustl(str_line))//'.*' )
     call mpi_file_open(mpi_comm_self, &
          trim(result_path)//z_s(1:len_trim(z_s))//'/RG/0.003/'//trim(adjustl(str_rank))//'/sout.'//trim(adjustl(str_line)), &
          MPI_MODE_WRONLY + MPI_MODE_CREATE, &
          mpi_info_null,fh_record(19),ierr)
     call mpi_file_open(mpi_comm_self, &
          trim(result_path)//z_s(1:len_trim(z_s))//'/RG/0.003/'//trim(adjustl(str_rank))//'/sins.'//trim(adjustl(str_line))//'.0', &
          MPI_MODE_WRONLY + MPI_MODE_CREATE, &
          mpi_info_null,fh_record(20),ierr)
     call mpi_file_open(mpi_comm_self, &
          trim(result_path)//z_s(1:len_trim(z_s))//'/RG/0.003/'//trim(adjustl(str_rank))//'/sins.'//trim(adjustl(str_line))//'.1', &
          MPI_MODE_WRONLY + MPI_MODE_CREATE, &
          mpi_info_null,fh_record(21),ierr)
#endif

     curHalo = headofline(i)
     do while(curHalo /= 0)
        curHaloid = haloid(curHalo)
        nu_dist = d_to_nu(distorted_dist(curhalo))
        nu_undist = d_to_nu(undistorted_dist(curhalo))
        M0 = convert_mass2physical(real(halodata(4,curHaloid),8))/M_sol
        impact_param = convert_length2physical(real(toline(curHalo),8),z)
        
        if(M0 > 1.e5 .and. M0 < 1.e8) then

           ! call tau_cal(M0,z,impact_param,n,r,rho,radius,delta_nu,tau,area_tau)
           ! absorp = 1.d0 - exp(-1*tau)
           ! extend_absorp =  1.d0 - exp(-1*area_tau)
           ! if(rank==0) then
           !    if(impact_param <= radius) then
           !       print*, 'actual'

           !       print*, 'tau',tau
           !       print*,'area tau',area_tau
           !       print*, 'delta nu',delta_nu
           !    end if
           ! end if

           mass_index = int(M0/(1.d8/n_cache)) 
           radius  = (3.*M0*M_sol/(4.*pi*Delta_c*rho_crit_z))**(1./3.)
           r0 = radius/zeta_t
           do r_index=1,n
              if(impact_param/r0 < r(r_index)) then
                 goto 227
              endif
           end do
227        continue
           tau = tau_cache(r_index-1,mass_index) + ( impact_param/r0 - r(r_index-1) )/(r(r_index)-r(r_index-1)) * (tau_cache(r_index,mass_index) - tau_cache(r_index-1,mass_index))
           absorp = 1.d0 - exp(-1.*tau)
           area_tau = areatau_cache(mass_index) + (M0-M0/(1.d8/n_cache)*mass_index)/(1.d8/n_cache)*(areatau_cache(mass_index+1)-areatau_cache(mass_index))
           extend_absorp =  1.d0 - exp(-1*area_tau)
           delta_nu = delta_nu_cache(mass_index) + (M0-M0/(1.d8/n_cache)*mass_index)/(1.d8/n_cache)*(delta_nu_cache(mass_index+1)-delta_nu_cache(mass_index))
           ! if(rank==0) then

           !    print*,'cache'

           !    print*, 'tau',tau
           !    print*,'area tau',areatau_cache(mass_index)
           !    print*, 'delta nu',delta_nu

           ! end if
           

           !point source case
           if(impact_param < radius) then
              ! not include z-distortion
#ifdef DEBUG
              print*, int(i), real(nu_undist), real(absorp),real(gaussian_sd)
              print*, int(i), real(nu_dist), real(absorp),real(gaussian_sd)
#endif
              this_absorp = absorp

#ifdef RR
              !source far away, pass random lines
              call MPI_FILE_WRITE(fh_record(1), M0, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
              call MPI_FILE_WRITE(fh_record(1), impact_param/Mpc, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
              call MPI_FILE_WRITE(fh_record(1), nu_dist, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
              call MPI_FILE_WRITE(fh_record(1), nu_undist, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
              call MPI_FILE_WRITE(fh_record(1), this_absorp, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
              call MPI_FILE_WRITE(fh_record(1), delta_nu, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
#else
              !source far away, los pass through cluster
              
              call MPI_FILE_WRITE(fh_record(1), M0, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
              call MPI_FILE_WRITE(fh_record(1), impact_param/Mpc, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
              call MPI_FILE_WRITE(fh_record(1), nu_dist, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
              call MPI_FILE_WRITE(fh_record(1), nu_undist, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
              call MPI_FILE_WRITE(fh_record(1), this_absorp, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
              call MPI_FILE_WRITE(fh_record(1), delta_nu, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)

              !source in cluster
              if(online(curHalo) > line_length_factor/2.*Boxsize ) then
                 call MPI_FILE_WRITE(fh_record(2), M0, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
                 call MPI_FILE_WRITE(fh_record(2), impact_param/Mpc, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
                 call MPI_FILE_WRITE(fh_record(2), d_to_nu(rev_distorted_dist(curhalo)), 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
                 call MPI_FILE_WRITE(fh_record(2), d_to_nu(rev_undistorted_dist(curhalo)), 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
                 call MPI_FILE_WRITE(fh_record(2), this_absorp, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
                 call MPI_FILE_WRITE(fh_record(2), delta_nu, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
              else
                 
                 call MPI_FILE_WRITE(fh_record(3), M0, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
                 call MPI_FILE_WRITE(fh_record(3), impact_param/Mpc, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
                 call MPI_FILE_WRITE(fh_record(3), nu_dist, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
                 call MPI_FILE_WRITE(fh_record(3), nu_undist, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
                 call MPI_FILE_WRITE(fh_record(3), this_absorp, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
                 call MPI_FILE_WRITE(fh_record(3), delta_nu, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
              endif
              
#endif
           end if !point source

           !extended source case
           !0.4   

           source_diameter = 0.4 
           source_radius = convert_length2physical(source_diameter/co_boxwidth*boxsize/2 ,z)
           if(impact_param < (radius + source_radius)) then
              overlap_index = min ( int (( (radius + source_radius) - impact_param)/radius), 10)
              theta = acos( 1.- 0.2 * overlap_index)
              block_area = radius**2/2.d0 *(2.*theta -sin(2*theta))
              block_ratio = block_area/(pi*source_radius**2.)*spherepart(overlap_index)
              this_absorp = extend_absorp*block_ratio
              ! if(rank == 0) then !
              !    print*,"index",overlap_index
              !    print*,"area index",block_area/(pi*radius**2.)
              !    print*,"block area",block_area/(pi*source_radius**2.)
              !    print*, "block factor",spherepart(overlap_index)
              !    print*,'block ratio', block_ratio
              !    print*,"absorp",extend_absorp*block_ratio
              ! end if
#ifdef RR
              !source far away, pass random lines
              call MPI_FILE_WRITE(fh_record(4), M0, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
              call MPI_FILE_WRITE(fh_record(4), impact_param/Mpc, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
              call MPI_FILE_WRITE(fh_record(4), nu_dist, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
              call MPI_FILE_WRITE(fh_record(4), nu_undist, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
              call MPI_FILE_WRITE(fh_record(4), this_absorp, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
              call MPI_FILE_WRITE(fh_record(4), delta_nu, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
#else
              !source far away, los pass through cluster

              call MPI_FILE_WRITE(fh_record(4), M0, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
              call MPI_FILE_WRITE(fh_record(4), impact_param/Mpc, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
              call MPI_FILE_WRITE(fh_record(4), nu_dist, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
              call MPI_FILE_WRITE(fh_record(4), nu_undist, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
              call MPI_FILE_WRITE(fh_record(4), this_absorp, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
              call MPI_FILE_WRITE(fh_record(4), delta_nu, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)

              !source in cluster
              if(online(curHalo) > line_length_factor/2.*Boxsize ) then
                 call MPI_FILE_WRITE(fh_record(5), M0, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
                 call MPI_FILE_WRITE(fh_record(5), impact_param/Mpc, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
                 call MPI_FILE_WRITE(fh_record(5), d_to_nu(rev_distorted_dist(curhalo)), 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
                 call MPI_FILE_WRITE(fh_record(5), d_to_nu(rev_undistorted_dist(curhalo)), 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
                 call MPI_FILE_WRITE(fh_record(5), this_absorp, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
                 call MPI_FILE_WRITE(fh_record(5), delta_nu, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
              else

                 call MPI_FILE_WRITE(fh_record(6), M0, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
                 call MPI_FILE_WRITE(fh_record(6), impact_param/Mpc, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
                 call MPI_FILE_WRITE(fh_record(6), nu_dist, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
                 call MPI_FILE_WRITE(fh_record(6), nu_undist, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
                 call MPI_FILE_WRITE(fh_record(6), this_absorp, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
                 call MPI_FILE_WRITE(fh_record(6), delta_nu, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
              endif

#endif
           end if ! extend source 0.4


           !extended source case
           !0.2
           source_diameter = 0.2 
           source_radius = convert_length2physical(source_diameter/co_boxwidth*boxsize/2 ,z)
           if(impact_param < (radius + source_radius)) then
              overlap_index = min ( int (( (radius + source_radius) - impact_param)/radius), 10)
              theta = acos( 1.- 0.2 * overlap_index)
              block_area = radius**2/2.d0 *(2.*theta -sin(2*theta))
              block_ratio = block_area/(pi*source_radius**2.)*spherepart(overlap_index)
              this_absorp = extend_absorp*block_ratio
              ! if(rank == 0) then
              !    print*,"index",overlap_index
              !    print*,"area index",block_area/(pi*radius**2.)
              !    print*,"block area",block_area/(pi*source_radius**2.)
              !    print*, "block factor",spherepart(overlap_index)
              !    print*,'block ratio', block_ratio
              !    print*,"absorp",extend_absorp*block_ratio
              ! end if
#ifdef RR
              !source far away, pass random lines
              call MPI_FILE_WRITE(fh_record(7), M0, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
              call MPI_FILE_WRITE(fh_record(7), impact_param/Mpc, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
              call MPI_FILE_WRITE(fh_record(7), nu_dist, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
              call MPI_FILE_WRITE(fh_record(7), nu_undist, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
              call MPI_FILE_WRITE(fh_record(7), this_absorp, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
              call MPI_FILE_WRITE(fh_record(7), delta_nu, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
#else
              !source far away, los pass through cluster

              call MPI_FILE_WRITE(fh_record(7), M0, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
              call MPI_FILE_WRITE(fh_record(7), impact_param/Mpc, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
              call MPI_FILE_WRITE(fh_record(7), nu_dist, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
              call MPI_FILE_WRITE(fh_record(7), nu_undist, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
              call MPI_FILE_WRITE(fh_record(7), this_absorp, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
              call MPI_FILE_WRITE(fh_record(7), delta_nu, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)

              !source in cluster
              if(online(curHalo) > line_length_factor/2.*Boxsize ) then
                 call MPI_FILE_WRITE(fh_record(8), M0, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
                 call MPI_FILE_WRITE(fh_record(8), impact_param/Mpc, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
                 call MPI_FILE_WRITE(fh_record(8), d_to_nu(rev_distorted_dist(curhalo)), 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
                 call MPI_FILE_WRITE(fh_record(8), d_to_nu(rev_undistorted_dist(curhalo)), 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
                 call MPI_FILE_WRITE(fh_record(8), this_absorp, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
                 call MPI_FILE_WRITE(fh_record(8), delta_nu, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
              else

                 call MPI_FILE_WRITE(fh_record(9), M0, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
                 call MPI_FILE_WRITE(fh_record(9), impact_param/Mpc, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
                 call MPI_FILE_WRITE(fh_record(9), nu_dist, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
                 call MPI_FILE_WRITE(fh_record(9), nu_undist, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
                 call MPI_FILE_WRITE(fh_record(9), this_absorp, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
                 call MPI_FILE_WRITE(fh_record(9), delta_nu, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
              endif

#endif
           end if ! extend source 0.2


           !extended source case
           !0.1
           source_diameter = 0.1
           source_radius = convert_length2physical(source_diameter/co_boxwidth*boxsize/2 ,z)
           if(impact_param < (radius + source_radius)) then
              overlap_index = min ( int (( (radius + source_radius) - impact_param)/radius), 10)
              theta = acos( 1.- 0.2 * overlap_index)
              block_area = radius**2/2.d0 *(2.*theta -sin(2*theta))
              block_ratio = block_area/(pi*source_radius**2.)*spherepart(overlap_index)
              this_absorp = extend_absorp*block_ratio
              ! if(rank == 0) then
              !    print*,"index",overlap_index
              !    print*,"area index",block_area/(pi*radius**2.)
              !    print*,"block area",block_area/(pi*source_radius**2.)
              !    print*, "block factor",spherepart(overlap_index)
              !    print*,'block ratio', block_ratio
              !    print*,"absorp",extend_absorp*block_ratio
              ! end if
#ifdef RR
              !source far away, pass random lines
              call MPI_FILE_WRITE(fh_record(10), M0, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
              call MPI_FILE_WRITE(fh_record(10), impact_param/Mpc, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
              call MPI_FILE_WRITE(fh_record(10), nu_dist, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
              call MPI_FILE_WRITE(fh_record(10), nu_undist, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
              call MPI_FILE_WRITE(fh_record(10), this_absorp, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
              call MPI_FILE_WRITE(fh_record(10), delta_nu, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
#else
              !source far away, los pass through cluster

              call MPI_FILE_WRITE(fh_record(10), M0, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
              call MPI_FILE_WRITE(fh_record(10), impact_param/Mpc, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
              call MPI_FILE_WRITE(fh_record(10), nu_dist, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
              call MPI_FILE_WRITE(fh_record(10), nu_undist, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
              call MPI_FILE_WRITE(fh_record(10), this_absorp, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
              call MPI_FILE_WRITE(fh_record(10), delta_nu, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)

              !source in cluster
              if(online(curHalo) > line_length_factor/2.*Boxsize ) then
                 call MPI_FILE_WRITE(fh_record(11), M0, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
                 call MPI_FILE_WRITE(fh_record(11), impact_param/Mpc, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
                 call MPI_FILE_WRITE(fh_record(11), d_to_nu(rev_distorted_dist(curhalo)), 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
                 call MPI_FILE_WRITE(fh_record(11), d_to_nu(rev_undistorted_dist(curhalo)), 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
                 call MPI_FILE_WRITE(fh_record(11), this_absorp, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
                 call MPI_FILE_WRITE(fh_record(11), delta_nu, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
              else

                 call MPI_FILE_WRITE(fh_record(12), M0, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
                 call MPI_FILE_WRITE(fh_record(12), impact_param/Mpc, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
                 call MPI_FILE_WRITE(fh_record(12), nu_dist, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
                 call MPI_FILE_WRITE(fh_record(12), nu_undist, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
                 call MPI_FILE_WRITE(fh_record(12), this_absorp, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
                 call MPI_FILE_WRITE(fh_record(12), delta_nu, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
              endif

#endif
           end if ! extend source 0.1

           !extended source case
           !0.05
           source_diameter = 0.05
           source_radius = convert_length2physical(source_diameter/co_boxwidth*boxsize/2 ,z)
           if(impact_param < (radius + source_radius)) then
              overlap_index = min ( int (( (radius + source_radius) - impact_param)/radius), 10)
              theta = acos( 1.- 0.2 * overlap_index)
              block_area = radius**2/2.d0 *(2.*theta -sin(2*theta))
              block_ratio = block_area/(pi*source_radius**2.)*spherepart(overlap_index)
              this_absorp = extend_absorp*block_ratio
              ! if(rank == 0) then
              !    print*,"index",overlap_index
              !    print*,"area index",block_area/(pi*radius**2.)
              !    print*,"block area",block_area/(pi*source_radius**2.)
              !    print*, "block factor",spherepart(overlap_index)
              !    print*,'block ratio', block_ratio
              !    print*,"absorp",extend_absorp*block_ratio
              ! end if
#ifdef RR
              !source far away, pass random lines
              call MPI_FILE_WRITE(fh_record(13), M0, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
              call MPI_FILE_WRITE(fh_record(13), impact_param/Mpc, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
              call MPI_FILE_WRITE(fh_record(13), nu_dist, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
              call MPI_FILE_WRITE(fh_record(13), nu_undist, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
              call MPI_FILE_WRITE(fh_record(13), this_absorp, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
              call MPI_FILE_WRITE(fh_record(13), delta_nu, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
#else
              !source far away, los pass through cluster

              call MPI_FILE_WRITE(fh_record(13), M0, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
              call MPI_FILE_WRITE(fh_record(13), impact_param/Mpc, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
              call MPI_FILE_WRITE(fh_record(13), nu_dist, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
              call MPI_FILE_WRITE(fh_record(13), nu_undist, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
              call MPI_FILE_WRITE(fh_record(13), this_absorp, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
              call MPI_FILE_WRITE(fh_record(13), delta_nu, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)

              !source in cluster
              if(online(curHalo) > line_length_factor/2.*Boxsize ) then
                 call MPI_FILE_WRITE(fh_record(14), M0, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
                 call MPI_FILE_WRITE(fh_record(14), impact_param/Mpc, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
                 call MPI_FILE_WRITE(fh_record(14), d_to_nu(rev_distorted_dist(curhalo)), 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
                 call MPI_FILE_WRITE(fh_record(14), d_to_nu(rev_undistorted_dist(curhalo)), 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
                 call MPI_FILE_WRITE(fh_record(14), this_absorp, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
                 call MPI_FILE_WRITE(fh_record(14), delta_nu, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
              else

                 call MPI_FILE_WRITE(fh_record(15), M0, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
                 call MPI_FILE_WRITE(fh_record(15), impact_param/Mpc, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
                 call MPI_FILE_WRITE(fh_record(15), nu_dist, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
                 call MPI_FILE_WRITE(fh_record(15), nu_undist, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
                 call MPI_FILE_WRITE(fh_record(15), this_absorp, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
                 call MPI_FILE_WRITE(fh_record(15), delta_nu, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
              endif

#endif
           end if ! extend source 0.05

           !extended source case
           !0.02
           source_diameter = 0.03
           source_radius = convert_length2physical(source_diameter/co_boxwidth*boxsize/2 ,z)
           if(impact_param < (radius + source_radius)) then
              overlap_index = min ( int (( (radius + source_radius) - impact_param)/radius), 10)
              theta = acos( 1.- 0.2 * overlap_index)
              block_area = radius**2/2.d0 *(2.*theta -sin(2*theta))
              block_ratio = block_area/(pi*source_radius**2.)*spherepart(overlap_index)
              this_absorp = extend_absorp*block_ratio
              ! if(rank == 0) then
              !    print*,"index",overlap_index
              !    print*,"area index",block_area/(pi*radius**2.)
              !    print*,"block area",block_area/(pi*source_radius**2.)
              !    print*, "block factor",spherepart(overlap_index)
              !    print*,'block ratio', block_ratio
              !    print*,"absorp",extend_absorp*block_ratio
              ! end if
#ifdef RR
              !source far away, pass random lines
              call MPI_FILE_WRITE(fh_record(16), M0, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
              call MPI_FILE_WRITE(fh_record(16), impact_param/Mpc, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
              call MPI_FILE_WRITE(fh_record(16), nu_dist, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
              call MPI_FILE_WRITE(fh_record(16), nu_undist, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
              call MPI_FILE_WRITE(fh_record(16), this_absorp, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
              call MPI_FILE_WRITE(fh_record(16), delta_nu, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
#else
              !source far away, los pass through cluster

              call MPI_FILE_WRITE(fh_record(16), M0, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
              call MPI_FILE_WRITE(fh_record(16), impact_param/Mpc, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
              call MPI_FILE_WRITE(fh_record(16), nu_dist, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
              call MPI_FILE_WRITE(fh_record(16), nu_undist, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
              call MPI_FILE_WRITE(fh_record(16), this_absorp, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
              call MPI_FILE_WRITE(fh_record(16), delta_nu, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)

              !source in cluster
              if(online(curHalo) > line_length_factor/2.*Boxsize ) then
                 call MPI_FILE_WRITE(fh_record(17), M0, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
                 call MPI_FILE_WRITE(fh_record(17), impact_param/Mpc, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
                 call MPI_FILE_WRITE(fh_record(17), d_to_nu(rev_distorted_dist(curhalo)), 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
                 call MPI_FILE_WRITE(fh_record(17), d_to_nu(rev_undistorted_dist(curhalo)), 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
                 call MPI_FILE_WRITE(fh_record(17), this_absorp, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
                 call MPI_FILE_WRITE(fh_record(17), delta_nu, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
              else

                 call MPI_FILE_WRITE(fh_record(18), M0, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
                 call MPI_FILE_WRITE(fh_record(18), impact_param/Mpc, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
                 call MPI_FILE_WRITE(fh_record(18), nu_dist, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
                 call MPI_FILE_WRITE(fh_record(18), nu_undist, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
                 call MPI_FILE_WRITE(fh_record(18), this_absorp, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
                 call MPI_FILE_WRITE(fh_record(18), delta_nu, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
              endif

#endif
           end if ! extend source 0.020
           
           !extensource minihalo 0.003
           source_diameter = 0.003
           source_radius = convert_length2physical(source_diameter/co_boxwidth*boxsize/2 ,z)
           if(impact_param < (radius + source_radius)) then
              tmp_tau = 0.d0
              !find average
              do k=-10,10
                 point_distance = abs(impact_param + k*0.1*source_radius)
                 do r_index=1,n
                    if(point_distance/r0 < r(r_index)) then
                       goto 228
                    endif
                 end do
228              continue
                 tmp_tau = tmp_tau + tau_cache(r_index-1,mass_index) + (point_distance/r0 - r(r_index-1) )/(r(r_index)-r(r_index-1)) * (tau_cache(r_index,mass_index) - tau_cache(r_index-1,mass_index))
                             
              end do
              tau = tmp_tau/21.
              this_absorp = 1.d0 - exp(-1.*tau)
              ! if(rank == 0) then
              !    print*,"index",overlap_index
              !    print*,"area index",block_area/(pi*radius**2.)
              !    print*,"block area",block_area/(pi*source_radius**2.)
              !    print*, "block factor",spherepart(overlap_index)
              !    print*,'block ratio', block_ratio
              !    print*,"absorp",extend_absorp*block_ratio
              ! end if
#ifdef RR
              !source far away, pass random lines
              call MPI_FILE_WRITE(fh_record(19), M0, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
              call MPI_FILE_WRITE(fh_record(19), impact_param/Mpc, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
              call MPI_FILE_WRITE(fh_record(19), nu_dist, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
              call MPI_FILE_WRITE(fh_record(19), nu_undist, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
              call MPI_FILE_WRITE(fh_record(19), this_absorp, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
              call MPI_FILE_WRITE(fh_record(19), delta_nu, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
#else
              !source far away, los pass through cluster

              call MPI_FILE_WRITE(fh_record(19), M0, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
              call MPI_FILE_WRITE(fh_record(19), impact_param/Mpc, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
              call MPI_FILE_WRITE(fh_record(19), nu_dist, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
              call MPI_FILE_WRITE(fh_record(19), nu_undist, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
              call MPI_FILE_WRITE(fh_record(19), this_absorp, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
              call MPI_FILE_WRITE(fh_record(19), delta_nu, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)

              !source in cluster
              if(online(curHalo) > line_length_factor/2.*Boxsize ) then
                 call MPI_FILE_WRITE(fh_record(20), M0, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
                 call MPI_FILE_WRITE(fh_record(20), impact_param/Mpc, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
                 call MPI_FILE_WRITE(fh_record(20), d_to_nu(rev_distorted_dist(curhalo)), 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
                 call MPI_FILE_WRITE(fh_record(20), d_to_nu(rev_undistorted_dist(curhalo)), 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
                 call MPI_FILE_WRITE(fh_record(20), this_absorp, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
                 call MPI_FILE_WRITE(fh_record(20), delta_nu, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
              else

                 call MPI_FILE_WRITE(fh_record(21), M0, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
                 call MPI_FILE_WRITE(fh_record(21), impact_param/Mpc, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
                 call MPI_FILE_WRITE(fh_record(21), nu_dist, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
                 call MPI_FILE_WRITE(fh_record(21), nu_undist, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
                 call MPI_FILE_WRITE(fh_record(21), this_absorp, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
                 call MPI_FILE_WRITE(fh_record(21), delta_nu, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
              endif

#endif
           end if


        end if

        !gaussian_sd = sqrt(2.)*sigma_V/c*nu_undist
        !if(rank==0) write(53,*) int(i), real(nu_undist), real(absorp),real(extend_absorp),real(gaussian_sd)
        
        !@ next halo
        curHalo = linelinkedlist(curHalo)
     end do
     print*, 'rank',rank, 'line',i,omp_get_wtime() - std_cputime, 's' 
     !call system("free")
     do k=1,21
        call MPI_FILE_CLOSE(fh_record(k), ierr)
     end do
     print*, 'close files'
  end do
  !@!$omp end do
  !@!$omp end parallel

  print*,'rank',rank,'complete all lines'
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  deallocate(halodata,online,toline, &
       lineid,haloid,distorted_dist,rev_distorted_dist, &
       comov_dist,undistorted_dist,rev_undistorted_dist, &
       linelinkedlist,headofline)
  deallocate(tau_cache,delta_nu_cache,areatau_cache)
call MPI_BARRIER(MPI_COMM_WORLD,ierr)
#ifdef RR
end subroutine makeobservedlines_rr
#else
end subroutine makeobservedlines_rg
#endif



#ifndef GETCHECKLINE
#define GETCHECKLINE
function get_checkline(z)
  use mpi
  use mpitools
  use common_vars
  implicit none
  logical :: file_e
  real(kind=8) :: z
  integer :: get_checkline
  character(len=100) :: z_s,str_rank
  
  write(z_s,'(f10.3)') z
  z_s = adjustl(z_s)
  write(str_rank,'(i5)') rank
  str_rank = adjustl(str_rank)
#ifdef RR
  inquire(file=trim(result_path)//z_s(1:len_trim(z_s))//'/RR/status.'//str_rank(1:len_trim(str_rank)), exist=file_e )
#else
  inquire(file=trim(result_path)//z_s(1:len_trim(z_s))//'/RG/status.'//str_rank(1:len_trim(str_rank)), exist=file_e )
#endif
  if(file_e == .False.) then
#ifdef RR
     print*, "Status file "//trim(result_path)//z_s(1:len_trim(z_s))//'/RR/status.'//str_rank(1:len_trim(str_rank))//" does not exist"
#else
     print*, "Status file "//trim(result_path)//z_s(1:len_trim(z_s))//'/RG/status.'//str_rank(1:len_trim(str_rank))//" does not exist"
#endif
     get_checkline = first_l
  else
     print*, "Status file ",rank,"exists open the file"
#ifdef RR
     open(39,file=trim(result_path)//z_s(1:len_trim(z_s))//'/RR/status.'//str_rank(1:len_trim(str_rank)), status='old')
#else
     open(39,file=trim(result_path)//z_s(1:len_trim(z_s))//'/RG/status.'//str_rank(1:len_trim(str_rank)), status='old')
#endif
    read(39,*) get_checkline
    close(39)
  end if
  get_checkline = max(get_checkline,first_l)
  return
end function get_checkline
#endif
