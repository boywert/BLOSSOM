subroutine call_subcell_corr(z)
  use common_vars
  implicit none
  integer :: i
  real(kind=8) :: z
  call read_subcelllist()
  do i=1,SUBCELLLIST_MAXSIZE
     if(subcell_perdim_list(i) > 0) then
        call subcell_corr(z,subcell_perdim_list(i))
     end if
  end do
end subroutine call_subcell_corr

subroutine subcell_corr(z, GridLines)
  use omp_lib
  use vectortools
  use arraytools
  use datatools
  use ifport
  use mpi
  use mpitools
  use io_tools
  use conversiontools
  use common_vars
  use utilities_serial
  implicit none   

  integer (kind = 4)   :: n_threads, omp_thread
  integer(kind=4)   :: omp
  integer(kind=4) :: GridLines 
  integer(kind=4):: PrGridLines 
  real(kind=4) :: GridSize

  real(kind=4) :: redshift_index
  real(kind=4) :: overden_min,overden_max,overden_binsize
  real(kind=8) :: sum_overden
  integer(kind=4) :: overden_block
  integer(kind=4) :: mass_block,NumBlock,curHalo,unitid,totalhaloline !curHalo => 8byte
 

  real(kind=4) :: shiftdistance(3),pos(3),block_dummy(3)

  integer(kind=4) :: i,j,k,l,m,n,xStart,xStop,o,linenum,hitnum,halonumber
  integer(kind=4) :: length,curBox,box1(27),box2(27),mixbox(54)
 
  integer(kind=4) :: timearray(8),command
  real(kind=4),allocatable :: halodata(:,:), positions(:,:),radius(:),mass(:)
  integer(kind=4), allocatable :: headofchain(:,:),linkedlist(:),hitperthread(:),missperthread(:)
  real(kind=4) ,allocatable :: local_positions(:,:),inner_positions(:,:)
  real(kind=4) :: density_centre(1:overden_nbin,1:3), density_halfbin(1:overden_nbin,1:3)
  character(len=100) :: str_rank,filename
  character(len=20) :: z_s
  integer(kind=mpi_offset_kind) :: filesize
  integer(kind=4) :: n_point,totalpoint,tag,n_point_index,den_block
  integer(kind=4),allocatable :: n_point_arr(:)
  real(kind=4), allocatable :: halo_in(:,:),subcell_overden(:)
  integer(4), dimension(mpi_status_size) :: status
  integer(kind=4), allocatable :: file_read_per_node(:),n_point_file(:)
  logical :: file_e, element_flag(1:17)
  real(kind=8) :: z,d0,error_delta
  integer(kind=4) :: n_elements, haloperline,count_halo,boydID,status_checkfiles
  integer(kind=4), allocatable :: count_haloincell(:,:), count_cellbybin(:,:)
  real(kind=4) :: xi_out(1:n_bin_xi),local_overden
  real(kind=4) :: rad(0:n_bin_xi),drlog,r1,r2,dv
  real(kind=8), allocatable :: xi_mean(:,:,:),xi_mean_sq(:,:,:), final_xi_mean(:,:,:), final_xi_mean_sq(:,:,:)
  real(kind=4), allocatable :: rho_mean(:,:), delta_cell(:,:)
  integer(kind=4),allocatable :: count_cell_xi(:,:), cal_cell_per_node(:), density_block(:,:)
  real(kind=4) :: stop_time, start_time
  character(len=20) :: grid_s, boxwidth_s, mass_s, delta_s
  integer(kind=8) :: count_for_delta(1:n_bin_xi)

  if(rank == 0) call cpu_time(start_time)
  if(rank == 0) print*, "subgrind",z,gridlines

  status_checkfiles = checkfiles(real(z,4))

  element_flag(1:17) = .FALSE.
  element_flag(1:3) = .TRUE. ! 1:3 => positions
  element_flag(15) = .TRUE. ! 14:15 =>   mass
  n_elements = get_n_elements(element_flag)

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





  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  write(z_s,'(f10.3)') z
  z_s = adjustl(z_s)

  write(str_rank,'(i10)') rank
  omp_thread = omp_get_max_threads()

  write(str_rank,'(i10)') rank
  str_rank = adjustl(str_rank)


  write(grid_s,'(i5)') GridLines
  grid_s = adjustl(grid_s)
  write(boxwidth_s,'(i5)') int(co_boxwidth)
  boxwidth_s = adjustl(boxwidth_s)


  if (rank==0) then
     print*, 'finish reading'
     call cpu_time(stop_time)
     print*, "Use",stop_time-start_time,"seconds to read data"
     print*,'MPI node:',nodes_returned
     print*,'OMP node:',omp_thread

     write(*,*) 'Start using Node 0 to read and transfer the essentials to others'


     !allocate arrays

     allocate(positions(1:3,1:halonumber))

     allocate(mass(1:halonumber))
     positions(1:3,1:halonumber) = halodata(1:3,1:halonumber)

     mass(1:halonumber) = halodata(4,1:halonumber)
     deallocate(halodata)
  
     !GridSize = Boxsize/6.

     Gridsize = BoxSize/GridLines

     print*, 'Gridsize: ',Gridsize
     print*, 'GridLines: ',GridLines
     print*, 'Subvolume:',0,'-',GridLines**3-1

     allocate(local_positions(1:3,0:GridLines**3-1))

     do i = 0,Gridlines-1
        do j = 0,Gridlines-1
           do k = 0,Gridlines-1
              NumBlock = k*GridLines**2 + j*GridLines + i
              !print*,i,j,k,NumBlock
              local_positions(1:3,NumBlock) = (/ i,j,k /) * Gridsize
           end do
        end do
     end do
     
     
 
     !overden_min = -10.05
     !overden_max =  10.05
     !overden_binsize = 0.10
     !overden_nbin = ceiling((overden_max - overden_min)/overden_binsize)
     call read_data_bin_info(GridLines, density_centre, density_halfbin)

     allocate(subcell_overden(0:GridLines**3-1))
     call read_subcell_overden(real(z,4),GridLines,subcell_overden)

     allocate(headofchain(1:3,0:Gridlines**3-1))
     allocate(count_haloincell(1:3,0:Gridlines**3-1))
     allocate(linkedlist(1:halonumber))


     !set the initials for linked listing
     headofchain(1:3,0:Gridlines**3-1) = 0
     count_haloincell(1:3,0:Gridlines**3-1) = 0
     linkedlist(1:halonumber) = 0

     !put halos in small cells
     write(*,*) 'Putting halos in cells'
     do i=1,halonumber 

        mass(i) = convert_mass2physical(real(mass(i),8))/M_sol
        mass_block = return_massbin_number(mass(i))
        !print*, mass(i), mass_block

        block_dummy(1) = floor(positions(1,i)/GridSize)
        if(block_dummy(1) == GridLines) block_dummy(1) = GridLines-1
        block_dummy(2) = floor(positions(2,i)/GridSize)
        if(block_dummy(2) == GridLines) block_dummy(2) = GridLines-1
        block_dummy(3) = floor(positions(3,i)/GridSize)
        if(block_dummy(3) == GridLines) block_dummy(3) = GridLines-1

        NumBlock = block_dummy(3)*GridLines**2 + &
             block_dummy(2)*GridLines + &
             block_dummy(1)
#ifdef DEBUG
        if(NumBlock < 0 .or. NumBlock > Gridlines**3-1) then
           print*,'Error: subvolume',NumBlock,'Invalid'
           print*,'Order:',i
           print*, 'Coordinate:'
           print*, floor(positions(1,i)/GridSize),positions(1,i)
           print*, floor(positions(2,i)/GridSize),positions(2,i)
           print*, floor(positions(3,i)/GridSize),positions(3,i)
           !goto 145
           call abort
        endif
#endif
        count_haloincell(mass_block,Numblock) =  count_haloincell(mass_block,Numblock) +1
        linkedlist(i) = headofchain(mass_block,Numblock)
        headofchain(mass_block,Numblock) = i
        
     end do
     deallocate(mass)
     
     write(*,*) 'Making linkedlist'
     do i=0,GridLines**3-1
        do j=1,3
           curHalo = headofchain(j,i)
           do while(curHalo /= 0) 
#ifdef DEBUG
              if(curHalo < 0 .or. curHalo > halonumber) then
                 print*, 'PRE-check: before bcast'
                 print*,'RANK',rank,'Error: halo',curHalo,'Invalid'
              endif
#endif
              curHalo = linkedlist(curhalo)
           end do

        end do
     end do
     write(*,*) 'Broadcasting ......'
     call cpu_time(start_time)
  endif



  call mpi_barrier(mpi_comm_world,ierr)
  !broadcast halonumber
  call mpi_bcast(halonumber,1,mpi_integer,0,mpi_comm_world,ierr)
  call mpi_barrier(mpi_comm_world,ierr)
  call mpi_bcast(GridLines,1,mpi_integer,0,mpi_comm_world,ierr)
  call mpi_barrier(mpi_comm_world,ierr)
  call mpi_bcast(GridSize,1,mpi_real,0,mpi_comm_world,ierr)
  call mpi_barrier(mpi_comm_world,ierr)

  call mpi_bcast(density_centre,overden_nbin*3,mpi_real,0,mpi_comm_world,ierr)
  call mpi_barrier(mpi_comm_world,ierr)
  call mpi_bcast(density_halfbin,overden_nbin*3,mpi_real,0,mpi_comm_world,ierr)
  call mpi_barrier(mpi_comm_world,ierr)
  !call mpi_bcast(overden_nbin,1,mpi_integer,0,mpi_comm_world,ierr)
  !call mpi_barrier(mpi_comm_world,ierr)
  !call mpi_bcast(overden_binsize,1,mpi_real,0,mpi_comm_world,ierr)
  !call mpi_barrier(mpi_comm_world,ierr)
  !call mpi_bcast(overden_max,1,mpi_real,0,mpi_comm_world,ierr)
  !call mpi_barrier(mpi_comm_world,ierr)
  !call mpi_bcast(overden_min,1,mpi_real,0,mpi_comm_world,ierr)
  !call mpi_barrier(mpi_comm_world,ierr)
  if(rank==1) print*,density_centre
  if(rank==1) print*, ''
  if(rank==1) print*,density_halfbin
  !print*, "end"
  !call abort
  !allocate array in other nodes
  if (rank /= 0) then
     allocate(positions(1:3,1:halonumber))
     !allocate(mass(1:halonumber))
     allocate(headofchain(1:3,0:Gridlines**3-1))
     allocate(count_haloincell(1:3,0:Gridlines**3-1))
     allocate(linkedlist(1:halonumber))
     allocate(local_positions(1:3,0:GridLines**3-1))
     allocate(subcell_overden(0:GridLines**3-1))
  endif
  call mpi_barrier(mpi_comm_world,ierr)

  !broadcast variables
  call mpi_bcast(positions,halonumber*3,mpi_real,0,mpi_comm_world,ierr)
  if (ierr /= mpi_success) call mpi_abort(mpi_comm_world,ierr,ierr)

  call mpi_bcast(local_positions,3*GridLines**3,mpi_real,0,mpi_comm_world,ierr)
  if (ierr /= mpi_success) call mpi_abort(mpi_comm_world,ierr,ierr)

  !call mpi_bcast(mass,halonumber,mpi_real,0,mpi_comm_world,ierr)
  !if (ierr /= mpi_success) call mpi_abort(mpi_comm_world,ierr,ierr)


  call mpi_bcast(linkedlist,halonumber,mpi_integer,0,mpi_comm_world,ierr)
  if (ierr /= mpi_success) call mpi_abort(mpi_comm_world,ierr,ierr)

  call mpi_bcast(headofchain,3*Gridlines**3,mpi_integer,0,mpi_comm_world,ierr)
  if (ierr /= mpi_success) call mpi_abort(mpi_comm_world,ierr,ierr)

  call mpi_bcast(count_haloincell,3*Gridlines**3,mpi_integer,0,mpi_comm_world,ierr)
  if (ierr /= mpi_success) call mpi_abort(mpi_comm_world,ierr,ierr)

  call mpi_bcast(Subcell_overden,Gridlines**3,mpi_real,0,mpi_comm_world,ierr)
  if (ierr /= mpi_success) call mpi_abort(mpi_comm_world,ierr,ierr)

  call mpi_barrier(mpi_comm_world,ierr)
  
  
#ifdef DEBUG
  do i=0,GridLines**3-1
     do j=1,3
        curHalo = headofchain(j,i)
        do while(curHalo /= 0) 
           if(curHalo < 0 .or. curHalo > halonumber) then
              print*, 'PRE-check: after bcast'
              print*,'RANK',rank,'Error: halo',curHalo,'Invalid'
              call abort
           endif
           curHalo = linkedlist(curhalo)
        end do
     end do
  end do
#endif

  call mpi_barrier(mpi_comm_world, ierr)
  
  allocate(xi_mean(1:n_bin_xi,1:overden_nbin,1:3))
  allocate(xi_mean_sq(1:n_bin_xi,1:overden_nbin,1:3))
  allocate(count_cell_xi(1:overden_nbin,1:3))
  allocate(rho_mean(1:overden_nbin,1:3))
  allocate(density_block(0:GridLines**3-1,1:3))

  !GridLines = 4

  if(rank==0) then
     allocate(final_xi_mean(1:n_bin_xi,1:overden_nbin,1:3))
     final_xi_mean(1:n_bin_xi,1:overden_nbin,1:3) = 0.
     allocate(final_xi_mean_sq(1:n_bin_xi,1:overden_nbin,1:3))
     final_xi_mean_sq(1:n_bin_xi,1:overden_nbin,1:3) = 0.
     call cpu_time(stop_time)
     print*, "Use", stop_time-start_time,"seconds to distribute data"
     call cpu_time(start_time)
  end if
  xi_mean(1:n_bin_xi,1:overden_nbin,1:3) = 0.
  xi_mean_sq(1:n_bin_xi,1:overden_nbin,1:3) = 0.
  count_cell_xi(1:overden_nbin,1:3) = 0
  rho_mean(1:overden_nbin,1:3) = 0.

  do i = 0,GridLines**3-1
     local_overden = Subcell_overden(i)    
     do m=1,3
        density_block(i,m) = return_denbin_number(local_overden,m,density_centre, density_halfbin)
        den_block = density_block(i,m)
        if(den_block >= 1 .and. den_block <= overden_nbin) then
           if(count_haloincell(m,i) >= min_particle_cell) then
              count_cell_xi(den_block,m) = count_cell_xi(den_block,m) + 1
              rho_mean(den_block,m) = rho_mean(den_block,m) + real(count_haloincell(m,i))
           end if
        end if
     end do
  end do
  rho_mean = rho_mean/count_cell_xi

  allocate(cal_cell_per_node(0:nodes_returned))
  cal_cell_per_node(0:nodes_returned) = 0


  if(nodes_returned < GridLines**3) then
     cal_cell_per_node(0:mod(GridLines**3,nodes_returned)-1) = GridLines**3/nodes_returned + 1
     cal_cell_per_node(mod(GridLines**3,nodes_returned):nodes_returned-1) = GridLines**3/nodes_returned
  else 
     cal_cell_per_node(0:GridLines**3-1) = 1
  endif

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  if(rank==0) print*, "Start calculate xi"
  if(rank < GridLines**3) then
     !print*, "rank",rank,"start calculate xi for subcell"
     do k = 1,cal_cell_per_node(rank)
        i = (k-1)*nodes_returned + rank
        local_overden = Subcell_overden(i)
        den_block = ceiling((local_overden - overden_min)/overden_binsize)
        !print*, "den",local_overden,'block',den_block

        do m=1,3
           den_block = density_block(i,m)
           if(den_block >= 1 .and. den_block <= overden_nbin) then
              j = 0
              allocate(inner_positions(1:3,1:count_haloincell(m,i)))
              curhalo = headofchain(m,i)
              do while(curhalo /= 0)
                 j = j + 1
                 inner_positions(1:3,j) = positions(1:3,curhalo) 
                 inner_positions(1:3,j) = inner_positions(1:3,j) - local_positions(1:3,i)
                 inner_positions(1:3,j) = inner_positions(1:3,j) / GridSize
                 curhalo = linkedlist(curhalo)
              end do

              if(j /= count_haloincell(m,i)) then
                 print*, 'number in cell does not match'
                 call abort
              end if
              if(j >= min_particle_cell) then
                 !print*, "cell",i,"massbin=",m,"N=",j

                 if(j == 0 .or. j == 1) then
                    xi_out(1:n_bin_xi) = -1.
                 else
                    call correlation_function_unitcube(rho_mean(den_block,m), j, inner_positions, xi_out, count_for_delta)
                 end if
                 !if(den_block ==2) 
                 !print*,xi_out
                 xi_mean(1:n_bin_xi,den_block,m) = xi_mean(1:n_bin_xi,den_block,m) + real(xi_out(1:n_bin_xi),8)
                 xi_mean_sq(1:n_bin_xi,den_block,m) = xi_mean_sq(1:n_bin_xi,den_block,m) + real(xi_out(1:n_bin_xi),8)**2.
                 !count_cell_xi(0:n_bin_xi,den_block,m) = count_cell_xi(0:n_bin_xi,den_block,m) + 1
              end if
              deallocate(inner_positions)
           end if
        end do

     end do
  end if
  call mpi_barrier(mpi_comm_world, ierr)
  deallocate(positions,local_positions,headofchain,linkedlist,subcell_overden,count_haloincell)
  call mpi_barrier(mpi_comm_world, ierr)

  call mpi_reduce(xi_mean, &
       final_xi_mean,n_bin_xi*overden_nbin*3, &
       mpi_real8,mpi_sum,0,mpi_comm_world,ierr)
  call mpi_reduce(xi_mean_sq, &
       final_xi_mean_sq,n_bin_xi*overden_nbin*3, &
       mpi_real8,mpi_sum,0,mpi_comm_world,ierr)

  if(rank==0) then
     do m=1,3
        do i=1,overden_nbin
           !print*, "m=",m, "den=",(real(i-1)+0.5)*overden_binsize + overden_min
           !print*, "cell",count_cell_xi(i,m), "rho:", rho_mean(i,m)
           !print*, final_xi_mean(0:n_bin_xi,i,m) /count_cell_xi(i,m)
           final_xi_mean(1:n_bin_xi,i,m) = final_xi_mean(1:n_bin_xi,i,m) /real(count_cell_xi(i,m))
           final_xi_mean_sq(1:n_bin_xi,i,m) = final_xi_mean_sq(1:n_bin_xi,i,m) /real(count_cell_xi(i,m))
        end do
     end do

     allocate(delta_cell(1:overden_nbin,1:3))
     !drlog = (log10(max_radius_xi)-log10(min_radius_xi))/n_bin_xi
     !do i=0,n_bin_xi
     !   rad(i)=10.**(log10(min_radius_xi) + i*drlog)
     !end do
     delta_cell(1:overden_nbin,1:3) = 0.
     command = system("mkdir -p "//trim(result_path)//'/'//z_s(1:len_trim(z_s)))
     do m=1,3
        write(mass_s,'(i5)') m
        mass_s = adjustl(mass_s)
        command = system("rm -f "//trim(result_path)//'/'//z_s(1:len_trim(z_s))//'_Delta_'//boxwidth_s(1:len_trim(boxwidth_s))//'_'//grid_s(1:len_trim(grid_s))//'_'//mass_s(1:len_trim(mass_s))//'.Delta')
        open(70,file=trim(result_path)//'/'//z_s(1:len_trim(z_s))//'_Delta_'//boxwidth_s(1:len_trim(boxwidth_s))//'_'//grid_s(1:len_trim(grid_s))//'_'//mass_s(1:len_trim(mass_s))//".Delta",STATUS = 'NEW')
                
        print*, "m = ",m
        do i=1,overden_nbin
           write(delta_s,'(f10.3)') real(density_centre(i,m),8)
           delta_s = adjustl(delta_s)
           if(count_cell_xi(i,m) > 0) then
              command = system("rm -f "//trim(result_path)//'/'//z_s(1:len_trim(z_s))//'_xi_'//boxwidth_s(1:len_trim(boxwidth_s))//'_'//grid_s(1:len_trim(grid_s))//'_'//mass_s(1:len_trim(mass_s))//'delta_'//delta_s(1:len_trim(delta_s))//'.xi' )
              open(71,file=trim(result_path)//'/'//z_s(1:len_trim(z_s))//'_xi_'//boxwidth_s(1:len_trim(boxwidth_s))//'_'//grid_s(1:len_trim(grid_s))//'_'//mass_s(1:len_trim(mass_s))//'delta_'//delta_s(1:len_trim(delta_s))//'.xi' ,STATUS = 'NEW')
              do j=1,n_bin_xi
                 write(71,"(3F12.4,I18)") (real(j)-.5)*(max_radius_xi)/n_bin_xi*Gridsize,final_xi_mean(j,i,m), &
                      sqrt((final_xi_mean_sq(j,i,m)-final_xi_mean(j,i,m)**2.)/(real(count_cell_xi(i,m)) -1.)), &
                      count_for_delta(j)
              end do
              close(71)
           end if
           error_delta = 0.
           do j=1,n_bin_xi
              delta_cell(i,m) = delta_cell(i,m) + final_xi_mean(j,i,m)*real(rho_mean(i,m)**2.*count_for_delta(j)*(1./n_smallcell)**6.)
              error_delta = error_delta + &
                   (real(rho_mean(i,m)**2.*count_for_delta(j)*(1./n_smallcell)**6.)*sqrt((final_xi_mean_sq(j,i,m)-final_xi_mean(j,i,m)**2.)/(real(count_cell_xi(i,m)) -1.)))**2.
           end do
           error_delta = sqrt(error_delta)
           if(count_cell_xi(i,m) > 0) print*, density_centre(i,m),delta_cell(i,m),error_delta,rho_mean(i,m),count_cell_xi(i,m)
           if(count_cell_xi(i,m) > 0) write(70,"(4F12.4,I10)") density_centre(i,m),delta_cell(i,m),error_delta,rho_mean(i,m),count_cell_xi(i,m)
        end do
        close(70)
     end do
     deallocate(delta_cell)
     call cpu_time(stop_time)
     print*,"Use",stop_time-start_time,"seconds to calculate Delta"
  end if
  call mpi_barrier(mpi_comm_world, ierr)
  deallocate(xi_mean,count_cell_xi,rho_mean)
  !print*, "Finish" 
  return

end subroutine subcell_corr




subroutine read_subcell_overden(z,GridLines,subcell_overden)
  use common_vars
  implicit none
  integer :: i,j,k,m,GridLines,GridLines_alt,den_block,devider
  real(kind=4) :: subcell_overden(0:GridLines**3-1)
  real (kind=4) :: z
  character(len=20) :: z_s , GridLines_s
  character(len=200) :: filename
  real(kind=4) :: dummy,den
  real(kind=8) :: totalden,meanden
  LOGICAL :: file_exists

  totalden = 0.
  subcell_overden(0:GridLines**3-1) = 0.
  write(z_s,'(f10.3)') z
  z_s = adjustl(z_s)

  write(GridLines_s,'(I5)') GridLines
  GridLines_s = adjustl(GridLines_s)

  filename =trim(subcellden_path)//z_s(1:len_trim(z_s))//'halo_num'//trim(GridLines_s)//'.bin'
  
  INQUIRE(FILE=trim(filename), EXIST=file_exists)
  if(file_exists) then
     open(28, file=trim(filename), status='old')
     do i=0, GridLines-1
        do j=0,GridLines-1
           do k=0,GridLines-1
              read(28,*) den,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy
              den_block = i+ j*GridLines + k*GridLines**2 
              subcell_overden(den_block) = den     
              totalden = totalden + real(den,8)
           end do
        end do
     end do
     close(28)
     meanden = totalden/(GridLines**3)
     subcell_overden(0:GridLines**3-1) = subcell_overden(0:GridLines**3-1)/(meanden) -1.
     do i=0, GridLines-1
        do j=0,GridLines-1
           do k=0,GridLines-1
              den_block = i+ j*GridLines + k*GridLines**2 
              !print*,subcell_overden(den_block)
           end do
        end do
     end do
     return
  else
     GridLines_alt = 256
     write(GridLines_s,'(I5)') GridLines_alt
     GridLines_s = adjustl(GridLines_s)     
     filename =trim(subcellden_path)//z_s(1:len_trim(z_s))//'halo_num'//trim(GridLines_s)//'.bin'
     INQUIRE(FILE=trim(filename), EXIST=file_exists)
     if(file_exists == .FALSE.) then
        print* , "cannot open file",filename
        call abort
     end if
     open(28, file=trim(filename), status='old')
     do i=0, GridLines_alt-1
        do j=0,GridLines_alt-1
           do k=0,GridLines_alt-1
              read(28,*) den,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy
              den_block = i/(GridLines_alt/GridLines)+ j/(GridLines_alt/GridLines)*GridLines + k/(GridLines_alt/GridLines)*GridLines**2 

              subcell_overden(den_block) = subcell_overden(den_block) + den     
              totalden = totalden + real(den,8)
           end do
        end do
     end do
     close(28)
     meanden = totalden/(GridLines**3)
     do i=0, GridLines**3-1
        print*,subcell_overden(i)/meanden
     end do
     subcell_overden(0:GridLines**3-1) = subcell_overden(0:GridLines**3-1)/(meanden) -1.
     return
  end if
  
end subroutine read_subcell_overden

subroutine read_data_bin_info(GridLines, density_centre, density_halfbin)
  use common_vars
  implicit none
  integer :: GridLines,i
  integer :: mass_block, density_block, numden
  character(len=20) :: grid_s, boxwidth_s, mass_s
  logical :: file_e
  real :: density_centre(1:overden_nbin,1:3), density_halfbin(1:overden_nbin,1:3)
  write(grid_s,'(i5)') GridLines
  grid_s = adjustl(grid_s)
  write(boxwidth_s,'(i5)') int(co_boxwidth)
  boxwidth_s = adjustl(boxwidth_s)
  do mass_block = 1,3
     write(mass_s,'(i5)') mass_block
     mass_s = adjustl(mass_s)
     print*,trim(databin_path)//'delta_'//boxwidth_s(1:len_trim(boxwidth_s))//'_'//grid_s(1:len_trim(grid_s))//'_'//mass_s(1:len_trim(mass_s))
     inquire( file=trim(databin_path)//'delta_'//boxwidth_s(1:len_trim(boxwidth_s))//'_'//grid_s(1:len_trim(grid_s))//'_'//mass_s(1:len_trim(mass_s)), exist=file_e )
     if(file_e == .TRUE.) then
        open(28, file=trim(databin_path)//'delta_'//boxwidth_s(1:len_trim(boxwidth_s))//'_'//grid_s(1:len_trim(grid_s))//'_'//mass_s(1:len_trim(mass_s)), status='old')
        read(28, fmt=*) numden
        do i= 1, numden
           read(28, fmt=*) density_centre(i,mass_block),density_halfbin(i,mass_block)
           !print*,density_centre(i,mass_block),density_halfbin(i,mass_block)
        end do
        do i=numden+1,overden_nbin
           density_centre(i,mass_block) = 1.e30
           density_halfbin(i,mass_block) = 0.
        end do
     else
        density_centre(1:overden_nbin,mass_block) = 1.e30
        density_halfbin(1:overden_nbin,mass_block) = 0.
     end if
  end do
  print*, density_centre
  print*, density_halfbin
  return
end subroutine read_data_bin_info

function return_denbin_number(overden,m,density_centre, density_halfbin )
  use common_vars
  implicit none
  integer :: i,m
  real :: density_centre(1:overden_nbin,1:3), density_halfbin(1:overden_nbin,1:3)
  real :: overden
  integer :: return_denbin_number
  return_denbin_number = -1
  do i=1,overden_nbin
     if(overden > density_centre(i,m)-density_halfbin(i,m) .and.  overden < density_centre(i,m)+density_halfbin(i,m)) then
        return_denbin_number = i
     end if
  end do
  return
end function return_denbin_number

function return_massbin_number(mass)
  implicit none
  integer(kind=4) :: return_massbin_number
  real(kind=4) :: mass
  if (mass .gt. 1.e5 .and. mass .le. 1.e8) then
     return_massbin_number = 1
  else if (mass .gt. 1.8 .and. mass .le. 1.e9 ) then
     return_massbin_number = 2
  else if (mass .gt. 1.e9 .and. mass .le. 1.5e12) then
     return_massbin_number = 3
  end if
  return
end function return_massbin_number

subroutine read_subcelllist()
  use common_vars
  implicit none
  integer :: i
  i = 1
  subcell_perdim_list(:) = 0
  open(unit=20,file=trim(adjustl(subcell_perdim_file)),status="old")
  do
     read(20,'(I5)',end=201) subcell_perdim_list(i)
     i = i+1
     if(i > SUBCELLLIST_MAXSIZE) call abort("FATAL: zlist is too large")
  end do
201 close(20)
end subroutine read_subcelllist
