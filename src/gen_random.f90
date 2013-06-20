subroutine gen_random(z)
  use datatools
  use ifport
  use mpi
  use mpitools
  use io_tools
  implicit none   
  !integer(kind=4) , parameter :: line_length_factor = 2
  integer (kind = 4)   :: n_threads, omp_thread

  integer(kind=4)   :: omp,fh_readhalo,fh_file

  real(kind=4) :: centre_point(1:3)

  integer(kind=4) :: GridLines 
  integer(kind=4):: PrGridLines 
  real(kind=4) :: GridSize

  real(kind=4) :: redshift_index 

  integer(kind=4) :: NumBlock,curHalo,unitid,totalhaloline !curHalo => 8byte

  real(kind=4) :: crossBlock(3),shiftdistance(3),pos(3),block_dummy(3)
  real(kind=4) :: baserandom1(3),baserandom2(3),baserandom3(3)
  integer(kind=4) :: i,j,k,l,m,n,xStart,xStop,o,linenum,hitnum
  integer(kind=4) :: length,curBox,box1(27),box2(27),mixbox(54)
  integer(kind=4) :: mixbox_x(54),mixbox_y(54), &
       mixbox_z(54),shiftcell(3),abscell(3)
  integer(kind=4) :: timearray(8),command
  real(kind=4),allocatable :: positions(:,:),radius(:),mass(:),spin(:,:)
  integer(kind=4), allocatable :: headofchain(:),linkedlist(:),hitperthread(:),missperthread(:)


  character(len=100) :: str_rank,str_omp,filename
  character(len=20) :: z_s
  integer(kind=mpi_offset_kind) :: filesize
  integer(kind=4) :: n_point,totalpoint,tag,n_point_index
  integer(kind=4),allocatable :: n_point_arr(:)
  real(kind=4), allocatable :: halo_in(:,:)
  integer(4), dimension(mpi_status_size) :: status
  real(kind=4) :: max_radius, dummy
  integer(kind=4), allocatable :: file_read_per_node(:),n_point_file(:)
  logical :: file_e, element_flag(1:17)
  real(kind=8) :: z,d0
  integer(kind=4) :: n_elements, haloperline,count_halo,boydID,status_checkfiles

  element_flag(1:17) = .TRUE.
  !element_flag(1:3) = .TRUE. ! 1:3 => positions
  !element_flag(7:9) = .TRUE. ! 7:9 => velocities
  !element_flag(14:15) = .TRUE. ! 14:15 =>  radius: mass
  dummy = 0.

  n_elements = get_n_elements(element_flag)

  !halonumber = 10

  write(z_s,'(f10.3)') z
  z_s = adjustl(z_s)

  omp_thread = omp_get_max_threads()

  write(str_rank,'(i10)') rank
  str_rank = adjustl(str_rank)


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
#ifdef DEBUG
        print*, "#DEBUG: halonumber in node 0",halonumber
#endif
     end if
     call MPI_BARRIER(MPI_COMM_WORLD,ierr)
     call mpi_bcast(halonumber,1,mpi_integer,0,mpi_comm_world,ierr)
     call MPI_BARRIER(MPI_COMM_WORLD,ierr)
     !print*, "total halos",halonumber
#ifdef DEBUG
     print*, "#DEBUG: halonumber in node ",rank,halonumber
#endif
     if(rank == 0) allocate(halodata(1:n_elements,1:halonumber))
     call mpi_read_halo(real(z,4),element_flag,n_elements)
  else
     if(rank==0) print*, "No file to read"
     call abort
  end if
  !## End read data to node 0
  if(rank==0) then
     print*, "Start gridding"
     GridLines = file_dimension
     GridSize = Boxsize/GridLines
     PrGridLines = GridLines*(2*line_length_factor+1)
     print*, 'Gridsize: ',Gridsize
     print*, 'GridLines: ',GridLines

     print*, 'Subvolume:',0,'-',GridLines**3-1

     allocate(headofchain(0:Gridlines**3-1))
     allocate(linkedlist(1:halonumber))


     !set the initials for linked listing
     headofchain(0:Gridlines**3-1) = 0
     linkedlist(1:halonumber) = 0


     !put halos in small cells
     write(*,*) 'Putting halos in cells'
     do i=1,halonumber 
        halodata(1:3,i) =  (/ rand(0), rand(0), rand(0) /) * Boxsize
        block_dummy(1) = floor(halodata(1,i)/GridSize)
        if(block_dummy(1) == GridLines) block_dummy(1) = GridLines-1
        block_dummy(2) = floor(halodata(2,i)/GridSize)
        if(block_dummy(2) == GridLines) block_dummy(2) = GridLines-1
        block_dummy(3) = floor(halodata(3,i)/GridSize)
        if(block_dummy(3) == GridLines) block_dummy(3) = GridLines-1

        NumBlock = block_dummy(1)*GridLines**2 + &
             block_dummy(2)*GridLines + &
             block_dummy(3)


        linkedlist(i) = headofchain(Numblock)
        headofchain(Numblock) = i
        
     end do

     !write files
     print*, "Writing outputs"
     command = system('mkdir -p '//trim(adjustl(randomhalo_path)))
     command = system('rm -f '//trim(randomhalo_path)//'/*')
     do i=0,Gridlines**3-1
        write(str_rank,'(i10)') i
        str_rank = adjustl(str_rank)
        call mpi_file_open(mpi_comm_self, &
             trim(randomhalo_path)//z_s(1:len_trim(z_s))//'halo'//str_rank(1:len_trim(str_rank))//'.dat', &
             MPI_MODE_WRONLY + MPI_MODE_CREATE, &
             mpi_info_null,fh_file,ierr)
        call MPI_FILE_WRITE(fh_file, dummy, 1, MPI_REAL, MPI_STATUS_IGNORE, ierr)
        curHalo = headofchain(i)
        do while(curHalo /= 0) 
           call MPI_FILE_WRITE(fh_file, halodata(1:17,curHalo), 17, MPI_REAL, MPI_STATUS_IGNORE, ierr)
           curHalo = linkedlist(curHalo)
        end do
        call MPI_FILE_CLOSE(fh_file, ierr)
     end do

     deallocate(halodata,headofchain,linkedlist)
  endif
end subroutine gen_random
