subroutine gen_random(z)
  use omp_lib
  use vectortools
  use arraytools
  use datatools
  use ifport
  use mpi
  use mpitools
  use io_tools
  implicit none   
  !integer(kind=4) , parameter :: line_length_factor = 2
  integer (kind = 4)   :: n_threads, omp_thread

  integer(kind=4)   :: omp,fh_readhalo

  real(kind=4) :: centre_point(1:3)

  integer(kind=4) :: GridLines 
  integer(kind=4):: PrGridLines 
  real(kind=4) :: GridSize

  real(kind=4) :: redshift_index 

  integer(kind=4) :: NumBlock,curHalo,unitid,totalhaloline !curHalo => 8byte

  real(kind=4) :: crossBlock(3),shiftdistance(3),pos(3),block_dummy(3)
  real(kind=4) :: baserandom1(3),baserandom2(3),baserandom3(3)
  integer(kind=4) :: i,j,k,l,m,n,xStart,xStop,o,linenum,hitnum,halonumber
  integer(kind=4) :: length,curBox,box1(27),box2(27),mixbox(54)
  integer(kind=4) :: mixbox_x(54),mixbox_y(54), &
       mixbox_z(54),shiftcell(3),abscell(3)
  integer(kind=4) :: timearray(8),command
  real(kind=4),allocatable :: halodata(:,:), positions(:,:),radius(:),mass(:),spin(:,:)
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
  integer(kind=4) :: n_elements, haloperline,count_halo,boydID

  element_flag(1:17) = .TRUE.
  !element_flag(1:3) = .TRUE. ! 1:3 => positions
  !element_flag(7:9) = .TRUE. ! 7:9 => velocities
  !element_flag(14:15) = .TRUE. ! 14:15 =>  radius: mass
  if(nodes_returned > 1) call abort()

  halonumber = get_halo_number(real(z,4))
  n_elements = get_n_elements(element_flag)

  !halonumber = 10
  allocate(halodata(1:n_elements,1:halonumber))
  call mpi_read_halo(real(z,4),element_flag,n_elements,halonumber, halodata)
  
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  write(z_s,'(f10.3)') z
  z_s = adjustl(z_s)

  write(str_rank,'(i10)') rank
  omp_thread = omp_get_max_threads()

  write(str_rank,'(i10)') rank
  str_rank = adjustl(str_rank)

  if (rank==0) print*, 'finish reading'


  print*,'MPI node:',nodes_returned
  print*,'OMP node:',omp_thread

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

  deallocate(halodata,headofchain,linkedlist)
 
end subroutine gen_random
