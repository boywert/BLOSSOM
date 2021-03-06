#ifdef RR
subroutine gen_los_rr(z)
#else
subroutine gen_los_rg(z)
#endif
  use omp_lib
  use vectortools
  use arraytools
  use datatools
  use ifport
  use mpi
  use mpitools
  use io_tools
  implicit none   

  integer (kind = 4)   :: n_threads, omp_thread

  integer(kind=4)   :: omp
  integer(kind=4)   :: fh_haloid,fh_lineid,fh_online,fh_toline,fh_logs,fh_centreref,fh_readhalo,fh_direction,fh_hitpoint,fh_startpoint
  integer(kind=4)   :: fh_haloid_disk,fh_lineid_disk,fh_online_disk,fh_toline_disk,fh_direction_disk,fh_hitpoint_disk
  real(kind=4) :: centre_point(1:3)

  integer(kind=4) :: GridLines 
  integer(kind=4):: PrGridLines 
  real(kind=4) :: GridSize 

  real(kind=4) :: redshift_index 

  integer(kind=4) :: NumBlock,curHalo,unitid,totalhaloline !curHalo => 8byte
  real(kind=4) :: InitPoint(3), TargetPoint(3), Direction(3)
  real(kind=4) :: FinishPoint(3),maxradius,r4_uniform_01

  real(kind=4) :: crossBlock(3),shiftdistance(3),pos(3),block_dummy(3)
  real(kind=4) :: baserandom1(3),baserandom2(3),baserandom3(3)
  integer(kind=4) :: i,j,k,l,m,n,xStart,xStop,o,linenum,hitnum,halonumber
  integer(kind=4) :: length,curBox,box1(27),box2(27),mixbox(54)
  integer(kind=4) :: mixbox_x(54),mixbox_y(54), &
       mixbox_z(54),shiftcell(3),abscell(3)
  integer(kind=4) :: timearray(8),command
  real(kind=4),allocatable :: halodata(:,:), positions(:,:),radius(:),mass(:),spin(:,:)
  integer(kind=4), allocatable :: headofchain(:),linkedlist(:),hitperthread(:),missperthread(:)
  real(kind=4) :: distanceonline,distancetoline,initarray(max_l,3)
  real(kind=4) :: MassivePos(1:3,1:MassiveNumber)
  integer(kind=4), allocatable :: totalcellx(:),&
       totalcelly(:), &
       totalcellz(:),&
       totalcell_dummy(:)
  character(len=100) :: str_rank,filename
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

  element_flag(1:17) = .FALSE.
  element_flag(1:3) = .TRUE. ! 1:3 => positions
  !element_flag(7:9) = .TRUE. ! 7:9 => velocities
  element_flag(14:15) = .TRUE. ! 14:15 =>  radius: mass


  halonumber = get_halo_number(real(z,4))
  n_elements = get_n_elements(element_flag)

  !halonumber = 10
  allocate(halodata(1:n_elements,1:halonumber))
  call mpi_read_halo(real(z,4),element_flag,n_elements,halonumber, halodata)
  
  !halodata(1:3,1:9) = 0.1
  !halodata(1:3,10) = Boxsize-0.1
  !halodata(4,1:10) = 1.0
  !halodata(5,1:9) = 0.1
  !halodata(5,10) = 1.0
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  write(z_s,'(f10.3)') z
  z_s = adjustl(z_s)

  write(str_rank,'(i10)') rank
  omp_thread = omp_get_max_threads()

  write(str_rank,'(i10)') rank
  str_rank = adjustl(str_rank)
  print*, 'finish reading'


  if (rank==0) then

     print*,'MPI node:',nodes_returned
     print*,'OMP node:',omp_thread

     write(*,*) 'Start using Node 0 to read and transfer the essentials to others'

     !count total halo from shell
     !halonumber = getlinenumber(halo_path, 30)

     !allocate arrays

     allocate(positions(1:3,1:halonumber))
     !allocate(velocity(1:halonumber,1:3))
     allocate(radius(1:halonumber))
     allocate(mass(1:halonumber))
     positions(1:3,1:halonumber) = halodata(1:3,1:halonumber)
     !velocity(1:halonumber,1:3) = halodata(1:halonumber,4:6)
     radius(1:halonumber) = halodata(4,1:halonumber)
     mass(1:halonumber) = halodata(5,1:halonumber)
     deallocate(halodata)
     !find top 100 massive halos

     call findmassive(positions, mass, halonumber, MassiveNumber, MassivePos)
     deallocate(mass)
     max_radius = 0.
     do i=1,halonumber
        if(radius(i) > max_radius) max_radius = radius(i)
     end do
     print*, 'max radius', max_radius
     do i=4*ceiling(max_radius), ceiling(Boxsize)
        if(mod(int(Boxsize),i) == 0  ) then
           goto 151
        end if
     end do
     print* , 'Cannot find Gridsize' 
     call abort
151  GridSize = real(i)
     !GridSize = Boxsize/6.
     GridLines = int(BoxSize/GridSize)
     PrGridLines = GridLines*3
     print*, 'Gridsize: ',Gridsize
     print*, 'GridLines: ',GridLines

     print*, 'Subvolume:',0,'-',GridLines**3-1
     print*, 'All Subvolume:',0,'-',PrGridLines**3-1 

     allocate(headofchain(0:Gridlines**3-1))
     allocate(linkedlist(1:halonumber))


     !set the initials for linked listing
     headofchain(0:Gridlines**3-1) = 0
     linkedlist(1:halonumber) = 0

     do i=1,MassiveNumber
        block_dummy(1) = floor(MassivePos(1,i)/GridSize)
        if(block_dummy(1) == GridLines) block_dummy(1) = GridLines-1
        block_dummy(2) = floor(MassivePos(2,i)/GridSize)
        if(block_dummy(2) == GridLines) block_dummy(2) = GridLines-1
        block_dummy(3) = floor(Massivepos(3,i)/GridSize)
        if(block_dummy(3) == GridLines) block_dummy(3) = GridLines-1

        NumBlock = block_dummy(1)*GridLines**2 + &
             block_dummy(2)*GridLines + &
             block_dummy(3)
        print*, "massive", i, NumBlock
     end do
     !put halos in small cells
     write(*,*) 'Putting halos in cells'
     do i=1,halonumber 
        block_dummy(1) = floor(positions(1,i)/GridSize)
        if(block_dummy(1) == GridLines) block_dummy(1) = GridLines-1
        block_dummy(2) = floor(positions(2,i)/GridSize)
        if(block_dummy(2) == GridLines) block_dummy(2) = GridLines-1
        block_dummy(3) = floor(positions(3,i)/GridSize)
        if(block_dummy(3) == GridLines) block_dummy(3) = GridLines-1

        NumBlock = block_dummy(1)*GridLines**2 + &
             block_dummy(2)*GridLines + &
             block_dummy(3)
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
        
        linkedlist(i) = headofchain(Numblock)
        headofchain(Numblock) = i
        if(positions(1,i) == Massivepos(1,1) .and. positions(2,i) == Massivepos(2,1) .and. positions(3,i) == Massivepos(3,1)) then
           print*, "found init", NumBlock,i
           boydID = i
        end if
        !145 continue
     end do



     write(*,*) 'Making linkedlist'
     do i=0,GridLines**3-1
        curHalo = headofchain(i)
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
     write(*,*) 'Broadcasting ......'
  endif
  call mpi_barrier(mpi_comm_world,ierr)
  !broadcast halonumber
  call mpi_bcast(halonumber,1,mpi_integer,0,mpi_comm_world,ierr)
  call mpi_barrier(mpi_comm_world,ierr)
  call mpi_bcast(GridLines,1,mpi_integer,0,mpi_comm_world,ierr)
  call mpi_bcast(PrGridLines,1,mpi_integer,0,mpi_comm_world,ierr)
  call mpi_bcast(GridSize,1,mpi_real,0,mpi_comm_world,ierr)
  call mpi_barrier(mpi_comm_world,ierr)


  !allocate array in other nodes
  if (rank /= 0) then
     allocate(positions(1:3,1:halonumber))
     allocate(radius(1:halonumber))
     allocate(headofchain(0:Gridlines**3-1))
     allocate(linkedlist(1:halonumber))
  endif
  call mpi_barrier(mpi_comm_world,ierr)
  !broadcast variables
  call mpi_bcast(positions,halonumber*3,mpi_real,0,mpi_comm_world,ierr)
  if (ierr /= mpi_success) call mpi_abort(mpi_comm_world,ierr,ierr)

  call mpi_bcast(radius,halonumber,mpi_real,0,mpi_comm_world,ierr)
  if (ierr /= mpi_success) call mpi_abort(mpi_comm_world,ierr,ierr)

  call mpi_bcast(MassivePos,MassiveNumber*3,mpi_real,0,mpi_comm_world,ierr)
  if (ierr /= mpi_success) call mpi_abort(mpi_comm_world,ierr,ierr)

  call mpi_bcast(linkedlist,halonumber,mpi_integer,0,mpi_comm_world,ierr)
  if (ierr /= mpi_success) call mpi_abort(mpi_comm_world,ierr,ierr)

  call mpi_bcast(headofchain,Gridlines**3,mpi_integer,0,mpi_comm_world,ierr)
  if (ierr /= mpi_success) call mpi_abort(mpi_comm_world,ierr,ierr)

  call mpi_barrier(mpi_comm_world,ierr)

#ifdef DEBUG
  do i=0,GridLines**3-1
     curHalo = headofchain(i)
     do while(curHalo /= 0) 
        if(curHalo < 0 .or. curHalo > halonumber) then
           print*, 'PRE-check: after bcast'
           print*,'RANK',rank,'Error: halo',curHalo,'Invalid'
           call abort
        endif
        curHalo = linkedlist(curhalo)
     end do
  end do
#endif

  totalhaloline = 0

  if(rank == 0) then
     print*, 'creating output folders....' 

#ifdef RR
     command = system('mkdir -p '//trim(los_path)//z_s(1:len_trim(z_s))//'/RR/'//'LINEID/')
     command = system('rm -f '//trim(los_path)//z_s(1:len_trim(z_s))//'/RR/'//'LINEID/*')
     command = system('mkdir -p '//trim(los_path)//z_s(1:len_trim(z_s))//'/RR/'//'HALOID/')
     command = system('rm -f '//trim(los_path)//z_s(1:len_trim(z_s))//'/RR/'//'HALOID/*')
     command = system('mkdir -p '//trim(los_path)//z_s(1:len_trim(z_s))//'/RR/'//'ONLINE/')
     command = system('rm -f '//trim(los_path)//z_s(1:len_trim(z_s))//'/RR/'//'ONLINE/*')
     command = system('mkdir -p '//trim(los_path)//z_s(1:len_trim(z_s))//'/RR/'//'TOLINE/')
     command = system('rm -f '//trim(los_path)//z_s(1:len_trim(z_s))//'/RR/'//'TOLINE/*')
     command = system('mkdir -p '//trim(los_path)//z_s(1:len_trim(z_s))//'/RR/'//'LOGS/')
     command = system('rm -f '//trim(los_path)//z_s(1:len_trim(z_s))//'/RR/'//'LOGS/*')
     command = system('mkdir -p '//trim(los_path)//z_s(1:len_trim(z_s))//'/RR/'//'DIRECTION/')
     command = system('rm -f '//trim(los_path)//z_s(1:len_trim(z_s))//'/RR/'//'DIRECTION/*')
     command = system('mkdir -p '//trim(los_path)//z_s(1:len_trim(z_s))//'/RR/'//'HITPOINT/')
     command = system('rm -f '//trim(los_path)//z_s(1:len_trim(z_s))//'/RR/'//'HITPOINT/*')
     command = system('mkdir -p '//trim(los_path)//z_s(1:len_trim(z_s))//'/RR/'//'STARTPOINT/')
     command = system('rm -f '//trim(los_path)//z_s(1:len_trim(z_s))//'/RR/'//'STARTPOINT/*')
#else
     command = system('mkdir -p '//trim(los_path)//z_s(1:len_trim(z_s))//'/RG/'//'LINEID/')
     command = system('rm -f '//trim(los_path)//z_s(1:len_trim(z_s))//'/RG/'//'LINEID/*')
     command = system('mkdir -p '//trim(los_path)//z_s(1:len_trim(z_s))//'/RG/'//'HALOID/')
     command = system('rm -f '//trim(los_path)//z_s(1:len_trim(z_s))//'/RG/'//'HALOID/*')
     command = system('mkdir -p '//trim(los_path)//z_s(1:len_trim(z_s))//'/RG/'//'ONLINE/')
     command = system('rm -f '//trim(los_path)//z_s(1:len_trim(z_s))//'/RG/'//'ONLINE/*')
     command = system('mkdir -p '//trim(los_path)//z_s(1:len_trim(z_s))//'/RG/'//'TOLINE/')
     command = system('rm -f '//trim(los_path)//z_s(1:len_trim(z_s))//'/RG/'//'TOLINE/*')
     command = system('mkdir -p '//trim(los_path)//z_s(1:len_trim(z_s))//'/RG/'//'LOGS/')
     command = system('rm -f '//trim(los_path)//z_s(1:len_trim(z_s))//'/RG/'//'LOGS/*')
     command = system('mkdir -p '//trim(los_path)//z_s(1:len_trim(z_s))//'/RG/'//'DIRECTION/')
     command = system('rm -f '//trim(los_path)//z_s(1:len_trim(z_s))//'/RG/'//'DIRECTION/*')
     command = system('mkdir -p '//trim(los_path)//z_s(1:len_trim(z_s))//'/RG/'//'HITPOINT/')
     command = system('rm -f '//trim(los_path)//z_s(1:len_trim(z_s))//'/RG/'//'HITPOINT/*')
     command = system('mkdir -p '//trim(los_path)//z_s(1:len_trim(z_s))//'/RG/'//'STARTPOINT/')
     command = system('rm -f '//trim(los_path)//z_s(1:len_trim(z_s))//'/RG/'//'STARTPOINT/*')
#endif

  endif


  call mpi_barrier(mpi_comm_world, ierr)

#ifdef RR
  call mpi_file_open(mpi_comm_self, &
       trim(los_path)//z_s(1:len_trim(z_s))//'/RR/'//'LINEID/'//trim(adjustl(str_rank)), &
       MPI_MODE_WRONLY + MPI_MODE_CREATE, &
       mpi_info_null,fh_lineid,ierr)
  call mpi_file_open(mpi_comm_self, &
       trim(los_path)//z_s(1:len_trim(z_s))//'/RR/'//'HALOID/'//trim(adjustl(str_rank)), &
       MPI_MODE_WRONLY + MPI_MODE_CREATE, &
       mpi_info_null,fh_haloid,ierr)
  call mpi_file_open(mpi_comm_self, &
       trim(los_path)//z_s(1:len_trim(z_s))//'/RR/'//'ONLINE/'//trim(adjustl(str_rank)), &
       MPI_MODE_WRONLY + MPI_MODE_CREATE, &
       mpi_info_null,fh_online,ierr)
  call mpi_file_open(mpi_comm_self, &
       trim(los_path)//z_s(1:len_trim(z_s))//'/RR/'//'TOLINE/'//trim(adjustl(str_rank)), &
       MPI_MODE_WRONLY + MPI_MODE_CREATE, &
       mpi_info_null,fh_toline,ierr)
  call mpi_file_open(mpi_comm_self, &
       trim(los_path)//z_s(1:len_trim(z_s))//'/RR/'//'DIRECTION/'//trim(adjustl(str_rank)), &
       MPI_MODE_WRONLY + MPI_MODE_CREATE, &
       mpi_info_null,fh_direction,ierr)
  call mpi_file_open(mpi_comm_self, &
       trim(los_path)//z_s(1:len_trim(z_s))//'/RR/'//'HITPOINT/'//trim(adjustl(str_rank)), &
       MPI_MODE_WRONLY + MPI_MODE_CREATE, &
       mpi_info_null,fh_hitpoint,ierr)
  call mpi_file_open(mpi_comm_self, &
       trim(los_path)//z_s(1:len_trim(z_s))//'/RR/'//'STARTPOINT/'//trim(adjustl(str_rank)), &
       MPI_MODE_WRONLY + MPI_MODE_CREATE, &
       mpi_info_null,fh_startpoint,ierr)
#else
  call mpi_file_open(mpi_comm_self, &
       trim(los_path)//z_s(1:len_trim(z_s))//'/RG/'//'LINEID/'//trim(adjustl(str_rank)), &
       MPI_MODE_WRONLY + MPI_MODE_CREATE, &
       mpi_info_null,fh_lineid,ierr)
  call mpi_file_open(mpi_comm_self, &
       trim(los_path)//z_s(1:len_trim(z_s))//'/RG/'//'HALOID/'//trim(adjustl(str_rank)), &
       MPI_MODE_WRONLY + MPI_MODE_CREATE, &
       mpi_info_null,fh_haloid,ierr)
  call mpi_file_open(mpi_comm_self, &
       trim(los_path)//z_s(1:len_trim(z_s))//'/RG/'//'ONLINE/'//trim(adjustl(str_rank)), &
       MPI_MODE_WRONLY + MPI_MODE_CREATE, &
       mpi_info_null,fh_online,ierr)
  call mpi_file_open(mpi_comm_self, &
       trim(los_path)//z_s(1:len_trim(z_s))//'/RG/'//'TOLINE/'//trim(adjustl(str_rank)), &
       MPI_MODE_WRONLY + MPI_MODE_CREATE, &
       mpi_info_null,fh_toline,ierr)
  call mpi_file_open(mpi_comm_self, &
       trim(los_path)//z_s(1:len_trim(z_s))//'/RG/'//'DIRECTION/'//trim(adjustl(str_rank)), &
       MPI_MODE_WRONLY + MPI_MODE_CREATE, &
       mpi_info_null,fh_direction,ierr)
  call mpi_file_open(mpi_comm_self, &
       trim(los_path)//z_s(1:len_trim(z_s))//'/RG/'//'HITPOINT/'//trim(adjustl(str_rank)), &
       MPI_MODE_WRONLY + MPI_MODE_CREATE, &
       mpi_info_null,fh_hitpoint,ierr)
  call mpi_file_open(mpi_comm_self, &
       trim(los_path)//z_s(1:len_trim(z_s))//'/RG/'//'STARTPOINT/'//trim(adjustl(str_rank)), &
       MPI_MODE_WRONLY + MPI_MODE_CREATE, &
       mpi_info_null,fh_startpoint,ierr)
#endif  

  fh_logs = 50
  fh_centreref = 51
  !open (unit=fh_logs,file=trim(los_path)//z_s(1:len_trim(z_s))//'/'//'LOGS/'//trim(adjustl(str_rank)))
  !open (unit=fh_centreref,file=trim(los_path)//z_s(1:len_trim(z_s))//'/'//'CENTREREF/'//trim(adjustl(str_rank)))
  call mpi_barrier(mpi_comm_world, ierr)


  

  
  call random_seed

  if(rank==0) print*, 'Start LOS finder'

  !$omp parallel private(i,j,k,l,o,InitPoint,TargetPoint,Direction,FinishPoint,xStart,xStop,totalcellx,totalcelly,totalcellz,totalcell_dummy,box1,box2,crossBlock,length,shiftcell,shiftdistance,Numblock,curHalo,pos,distancetoline,distanceonline,unitid,timearray,baserandom3,baserandom2,baserandom1,m,haloperline) &
  !$omp shared(halodata,initarray,linenum, hitnum)
  call date_and_time(values=timearray)
  call srand(rank*omp_get_thread_num()+timearray(5)+timearray(6)+timearray(7)+timearray(8)+timearray(1)*(rank+1)+1)
  !call srand(rank + omp_get_thread_num() + 2)

  !$omp do
  do l=1,max_l

     allocate(totalcellx(1:54*(GridLines+1)))
     allocate(totalcelly(1:54*(GridLines+1)))
     allocate(totalcellz(1:54*(GridLines+1)))
     allocate(totalcell_dummy(1:54*(GridLines+1)*3))

     totalcellx(1:size(totalcellx)) = -1
     totalcelly(1:size(totalcelly)) = -1
     totalcellz(1:size(totalcellz)) = -1
     totalcell_dummy(1:size(totalcell_dummy)) = -1

#ifdef RR
     !print*, "RR mode"
     InitPoint = (/ rand(0), rand(0), rand(0) /)
     TargetPoint = (/ rand(0), rand(0), rand(0) /)
     InitPoint(1:3) = InitPoint(1:3)*BoxSize
     TargetPoint(1:3) = TargetPoint(1:3)*BoxSize
     !write(fh_logs,*) hitperthread(omp_get_thread_num()),omp_get_thread_num(),TargetPoint
#else
     !print*, "RG mode"
     TargetPoint = (/ rand(0), rand(0), rand(0) /)
     !TargetPoint = (/ 0.5, 0.5, 0.5 /)
     TargetPoint(1:3) = TargetPoint(1:3)*BoxSize
     m = int(rand(0)*MassiveNumber)+1
     InitPoint(1:3) = MassivePos(1:3,m)
     !TargetPoint = (/ initpoint(1), initpoint(1), initpoint(3)+BoxSize/2. /)
     !TargetPoint(1:3) = TargetPoint(1:3)*BoxSize
     !initpoint = (/0.51,0.51, 2000. /)
     !write(fh_logs,*) hitperthread(omp_get_thread_num()),omp_get_thread_num(),TargetPoint
#endif
     !write(*,*) InitPoint, omp_get_thread_num(),rank
     !write(fh_logs,*)  TargetPoint
     initPoint = (/ BoxSize, BoxSize, BoxSize/) + initPoint
     TargetPoint = (/ BoxSize, BoxSize, BoxSize/) + TargetPoint

     !print*,'ori',  floor((InitPoint(1)-BoxSize)/GridSize)*GridLines**2 + &
     !     floor((InitPoint(2)-BoxSize)/GridSize)*GridLines + &
     !     floor((InitPoint(3)-BoxSize)/GridSize)
     !print*,'update',  floor(InitPoint(1)/GridSize)*PrGridLines**2 + &
     !     floor(InitPoint(2)/GridSize)*PrGridLines + &
     !     floor(InitPoint(3)/GridSize)


     Direction(1:3) =(TargetPoint(1:3)-InitPoint(1:3))/vectorabs(TargetPoint - InitPoint)
     !initpoint(1:3) = InitPoint(1:3) - BoxSize * Direction(1:3) / 2.
     FinishPoint(1:3) = InitPoint(1:3) + BoxSize * Direction(1:3)
     !if(rank == 0 .and. omp_get_thread_num() == 0) then
     !                 print*, 'finishppoint',FinishPoint
     !                 print*, 'target',TargetPoint
     !                 print*, 'init',InitPoint
     !                 print*, 'direction',Direction
     !endif
     haloperline = 0
     centre_point = (/ BoxSize*3./2., BoxSize*3./2., BoxSize*3./2. /)
!!$omp critical 
     !write(fh_centreref,*)  int((l+(m-1)*max_l),4),linedistance(initPoint,FinishPoint,centre_point), finddistance(initPoint,FinishPoint,centre_point)
!!$omp end critical
     !print*, m,finddistance(initPoint,FinishPoint,MassivePos(1:3,m)+(/ BoxSize, BoxSize, BoxSize/)),linedistance(initPoint,FinishPoint,MassivePos(1:3,m)+(/ BoxSize, BoxSize, BoxSize/))


     !increase in x-direction     
     if(Direction(1) >= 0) then
        xStart = ceiling(InitPoint(1)/(GridSize))
        xStop = floor(FinishPoint(1)/(GridSize))
     else 
        xStart = ceiling(FinishPoint(1)/(GridSize))
        xStop = floor(InitPoint(1)/(GridSize))
     end if

     !    if(rank == 0 .and. omp_get_thread_num() == 0) then
     !                write(*,*),'x', xStart, xStop
     !    endif
     !print*, InitPoint(1)/(GridSize), FinishPoint(1)/(GridSize)

     !write(*,*), xStart, xStop,iabs(xStop-xStart)+1
     if(int(InitPoint(1)/(GridSize)) == int(FinishPoint(1)/(GridSize))) then
        box1(1) = floor(InitPoint(1)/GridSize)*PrGridLines**2 + &
             floor(InitPoint(2)/GridSize)*PrGridLines + &
             floor(InitPoint(3)/GridSize)

        box1(2) = box1(1)+1
        box1(8) = box1(2)+PrGridLines
        box1(9) = box1(2)-PrGridLines
        box1(10) = box1(2)+PrGridLines**2
        box1(11) = box1(2)-PrGridLines**2
        box1(3) = box1(1)-1
        box1(12) = box1(3)+PrGridLines
        box1(13) = box1(3)-PrGridLines
        box1(14) = box1(3)+PrGridLines**2
        box1(15) = box1(3)-PrGridLines**2
        box1(4) = box1(1)+PrGridLines
        box1(16) = box1(4)+PrGridLines**2
        box1(17) = box1(4)-PrGridLines**2
        box1(5) = box1(1)-PrGridLines
        box1(18) = box1(5)+PrGridLines**2
        box1(19) = box1(5)-PrGridLines**2
        box1(6) = box1(1)+PrGridLines**2
        box1(7) = box1(1)-PrGridLines**2
        box1(20) = box1(8)+PrGridLines**2
        box1(21) = box1(8)-PrGridLines**2
        box1(22) = box1(9)+PrGridLines**2
        box1(23) = box1(9)-PrGridLines**2
        box1(24) = box1(12)+PrGridLines**2
        box1(25) = box1(12)-PrGridLines**2
        box1(26) = box1(13)+PrGridLines**2
        box1(27) = box1(13)-PrGridLines**2
        totalcellx(1:27) = box1(1:27)
     else
        j=0

        do i = xStart,xStop,1
     
           crossBlock(1) = real(i)*(GridSize)
           if(Direction(1) /= 0.) then

              length = (crossBlock(1) - InitPoint(1))/Direction(1)
           else 
              length = 0.
           endif
           crossBlock(2:3) = InitPoint(2:3) + length*Direction(2:3)
#ifdef DEBUG   
           if(crossBlock(1)<0 .or. crossBlock(2)<0 .or. crossBlock(3) < 0) then
              print*,'x',crossBlock
              print*,'direction',Direction
              print*,'length',length
              print*,'start',xstart,'stop',xstop
              call abort
           endif
#endif

           box1(1) = floor(crossBlock(1)/GridSize)*PrGridLines**2 + &
                floor(crossBlock(2)/GridSize)*PrGridLines + &
                floor(crossBlock(3)/GridSize)

           box1(2) = box1(1)+1
           box1(8) = box1(2)+PrGridLines
           box1(9) = box1(2)-PrGridLines
           box1(10) = box1(2)+PrGridLines**2
           box1(11) = box1(2)-PrGridLines**2
           box1(3) = box1(1)-1
           box1(12) = box1(3)+PrGridLines
           box1(13) = box1(3)-PrGridLines
           box1(14) = box1(3)+PrGridLines**2
           box1(15) = box1(3)-PrGridLines**2
           box1(4) = box1(1)+PrGridLines
           box1(16) = box1(4)+PrGridLines**2
           box1(17) = box1(4)-PrGridLines**2
           box1(5) = box1(1)-PrGridLines
           box1(18) = box1(5)+PrGridLines**2
           box1(19) = box1(5)-PrGridLines**2
           box1(6) = box1(1)+PrGridLines**2
           box1(7) = box1(1)-PrGridLines**2
           box1(20) = box1(8)+PrGridLines**2
           box1(21) = box1(8)-PrGridLines**2
           box1(22) = box1(9)+PrGridLines**2
           box1(23) = box1(9)-PrGridLines**2
           box1(24) = box1(12)+PrGridLines**2
           box1(25) = box1(12)-PrGridLines**2
           box1(26) = box1(13)+PrGridLines**2
           box1(27) = box1(13)-PrGridLines**2

           box2(1) =(floor(crossBlock(1)/GridSize)+1)*PrGridLines**2 + &
                floor(crossBlock(2)/GridSize)*PrGridLines + &
                floor(crossBlock(3)/GridSize)

           box2(2) = box2(1)+1
           box2(8) = box2(2)+PrGridLines
           box2(9) = box2(2)-PrGridLines
           box2(10) = box2(2)+PrGridLines**2
           box2(11) = box2(2)-PrGridLines**2
           box2(3) = box2(1)-1
           box2(12) = box2(3)+PrGridLines
           box2(13) = box2(3)-PrGridLines
           box2(14) = box2(3)+PrGridLines**2
           box2(15) = box2(3)-PrGridLines**2
           box2(4) = box2(1)+PrGridLines
           box2(16) = box2(4)+PrGridLines**2
           box2(17) = box2(4)-PrGridLines**2
           box2(5) = box2(1)-PrGridLines
           box2(18) = box2(5)+PrGridLines**2
           box2(19) = box2(5)-PrGridLines**2
           box2(6) = box2(1)+PrGridLines**2
           box2(7) = box2(1)-PrGridLines**2
           box2(20) = box2(8)+PrGridLines**2
           box2(21) = box2(8)-PrGridLines**2
           box2(22) = box2(9)+PrGridLines**2
           box2(23) = box2(9)-PrGridLines**2
           box2(24) = box2(12)+PrGridLines**2
           box2(25) = box2(12)-PrGridLines**2
           box2(26) = box2(13)+PrGridLines**2
           box2(27) = box2(13)-PrGridLines**2

           totalcellx((j*54)+1:(j*54)+27) = box1(1:27)
           totalcellx((j*54)+28:(j*54)+54) = box2(1:27)
           j = j+1
        end do
     end if


     totalcellx = remove_dups(totalcellx,size(totalcellx))

     !print*,rank, omp_get_thread_num(),'calculate in y-direction'
     if(Direction(2) >= 0) then
        xStart = ceiling(InitPoint(2)/(GridSize))
        xStop = floor(FinishPoint(2)/(GridSize))
     else 
        xStart = ceiling(FinishPoint(2)/(GridSize))
        xStop = floor(InitPoint(2)/(GridSize))
     end if
     !if(rank == 0 .and. omp_get_thread_num() == 0) then
     !                  write(*,*),'y', xStart, xStop
     !endif
     !write(*,*), xStart, xStop,iabs(xStop-xStart)+1
     
     if(int(InitPoint(1)/(GridSize)) == int(FinishPoint(1)/(GridSize))) then
        box1(1) = floor(InitPoint(1)/GridSize)*PrGridLines**2 + &
             floor(InitPoint(2)/GridSize)*PrGridLines + &
             floor(InitPoint(3)/GridSize)

        box1(2) = box1(1)+1
        box1(8) = box1(2)+PrGridLines
        box1(9) = box1(2)-PrGridLines
        box1(10) = box1(2)+PrGridLines**2
        box1(11) = box1(2)-PrGridLines**2
        box1(3) = box1(1)-1
        box1(12) = box1(3)+PrGridLines
        box1(13) = box1(3)-PrGridLines
        box1(14) = box1(3)+PrGridLines**2
        box1(15) = box1(3)-PrGridLines**2
        box1(4) = box1(1)+PrGridLines
        box1(16) = box1(4)+PrGridLines**2
        box1(17) = box1(4)-PrGridLines**2
        box1(5) = box1(1)-PrGridLines
        box1(18) = box1(5)+PrGridLines**2
        box1(19) = box1(5)-PrGridLines**2
        box1(6) = box1(1)+PrGridLines**2
        box1(7) = box1(1)-PrGridLines**2
        box1(20) = box1(8)+PrGridLines**2
        box1(21) = box1(8)-PrGridLines**2
        box1(22) = box1(9)+PrGridLines**2
        box1(23) = box1(9)-PrGridLines**2
        box1(24) = box1(12)+PrGridLines**2
        box1(25) = box1(12)-PrGridLines**2
        box1(26) = box1(13)+PrGridLines**2
        box1(27) = box1(13)-PrGridLines**2
        totalcelly(1:27) = box1(1:27)
     else


        j=0
        do i = xStart,xStop,1
           crossBlock(2) = real(i)*(GridSize)
           if (Direction(2) /= 0. ) then
              length = (crossBlock(2) - InitPoint(2))/Direction(2)
           else
              length = 0.
           endif
           crossBlock(1) = InitPoint(1) + length*Direction(1)
           crossBlock(3) = InitPoint(3) + length*Direction(3)
#ifdef DEBUG
           if(crossBlock(1)<0 .or. crossBlock(2)<0 .or. crossBlock(3) < 0) then
              print*,'y',crossBlock
              print*,'direction',Direction
              print*,'length',length
              print*,'start',xstart,'stop',xstop
              call abort
           endif
#endif


           box1(1) = floor(crossBlock(1)/GridSize)*PrGridLines**2 + &
                floor(crossBlock(2)/GridSize)*PrGridLines + &
                floor(crossBlock(3)/GridSize)

           box1(2) = box1(1)+1
           box1(8) = box1(2)+PrGridLines
           box1(9) = box1(2)-PrGridLines
           box1(10) = box1(2)+PrGridLines**2
           box1(11) = box1(2)-PrGridLines**2
           box1(3) = box1(1)-1
           box1(12) = box1(3)+PrGridLines
           box1(13) = box1(3)-PrGridLines
           box1(14) = box1(3)+PrGridLines**2
           box1(15) = box1(3)-PrGridLines**2
           box1(4) = box1(1)+PrGridLines
           box1(16) = box1(4)+PrGridLines**2
           box1(17) = box1(4)-PrGridLines**2
           box1(5) = box1(1)-PrGridLines
           box1(18) = box1(5)+PrGridLines**2
           box1(19) = box1(5)-PrGridLines**2
           box1(6) = box1(1)+PrGridLines**2
           box1(7) = box1(1)-PrGridLines**2
           box1(20) = box1(8)+PrGridLines**2
           box1(21) = box1(8)-PrGridLines**2
           box1(22) = box1(9)+PrGridLines**2
           box1(23) = box1(9)-PrGridLines**2
           box1(24) = box1(12)+PrGridLines**2
           box1(25) = box1(12)-PrGridLines**2
           box1(26) = box1(13)+PrGridLines**2
           box1(27) = box1(13)-PrGridLines**2
           box2(1) = floor(crossBlock(1)/GridSize)*PrGridLines**2 + &
                (floor(crossBlock(2)/GridSize)+1)*PrGridLines + &
                floor(crossBlock(3)/GridSize)

           box2(2) = box2(1)+1
           box2(8) = box2(2)+PrGridLines
           box2(9) = box2(2)-PrGridLines
           box2(10) = box2(2)+PrGridLines**2
           box2(11) = box2(2)-PrGridLines**2
           box2(3) = box2(1)-1
           box2(12) = box2(3)+PrGridLines
           box2(13) = box2(3)-PrGridLines
           box2(14) = box2(3)+PrGridLines**2
           box2(15) = box2(3)-PrGridLines**2
           box2(4) = box2(1)+PrGridLines
           box2(16) = box2(4)+PrGridLines**2
           box2(17) = box2(4)-PrGridLines**2
           box2(5) = box2(1)-PrGridLines
           box2(18) = box2(5)+PrGridLines**2
           box2(19) = box2(5)-PrGridLines**2
           box2(6) = box2(1)+PrGridLines**2
           box2(7) = box2(1)-PrGridLines**2
           box2(20) = box2(8)+PrGridLines**2
           box2(21) = box2(8)-PrGridLines**2
           box2(22) = box2(9)+PrGridLines**2
           box2(23) = box2(9)-PrGridLines**2
           box2(24) = box2(12)+PrGridLines**2
           box2(25) = box2(12)-PrGridLines**2
           box2(26) = box2(13)+PrGridLines**2
           box2(27) = box2(13)-PrGridLines**2
           !        mixbox_y = concat_array(box1,box2,27)
           totalcelly((j*54)+1:(j*54)+27) = box1(1:27)
           totalcelly((j*54)+28:(j*54)+54) = box2(1:27)

           j = j+1
        end do

     end if


     totalcelly = remove_dups(totalcelly,size(totalcelly))

     if(Direction(3) >= 0) then
        xStart = ceiling(InitPoint(3)/(GridSize))
        xStop = floor(FinishPoint(3)/(GridSize))
     else 
        xStart = ceiling(FinishPoint(3)/(GridSize))
        xStop = floor(InitPoint(3)/(GridSize))
     end if

     !write(*,*), xStart, xStop,iabs(xStop-xStart)+1

     !if(rank == 0 .and. omp_get_thread_num() == 0) then
     !                  write(*,*),'z', xStart, xStop
     !          endif

     if(int(InitPoint(1)/(GridSize)) == int(FinishPoint(1)/(GridSize))) then
        box1(1) = floor(InitPoint(1)/GridSize)*PrGridLines**2 + &
             floor(InitPoint(2)/GridSize)*PrGridLines + &
             floor(InitPoint(3)/GridSize)

        box1(2) = box1(1)+1
        box1(8) = box1(2)+PrGridLines
        box1(9) = box1(2)-PrGridLines
        box1(10) = box1(2)+PrGridLines**2
        box1(11) = box1(2)-PrGridLines**2
        box1(3) = box1(1)-1
        box1(12) = box1(3)+PrGridLines
        box1(13) = box1(3)-PrGridLines
        box1(14) = box1(3)+PrGridLines**2
        box1(15) = box1(3)-PrGridLines**2
        box1(4) = box1(1)+PrGridLines
        box1(16) = box1(4)+PrGridLines**2
        box1(17) = box1(4)-PrGridLines**2
        box1(5) = box1(1)-PrGridLines
        box1(18) = box1(5)+PrGridLines**2
        box1(19) = box1(5)-PrGridLines**2
        box1(6) = box1(1)+PrGridLines**2
        box1(7) = box1(1)-PrGridLines**2
        box1(20) = box1(8)+PrGridLines**2
        box1(21) = box1(8)-PrGridLines**2
        box1(22) = box1(9)+PrGridLines**2
        box1(23) = box1(9)-PrGridLines**2
        box1(24) = box1(12)+PrGridLines**2
        box1(25) = box1(12)-PrGridLines**2
        box1(26) = box1(13)+PrGridLines**2
        box1(27) = box1(13)-PrGridLines**2
        totalcellz(1:27) = box1(1:27)
     else


        j=0

        do i = xStart,xStop,1
           crossBlock(3) = real(i)*(GridSize)
           if(Direction(3) /= 0.) then
              length = (crossBlock(3) - InitPoint(3))/Direction(3)
           else
              length = 0.
           endif
           crossBlock(1:2) = InitPoint(1:2)+length*Direction(1:2)
#ifdef DEBUG
           if(crossBlock(1)<0 .or. crossBlock(2)<0 .or. crossBlock(3) < 0) then
              print*,'z',crossBlock
              print*,'direction',Direction
              print*,'length',length
              print*,'start',xstart,'stop',xstop
              call abort
           endif
#endif


           !write(*,*), floor(crossBlock(1)/GridSize), floor(crossBlock(2)/GridSize),floor(crossBlock(3)/GridSize)
           box1(1) = floor(crossBlock(1)/GridSize)*PrGridLines**2 + &
                floor(crossBlock(2)/GridSize)*PrGridLines + &
                floor(crossBlock(3)/GridSize)

           box1(2) = box1(1)+1
           box1(8) = box1(2)+PrGridLines
           box1(9) = box1(2)-PrGridLines
           box1(10) = box1(2)+PrGridLines**2
           box1(11) = box1(2)-PrGridLines**2
           box1(3) = box1(1)-1
           box1(12) = box1(3)+PrGridLines
           box1(13) = box1(3)-PrGridLines
           box1(14) = box1(3)+PrGridLines**2
           box1(15) = box1(3)-PrGridLines**2
           box1(4) = box1(1)+PrGridLines
           box1(16) = box1(4)+PrGridLines**2
           box1(17) = box1(4)-PrGridLines**2
           box1(5) = box1(1)-PrGridLines
           box1(18) = box1(5)+PrGridLines**2
           box1(19) = box1(5)-PrGridLines**2
           box1(6) = box1(1)+PrGridLines**2
           box1(7) = box1(1)-PrGridLines**2
           box1(20) = box1(8)+PrGridLines**2
           box1(21) = box1(8)-PrGridLines**2
           box1(22) = box1(9)+PrGridLines**2
           box1(23) = box1(9)-PrGridLines**2
           box1(24) = box1(12)+PrGridLines**2
           box1(25) = box1(12)-PrGridLines**2
           box1(26) = box1(13)+PrGridLines**2
           box1(27) = box1(13)-PrGridLines**2
           box2(1) = floor(crossBlock(1)/GridSize)*PrGridLines**2 + &
                floor(crossBlock(2)/GridSize)*PrGridLines + &
                floor(crossBlock(3)/GridSize)+1

           box2(2) = box2(1)+1
           box2(8) = box2(2)+PrGridLines
           box2(9) = box2(2)-PrGridLines
           box2(10) = box2(2)+PrGridLines**2
           box2(11) = box2(2)-PrGridLines**2
           box2(3) = box2(1)-1
           box2(12) = box2(3)+PrGridLines
           box2(13) = box2(3)-PrGridLines
           box2(14) = box2(3)+PrGridLines**2
           box2(15) = box2(3)-PrGridLines**2
           box2(4) = box2(1)+PrGridLines
           box2(16) = box2(4)+PrGridLines**2
           box2(17) = box2(4)-PrGridLines**2
           box2(5) = box2(1)-PrGridLines
           box2(18) = box2(5)+PrGridLines**2
           box2(19) = box2(5)-PrGridLines**2
           box2(6) = box2(1)+PrGridLines**2
           box2(7) = box2(1)-PrGridLines**2
           box2(20) = box2(8)+PrGridLines**2
           box2(21) = box2(8)-PrGridLines**2
           box2(22) = box2(9)+PrGridLines**2
           box2(23) = box2(9)-PrGridLines**2
           box2(24) = box2(12)+PrGridLines**2
           box2(25) = box2(12)-PrGridLines**2
           box2(26) = box2(13)+PrGridLines**2
           box2(27) = box2(13)-PrGridLines**2
           !  mixbox_z = concat_array(box1,box2,27)
           totalcellz((j*54)+1:(j*54)+27) = box1(1:27)
           totalcellz((j*54)+28:(j*54)+54) = box2(1:27)
           j = j+1
        end do

     end if
     totalcellz = remove_dups(totalcellz,size(totalcellz))
     do i=1,54
        if(totalcellz(i) < -1) then
           print*, crossBlock
           print*, Direction
           print*,totalcellz
           call abort
        end if
     end do

     !if(rank == 0 .and. omp_get_thread_num() == 0) then
     !        write(*,*),'merging x-y-z to dummy'
     !endif


     totalcell_dummy(1:size(totalcellx)) = totalcellx(1:size(totalcellx))
     totalcell_dummy(size(totalcellx)+1:size(totalcellx)+size(totalcelly)) = totalcelly(1:size(totalcelly))
     totalcell_dummy(size(totalcellx)+size(totalcelly)+1:size(totalcellx)+size(totalcelly)+size(totalcellz)) = totalcellz(1:size(totalcellz))
     totalcell_dummy = remove_dups(totalcell_dummy,size(totalcell_dummy))
     count_halo = 0
     !print*, totalcell_dummy
     i=1
     deallocate(totalcellx,totalcelly,totalcellz)
     !print*,totalcell_dummy(1:20)
     !if(totalcell_dummy(i) == -1) write(fh_logs,*) 'start dump',totalcell_dummy
     do while (totalcell_dummy(i) /= -1)
        !print*, totalcell_dummy(i)
#ifdef DEBUG
        if(totalcell_dummy(i) < -1) then
           print*, 'Error cell',totalcell_dummy,'is invalid'
           call abort
        endif
#endif 
        !write(fh_logs,*) 'cell',totalcell_dummy(i)
        abscell(1) = totalcell_dummy(i)/PrGridlines**2
        abscell(2) = (totalcell_dummy(i)-abscell(1)*PrGridlines**2)/(PrGridlines)
        abscell(3) = totalcell_dummy(i)-abscell(1)*PrGridlines**2-abscell(2)*PrGridlines

        shiftcell(1:3) = mod(abscell(1:3),GridLines)

        Numblock = shiftcell(1)*GridLines**2 + shiftcell(2)*GridLines+shiftcell(3)
        !print*, "block", totalcell_dummy(i), Numblock
#ifdef DEBUG
        if(Numblock < 0 .or. Numblock > GridLines**3-1) then
           print*,'Error: Block',Numblock,'Invalid'
           call abort
        endif
#endif 
        curHalo = headofchain(Numblock)

        shiftdistance(1:3) = BoxSize* real(abscell(1:3)/GridLines)

        !if(rank == 0 .and. omp_get_thread_num() == 0) then
        !           print*, abscell
        !           write(*,*),'numblock',Numblock,'curhalo',curHalo
        !   endif


        !curHalo = 100
        !print*, m,finddistance(initPoint,FinishPoint,MassivePos(1:3,m)),linedistance(initPoint,FinishPoint,MassivePos(1:3,m))

        do while (curHalo /= 0)

#ifdef DEBUG   
           if(curHalo > halonumber .or. curhalo < 0) then
              print*, 'Error: Halo ',curhalo,'is out of range'
              print*,'Dump',totalcell_dummy(:)
              call abort
           endif
#endif
           count_halo = count_halo + 1
           pos(1:3) = positions(1:3,curHalo) + shiftdistance(1:3) 
           !print*, rank, omp_get_thread_num(),positions(curHalo,1:3) + shiftdistance(1:3) 
           distancetoline = finddistance(initPoint,FinishPoint,pos)
           distanceonline = linedistance(initPoint,FinishPoint,pos)
           !if(pos(1) == InitPoint(1) .and. pos(2) == InitPoint(2) .and. pos(3) == InitPoint(3) ) then
           !if(distancetoline .lt. radius(curHalo)) then
           !   print*,l,m,curHalo,distanceonline,distancetoline
           !endif
!           if(curHalo == boydID) then
!              print*, "hasve minit halo"
!              print*, 'Block', Numblock
!              print*, 'GlobalBlock', totalcell_dummy(i)
!              print*, "abscell", abscell
!              print*, "shiftcell",abscell(1:3)/GridLines
!              print*, curHalo, distancetoline, distanceonline
!           end if

           if(distancetoline < radius(curHalo) .and. distanceonline <= BoxSize .and. distanceonline >= 0.) then
              haloperline = haloperline + 1

              !$omp critical
              call MPI_FILE_WRITE(fh_haloid, int(curHalo,4), 1, MPI_INTEGER, MPI_STATUS_IGNORE, ierr) 
              call MPI_FILE_WRITE(fh_online, distanceonline, 1, MPI_REAL, MPI_STATUS_IGNORE, ierr) 
              call MPI_FILE_WRITE(fh_lineid, int((l),4), 1, MPI_INTEGER, MPI_STATUS_IGNORE, ierr)
              call MPI_FILE_WRITE(fh_toline, distancetoline, 1, MPI_REAL, MPI_STATUS_IGNORE, ierr)
              call MPI_FILE_WRITE(fh_direction, direction, 3, MPI_REAL, MPI_STATUS_IGNORE, ierr)
              call MPI_FILE_WRITE(fh_startpoint, real(initpoint - (/ BoxSize, BoxSize, BoxSize/),4) , 3, MPI_REAL, MPI_STATUS_IGNORE, ierr)
              call MPI_FILE_WRITE(fh_hitpoint, direction(1:3)*(distanceonline)-(pos(1:3)-initpoint(1:3)), 3, MPI_REAL, MPI_STATUS_IGNORE, ierr)
              !$omp end critical

           endif

           curHalo = linkedlist(curHalo)
        end do
        i = i+1
     end do
152  continue
     deallocate(totalcell_dummy)
     !if(haloperline<1) 
     !print*, "total halo", halonumber, count_halo
     !print*, haloperline, "halos"
  end do
  !$omp end do
  !$omp end parallel


  call mpi_barrier(mpi_comm_world,ierr)        


  call MPI_FILE_CLOSE(fh_toline, ierr)
  call MPI_FILE_CLOSE(fh_online, ierr)
  call MPI_FILE_CLOSE(fh_lineid, ierr)
  call MPI_FILE_CLOSE(fh_haloid, ierr)
  call MPI_FILE_CLOSE(fh_hitpoint, ierr)
  call MPI_FILE_CLOSE(fh_direction, ierr)
  call MPI_FILE_CLOSE(fh_startpoint, ierr)

#ifdef separatedisk
  call MPI_FILE_CLOSE(fh_toline_disk, ierr)
  call MPI_FILE_CLOSE(fh_online_disk, ierr)
  call MPI_FILE_CLOSE(fh_lineid_disk, ierr)
  call MPI_FILE_CLOSE(fh_haloid_disk, ierr)
  call MPI_FILE_CLOSE(fh_hitpoint_disk, ierr)
  call MPI_FILE_CLOSE(fh_direction_disk, ierr)
#endif
  !close(fh_logs)
  !close(fh_centreref)

  deallocate(positions,radius,headofchain,linkedlist)
  return
#ifdef RR
end subroutine gen_los_rr
#else
end subroutine gen_los_rg
#endif

