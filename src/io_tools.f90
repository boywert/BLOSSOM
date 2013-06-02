module io_tools
  use mpitools
  use omp_lib
  use common_vars
  real(kind=4), allocatable :: tmpfloat(:,:)
contains
  integer function get_halo_number(z)
    implicit none
    character(len=20) :: z_s,str_i
    real(kind=4) :: z
    integer(kind=4) :: i,j,n
    integer :: buff(13)
    integer :: status
    write(z_s,'(f10.3)') z
    z_s = adjustl(z_s)
    get_halo_number = 0
    do i=0, total_file-1
       write(str_i,'(i10)') i
       str_i = adjustl(str_i)
       n =  (get_file_size(trim(halo_path)//z_s(1:len_trim(z_s))//'halo'//str_i(1:len_trim(str_i))//'.dat')-4)/4/17
       if(n .eq. 0) then
          get_halo_number = 0
          return
       end if
       get_halo_number = get_halo_number + n
    end do
  end function get_halo_number

  function ascii_count_line(z)
    implicit none
    character (len=512) :: filename
    character (len=50) :: tmpfile,tmpstr,z_s
    integer  :: i,j,ascii_count_line,readnumber
    integer :: tmpint,seed,time_array(8)
    real(kind=4) :: z

    write(z_s,'(f10.3)') z
    z_s = adjustl(z_s)
    filename =trim(halo_path)//z_s(1:len_trim(z_s))//'halo.dat'
    readnumber = 30
    CALL date_and_time(values=time_array)
    tmpint = time_array(8)*readnumber
    write(tmpstr, '(i10)') tmpint
    !write(*,*),tmpint,tmpstr,seed
    tmpfile = '/tmp/'//trim(adjustl(tmpstr))//'countline'
    call system('rm -f '//trim(adjustl(tmpfile)))
    call system('wc '//trim(filename)//' > '//trim(adjustl(tmpfile)))
    open (readnumber,file=trim(adjustl(tmpfile)), status='old')
    read(readnumber,*) ascii_count_line
    close(readnumber)  
    call system('rm -f '//trim(adjustl(tmpfile)))
  end function ascii_count_line

  integer function get_file_size(filename)
    implicit none
    character(len=100) :: filename
    integer :: buff(13)
    filename = adjustl(filename)
    call stat(filename(1:len_trim(filename)),buff)
    get_file_size = buff(8)
  end function get_file_size


  subroutine ascii_read_halo(z,element_flag,n_elements,halonumber)
    implicit none
    integer(kind=4) :: halonumber
    integer  :: i,j,k,l,n_elements
    character (len=512) :: filename
    logical :: element_flag(1:17)  
    real (kind=4):: dataout(1:17)
    character(len=20) :: z_s
    real :: z

    write(z_s,'(f10.3)') z
    z_s = adjustl(z_s)

    allocate(tmpfloat(1:n_elements,1:maxhalos_persnap))

    filename =trim(halo_path)//z_s(1:len_trim(z_s))//'halo.dat'
    filename = adjustl(filename)
    print*, "STart to read ASCII",filename
    open(28, file=trim(filename), status='old')

    do i=1, maxhalos_persnap
       read(28,fmt = *,end=131) (dataout(j), j=1,17)
       l=1
       do k=1,17
          if(element_flag(k)) then
             tmpfloat(l,i) = dataout(k)
             l = l+1
          endif
       end do
    end do

131 close (28)
    halonumber = i-1
    return  
  end subroutine ascii_read_halo


  function checkfiles(z)
    implicit none
    real(kind=4) :: z
    integer :: status, checkfiles,i
    logical :: file_e
    character(len=20) :: z_s
    character(len=100) :: str_rank
    write(z_s,'(f10.3)') z
    z_s = adjustl(z_s)

    if(rank==0) then
       status = 2 ! set use mpi read as default
       do i=0,total_file-1
          write(str_rank,'(i10)') i
          str_rank = adjustl(str_rank)
          inquire( file=trim(halo_path)//z_s(1:len_trim(z_s))//'halo'//str_rank(1:len_trim(str_rank))//'.dat', exist=file_e )
          if(file_e == .FALSE.) then
             print*, 'file', i, 'does not exist'
             status = 0  ! set to have none file
             goto 129
          endif
       end do
       goto 130
       
129    inquire( file=trim(halo_path)//z_s(1:len_trim(z_s))//'halo.dat', exist=file_e )
       
       if(file_e == .FALSE.) then
          print*, 'Full catalogue'//trim(halo_path)//z_s(1:len_trim(z_s))//'halo.dat '//'does not exist.... aborting'
          call abort
       else
          status = 1   ! set to use full halo catalogue ASCII
       endif
130 continue 
    endif
    call mpi_bcast(status,1,mpi_integer,0,mpi_comm_world,ierr)
    checkfiles = status
  end function checkfiles

  subroutine mpi_read_halo(z,element_flag,n_elements,halo_number, halodata)
    use mpi
    use mpitools
    implicit none   
    logical :: element_flag(1:17)    
    integer(kind=4)   :: fh_readhalo
    integer(kind=4) :: n_elements

    integer(kind=4) :: i,j,k,l,m,n,halo_number
    real(kind=4) :: z

    real(kind=4) :: halodata(1:n_elements,1:halo_number)

    character(len=100) :: str_rank,filename
    character(len=20) :: z_s
    integer(kind=4) :: n_point,totalpoint,tag,n_point_index
    integer(kind=4),allocatable :: n_point_arr(:)
    real(kind=4), allocatable :: halo_in(:,:),halo_cache(:,:)
    integer(4), dimension(mpi_status_size) :: status
    real(kind=4) :: halo_dummy(17)
    integer(kind=4), allocatable :: file_read_per_node(:),n_point_file(:)
    logical :: file_e
    !! Initialise MPI

    write(z_s,'(f10.3)') z
    z_s = adjustl(z_s)


    if(rank==0) then
       do i=0,total_file-1
          write(str_rank,'(i10)') i
          str_rank = adjustl(str_rank)
          inquire( file=trim(halo_path)//z_s(1:len_trim(z_s))//'halo'//str_rank(1:len_trim(str_rank))//'.dat', exist=file_e )
          if(file_e == .FALSE.) then
             print*, 'file', i, 'does not exist.... aborting'
             call abort
          endif
       end do
    endif

    if(rank == 0) allocate(n_point_arr(0:nodes_returned-1))
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    allocate(file_read_per_node(0:nodes_returned-1))
    file_read_per_node = 0

    if(nodes_returned < total_file) then
       file_read_per_node(0:mod(total_file,nodes_returned)-1) = total_file/nodes_returned + 1
       file_read_per_node(mod(total_file,nodes_returned):nodes_returned-1) = total_file/nodes_returned
    else 
       file_read_per_node(0:total_file-1) = 1
    endif
    do i=0,nodes_returned-1
       !if(rank == 0) print*,i,'files/node',file_read_per_node(i)
    end do
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    if(rank < total_file) then
       n_point = 0
       allocate(n_point_file(1:file_read_per_node(rank)))
       do j=1,file_read_per_node(rank)
          write(str_rank,'(i10)') (j-1)*nodes_returned + rank
          str_rank = adjustl(str_rank)  
          n_point_file(j) = (get_file_size(trim(halo_path)//z_s(1:len_trim(z_s))//'halo'//str_rank(1:len_trim(str_rank))//'.dat')-4)/4/17
          n_point = n_point + n_point_file(j)
       end do
       !read only positio(3), radius, mass
       allocate(halo_in(1:n_elements,1:n_point))
       n_point_index = 1
       !print*,'rank',rank,'starts to read file'
       do j=1,file_read_per_node(rank)
          write(str_rank,'(i10)') (j-1)*nodes_returned + rank
          str_rank = adjustl(str_rank)  
          fh_readhalo = 52
          open (unit=fh_readhalo, &
               file=trim(halo_path)//z_s(1:len_trim(z_s))//'halo'//str_rank(1:len_trim(str_rank))//'.dat', &
               status='old', &
               form='binary')
          read(unit=fh_readhalo) n
          allocate(halo_cache(1:17,1:n_point_file(j)))
          read(unit=fh_readhalo) halo_cache(1:17,1:n_point_file(j))
          close(fh_readhalo)
          l=1
          do k=1,17
             if(element_flag(k)) then
                halo_in(l,n_point_index:n_point_index+n_point_file(j)-1) = halo_cache(k,1:n_point_file(j))
                l = l+1
             endif
          end do
          deallocate(halo_cache)
          !do i=n_point_index,n_point_index+n_point_file(j)-1
          !   read(unit=fh_readhalo) halo_dummy(1:17)
          !   l=1
          !   do k=1,17
          !      if(element_flag(k)) then
          !         halo_in(l,i) = halo_dummy(k)
          !         l = l+1
          !      endif
          !   end do

          !enddo
          n_point_index = n_point_index + n_point_file(j)
          !close(fh_readhalo)
       end do

       deallocate(n_point_file)
    endif
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    deallocate(file_read_per_node)
    if(rank==0) print*, 'Finish reading data'
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    call mpi_gather(n_point,1,mpi_integer,n_point_arr,1,mpi_integer,0, &
         mpi_comm_world,ierr)
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)

    if(rank ==0) then

       halo_number = 0
       do i=0, min(total_file,nodes_returned)-1
          halo_number = halo_number + n_point_arr(i)
       enddo
       totalpoint = 1 
       halodata(1:n_elements,1:n_point) = halo_in(1:n_elements,1:n_point)
       totalpoint = totalpoint + n_point
    endif
    if(rank==0) print*, 'Start transfering to node 0'

    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    
    do i=1,min(total_file,nodes_returned)-1
       tag=i

       if (rank == 0) then
          !print*, '.. recieving',n_point_arr(i),' point from',i
          call mpi_recv(halodata(1:n_elements,totalpoint:totalpoint+n_point_arr(i)-1),n_elements*n_point_arr(i),mpi_real, &
               i,tag,mpi_comm_world,status,ierr)
          !print*,halodata(totalpoint+n_point_arr(i)-1,1:n_elements)
       elseif (rank == i) then
          !print*, '.. transfering',n_point,' points from',i
          !print*,halo_in(n_point,1:n_elements)
          call mpi_send(halo_in(1:n_elements,1:n_point),n_elements*n_point,mpi_real,0,tag,mpi_comm_world,ierr)
       endif
       call mpi_barrier(mpi_comm_world,ierr)
       if (rank ==0)  totalpoint = totalpoint + n_point_arr(i)
       call mpi_barrier(mpi_comm_world,ierr)
    enddo
    if(rank==0) print*, 'Finish transfering'


    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    if(rank < total_file ) deallocate(halo_in)
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)

    return
  end subroutine mpi_read_halo

  integer function get_n_elements(element_flag)
    implicit none
    logical :: element_flag(17)
    integer :: i
    get_n_elements = 0
    do i=1,17
       if(element_flag(i)) get_n_elements = get_n_elements + 1
    end do
  end function get_n_elements

end module io_tools
