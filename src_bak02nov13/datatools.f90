module datatools
  implicit none
contains
  ! subroutine findmassive (pos_in, mass, row, n, pos)
  !   use conversiontools
  !   use 
  !   implicit none
  !   integer(kind=4) :: row,n
  !   real (kind=4) :: pos_in(1:3,1:row),mass(1:row)
  !   real (kind=4) :: pos(1:3,1:n)
  !   integer (kind=4)  :: i,j, reforder
  !   real (kind=4) :: refmass
  !   do j=1,n
  !      refmass = 0.
  !      reforder = 0
  !      do i=1,row
  !         if(refmass < mass(i)) then
  !            refmass = mass(i)
  !            pos(1:3,j) = pos_in(1:3,i)
  !            reforder = i
  !         end if
  !      end do
  !      print*, "massive", j, convert_mass2physical(real(mass(reforder),8))/M_sol
  !      mass(reforder) = 1.0e-32
  !   end do
  !   return
  ! end subroutine findmassive
  subroutine readhalo(readnumber, filename, row, dataout)
    implicit none
    integer(kind=4), intent(in) :: row
    integer, intent(in) :: readnumber
    character (len=*) :: filename
    real (kind=4):: dataout(1:row, 1:17)
    integer :: i,j         
    open(readnumber, file=trim(filename), status='old')
    do i=1, row
       read(readnumber, *) (dataout(i,j), j=1,17)
    end do

    close (readnumber)
    return  
  end subroutine readhalo
  function getlinenumber(filename,readnumber)
    implicit none
    character (len=*) :: filename
    character (len=50) :: tmpfile,tmpstr
    integer(kind=4)  :: i,j,getlinenumber,readnumber
    integer :: tmpint,seed,time_array(8)
    CALL date_and_time(values=time_array)
    tmpint = time_array(8)*readnumber
    write(tmpstr, '(i10)') tmpint
    !write(*,*),tmpint,tmpstr,seed
    tmpfile = '/tmp/'//trim(adjustl(tmpstr))//'countline'
    call system('rm -f '//trim(adjustl(tmpfile)))
    call system('wc '//trim(filename)//' > '//trim(adjustl(tmpfile)))
    open (readnumber,file=trim(adjustl(tmpfile)),    status='old')
    read(readnumber,*) getlinenumber
    close(readnumber)  
    call system('rm -f '//trim(adjustl(tmpfile)))
  end function getlinenumber
  SUBROUTINE init_random_seed()
    INTEGER :: i, n, clock
    INTEGER, DIMENSION(:), ALLOCATABLE :: seed

    CALL RANDOM_SEED(size = n)
    ALLOCATE(seed(n))

    CALL SYSTEM_CLOCK(COUNT=clock)

    seed = clock + 37 * (/ (i - 1, i = 1, n) /)
    CALL RANDOM_SEED(PUT = seed)

  END SUBROUTINE init_random_seed

end module datatools
