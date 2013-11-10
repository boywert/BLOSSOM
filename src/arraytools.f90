module arraytools
  implicit none
contains 
  function remove_dups(input,N, PrGridLines)
    implicit none
    integer(kind=4) :: N
    integer(kind=4) :: i,j,k,l,PrGridLines
    integer(kind=4) :: remove_dups(N),input(N)

    input = remove_excess(input,N,PrGridLines)
    k = 1
    remove_dups(1:N) = -1

    remove_dups(1) = input(1)
    outer: do i=2,N
       do j=1,k
          if(remove_dups(j)==input(i)) then
             cycle outer

          end if
       end do
       k = k+1
       remove_dups(k) = input(i) 
    end do outer
    do i=1,N
       if(remove_dups(i) == -1) then
          remove_dups(i:N-1) = remove_dups(i+1:N)  
          goto 101
       endif
    enddo
101 continue
  end function remove_dups
  
  function remove_excess(input,N,PrGridLines)
    implicit none
    integer(kind=4) :: N
    integer(kind=4) :: i,PrGridLines
    integer(kind=4) :: remove_excess(N),input(N)
    do i=1,N
       if(input(i) < -1 .or. input(i) > PrGridLines**3-1) then
          input(i) = -1
       end if
    end do
    remove_excess(1:N) = input(1:N)
  end function remove_excess
  
  ! interpolate arrays, x_in,x_out must be increasing series
  subroutine array_intrpol(x_in,y_in,N_in,x_out,y_out,N_out)
    implicit none
    integer :: N_in, N_out
    real(kind=8) :: x_in(1:N_in),y_in(1:N_in),x_out(1:N_out),y_out(1:N_out)
    real(kind=8) :: ref_x
    integer :: i,j
   
    if(x_out(1) < x_in(1)) then
       print*,"x_out is out of range"
       call exit()
    end if
    if(x_out(N_out) > x_in(N_in)) then
       print*,"x_out is out of range"
       call exit()
    end if

    j = 1
    do i=1,N_out
       ref_x = x_out(i)
       do while (ref_x .ge. x_in(j))
          j = j+1
       end do
       y_out(i) = y_in(j-1) + (y_in(j)-y_in(j-1))/(x_in(j) - x_in(j-1))*(ref_x - x_in(j-1))
    end do
    return
  end subroutine array_intrpol
end module arraytools
