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
end module arraytools
