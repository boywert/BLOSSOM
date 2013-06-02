program test
  implicit none
  character(len=1024) :: filename
  
  filename = "inputs/Config"
  call getarg(1,filename)
  filename = adjustl(filename)
  print*, iargc(), trim(filename)

end program test

