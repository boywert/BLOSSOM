module read_parameters
  use mpitools
  use omp_lib
  use common_vars
  implicit none
  integer,parameter :: max_tags = 300
  integer :: n
  character(len=STRINGLEN) :: tag(max_tags)
  integer :: var_type(max_tags) !string=0, int4=1, int8=2,real4=3,real8=4, bool=5
  logical :: checked(max_tags)

  type container
     integer(kind=4), pointer :: int4
     integer(kind=8), pointer :: int8
     real(kind=4), pointer :: real4
     real(kind=8), pointer :: real8
     character(len=STRINGLEN),pointer :: str
     logical, pointer :: bool 
  end type container

  type(container), allocatable :: var_pointer(:)

contains
  subroutine allocate_pointers()
    implicit none
    allocate(var_pointer(1:max_tags))
  end subroutine allocate_pointers

  subroutine read_zlist()
    implicit none
    integer :: i
    i = 1
    zlist(:) = -1.0
    open(unit=20,file=trim(adjustl(zlist_file)),status="old")
    do
       read(20,'(e12.4)',end=21) zlist(i)
       i = i+1
       if(i > ZLIST_MAXSIZE) call abort("FATAL: zlist is too large")
    end do
21  close(20)
  end subroutine read_zlist

  subroutine read_config (filename)
    implicit none
    character(len=STRINGLEN) :: filename
    character(len=STRINGLEN) :: buff
    call allocate_pointers
    call list_toread_vars

    checked(:) = .FALSE.
    
    if(rank ==0) then
       print*, ""
       print*, "#######################################"
       print*, "READ PARAMETERS FROM CONFIGURATION FILE"
       print*, "#######################################"
       print*, ""
       print*, "READ CONFIG FILE"//achar(9)// trim(adjustl(filename))
       print*, ""
    end if
    open(unit=10,file=trim(adjustl(filename)),status="old")
    do
       read(10,12,end=11)  buff
       call check_vars(buff)
    end do
12  format(A256)
11  close(10)
    
    call calculate_variables
  end subroutine read_config

  subroutine list_toread_vars()
    implicit none
    n = 1

    tag(n) = 'zlist_file'
    var_type(n) = 0
    var_pointer(n)%str => zlist_file

    n=n+1
    tag(n) = 'den_profile_file'
    var_type(n) = 0
    var_pointer(n)%str => den_profile_file

    n=n+1
    tag(n) = 'subcell_perdim_file'
    var_type(n) = 0
    var_pointer(n)%str => subcell_perdim_file

    n=n+1
    tag(n) = 'Omega_0'
    var_type(n) = 4
    var_pointer(n)%real8 => Omega_0

    n=n+1
    tag(n) = 'Lambda0'
    var_type(n) = 4
    var_pointer(n)%real8 => Lambda0

    n=n+1
    tag(n) = 'h'
    var_type(n) = 4
    var_pointer(n)%real8 => h

    n=n+1
    tag(n) = 'rho_crit_0'
    var_type(n) = 4
    var_pointer(n)%real8 => rho_crit_0

    n=n+1
    tag(n) = 'Omega_b'
    var_type(n) = 4
    var_pointer(n)%real8 => Omega_b

    n=n+1
    tag(n) = 'file_dimension'
    var_type(n) = 1
    var_pointer(n)%int4 => file_dimension

    n=n+1
    tag(n) = 'line_length_factor'
    var_type(n) = 1
    var_pointer(n)%int4 => line_length_factor

    n=n+1
    tag(n) = 'particle_per_dim'
    var_type(n) = 2
    var_pointer(n)%int8 => particle_per_dim

    n=n+1
    tag(n) = 'co_boxwidth'
    var_type(n) = 4
    var_pointer(n)%real8 => co_boxwidth

    n=n+1
    tag(n) = 'MaxSourceSize'
    var_type(n) = 4
    var_pointer(n)%real8 => MaxSourceSize
   
    n=n+1
    tag(n) = 'los_path'
    var_type(n) = 0
    var_pointer(n)%str => los_path

    n=n+1
    tag(n) = 'halo_path'
    var_type(n) = 0
    var_pointer(n)%str => halo_path

    n=n+1
    tag(n) = 'subcellden_path'
    var_type(n) = 0
    var_pointer(n)%str => subcellden_path

    n=n+1
    tag(n) = 'result_path'
    var_type(n) = 0
    var_pointer(n)%str => result_path

    n=n+1
    tag(n) = 'databin_path'
    var_type(n) = 0
    var_pointer(n)%str => databin_path


    n=n+1
    tag(n) = 'randomhalo_path'
    var_type(n) = 0
    var_pointer(n)%str => randomhalo_path


    n=n+1
    tag(n) = 'rr'
    var_type(n) = 5
    var_pointer(n)%bool => rr

    n=n+1
    tag(n) = 'MassiveNumber'
    var_type(n) = 1
    var_pointer(n)%int4 => MassiveNumber

    n=n+1
    tag(n) = 'Max_l'
    var_type(n) = 1
    var_pointer(n)%int4 => Max_l

    n=n+1
    tag(n) = 'first_l'
    var_type(n) = 1
    var_pointer(n)%int4 => first_l


    n=n+1
    tag(n) = 'last_l'
    var_type(n) = 1
    var_pointer(n)%int4 => last_l

    n=n+1
    tag(n) = 'obsfreqresolution'
    var_type(n) = 4
    var_pointer(n)%real8 => obsfreqresolution

    n=n+1
    tag(n) = 'maxfreqresolution'
    var_type(n) = 4
    var_pointer(n)%real8 => maxfreqresolution

    n=n+1
    tag(n) = 'Mttil'
    var_type(n) = 4
    var_pointer(n)%real8 => Mttil

    n=n+1
    tag(n) = 'etaTIS'
    var_type(n) = 4
    var_pointer(n)%real8 => etaTIS

    n=n+1
    tag(n) = 'etaSUS'
    var_type(n) = 4
    var_pointer(n)%real8 => etaSUS

    n=n+1
    tag(n) = 'zeta_t'
    var_type(n) = 4
    var_pointer(n)%real8 => zeta_t

    n=n+1
    tag(n) = 'A10'
    var_type(n) = 4
    var_pointer(n)%real8 => A10

    n=n+1
    tag(n) = 'nu0'
    var_type(n) = 4
    var_pointer(n)%real8 => nu0

    n=n+1
    tag(n) = 'T_star'
    var_type(n) = 4
    var_pointer(n)%real8 => T_star

    n=n+1
    tag(n) = 'grav_const'
    var_type(n) = 4
    var_pointer(n)%real8 => grav_const

    n=n+1
    tag(n) = 'AMASS'
    var_type(n) = 4
    var_pointer(n)%real8 => AMASS

    n=n+1
    tag(n) = 'boltzk'
    var_type(n) = 4
    var_pointer(n)%real8 => boltzk

    n=n+1
    tag(n) = 'nH'
    var_type(n) = 4
    var_pointer(n)%real8 => nH

    n=n+1
    tag(n) = 'nHe'
    var_type(n) = 4
    var_pointer(n)%real8 => nHe

    n=n+1
    tag(n) = 'mu'
    var_type(n) = 4
    var_pointer(n)%real8 => mu

    n=n+1
    tag(n) = 'mu_H'
    var_type(n) = 4
    var_pointer(n)%real8 => mu_H

    n=n+1
    tag(n) = 'b_T'
    var_type(n) = 4
    var_pointer(n)%real8 => b_T

    n=n+1
    tag(n) = 'm_tis'
    var_type(n) = 1
    var_pointer(n)%int4 => m_tis

    n=n+1
    tag(n) = 'max_size'
    var_type(n) = 1
    var_pointer(n)%int4 => max_size

    n=n+1
    tag(n) = 'c'
    var_type(n) = 4
    var_pointer(n)%real8 => c

    n=n+1
    tag(n) = 'Mpc'
    var_type(n) = 4
    var_pointer(n)%real8 => Mpc

    n=n+1
    tag(n) = 'M_sol'
    var_type(n) = 4
    var_pointer(n)%real8 => M_sol




  end subroutine list_toread_vars

  subroutine calculate_variables
    implicit none
    integer :: i
    !@ Cosmological params
    rho_crit_0 = rho_crit_0 * h**2.d0
    Omega_b = Omega_b/ h**2.d0
    H0 = h*1.d7/Mpc   

    !@ CubeP3M
    total_file = file_dimension**3
    BoxSize = 2.*real(particle_per_dim)
    PrBoxSize = (2*line_length_factor+1)*BoxSize

    !@ unit convertion tools
    kpc = Mpc/1000.0
    
    
    dscale = omega_0*rho_crit_0 
    tscale = 2./(3.*sqrt(omega_0)*H0) ! s
    lscale = co_boxwidth/h/Boxsize*Mpc ! cm
    rho_matter = dscale
    M_box = rho_matter*(co_boxwidth*Mpc/h)**3. ! g
    M_particle = M_box/(Boxsize/2.)**3. ! g
    M_grid = M_box/(Boxsize)**3. !g
    
    !@ LOS line
    max_line = max_l !* MassiveNumber
    line_use_to_save = 0 !use 0 as default

    !@ TIS model
    
    !@ Files and Folders

    !@ OMP parameters
    omp_workthreads = omp_get_max_threads()
    
    if(rank==0) then
       print*, ""
       print*, "#######################################"
       print*, "COMPUTED PARAMETERS"
       print*, "#######################################"
       print*, ""

       print*, "dscale : ", dscale
       print*, "tscale : ", tscale
       print*, "lscale : ", lscale
       print*, "M_box : ", M_box
       print*, "M_particle : ", M_particle
       print*, "M_grid : ", M_grid

       print*, ""

       print*, "total_file : ", total_file
       print*, "Boxsize : ", Boxsize
       print*, "PrBoxSize : ", PrBoxSize

       print*, ""

       print*, "rho_crit_0 : ", rho_crit_0
       print*, "Omega_b : ", Omega_b
       print*, "H0 : ", H0
    end if
  end subroutine calculate_variables


  subroutine check_vars (buff)
    implicit none
    character(len=STRINGLEN) :: buff
    character(len=1) :: ext
    integer :: i,tag_length,iden


    buff = adjustl(buff)
    !@ remove everything after #
    iden = index(buff,"#")
    if(iden > 0) then
       buff(iden:STRINGLEN) = ""
       buff = adjustl(buff)
    end if

    if(len_trim(buff) == 0) return
    !@ remove tab and space in the front
    do i=1, STRINGLEN
       !@ variable starts with A-Z a-z only
       if(ichar(buff(i:i)) >= 65  .and. ichar(buff(i:i)) <= 122) goto 13
    end do
13  buff(1:i-1) = ""
    buff = adjustl(buff)

    do i=1,n
       tag(i) = adjustl(tag(i))
       tag_length = len_trim(tag(i))
       ext = adjustl(buff(tag_length+1:tag_length+1))

       !@ check for tab or space
       if(buff(1:tag_length) == trim(tag(i)) .and. (ichar(ext) == 9 .or. ichar(ext)==32)) then
          buff(1:tag_length) = ""
          buff = adjustl(buff)
          buff = removespace(buff)

          if(var_type(i) == 0) then
             var_pointer(i)%str = buff
             checked(i) = .TRUE.
             if(rank==0) print*, trim(tag(i)),achar(9),achar(9),trim(var_pointer(i)%str)
          else if (var_type(i) == 1) then
             read(buff, '(i20)') var_pointer(i)%int4
             checked(i) = .TRUE.
             if(rank==0) print*, trim(tag(i)),achar(9),achar(9),var_pointer(i)%int4
          else if (var_type(i) == 2) then
             read(buff, '(i20)') var_pointer(i)%int8
             checked(i) = .TRUE.
             if(rank==0) print*, trim(tag(i)),achar(9),achar(9),var_pointer(i)%int8
          else if (var_type(i) == 3) then
             read(buff, '(e15.7)') var_pointer(i)%real4
             checked(i) = .TRUE.
             if(rank==0) print*, trim(tag(i)),achar(9),achar(9),var_pointer(i)%real4
          else if (var_type(i) == 4) then
             read(buff, '(e15.7)') var_pointer(i)%real8
             checked(i) = .TRUE.
             if(rank==0) print*, trim(tag(i)),achar(9),achar(9),var_pointer(i)%real8
          else if (var_type(i) == 5) then
             if (trim(buff) == '0') then
                var_pointer(i)%bool = .FALSE.
             else
                var_pointer(i)%bool = .TRUE.
             end if
             checked(i) = .TRUE.
             if(rank==0) print*, trim(tag(i)),achar(9),achar(9),var_pointer(i)%bool
          else 
             return
          end if

          return
       end if
    end do
    return
  end subroutine check_vars

  function removespace(buff)
    implicit none
    character(len=STRINGLEN) :: buff,removespace
    integer :: i, len
    len = len_trim(adjustl(buff))
    do i=1,len
       if(ichar(buff(i:i)) == 9 .or. ichar(buff(i:i)) == 32) then
          buff(i:i) = ""
       endif
    end do
    removespace = trim(adjustl(buff))
  end function removespace
end module read_parameters
