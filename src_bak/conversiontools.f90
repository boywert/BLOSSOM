module conversiontools
  use common_vars
  implicit none

contains

  function z_to_d(z)
    implicit none
    real(kind=8) :: z_to_d, z
    z_to_d = c/H0*((1.d0+z)**2.-1.d0)/((1.d0+z)**2.+1.d0)
  end function z_to_d

  function d_to_z(d)
    implicit none
    real(kind=8) :: d_to_z,d
    d_to_z = sqrt(1.d0+H0*d/c)/sqrt(1.d0-H0*d/c) - 1.d0
  end function d_to_z
  
  function d_to_nu(d)
    implicit none
    real(kind=8) :: d, d_to_nu
    d_to_nu = nu0*sqrt(c-H0*d)/sqrt(c+H0*d)
  end function d_to_nu

  function nu_to_d(nu)
    implicit none
    real(kind=8) :: nu, nu_to_d
    nu_to_d = c/H0*(1-(nu/nu0)**2.)/(1+(nu/nu0)**2.)
  end function nu_to_d

  function convert_mass2physical(grid_mass)
    implicit none
    real(kind=8) :: grid_mass,convert_mass2physical
    convert_mass2physical = M_grid*grid_mass
  end function convert_mass2physical

  function convert_length2physical(length,z)
    implicit none
    real(kind=8) :: convert_length2physical,z,length
    convert_length2physical = length * lscale / (1.+z)
  end function convert_length2physical

  function convert_vel2physical(vel,z)
    implicit none
    real(kind=8) :: convert_vel2physical,vel,z
    convert_vel2physical = vel*lscale/tscale *(1.+z)
  end function convert_vel2physical

  subroutine run_unit_test(z)
    implicit none
    real(kind=8) :: z
    print*, 'Omega_0:' , Omega_0
    print*, 'Omega_b:' , Omega_b
    print*, 'h:', h
    print*, 'Lambda_0:', lambda0
    print*, 'H0:', H0, '1/s'
    print*, 'Boxsize:' , co_boxwidth, 'Mpc/h'
    print*, 'Boxgrid per dimension:', int(boxsize,8)
    print*, '1u velocity:',convert_vel2physical(1.d0,z), 'cm/s'
    print*, '1u length:', convert_length2physical(1.d0,z), 'cm'
    print*, '1u mass:', convert_mass2physical(1.d0)/M_sol,'M_sol'
  end subroutine run_unit_test
end module conversiontools
