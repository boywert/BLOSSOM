module common_vars
  use mpitools
  implicit none
  !########################################################################
  !@ SYSTEM params
  integer, parameter :: STRINGLEN = 256
  integer, parameter :: ZLIST_MAXSIZE = 500
  !########################################################################
  !@ Constant
  real(kind=8), parameter :: pi = 3.14159265359
  !########################################################################
  !@   TIS parameters
  real(kind=8), target :: Mttil 
  real(kind=8), target :: etaTIS
  real(kind=8), target :: etaSUS 
  real(kind=8), target :: zeta_t
  real(kind=8), target :: M_sol
  real(kind=8), target :: A10
  real(kind=8), target :: nu0
  real(kind=8), target :: T_star
  real(kind=8), target :: grav_const
  real(kind=8), target :: AMASS
  real(kind=8), target :: boltzk
  real(kind=8), target :: mu
  real(kind=8), target :: mu_H
  real(kind=8), target :: nH
  real(kind=8), target :: nHe
  integer, target :: m_tis 
  integer, target :: max_size 

  !#######################################################################
  !@ Cosmological parameters
  real(kind=8), target :: Omega_0 
  real(kind=8), target :: lambda0 
  real(kind=8), target :: c
  real(kind=8), target :: h 
  real(kind=8), target :: rho_crit_0

  real(kind=8), target :: Omega_b 
  real(kind=8), target :: Mpc 

  real(kind=8) :: H0 

  !#######################################################################
  !@ cubep3m parameter

  character(len=STRINGLEN), target :: los_path 
  logical , target :: rr
  character(len=STRINGLEN), target :: halo_path 
  character(len=STRINGLEN), target :: result_path 

  integer(kind=4), target :: file_dimension
  integer(kind=4) :: total_file 
  real(kind=4) :: BoxSize 
  real(kind=4) :: PrBoxSize
  real(kind=8), target :: co_boxwidth 
  integer(kind=4), target :: MassiveNumber 
  integer(kind=4), target ::  max_l
  integer(kind=8), target  :: particle_per_dim   

  !#######################################################################
  !@ Filenames

  character(len=STRINGLEN), target :: zlist_file
  character(len=STRINGLEN), target :: den_profile_file

  !#######################################################################
  !@ Unit conversion
  real(kind=8) :: dscale 
  real(kind=8) :: tscale 
  real(kind=8) :: lscale 
  real(kind=8) :: rho_matter 
  real(kind=8) :: M_box 
  real(kind=8) :: M_particle 
  real(kind=8) :: M_grid 

  !#######################################################################
  !@ for LOSs
  integer :: max_line
  integer :: line_use_to_save

  !#######################################################################
  !@ zlist
  real(kind=8) :: zlist(ZLIST_MAXSIZE)

#define common_var
end module common_vars


