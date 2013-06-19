module common_vars
  use mpitools
  implicit none
  !########################################################################
  !@ SYSTEM params
  integer, parameter :: STRINGLEN = 256
  integer, parameter :: ZLIST_MAXSIZE = 500
  integer, parameter :: SUBCELLLIST_MAXSIZE = 500
  integer, parameter :: maxhalos_persnap = 200000000
  !########################################################################
  !@ OMP
  integer :: omp_workthreads
  !########################################################################
  !@ Correlation function
  integer, parameter :: maxpoint_serial_xi = 500000
  real, parameter :: min_radius_xi = 0.0001
  real, parameter :: max_radius_xi = (3.)**(0.5)
  integer, parameter :: n_bin_xi = 10
  integer, parameter :: n_smallcell = 25
  !########################################################################
  !@ Overdensity
  integer, parameter :: overden_nbin = 25
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
  character(len=STRINGLEN), target :: subcellden_path
  character(len=STRINGLEN), target :: databin_path
  character(len=STRINGLEN), target :: randomhalo_path

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
  character(len=STRINGLEN), target :: subcell_perdim_file

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
  integer, target :: line_length_factor
  real(kind=8), target :: MaxSourceSize
  !#######################################################################
  !@ zlist
  real(kind=8) :: zlist(ZLIST_MAXSIZE)
  !@ subcell_perdim list
  integer(kind=4) :: subcell_perdim_list(SUBCELLLIST_MAXSIZE)
  integer, parameter :: min_particle_cell = 0

#define common_var
end module common_vars


