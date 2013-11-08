module runprocs
  use mod_subcell_corr
  implicit none
contains

!#include "massfunction.f90"

#include "gen_los.f90"
#include "makeobservedline.f90"

#define RR
#include "gen_los.f90"
#include "makeobservedline.f90"
#include "gen_random.f90"

end module runprocs
