module runprocs
  implicit none
contains

!#include "massfunction.f90"
#include "subcell_corr.f90"

#include "gen_los.f90"
!#include "makecorrelation.f90"
!#include "absorpstats.f90"
#include "makeobservedline.f90"

#define RR
#include "gen_los.f90"
!#include "makecorrelation.f90"
#include "makeobservedline.f90"
#include "gen_random.f90"

end module runprocs
