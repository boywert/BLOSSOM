module runprocs
  implicit none
contains

#include "massfunction.f90"

#include "gen_los.f90"
#include "makecorrelation.f90"
#include "absorpstats.f90"

#define RR
#include "gen_los.f90"
#include "makecorrelation.f90"

end module runprocs
