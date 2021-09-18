#define I_AM_MAIN

#include "MAPL_Generic.h"
#define _RC rc=status); _VERIFY(status

program CreateInterpWeights
   use MAPL
   use CreateInterpWeights_GridCompMod,      only: SetServices
   implicit none

!EOP

!EOC

   character(*), parameter :: IAM = __FILE__

   type (MAPL_Cap) :: cap
   type (MAPL_FlapCLI) :: cli
   type (MAPL_CapOptions) :: cap_options
   integer :: status

   cli = MAPL_FlapCLI(description = 'Standalone Interp Weights',&
                              authors      =  'W. Putman')
   cap_options = MAPL_CapOptions(cli)
   cap = MAPL_Cap('GCM', SetServices, cap_options = cap_options)
   call cap%run(_RC)

 end Program CreateInterpWeights

