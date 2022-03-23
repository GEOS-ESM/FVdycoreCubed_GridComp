#define I_AM_MAIN

#include "MAPL_Generic.h"

program StandAlone_DynAdvCore
   use MAPL
  use StandAlone_DynAdvCore_GridCompMod, only: SetServices
   use MPI

   implicit none

!EOP

!EOC

   character(*), parameter :: IAM = __FILE__

   type (MAPL_Cap) :: cap
   type (MAPL_CapOptions) :: cap_options
   integer :: status

   cap_options = FlapCLI( &
        description = 'FV Standalone DyAdvCore', &
        authors     = 'S.J. Lin, R. Rood, W. Putman')
   cap = MAPL_Cap('Standalone FV3 DynAdvCore', SetServices, cap_options = cap_options)
   call cap%run(_RC)

 end Program StandAlone_DynAdvCore

!EOC




