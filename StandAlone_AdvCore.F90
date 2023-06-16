#define I_AM_MAIN

#include "MAPL_Generic.h"

program StandAlone_AdvCore
   use MAPL
   use AdvCore_GridCompMod, only: SetServices
   use MPI

   implicit none

!EOP

!EOC

   character(*), parameter :: IAM = __FILE__

   type (MAPL_Cap) :: cap
   type (MAPL_FargparseCLI) :: cli
   type (MAPL_CapOptions) :: cap_options
   integer :: status

   cli = MAPL_FargparseCLI()

   cap_options = MAPL_CapOptions(cli)
   cap = MAPL_Cap('Standalone FV3 AdvCore', SetServices, cap_options=cap_options)
   call cap%run(_RC)

end program StandAlone_AdvCore

!EOC


