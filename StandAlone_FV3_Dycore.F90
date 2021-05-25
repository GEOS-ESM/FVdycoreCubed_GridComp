#define I_AM_MAIN

#include "MAPL_Generic.h"
#define _RC rc=status); _VERIFY(status

program StandAlone_FV3_Dycore
   use MAPL
   use FVdycoreCubed_GridComp,      only: SetServices
   implicit none

!EOP

!EOC

   character(*), parameter :: IAM = __FILE__

   type (MAPL_Cap) :: cap
   type (MAPL_FlapCLI) :: cli
   type (MAPL_CapOptions) :: cap_options
   integer :: status

   cli = MAPL_FlapCLI(description = 'FV Standalone Dycore',&
                              authors      =  'S.J. Lin, R. Rood, W. Putman')
   cap_options = MAPL_CapOptions(cli)
   cap = MAPL_Cap('GCM', SetServices, cap_options = cap_options)
   call cap%run(_RC)

 end Program StandAlone_FV3_Dycore

!EOC

