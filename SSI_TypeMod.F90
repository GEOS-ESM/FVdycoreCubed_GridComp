module SSI_TypeMod
! contains info to map local SSI arrays to coarse decomposition PETs
! Relies of Tom C.'s compact communicator

implicit none
private

type SSI_Type
   integer :: nnx  ! node_topology(1) from Tom's compact communicator
   integer :: nny  ! node_topology(2) from Tom's compact communicator
   integer :: nth_x ! number of threads in x-direction
   integer :: nth_y ! number of threads in y-direction
   integer :: pet_id_x ! node-local pet id in x-direction in a 2-D mapping of PETs on a node 
   integer :: pet_id_y ! node-local pet id in y-direction in a 2-D mapping of PETs on a node 
   integer :: npet_x  ! node-local num pets in x-direction
   integer :: npet_y  ! node-local num pets in y-direction
   integer :: is  
   integer :: js 
end type SSI_Type

public SSI_Type
end module SSI_TypeMod
