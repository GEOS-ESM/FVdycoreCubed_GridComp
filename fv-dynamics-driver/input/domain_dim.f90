module domain_dim_mod

  implicit none

  private
  public :: DomainDim_T

  type DomainDim_T
     integer :: npx, npy, npz
   contains
     procedure :: wr1te
  end type DomainDim_T

contains
  
  subroutine wr1te(self)

    ! Arguments
    class(DomainDim_T), intent(in) :: self

    ! Start
    print *, 'npx: ', self%npx
    print *, 'npy: ', self%npy
    print *, 'npz: ', self%npz

  end subroutine wr1te

end module domain_dim_mod
