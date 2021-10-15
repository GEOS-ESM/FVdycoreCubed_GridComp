module grid_bounds_mod

  implicit none

  private
  public :: GridBounds_T
  
  type GridBounds_T
     integer :: isd, ied
     integer :: jsd, jed
   contains
     procedure :: wr1te
  end type GridBounds_T

contains
  
  subroutine wr1te(self)

    ! Arguments
    class(GridBounds_T), intent(in) :: self

    ! Start
    print *, 'isd: ', self%isd
    print *, 'ied: ', self%ied
    print *, 'jsd: ', self%jsd
    print *, 'jed: ', self%jed

  end subroutine wr1te
  
end module grid_bounds_mod
