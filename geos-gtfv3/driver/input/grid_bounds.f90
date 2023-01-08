module grid_bounds_mod

  implicit none

  private
  public :: GridBounds_T
  
  type GridBounds_T
     integer :: is, ie, js, je
     integer :: isd, ied, jsd, jed
   contains
     procedure :: wr1te
  end type GridBounds_T

contains
  
  subroutine wr1te(self)

    ! Arguments
    class(GridBounds_T), intent(in) :: self

    ! Start
    print *, 'is/ie/js/je: ', self%is, self%ie, self%js, self%je
    print *, 'isd/ied/jsd/jed: ', self%isd, self%ied, self%jsd, self%jed

  end subroutine wr1te
  
end module grid_bounds_mod
