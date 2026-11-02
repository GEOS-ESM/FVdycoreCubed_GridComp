#define VERIFY_(status) if (status /= 0) error stop 'MAPL pointer failure'
#define _ASSERT(condition, message) if (.not. (condition)) error stop message
#define RETURN_(status) rc = status; return

module split_test
  use, intrinsic :: iso_fortran_env, only: REAL8 => real64
  implicit none
  integer, parameter :: ESMF_SUCCESS=0, FVPRC=kind(1.0)
  logical :: ADJUST_SPLITS=.true., FV_HYDROSTATIC=.false., SW_DYNAMICS=.false.
  real :: setup_base_acoustic_dt=50.0
  integer :: pointer_calls=0
  type T_FVDYCORE_VARS
  end type
  type T_FVDYCORE_GRID
  end type
  type ESMF_State
     real(REAL8), pointer :: values(:) => null()
  end type
  type flags
     integer :: k_split=2, n_split=6, q_split=0
  end type
  type atmosphere
     type(flags) :: flagstruct
  end type
  type(atmosphere) :: FV_Atm(1)
#include "split_state.inc"
contains
  subroutine MAPL_GetPointer(internal, ptr, name, rc)
    type(ESMF_State), intent(inout) :: internal
    real(REAL8), pointer :: ptr(:)
    character(len=*), intent(in) :: name
    integer, intent(out) :: rc
    integer :: i
    pointer_calls=pointer_calls+1
    select case(name)
    case('K_SPLIT')
       i=1
    case('N_SPLIT')
       i=2
    case('Q_SPLIT')
       i=3
    case('BASE_ACOUSTIC_DT')
       i=4
    case default
       error stop 'Unknown field'
    end select
    ptr => internal%values(i:i)
    rc=0
  end subroutine

  subroutine WRITE_PARALLEL(message)
    character(len=*), intent(in) :: message
    print *, message
  end subroutine

  subroutine mp_reduce_max(value)
    real, intent(inout) :: value
  end subroutine

  subroutine mp_bcst(value)
    integer, intent(inout) :: value
  end subroutine

  integer function mpp_pe()
    mpp_pe=0
  end function

  subroutine step(state, index)
    type(T_FVDYCORE_STATE), intent(inout) :: state
    integer, intent(in) :: index
    real :: u(1,2,1), v(2,1,1), w(1,1,1), dz(1,1,1), dx(1,1), dy(1,1)
    real, parameter :: vertical(8)=[20.0,10.0,7.0,30.0,1.0,20.0,10.0,7.0]
    u=real(index)
    v=0.0
    w=vertical(index)/real(state%DT)
    dz=1.0
    dx=600.0
    dy=600.0
    call adjust_fv3_splits(state,1,1,1,1,1,1,1,1,1,u,v,w,dz,dx,dy)
    call FV_SaveSplits(state)
  end subroutine

#include "split_routines.inc"
end module

program test_restart
  use split_test
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  implicit none
  type(T_FVDYCORE_STATE) :: state, continuous, resumed
  type(ESMF_State) :: internal, full_restart, segment_restart, empty_restart
  real(REAL8) :: history(4,8), checkpoint(4)
  integer :: rc, index, mode
  character(len=16) :: argument

  call get_command_argument(1,argument)
  read(argument,*) mode
  allocate(internal%values(4),full_restart%values(4),segment_restart%values(4))
  state%DT=600.0
  internal%values=-1.0_REAL8

  if (mode /= 0) then
     select case(mode)
     case(1)
        internal%values(1)=0.0_REAL8
     case(2)
        internal%values(2)=-2.0_REAL8
     case(3)
        internal%values(3)=1.5_REAL8
     case(4)
        internal%values(4)=0.0_REAL8
     case(5)
        internal%values(1)=ieee_value(0.0_REAL8,ieee_quiet_nan)
     case(6)
        internal%values(1)=real(huge(0),REAL8)+1.0_REAL8
     case(7)
        internal%values(4)=ieee_value(0.0_REAL8,ieee_quiet_nan)
     end select
     call FV_InitSplits(state,internal,rc)
     stop 0
  endif

  call FV_InitSplits(state,internal,rc)
  if (rc /= 0 .or. any(internal%values /= [2.0_REAL8,6.0_REAL8,0.0_REAL8,50.0_REAL8])) &
       error stop 'Missing fields did not use setup'

  internal%values=[4.0_REAL8,-1.0_REAL8,0.0_REAL8,-1.0_REAL8]
  call FV_InitSplits(state,internal,rc)
  if (any(internal%values /= [4.0_REAL8,6.0_REAL8,0.0_REAL8,50.0_REAL8])) &
       error stop 'Partial restart did not use original baseline'

  internal%values=[3.0_REAL8,0.0_REAL8,0.0_REAL8,25.0_REAL8]
  call FV_InitSplits(state,internal,rc)
  if (state%k_split /= 3 .or. state%n_split /= 0 .or. state%q_split /= 0) &
       error stop 'Valid zeros were not restored'
  if (state%base_acoustic_dt /= 25.0) error stop 'Saved baseline overwritten'

  ADJUST_SPLITS=.false.
  pointer_calls=0
  checkpoint=internal%values
  call FV_InitSplits(state,empty_restart,rc)
  call FV_SaveSplits(state)
  if (pointer_calls /= 0) error stop 'Disabled adjustment accessed restart fields'
  if (any(internal%values /= checkpoint)) error stop 'Disabled adjustment modified restart fields'
  if (state%k_split /= 2 .or. state%n_split /= 6 .or. state%q_split /= 0) &
       error stop 'Disabled adjustment changed configured splits'
  ADJUST_SPLITS=.true.

  full_restart%values=-1.0_REAL8
  continuous%DT=600.0
  call FV_InitSplits(continuous,full_restart,rc)
  do index=1,8
     call step(continuous,index)
     history(:,index)=full_restart%values
  enddo
  if (all(history(1:3,4) == [2.0_REAL8,6.0_REAL8,0.0_REAL8])) &
       error stop 'Start/stop test did not change splits before checkpoint'

  segment_restart%values=-1.0_REAL8
  resumed%DT=600.0
  call FV_InitSplits(resumed,segment_restart,rc)
  do index=1,4
     call step(resumed,index)
  enddo
  checkpoint=segment_restart%values
  deallocate(segment_restart%values)
  allocate(segment_restart%values(4))
  segment_restart%values=checkpoint
  setup_base_acoustic_dt=100.0
  call FV_InitSplits(resumed,segment_restart,rc)
  if (any(segment_restart%values /= checkpoint)) error stop 'Initialization adjusted restored splits'
  do index=5,8
     call step(resumed,index)
     if (any(segment_restart%values /= history(:,index))) error stop 'Start/stop split sequence differs'
  enddo
end program
