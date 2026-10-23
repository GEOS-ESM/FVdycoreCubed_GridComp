program test_fv3_splits
   use FV_StateMod, only: select_fv3_splits
   use fv_arrays_mod, only: FVPRC
   implicit none

   integer :: iteration, k_split, n_split, q_split, target_k_split

   ! A constant requirement between factors must not alternate between them.
   do iteration = 1, 20
      call check(5, 8, 0, 0.0_FVPRC, 24.0_FVPRC, 8, 5, 0, 6)
   enddo
   call check(5, 8, 0, 0.0_FVPRC, 20.0_FVPRC, 5, 8, 0, 5)
   call check(5, 8, 0, 0.0_FVPRC, 20.5_FVPRC, 8, 5, 0, 6)
   call check(5, 8, 0, 0.0_FVPRC, 32.0_FVPRC, 8, 5, 0, 8)
   call check(5, 8, 0, 0.0_FVPRC, 32.5_FVPRC, 10, 4, 0, 9)

   ! Increase immediately, including beyond the baseline work budget.
   call check(5, 8, 0, 0.0_FVPRC, 80.0_FVPRC, 20, 2, 0, 20)
   call check(5, 8, 0, 0.0_FVPRC, 160.0_FVPRC, 40, 1, 0, 40)
   call check(5, 8, 0, 0.0_FVPRC, 164.0_FVPRC, 41, 1, 0, 41)

   ! Relax directly back to baseline; no previous selection is an input.
   call check(5, 8, 0, 0.0_FVPRC, 0.0_FVPRC, 5, 8, 0, 5)
   call check(5, 8, 0, 30.0_FVPRC, 24.0_FVPRC, 8, 15, 0, 6)
   call check(5, 8, 1, 20.0_FVPRC, 24.0_FVPRC, 8, 10, 3, 6)
   call check(5, 8, 4, 0.0_FVPRC, 24.0_FVPRC, 8, 5, 4, 6)

   ! Preserve the core's automatic n_split and q_split sentinel modes.
   call check(5, 0, 0, 30.0_FVPRC, 24.0_FVPRC, 6, 0, 0, 6)
   call check(5, 0, 1, 30.0_FVPRC, 24.0_FVPRC, 6, 0, 5, 6)
   call check(1, 1, 0, 0.0_FVPRC, 0.0_FVPRC, 1, 1, 0, 1)

contains

   subroutine check(base_k, base_n, base_q, cfl_h, cfl_z, expected_k, expected_n, expected_q, expected_target)
      integer, intent(in) :: base_k, base_n, base_q, expected_k, expected_n, expected_q, expected_target
      real(FVPRC), intent(in) :: cfl_h, cfl_z

      call select_fv3_splits(base_k, base_n, base_q, cfl_h, cfl_z, &
                             k_split, n_split, q_split, target_k_split)
      if (any([k_split, n_split, q_split, target_k_split] /= &
              [expected_k, expected_n, expected_q, expected_target])) then
         print *, 'Actual k, n, q, target: ', k_split, n_split, q_split, target_k_split
         print *, 'Expected: ', expected_k, expected_n, expected_q, expected_target
         stop 1
      endif
      if (k_split < target_k_split) stop 2
      if (base_n > 0) then
         if (k_split * n_split < base_k * base_n) stop 3
         if (real(k_split * n_split, FVPRC) < 4.0_FVPRC * cfl_h) stop 4
      endif
   end subroutine check

end program test_fv3_splits
