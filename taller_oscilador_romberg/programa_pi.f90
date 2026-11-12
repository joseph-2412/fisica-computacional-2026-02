! ============================================================
! Programa:  programa_pi
! Propósito: Calcular pi = int_0^1 4/(1+x^2) dx usando modulo_romberg sin modificar.
! ============================================================
program programa_pi
  use modulo_romberg, only: dp, romberg
  implicit none
  real(dp), parameter :: pi_exacto = 3.14159265358979323846_dp
  real(dp) :: a, b, tol, res_pi, error_abs
  integer :: kmax, k_final, n_eval
  real(dp), allocatable :: S(:,:)

  a = 0.0_dp
  b = 1.0_dp
  tol = 1.0e-5_dp
  kmax = 10

  call romberg(2, a, b, tol, kmax, S, k_final, n_eval)
  res_pi = S(k_final, k_final)
  error_abs = abs(res_pi - pi_exacto)

  open(unit=12, file='resultado_pi.dat', status='replace')
  write(12, '(i5,1x,i5,1x,es24.16,1x,es24.16)') &
        k_final, n_eval, res_pi, error_abs
  close(12)

  deallocate(S)
end program programa_pi
