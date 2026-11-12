! ============================================================
! Programa:  programa_energia
! Propósito: Calcular la energía disipada E_dis en [0, 3] s mediante Romberg
!            y compararla con E(0) - E(T).
! ============================================================
program programa_energia
  use modulo_romberg, only: dp, romberg, gamma_param, omega_param
  implicit none
  real(dp) :: a, b, tol, e_dis, e_0, e_t, e_ref, dif, omega02, v0, vT, x0, xT
  integer :: kmax, k_final, n_eval, i, j
  real(dp), allocatable :: S(:,:)

  a = 0.0_dp
  b = 3.0_dp
  tol = 1.0e-6_dp
  kmax = 12

  call romberg(1, a, b, tol, kmax, S, k_final, n_eval)
  e_dis = S(k_final, k_final)

  omega02 = omega_param**2 + gamma_param**2
  x0 = exp(-gamma_param * a) * cos(omega_param * a)
  v0 = -exp(-gamma_param * a) * (gamma_param * cos(omega_param * a) + omega_param * sin(omega_param * a))
  e_0 = 0.5_dp * v0**2 + 0.5_dp * omega02 * x0**2

  xT = exp(-gamma_param * b) * cos(omega_param * b)
  vT = -exp(-gamma_param * b) * (gamma_param * cos(omega_param * b) + omega_param * sin(omega_param * b))
  e_t = 0.5_dp * vT**2 + 0.5_dp * omega02 * xT**2

  e_ref = e_0 - e_t
  dif = abs(e_dis - e_ref)

  ! Exportar tabla_energia.dat (filas k=0..k_final)
  open(unit=10, file='tabla_energia.dat', status='replace')
  do i = 0, k_final
     write(10, '(13(es24.16,1x))') (S(i, j), j=0, i)
  end do
  close(10)

  ! Exportar resultado_energia.dat
  open(unit=11, file='resultado_energia.dat', status='replace')
  write(11, '(i5,1x,i5,1x,es24.16,1x,es24.16,1x,es24.16)') &
        k_final, n_eval, e_dis, e_ref, dif
  close(11)

  deallocate(S)
end program programa_energia
