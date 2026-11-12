module modulo_romberg
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  real(dp), parameter :: gamma_param = 0.5_dp
  real(dp), parameter :: omega_param = 2.0_dp * 3.14159265358979323846_dp

contains

  real(dp) function integrando_energia(t)
    real(dp), intent(in) :: t
    real(dp) :: v
    v = -exp(-gamma_param * t) * (gamma_param * cos(omega_param * t) + omega_param * sin(omega_param * t))
    integrando_energia = 2.0_dp * gamma_param * (v**2)
  end function integrando_energia

  real(dp) function integrando_pi(x)
    real(dp), intent(in) :: x
    integrando_pi = 4.0_dp / (1.0_dp + x**2)
  end function integrando_pi

  real(dp) function f_eval(opcion, x)
    integer, intent(in) :: opcion
    real(dp), intent(in) :: x
    if (opcion == 1) then
       f_eval = integrando_energia(x)
    else
       f_eval = integrando_pi(x)
    end if
  end function f_eval

  subroutine romberg(opcion, a, b, tol, kmax, S, k_final, n_eval)
    integer, intent(in) :: opcion, kmax
    real(dp), intent(in) :: a, b, tol
    real(dp), allocatable, intent(out) :: S(:,:)
    integer, intent(out) :: k_final, n_eval
    
    integer :: k, j, i, puntos_nuevos
    real(dp) :: h, suma, factor

    allocate(S(0:kmax, 0:kmax))
    S = 0.0_dp

    h = b - a
    S(0,0) = 0.5_dp * h * (f_eval(opcion, a) + f_eval(opcion, b))
    n_eval = 2
    k_final = 0

    do k = 1, kmax
       h = h / 2.0_dp
       puntos_nuevos = 2**(k - 1)
       suma = 0.0_dp
       do i = 1, puntos_nuevos
          suma = suma + f_eval(opcion, a + real(2*i - 1, dp) * h)
       end do
       n_eval = n_eval + puntos_nuevos
       S(k,0) = 0.5_dp * S(k-1,0) + h * suma

       do j = 1, k
          factor = 4.0_dp**j
          S(k,j) = S(k,j-1) + (S(k,j-1) - S(k-1,j-1)) / (factor - 1.0_dp)
       end do

       if (abs(S(k,k) - S(k-1,k-1)) < tol) then
          k_final = k
          exit
       end if
    end do

    if (k_final == 0) k_final = kmax
  end subroutine romberg

end module modulo_romberg
