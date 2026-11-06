module modulo_nucleo
  implicit none

  ! Requisito R7: Doble precisión obligatoria
  integer, parameter :: dp = kind(1.0d0)

  ! Requisito R3: Constantes del modelo en MeV
  real(dp), parameter :: a_V = 15.8_dp
  real(dp), parameter :: a_S = 18.3_dp
  real(dp), parameter :: a_C = 0.71_dp
  real(dp), parameter :: a_A = 92.7_dp
  real(dp), parameter :: M_n = 939.57_dp
  real(dp), parameter :: M_H = 938.78_dp

contains

  ! Requisito R4: Funciones para los 4 términos
  real(dp) function termino_volumen(A)
    integer, intent(in) :: A
    if (A > 0) continue ! Evita warning por argumento sin usar
    termino_volumen = -a_V
  end function termino_volumen

  real(dp) function termino_superficie(A)
    integer, intent(in) :: A
    termino_superficie = a_S * (real(A, dp)**(-1.0_dp / 3.0_dp))
  end function termino_superficie

  real(dp) function termino_coulomb(Z, A)
    integer, intent(in) :: Z, A
    termino_coulomb = a_C * (real(Z, dp)**2) * (real(A, dp)**(-4.0_dp / 3.0_dp))
  end function termino_coulomb

  real(dp) function termino_asimetria(Z, A)
    integer, intent(in) :: Z, A
    real(dp) :: Z_r, A_r
    Z_r = real(Z, dp)
    A_r = real(A, dp)
    termino_asimetria = a_A * ((Z_r - 0.5_dp * A_r)**2) / (A_r**2)
  end function termino_asimetria

  ! Requisito R4: Energía de enlace por nucleón
  real(dp) function energia_enlace(Z, A)
    integer, intent(in) :: Z, A
    energia_enlace = termino_volumen(A) + termino_superficie(A) + &
                     termino_coulomb(Z, A) + termino_asimetria(Z, A)
  end function energia_enlace

  ! Requisito R4: Masa atómica
  real(dp) function masa_atomica(Z, A)
    integer, intent(in) :: Z, A
    real(dp) :: Z_r, A_r
    Z_r = real(Z, dp)
    A_r = real(A, dp)
    masa_atomica = Z_r * M_H + (A_r - Z_r) * M_n + A_r * energia_enlace(Z, A)
  end function masa_atomica

  ! Requisito R5: Subrutina para buscar Z de máximo enlace (E_B más negativo)
  subroutine Z_estable(A, Z_min, E_B_min)
    integer, intent(in)   :: A
    integer, intent(out)  :: Z_min
    real(dp), intent(out) :: E_B_min
    integer  :: Z_c
    real(dp) :: E_actual

    ! Inicialización adecuada con el primer valor posible Z = 1
    Z_min = 1
    E_B_min = energia_enlace(1, A)

    do Z_c = 2, A
       E_actual = energia_enlace(Z_c, A)
       if (E_actual < E_B_min) then ! Busca el mínimo algebraico (más negativo)
          E_B_min = E_actual
          Z_min = Z_c
       end if
    end do
  end subroutine Z_estable

end module modulo_nucleo
