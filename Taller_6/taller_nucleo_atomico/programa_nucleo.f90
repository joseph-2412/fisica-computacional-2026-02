! Taller 6 -- Programa principal (escriba aquí su programa)
program programa_nucleo
  use modulo_nucleo, only: dp, termino_volumen, termino_superficie, &
                           termino_coulomb, termino_asimetria, &
                           energia_enlace, masa_atomica, Z_estable
  implicit none

  integer :: i, u_nuc, u_val, u_max
  integer :: A, Z_opt, A_max, Z_max
  real(dp) :: E_B_opt, E_B_global_min

  integer, dimension(8) :: Z_p = (/ 94, 32, 72, 18, 26, 20, 92, 92 /)
  integer, dimension(8) :: A_p = (/ 244, 70, 179, 38, 56, 40, 235, 238 /)

  ! 1. Escribir nucleos.dat
  open(newunit=u_nuc, file='nucleos.dat', status='replace')
  do i = 1, 8
     write(u_nuc, '(i4, 1x, i4, 1x, 5(f16.8, 1x), f16.4)') &
          Z_p(i), A_p(i), termino_volumen(A_p(i)), termino_superficie(A_p(i)), &
          termino_coulomb(Z_p(i), A_p(i)), termino_asimetria(Z_p(i), A_p(i)), &
          energia_enlace(Z_p(i), A_p(i)), masa_atomica(Z_p(i), A_p(i))
  end do
  close(u_nuc)

  ! 2. Escribir valle.dat y buscar el máximo enlace global (E_B mínimo)
  open(newunit=u_val, file='valle.dat', status='replace')

  ! Inicialización con A = 4
  call Z_estable(4, Z_opt, E_B_opt)
  write(u_val, '(i4, 1x, i4, 1x, f16.8)') 4, Z_opt, E_B_opt
  E_B_global_min = E_B_opt
  A_max = 4
  Z_max = Z_opt

  do A = 5, 260
     call Z_estable(A, Z_opt, E_B_opt)
     write(u_val, '(i4, 1x, i4, 1x, f16.8)') A, Z_opt, E_B_opt
     if (E_B_opt < E_B_global_min) then ! Búsqueda correcta del núcleo más ligado
        E_B_global_min = E_B_opt
        A_max = A
        Z_max = Z_opt
     end if
  end do
  close(u_val)

  ! 3. Escribir maximo.dat
  open(newunit=u_max, file='maximo.dat', status='replace')
  write(u_max, '(i4, 1x, i4, 1x, f16.8)') A_max, Z_max, E_B_global_min
  close(u_max)

end program programa_nucleo
