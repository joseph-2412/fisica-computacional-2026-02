module integracion
  implicit none
contains

  subroutine trapecio(n, y, h, s)
    integer, intent(in)  :: n
    real(8), intent(in)  :: y(n), h
    real(8), intent(out) :: s

    if (n < 2) then
       s = 0.0d0
       return
    end if

    s = h * (0.5d0 * (y(1) + y(n)) + sum(y(2:n-1)))
  end subroutine trapecio

  subroutine simpson(n, y, h, s)
    integer, intent(in)  :: n
    real(8), intent(in)  :: y(n), h
    real(8), intent(out) :: s

    if (n < 3 .or. mod(n, 2) == 0) then
       s = 0.0d0
       return
    end if

    s = (h / 3.0d0) * (y(1) + y(n) + 4.0d0 * sum(y(2:n-1:2)) + 2.0d0 * sum(y(3:n-2:2)))
  end subroutine simpson

end module integracion
