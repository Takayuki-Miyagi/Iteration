program test
  use Iteration, only: Optimizer
  integer :: i, n =1000, ndim
  real(8), allocatable :: fk(:), xk(:)
  type(Optimizer) :: opt

  ndim = 1
  allocate(fk(ndim), xk(ndim))

  !call opt%init(n=ndim, method="broyden", a = 0.9d0, m = 10)
  !call opt%init(n=ndim, method="bfgs", a = 0.1d0, m = 10)
  !call opt%init(n=ndim, method="mbroyden", a = 0.7d0, m = 10)
  call opt%init(n=ndim, method="l-bfgs", a = 0.7d0, m = 10)
  call opt%set_init([1.d0,1.d0])
  do i = 1, n
    call func(xk, fk)
    call opt%GetNextInput(fk, i)
    xk = opt%xj%v
    call func(xk, fk)
    write(*,'(i4, 4f12.5)') i, xk(:), fk(:) - xk(:)
    if(opt%r < 1.d-8) exit
  end do
  call func(xk, fk)
  write(*,'(2f12.5)') fk(:)
  call opt%Fin()

contains
  subroutine func(x, f)
    real(8), intent(in) :: x(:)
    real(8), intent(out) :: f(:)

    f(:) = cos(x(:))
  end subroutine func

end program test
