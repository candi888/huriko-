module subprogs
    implicit none

contains
    subroutine setParameters(n, t_max_cnt, delta_t)
        integer, intent(out) ::n, t_max_cnt
        real(8), intent(out)::  delta_t

        delta_t = 1d-3
        write (*, *) "input n:"
        read (*, *) n
        t_max_cnt = 30000

    end subroutine setParameters

    subroutine setInitL(l, n)
        integer, intent(in) :: n
        real(8), intent(out) :: l(n)

        l(:) = 0.1d0

    end subroutine setInitL

    subroutine setInitTheta(theta, n)
        integer, intent(in) :: n
        real(8), intent(out) :: theta(n)
        real(8), parameter :: pi = acos(-1d0)

        theta(:) = pi/2 - 0.1

    end subroutine setInitTheta

    subroutine setInitOmega(omega, n)
        integer, intent(in) :: n
        real(8), intent(out) :: omega(n)

        omega(:) = 0d0

    end subroutine setInitOmega

    subroutine setInitM(m, n)
        integer, intent(in) :: n
        real(8), intent(out) :: m(n)

        m(:) = 1.0d0

    end subroutine setInitM

    subroutine setInitAccM(accm, m, n)
        integer, intent(in) :: n
        real(8), intent(in) :: m(n)
        real(8), intent(out) :: accm(n + 1)
        integer :: i

        accm(1) = 0d0
        do i = 1, n
            accm(i + 1) = accm(i) + m(i)
        end do

    end subroutine setInitAccM

    function CalcOmegadot(cur_omega, cur_theta, accm, l, delta_t, n, g) result(nxt_omega)
        integer, intent(in) ::n
        real(8), intent(in) :: delta_t, g, cur_omega(n), cur_theta(n), accm(n), l(n)
        real(8):: nxt_omega(n)
        real(8) :: A(n, n), B(n, n), tmp(n, n)
        integer :: i, j
        B(:, :) = 0d0

        do i = 1, n
            B(i, i) = -(accm(n + 1) - accm(i))*g*sin(cur_theta(i))
            do j = 1, n
                A(i, j) = (accm(n + 1) - accm(max(i, j)))*l(j)*cos(cur_theta(i) - cur_theta(j))
                B(i, j) = B(i, j) - (accm(n + 1) - accm(max(i, j)))* &
                &(l(j)*(cur_omega(j)**2)*sin(cur_theta(i) - cur_theta(j)))
            end do
        end do
        do i = 1, n
            tmp = matmul(calc_inv_A(A, n), B)
            nxt_omega(i) = sum(tmp(i, :))
        end do

    end function CalcOmegadot

    subroutine CalcByRungeKutta(omega, theta, accm, l, delta_t, g, n)
        integer, intent(in) :: n
        real(8), intent(in) :: accm(n), l(n), delta_t, g
        real(8) :: omega(n), theta(n)
        real(8) :: p1(n), p2(n), p3(n), p4(n), y1(n), y2(n), y3(n), y4(n)

        ! pがomega, yがtheta
        p1(:) = delta_t*CalcOmegadot(omega, theta, accm, l, delta_t, n, g)
        y1(:) = delta_t*omega(:)
        p2(:) = delta_t*CalcOmegadot(omega(:) + p1(:)/2, theta(:) + y1(:)/2, accm, l, delta_t, n, g)
        y2(:) = delta_t*(omega(:) + p1(:)/2)
        p3(:) = delta_t*CalcOmegadot(omega(:) + p2(:)/2, theta(:) + y2(:)/2, accm, l, delta_t, n, g)
        y3(:) = delta_t*(omega(:) + p2(:)/2)
        p4(:) = delta_t*CalcOmegadot(omega(:) + p3(:), theta(:) + y3(:), accm, l, delta_t, n, g)
        y4(:) = delta_t*(omega(:) + p3(:))

        theta(:) = theta(:) + (y1(:) + 2*y2(:) + 2*y3(:) + y4(:))/6
        omega(:) = omega(:) + (p1(:) + 2*p2(:) + 2*p3(:) + p4(:))/6
    end subroutine CalcByRungeKutta

    function calc_inv_A(a0, n) result(inv_a)
        ! 部分ピボット選択なし
        integer, intent(in) :: n
        real(8), intent(in) :: a0(n, n)
        real(8) :: inv_a(n, n)
        integer i, k
        real(8) ar, a(n, 2*n)

        a(:, 1:n) = a0(:, :)
        a(:, n + 1:2*n) = 0.0d0
        do i = 1, n
            a(i, i + n) = 1.0d0
        end do

        do k = 1, n
            if (a(k, k) == 0.0d0) stop "pivot =0"
            ar = 1.0d0/a(k, k)
            a(k, k) = 1.0d0
            a(k, k + 1:2*n) = ar*a(k, k + 1:2*n)
            do i = 1, n
                if (i == k) cycle
                a(i, k + 1:2*n) = a(i, k + 1:2*n) - a(i, k)*a(k, k + 1:2*n)
                a(i, k) = 0.0d0
            end do
        end do

        inv_a(:, :) = a(:, n + 1:2*n)
    end function calc_inv_A
end module subprogs

! 単位[m,s,kg]
program main
    use subprogs
    implicit none

    real(8), allocatable :: theta(:), omega(:), l(:), m(:), accm(:)
    real(8), allocatable :: A(:, :), B(:, :), x(:), tmp(:, :)
    real(8), parameter ::  g = 9.8d0
    integer ::n, t_max_cnt
    real(8) ::  delta_t
    integer :: i, j, it

    call setParameters(n, t_max_cnt, delta_t)
    allocate (theta(n), omega(n), l(n), m(n), accm(n + 1), tmp(n, n))
    allocate (A(n, n), x(n), B(n, n))

    call setInitL(l, n)
    call setInitM(m, n)
    call setInitAccM(accm, m, n)
    call setInitTheta(theta, n)
    call setInitOmega(omega, n)

    open (11, file="./data.dat", status="replace")
    write (11, *) "各ひもの長さ"
    write (11, *) l(:)
    write (11, *) "時間刻み"
    write (11, *) delta_t
    write (11, *) 0, theta(:)

    do it = 1, t_max_cnt

        ! B(:, :) = 0d0

        ! do i = 1, n
        !     B(i, i) = -(accm(n + 1) - accm(i))*g*sin(theta(i))
        !     do j = 1, n
        !         A(i, j) = (accm(n + 1) - accm(max(i, j)))*l(j)*cos(theta(i) - theta(j))
        !         B(i, j) = B(i, j) - (accm(n + 1) - accm(max(i, j)))*(l(j)*(omega(j)**2)*sin(theta(i) - theta(j)))
        !     end do
        ! end do
        ! do i = 1, n
        !     tmp = matmul(calc_inv_A(A, n), B)
        !     x(i) = sum(tmp(i, :))
        ! end do
        ! omega(:) = omega(:) + delta_t*x(:)
        ! theta(:) = theta(:) + delta_t*omega(:)

        call CalcByRungeKutta(omega, theta, accm, l, delta_t, g, n)

        write (11, *) it*delta_t, theta(:)
        if (isnan(theta(n))) stop "NAN"
    end do

    close (11)
    deallocate (theta, omega, l, m)
end program main
