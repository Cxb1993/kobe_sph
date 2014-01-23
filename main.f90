program sph
  use particle
  implicit none
  integer :: itr = 0

  call init

  call out_result(itr)

  do itr = 1, nloop
    write(*,*) 'itr = ', itr
    call set_rho_and_p()

    call set_force()

    call move_particle()

    call out_result(itr)
  enddo

  write(*,*) 'finished!'

contains
  subroutine init
    implicit none
    integer i, j
    integer root
    real(8) inv
    real(8) length
    real(8) rn

    !root = sqrt(dble(n_particle))
    root = dble(n_particle)**(1.0d0/3.0d0)

    inv = 1.0d0 / root
    length = area / root
    write(*,*) "root = ", root, " length = ",length

    !do j = 1, n_particle
    !  do i = 1, ndim
    !    if(i == 1) then
    !      !x(i,j) = 0.0d0
    !      x(i,j) = length * mod(j-1,root)
    !      v(i,j) = 0.0d0
    !    elseif(i == 2) then
    !      !x(i,j) = length * mod(j-1,root)
    !      !v(i,j) = inv * mod(j-1,root)
    !      !v(i,j) = 1000 * inv * mod(j-1,root)
    !      x(i,j) = length * mod(int((j-1) * inv),root)
    !      v(i,j) = 0.0d0
    !    else
    !      x(i,j) = length * int((j-1) * inv * inv)
    !      v(i,j) = 0.0d0
    !    endif
    !  enddo
    !enddo

    do j = 1, n_particle
      do i = 1, ndim
        call random_number(rn)
        if(i == 1) then
          x(i,j) = 0.5d0 * Lx * rn 
          v(i,j) = 0.0d0
        elseif(i == 2) then
          x(i,j) = Ly * rn 
          v(i,j) = 0.0d0
        else
          x(i,j) = 2.0d0 / 3.0d0 * Lz * rn 
          v(i,j) = 0.0d0
        endif
      enddo
    enddo

    x_min(1) = 0.0d0
    x_min(2) = 0.0d0
    x_min(3) = 0.0d0
    !x_max(1) = 0.0001d0
    !x_max(1) = 0.1d0
    x_max(1) = Lx
    !x_max(1) = 2.0d0*area
    !x_max(2) = 0.1d0
    x_max(2) = Ly
    !x_max(2) = 2.0d0*area
    !x_max(2) = 0.01d0
    !x_max(3) = 0.1d0
    !x_max(3) = 0.2d0
    x_max(3) = Lz
    !x_max(3) = 2.0d0*area
    !x_max(3) = 0.01d0

  end subroutine init

  subroutine set_rho_and_p()
    integer i,j

    do i = 1, n_particle
      rho(i) = 0.0d0
      do j = 1, n_particle
        if(i == j) cycle
        rho(i) = rho(i) + mass * weight_poly6(x(:,i),x(:,j),h) 
        !write(*,*) 'i = ', i, 'mass = ', mass, ' rho(i) = ', rho(i)
        p(i) = k * (rho(i) - rho0)
      enddo
      !write(*,*) 'i = ', i, ' rho(i) = ', rho(i)
    enddo
  end subroutine set_rho_and_p

  subroutine set_force()
    integer i,j,k
    real(8) weight_vis_, weight_press_(ndim)
    real(8) f_visc, f_press

    do i = 1, n_particle
      f(:,i) = 0.0d0
      do j = 1, n_particle
        if(i == j) cycle
        weight_vis_ = weight_vis(x(:,i),x(:,j),h)
        weight_press_ = weight_press(x(:,i),x(:,j),h)
        !write(*,*) 'weight_vis_ = ', weight_vis_, ' weight_press_ = ', weight_press_
        do k = 1, ndim
          f_visc = mass * (v(k,j) - v(k,i)) / rho(j) * weight_vis_ 
          f_press = mass * (p(i) + p(j)) / (2 * rho(j)) * weight_press_(k)
          f(k,i) = f(k,i) + mu * f_visc - f_press
          !if(isnan(f(k,i))) then
          !  f(k,i) = 0.0d0
          !  !write(*,*) '-eps < rho(j) = ', -eps < rho(j), ' rho(j) < eps = ', rho(j) < eps
          !  !write(*,*) 'mass = ', mass, ' v(k,j) = ', v(k,j), ' v(k,i) = ', v(k,i)
          !  !write(*,*) 'p(i) = ', p(i), ' p(j) = ', p(j), ' rho(j) = ', rho(j)
          !  !write(*,*) 'f_visc = ', f_visc, ' f_press = ', f_press, ' f = ', f(k,i)
          !  !stop
          !endif
        enddo
      enddo
      !write(*,*) 'i = ', i, ' rho(i) = ', rho(i)
    enddo
  end subroutine set_force

  function weight_poly6(x1,x2,h) result(weight)
    real(8), intent(in) :: x1(:)
    real(8), intent(in) :: x2(:)
    real(8), intent(in) :: h 
    real(8) r(ndim)
    real(8) :: weight
    real(8) alpha
    real(8) norm,norm2
    real(8), parameter :: pi = acos(-1.0d0)

    weight = 0.0d0

    r(1) = x2(1) - x1(1)
    r(2) = x2(2) - x1(2)
    r(3) = x2(3) - x1(3)

    norm2 = r(1)*r(1) + r(2)*r(2) + r(3)*r(3)
    norm = sqrt(norm2)

    if(norm <= h) then
      alpha = (315d0/(64d0*pi*(h**9)))
      weight = alpha * (h**2 - norm2)**3 
    endif

    !write(*,*) 'x1 = ', x1
    !write(*,*) 'x2 = ', x2
    !write(*,*) 'r = ', r
    !write(*,*) 'norm2 = ', norm2
    !write(*,*) 'norm = ', norm
    !write(*,*) 'pi = ', pi
    !write(*,*) 'h**9 = ', h**9
    !write(*,*) 'alpha = ', alpha
    !write(*,*) 'weight = ', weight
  end function weight_poly6   

  function weight_vis(x1,x2,h) result(weight)
    real(8), intent(in) :: x1(:)
    real(8), intent(in) :: x2(:)
    real(8), intent(in) :: h 
    real(8) r(ndim)
    real(8) :: weight
    real(8) alpha
    real(8) norm,norm2
    real(8), parameter :: pi = acos(-1.0d0)

    weight = 0.0d0

    r(1) = x2(1) - x1(1)
    r(2) = x2(2) - x1(2)
    r(3) = x2(3) - x1(3)

    norm2 = r(1)*r(1) + r(2)*r(2) + r(3)*r(3)
    norm = sqrt(norm2)

    if(norm <= h) then
      alpha = (45d0/(pi*(h**6)))
      weight = alpha * (h - norm) 
    endif

    !write(*,*) '*** weight_vis ***'
    !write(*,*) 'x1 = ', x1
    !write(*,*) 'x2 = ', x2
    !write(*,*) 'r = ', r
    !write(*,*) 'norm2 = ', norm2
    !write(*,*) 'norm = ', norm
    !write(*,*) 'pi = ', pi
    !write(*,*) 'h**6 = ', h**6
    !write(*,*) 'alpha = ', alpha
    !write(*,*) 'weight = ', weight
  end function weight_vis    

  function weight_press(x1,x2,h) result(weight)
    real(8), intent(in) :: x1(:)
    real(8), intent(in) :: x2(:)
    real(8), intent(in) :: h 
    real(8) r(ndim)
    real(8) :: weight(ndim)
    real(8) alpha
    real(8) norm,norm2
    real(8), parameter :: pi = acos(-1.0d0)
    integer k

    weight(:) = 0.0d0

    r(1) = x2(1) - x1(1)
    r(2) = x2(2) - x1(2)
    r(3) = x2(3) - x1(3)

    norm2 = r(1)*r(1) + r(2)*r(2) + r(3)*r(3)
    norm = sqrt(norm2)

    if(norm <= h) then
      alpha = -(45d0/(pi*(h**6)))
      do k = 1, ndim
        weight(k) = alpha * (h - norm)**2 * r(k) / norm
      enddo
    endif

    !write(*,*) '*** weight_press ***'
    !write(*,*) 'x1 = ', x1
    !write(*,*) 'x2 = ', x2
    !write(*,*) 'r = ', r
    !write(*,*) 'norm2 = ', norm2
    !write(*,*) 'norm = ', norm
    !write(*,*) 'pi = ', pi
    !write(*,*) 'h**6 = ', h**6
    !write(*,*) 'alpha = ', alpha
    !write(*,*) 'weight = ', weight
  end function weight_press    

  subroutine out_result(itr)
    integer, intent(in) :: itr
    integer i, j
    character filename*128

    write(*, '(f10.3,a,f10.3,a,f10.3,a,f10.3,a,f10.3,a,f10.3,a)') x(1,1), ',', x(2,1), ',', x(3,1), ',', v(1,1), ',', v(2,1), ',', v(3,1), ','

    !do i = 1, n_particle
    !    write(*, '(a,i4,a,3e15.3)') 'x(', i, ") = ", x(:,i)
    !enddo
    !do i = 1, n_particle
    !    write(*, '(a,i4,a,3e15.3)') 'v(', i, ") = ", v(:,i)
    !enddo
    !do i = 1, n_particle
    !    write(*, '(a,i,a,e10.3,a,e15.3,a,3e15.3)') 'i = ', i,' rho = ', rho(i), ' p = ', p(i), ' f = ', f(:,i)
    !enddo

    write(filename,'("data",i5.5)') itr
    open(11, file=filename, status='replace')
    do i = 1, n_particle
      write(11, '(f15.3,a,f15.3,a,f15.3,a)') x(1,i), ',', x(2,i), ',', x(3,i), ',' 
      !write(11, '(f15.3,a,f15.3,a,f15.3,a)') 5.0d0*x(1,i), ',', 5.0d0*x(2,i), ',', 5.0d0*x(3,i), ',' 
      !write(11, '(f15.3,a,f15.3,a,f15.3,a)') 10.0d0*x(1,i), ',', 10.0d0*x(2,i), ',', 10.0d0*x(3,i), ',' 
      !write(11, '(f15.3,a,f15.3,a,f15.3,a)') 50.0d0*x(1,i), ',', 50.0d0*x(2,i), ',', 50.0d0*x(3,i), ',' 
      !write(11, '(f15.3,a,f15.3,a,f15.3,a)') 500.0d0*x(1,i), ',', 500.0d0*x(2,i), ',', 500.0d0*x(3,i), ',' 
    enddo
    close(11)
  end subroutine out_result

  subroutine move_particle()
    integer i,j
    real(8) accel
    real(8) diff, adj

    do i = 1, n_particle
      do j = 1, ndim
        !write(*,*) 'prep0 rho(j) = ', rho(j), ' f(j,i) = ', f(j,i), ' f(j,i) / rho(i) = ', f(j,i) / rho(i)
        accel = f(j,i) / rho(i)

        !if(isnan(accel)) then
        !  accel = 0.0d0
        !endif

        !write(*,*) 'prep1 --- v = ', v(:,i), ' accel = ', accel
        ! boundary condition
        diff = radius - (x(j,i) - x_min(j)) 
        if(diff > eps) then
          !if(v(j,i) < 0.0d0) then
            if(i == 1) then
              if(j == 2) then
                write(*,'(a,f10.3,a,f10.3,a,f10.3)') 'radius = ', radius, ' x(j,i) = ', x(j,i), ' x_min(j) = ', x_min(j)
                write(*,'(a,f10.3)') 'diff = ', diff
              endif
            endif
            adj = ext_stiff * diff - ext_damp * v(j,i)
            if(i == 1) then
              if(j == 2) then
                write(*,*) 'accel_pre = ', accel, 'v(j,i) = ', v(j,i), ' adj = ', adj
              endif
            endif
            accel = accel + adj
            if(i == 1) then
              if(j == 2) then
                write(*,*) 'accel_mod = ', accel
              endif
            endif
          !endif
        endif

        !write(*,*) 'prep2 --- v = ', v(:,i), ' accel = ', accel
        diff = radius - (x_max(j) - x(j,i)) 
        if(diff > eps) then
          !if(v(j,i) > 0.0d0) then
            adj = ext_stiff * diff + ext_damp * v(j,i)
            accel = accel - adj
          !endif
        endif

        !write(*,*) 'prep3 --- v = ', v(:,i), ' accel = ', accel
        if(j == 3) then
          accel = accel - 9.8d0
        endif
        v(j,i) = v(j,i) + accel * dt
        x(j,i) = x(j,i) + v(j,i) * dt
      enddo
      !write(*,*) 'i = ', i, ' f(i) = ', f(:,i), ' rho = ', rho(i), ' dt = ', dt
      !write(*,*) 'v = ', v(:,i), ' x = ', x(:,i)
    enddo
  end subroutine move_particle
end program sph
