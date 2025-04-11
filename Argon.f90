module liquid
  implicit none
  integer :: n = 864, bc = 1, tr = 0 !bc: periodic boundary, tr: poor-mans algorithm, 1: on
  real, allocatable :: dens(:,:,:)
  !crucial parameters
  real, parameter :: T = 94.4, kb= 1.380649*10.**(-16.)
  real, parameter :: m = 39.95*1.67*10.**(-24), sigma = 3.4*10**(-8.)
  real, parameter :: ep = 120.*kb, rc = 2.25, L = 10.229, halfL = 10.229/2.
  !time parameters
  integer :: nt, ntrelax
  real, parameter :: tmin = 0.0, tmax = 100.0, trelax = 5.0 !in ps
  real ::  dt = 0.01
contains
 subroutine Initial()
    implicit none
    integer :: i, j, k, p
    real :: spacing, v_stddev, rand1, rand2, Tp, E
    allocate(dens(n,3,2))
    ntrelax = int(trelax/dt); nt = int((tmax - tmin)/dt)
    dt = dt*1.0e-12*(ep/m)**0.5*1./sigma
    ! Initialize densitions randomly within the box
    spacing = L / (int(n**(1.0/3.0)) + 1) ! Adjust spacing to avoid boundary
    p = 0
    do i = 0, int(n**(1.0/3.0))
    do j = 0, int(n**(1.0/3.0))
    do k = 0, int(n**(1.0/3.0))
      p = p + 1
      if (p > n) exit
      dens(p,1,1) = i * spacing
      dens(p,2,1) = j * spacing
      dens(p,3,1) = k * spacing
    end do
    if (p > n) exit
    end do
    if (p > n) exit
    end do
    v_stddev = sqrt(kb * T / m)
    do i = 1, n
      do j = 1, 3
      call random_number(rand1); call random_number(rand2)
      dens(i,j,2) = v_stddev * sqrt(-2.0 * log(rand1)) &
                    * cos(2.0 * 3.14159265359 * rand2)
      end do
    end do
    do i = 1, 3
     dens(:,i,2) = dens(:,i,2) - sum(dens(:,i,2))/n
    end do
    Tp = sum(dens(:,:,2)**2.) * m / (3.0 * n * kb)
    dens(:,:,2) = dens(:,:,2) * sqrt(T / Tp)
    dens(:,:,2) = dens(:,:,2) / (ep/m)**0.5 ! Convert to dimensionless units
  end subroutine

 subroutine calaccf(acc, dens)
    implicit none
  integer :: i, j, k
  real :: rr, rri, rri2, rri6, dr(3), coef(3)
    real, intent(out) :: acc(n,3)
  real, intent(in) :: dens(n,3)
    acc = 0.0
    do i = 1, n
      do j = i + 1, n
      dr(:) = dens(i,:) - dens(j,:)
if (dr(1) >  halfL)  dr(1) = dr(1) - L
if (dr(1) < -halfL)  dr(1) = dr(1) + L
if (dr(2) >  halfL)  dr(2) = dr(2) - L
if (dr(2) < -halfL)  dr(2) = dr(2) + L
if (dr(3) >  halfL)  dr(3) = dr(3) - L
if (dr(3) < -halfL)  dr(3) = dr(3) + L
      rr = sum(dr**2)
      if (rr < rc**2.) then
        rri2 = 1.0 / rr; rri6 = rri2 * rri2 * rri2
        coef = rri2 * (2.0 * rri6 * rri6 - rri6)
        acc(i,:) = acc(i,:) + coef*dr(:)
        acc(j,:) = acc(j,:) - coef*dr(:)
      end if
      end do
    end do
  acc = acc*24.0
  end subroutine


  subroutine verlexp()
    implicit none
    integer :: i, j
    real :: acc(n,3), temp(n,3,2), Tp, vhalf(n,3)
    call calaccf(acc, dens(:,:,1))
    !print*, maxval(acc), minval(acc)
    !$OMP PARALLEL DO PRIVATE(i,j) SCHEDULE(dynamic)
    do i = 1, n
      do j = 1, 3
        vhalf(i,j) = dens(i,j,2) + 0.5 * dt * acc(i,j)
        dens(i,j,1) = dens(i,j,1) + dt * vhalf(i,j)
      end do
    end do
    !$OMP END PARALLEL DO
    call calaccf(acc, dens(:,:,1))
    !$OMP PARALLEL DO PRIVATE(i,j) SCHEDULE(dynamic)
    do i = 1, n
      do j = 1, 3
        dens(i,j,2) = vhalf(i,j) + 0.5 * dt * acc(i,j)
      end do
    end do
    !$OMP END PARALLEL DO
    ! Apply periodic boundary conditions
    do i = 1, n
      do j = 1, 3
        dens(i,j,1) = mod(dens(i,j,1), L)
        if (dens(i,j,1) < 0.) dens(i,j,1) = dens(i,j,1) + L
      end do
    end do
    Tp = sum(dens(:,:,2)**2.) * m / (3. * kb * n) * (ep / m) ! Temperature in Kelvin
    if ((abs(Tp - T) > 5) .AND. (tr == 1)) then
      !$OMP PARALLEL DO PRIVATE(i,j) SCHEDULE(dynamic)
      do i = 1, n
        do j = 1, 3
          dens(i,j,2) = dens(i,j,2) * sqrt(T / Tp)
        end do
      end do
      !$OMP END PARALLEL DO
    end if
  end subroutine


  subroutine writetraject(ns)
    implicit none
    integer, intent(in) :: ns
    integer :: i
    if (ns == 0) then
      open(unit = 2, file = "Data/Output.lammpstrj", status = "replace", action = "write")
    else
      open(unit = 2, file = "Data/Output.lammpstrj", action = "write", position = "append")
    end if
    write(2,*) "ITEM: TIMESTEP"
    write(2,*) ns
    write(2,*) "ITEM: NUMBER OF ATOMS"
    write(2,*) n
    write(2,*) "ITEM: BOX BOUNDS pp pp pp"
    write(2,*) 0.0, L
    write(2,*) 0.0, L
    write(2,*) 0.0, L
    write(2,*) "ITEM: ATOMS id x y z vx vy vz"
    do i = 1, n
    write(2,*) i, dens(i,1,1), dens(i,2,1), dens(i,3,1), &
                dens(i,1,2)*(m/ep)**0.5, dens(i,2,2)*(m/ep)**0.5, dens(i,3,2)*(m/ep)**0.5
  end do
  close(2)
  end subroutine
  
  real function Potent(pos)
    implicit none
    integer :: i, j
    real :: rr, rri, rri2, rri6, dr(3), pot
    real, intent(in) :: pos(n,3)
    pot = 0.0
    !$OMP PARALLEL DO PRIVATE(i, j, dr, rr, rri2, rri6) REDUCTION(+:pot) SCHEDULE(dynamic)
    do i = 1, n
      do j = i + 1, n
      dr(:) = pos(i,:) - pos(j,:)
      if (dr(1) >  halfL)  dr(1) = dr(1) - L
      if (dr(1) < -halfL)  dr(1) = dr(1) + L
      if (dr(2) >  halfL)  dr(2) = dr(2) - L
      if (dr(2) < -halfL)  dr(2) = dr(2) + L
      if (dr(3) >  halfL)  dr(3) = dr(3) - L
      if (dr(3) < -halfL)  dr(3) = dr(3) + L
      rr = sum(dr**2)
      if (rr < rc**2.) then
        rri2 = 1.0 / rr; rri6 = rri2 * rri2 * rri2
        pot = pot + 4.0 * ep * (rri6 * rri6 - rri6)
      end if
      end do
    end do
    !$OMP END PARALLEL DO
    Potent = pot
  end function
  subroutine write_energy(ns)
    implicit none
    integer, intent(in) :: ns
    integer :: i
    real :: Ekin, Epot, E, Tp
    Ekin = sum(dens(:,:,2)**2.) * ep / 2.0
    Epot = Potent(dens(:,:,1))
    E = Ekin + Epot
    Tp = Ekin / (1.5 * n* kb)
    if (ns == 0) then
      open(unit = 2, file = "Data/Energy.txt", status = "replace", action = "write")
      write(2,*) "Time (ps)", "Ekin (1e12erg)", "Epot (1e12erg)", "E (1e12erg)", "T (K)"
    else
      open(unit = 2, file = "Data/Energy.txt", action = "write", position = "append")
    end if
    write(2,*) ns*dt*1e12*(m/ep)**0.5*sigma, Ekin*1e12, Epot*1e12, E*1e12, Tp*1e12
    close(2)
  end subroutine

end module
program main
 use liquid
 implicit none
 integer :: it
 real :: Tp
 call initial()
 do it = 0, ntrelax
    if (mod(it, 100) == 0) then
      print*, "Current time: ", it*dt*1e12*(m/ep)**0.5*sigma
      print*, "Temperature: ", sum(dens(:,:,2)**2.) * m / (3.0 * n * kb) * (ep / m)
    end if
    call writetraject(it)
    call write_energy(it)
    call verlexp()
 end do
 print*, "Relaxation finished"
 Tp = sum(dens(:,:,2)**2.) * m / (3.0 * n * kb) * (ep / m)
 dens(:,:,2) = dens(:,:,2) * sqrt(T / Tp)
 do it = ntrelax + 1, ntrelax + nt
    if (mod(it, 1000) == 0) then
      print*, "Current time: ", it*dt*1e12*(m/ep)**0.5*sigma
      print*, "Temperature: ", sum(dens(:,:,2)**2.) * m / (3.0 * n * kb) * (ep / m)
    end if
    call writetraject(it)
    call write_energy(it)
    call verlexp()
 end do
 deallocate(dens)
end program
