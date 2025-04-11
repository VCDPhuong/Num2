module liquid
  implicit none
  integer :: n = 864
  integer :: nt, ntrelax, bc = 1, tr = 0 !(bc = 1) for periodic boundary condition, (tr = 1) for poor-mans algorithm
  real, parameter :: T = 94.4 !K
  real :: dt = 0.01, tmin = 0., tmax = 10., trelax = 6. !ps
  real, parameter :: kb = 1.380649*10.**(-16.) ! erg/K (1 erg = 1 g·cm²/s²)
  real, parameter :: m = 39.95*1.67*10.**(-24) ! g (mass of Argon atom)
  real, parameter :: sigma = 3.4e-8 !cm
  real, allocatable :: pos(:,:,:) !index: particle, x, v, a
  real, parameter :: ep = 120.*kb, rc = 2.25, L = 10.229, halfL = 10.229/2. !unitless length in sigma
  real, allocatable :: initial_velocities(:,:)
contains
  subroutine cal()
    implicit none
    nt = int((tmax-tmin)/dt)
    ntrelax = int(trelax/dt)
    if (rc > L/2.) then
      print*, "Error: cut off radius is too large"
      print*, "nr must be less than nl/2"
      stop
    end if
    print*, "Relaxation time: ", trelax, "ps"
    print*, "Number of time steps: ", nt
    print*, "Number of relaxation time steps: ", ntrelax
    print*, "Box length: ", L*sigma*1e8, "angstrong"
    print*, "Cut off radius: ", rc*sigma*1e8, "angstrong"
    print*, "Number of particles: ", n
    print*, "Time step: ", dt, "ps"
    print*, "Total time simuation: ", tmax - tmin, "ps"
    print*, "Minimum energy: ", LJpotent(2.**(1./6.)), -ep, "erg" !2^(1./6.) = 1.122462048
    print *, "Dimensionless transform: Begin"
    dt = dt*1e-12 * (ep/m)**(1./2.)*(1./sigma) ! Convert to dimensionless units
  end subroutine

  subroutine Initial()
    implicit none
    integer :: i, j
    integer :: nx, ny, nz, ix, iy, iz, p
    real :: spacing
    real :: v_mean, v_stddev, rand1, rand2, Tp, E, K
    allocate(pos(n,3,2))
    ! Initialize positions randomly within the box
    spacing = L / (int(n**(1.0/3.0)) + 1) ! Adjust spacing to avoid boundary
    p = 0
    do iz = 0, int(n**(1.0/3.0))
    do iy = 0, int(n**(1.0/3.0))
    do ix = 0, int(n**(1.0/3.0))
      p = p + 1
      if (p > n) exit
      pos(p,1,1) = ix * spacing
      pos(p,2,1) = iy * spacing
      pos(p,3,1) = iz * spacing
    end do
    if (p > n) exit
    end do
    if (p > n) exit
    end do
    ! Initialize velocities using a Maxwell-Boltzmann distribution
    v_mean = 0.0
    v_stddev = sqrt(kb * T / m)
    do i = 1, n
      do j = 1, 3
      call random_number(rand1)
      call random_number(rand2)
      pos(i,j,2) = v_mean + v_stddev * sqrt(-2.0 * log(rand1)) * cos(2.0 * 3.14159265359 * rand2)
      end do
    end do
    print*, "Initial velocities generated using Maxwell-Boltzmann distribution: Done"
    ! Remove any net momentum
    do i = 1, 3
     pos(:,i,2) = pos(:,i,2) - sum(pos(:,i,2))/n
    end do
    Tp = sum(pos(:,:,2)**2.)*m/(3.*kb*n)
    print*, "Final value: ", Tp, "Incompare with T: ", T, "relative error", abs(Tp - T)/T
    print*, "Expect Kinetic Energy: ", 3./2.*n*kb*T, "erg"
    print*, "Initial Kinetic Energy: ", 1./2.*m*sum(pos(:,:,2)**2.), "erg"
    pos(:,:,2) = pos(:,:,2) / (ep/m)**0.5 ! Convert to dimensionless units
  end subroutine
  subroutine calaccf(acc, rhos)
    implicit none
    integer :: i, j
    real :: dx, dy, dz, rr
    real :: rrinverse2, rrinverse6, coef
    real, intent(out) :: acc(n,3)
    real, intent(in) :: rhos(n,3)
    acc = 0.0
    !$OMP PARALLEL DO PRIVATE(i,j,dx,dy,dz,rr,rrinverse2,rrinverse6,coef) REDUCTION(+:acc) SCHEDULE(dynamic)
    do i = 1, n
      do j = i + 1, n
      dx = rhos(i,1) - rhos(j,1)
      dy = rhos(i,2) - rhos(j,2)
      dz = rhos(i,3) - rhos(j,3)
      if (dx > halfL) dx = dx - L
      if (dx < -halfL) dx = dx + L
      if (dy > halfL) dy = dy - L
      if (dy < -halfL) dy = dy + L
      if (dz > halfL) dz = dz - L
      if (dz < -halfL) dz = dz + L
      rr = dx**2 + dy**2 + dz**2
      if (rr < rc**2) then
        rrinverse2 = 1.0 / rr
        rrinverse6 = rrinverse2**3
        coef = rrinverse2 * (2.0 * rrinverse6 * rrinverse6 - rrinverse6)
        acc(i,1) = acc(i,1) + coef * dx
        acc(i,2) = acc(i,2) + coef * dy
        acc(i,3) = acc(i,3) + coef * dz
        acc(j,1) = acc(j,1) - coef * dx
        acc(j,2) = acc(j,2) - coef * dy
        acc(j,3) = acc(j,3) - coef * dz
      end if
      end do
    end do
    !$OMP END PARALLEL DO
    acc = acc * 24.0
    !acc = 0.
  end subroutine
subroutine verlexp(it, pos)
  implicit none
  integer, intent(in) :: it
  integer :: i, j
  real, intent(inout) :: pos(n,3,2)
  real :: acc(n,3), temp(n,3,2), Tp, vhalf(n,3)
  call calaccf(acc, pos(:,:,1))
  !$OMP PARALLEL DO PRIVATE(i,j) SCHEDULE(dynamic)
  do i = 1, n
    do j = 1, 3
      vhalf(i,j) = pos(i,j,2) + 0.5 * dt * acc(i,j)
      pos(i,j,1) = pos(i,j,1) + dt * vhalf(i,j)
    end do
  end do
  !$OMP END PARALLEL DO
  call calaccf(acc, pos(:,:,1))
  !$OMP PARALLEL DO PRIVATE(i,j) SCHEDULE(dynamic)
  do i = 1, n
    do j = 1, 3
      pos(i,j,2) = vhalf(i,j) + 0.5 * dt * acc(i,j)
    end do
  end do
  !$OMP END PARALLEL DO
  ! Apply periodic boundary conditions
  do i = 1, n
    do j = 1, 3
      pos(i,j,1) = mod(pos(i,j,1), L)
      if (pos(i,j,1) < 0.) pos(i,j,1) = pos(i,j,1) + L
    end do
  end do
  Tp = sum(pos(:,:,2)**2.) * m / (3. * kb * n) * (ep / m) ! Temperature in Kelvin
  if ((abs(Tp - T) > 5) .AND. (tr == 1)) then
    !$OMP PARALLEL DO PRIVATE(i,j) SCHEDULE(dynamic)
    do i = 1, n
      do j = 1, 3
        pos(i,j,2) = pos(i,j,2) * sqrt(T / Tp)
      end do
    end do
    !$OMP END PARALLEL DO
  end if
end subroutine
subroutine writedata(ns)
  implicit none
  integer, intent(in) :: ns
  real :: Tp, E, P, rr, sigs, dx, dy, dz
  real :: rrinverse2, rrinverse6
  integer :: i, j
  if (ns == 0) then
   open(unit = 1, file = "Data/Argon.txt", status = "replace", action = "write")
  else
   open(unit = 1, file = "Data/Argon.txt", status = "old", action = "write", position = "append")
  end if
  Tp = sum(pos(:,:,2)**2.) * m/(3.*kb*n) * (ep/m) ! Temperature in Kelvin
  call calpotent(E, pos(:,:,1))
  write(1,'(F20.7,F20.7,F20.7,F20.7,F20.7)') tmin + ns*dt*sigma*(M/ep)**(1/2.)*1e12, Tp, &
                3./2.*n*kb*Tp, E
  close(1)
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
    write(2,*) i, pos(i,1,1), pos(i,2,1), pos(i,3,1), &
                pos(i,1,2)*(m/ep)**0.5, pos(i,2,2)*(m/ep)**0.5, pos(i,3,2)*(m/ep)**0.5
  end do
  close(2)
end subroutine
subroutine calpotent(E, posi)
  implicit none
  real :: posi(n,3), dx, dy, dz, rr, rrmin
  real :: E, rrinverse2, rrinverse6, rrinverse
  integer :: i, j
  E = 0.0
  rrmin = 0.0

  !$OMP PARALLEL DO PRIVATE(i, j, dx, dy, dz, rr, rrinverse2, rrinverse6) REDUCTION(+:E)
  do i = 1, n
    do j = i + 1, n
      dx = posi(i,1) - posi(j,1)
      dy = posi(i,2) - posi(j,2)
      dz = posi(i,3) - posi(j,3)
      if (dx > halfL) dx = dx - L
      if (dx < -halfL) dx = dx + L
      if (dy > halfL) dy = dy - L
      if (dy < -halfL) dy = dy + L
      if (dz > halfL) dz = dz - L
      if (dz < -halfL) dz = dz + L
      rr = sqrt(dx**2 + dy**2 + dz**2) ! Convert to cm
      if (rr == 0.) then
        !$OMP CRITICAL
        print*, "Error: Two particles are at the same position"
        print*, "Particle 1: ", i, "Position: ", posi(i,1), posi(i,2), posi(i,3)
        print*, "Particle 2: ", j, "Position: ", posi(j,1), posi(j,2), posi(j,3)
        stop
        !$OMP END CRITICAL
      end if
      if (rr < rc) then
        E = E + LJpotent(rr) ! Energy in erg
      end if
    end do
  end do
  !$OMP END PARALLEL DO
end subroutine
subroutine write_distant()
  implicit none
  integer :: i, j, bin
  real :: dx, dy, dz, rr, dr
  real, allocatable :: g(:)
  real :: volume, density, normalization
  integer, parameter :: nbins = 1000
  real, parameter :: rmax = 4.0

  ! Allocate and initialize the radial distribution function array
  allocate(g(nbins))
  g = 0.0

  ! Bin width
  dr = rmax / nbins

  ! Calculate the radial distribution function
  do i = 1, n
    do j = i + 1, n
      dx = pos(i,1,1) - pos(j,1,1)
      dy = pos(i,2,1) - pos(j,2,1)
      dz = pos(i,3,1) - pos(j,3,1)

      ! Apply periodic boundary conditions
      if (dx > halfL) dx = dx - L
      if (dx < -halfL) dx = dx + L
      if (dy > halfL) dy = dy - L
      if (dy < -halfL) dy = dy + L
      if (dz > halfL) dz = dz - L
      if (dz < -halfL) dz = dz + L

      rr = sqrt(dx**2 + dy**2 + dz**2)
      if (rr < rmax) then
        bin = int(rr / dr) + 1
        g(bin) = g(bin) + 2.0
      end if
    end do
  end do

  ! Normalize the radial distribution function
  volume = L**3
  density = n / volume
  do bin = 1, nbins
    rr = (bin - 0.5) * dr
    normalization = 4.0 / 3.0 * 3.14159265359 * ((rr + dr)**3 - rr**3) * density * n
    g(bin) = g(bin) / normalization
  end do

  ! Write the radial distribution function to a file
  open(unit = 6, file = "Data/RadialDistribution.txt", status = "replace", action = "write")
  write(6,*) "Distance (sigma)", "g(r)"

  do bin = 1, nbins
    rr = (bin - 0.5) * dr
    write(6, '(F20.7, F20.7)') rr, g(bin)
  end do

  close(6)
  deallocate(g)
end subroutine
subroutine write_velocity_autocorrelation(it)
  implicit none
  integer, intent(in) :: it
  real :: autocorrelation(3)
  integer :: i, j

  ! Allocate and initialize initial velocities on the first call
  if (it == 0) then
    allocate(initial_velocities(n,3))
    initial_velocities = pos(:,:,2)
    open(unit = 5, file = "Data/VelocityAutocorrelation.txt", status = "replace", action = "write")
    write(5,*) "Time (ps)", "VACF_x", "VACF_y", "VACF_z"
  else
    open(unit = 5, file = "Data/VelocityAutocorrelation.txt", action = "write", position = "append")
  end if

  ! Calculate velocity autocorrelation function
  autocorrelation = 0.0
  do i = 1, n
    do j = 1, 3
      autocorrelation(j) = autocorrelation(j) + initial_velocities(i,j) * pos(i,j,2)
    end do
  end do
  autocorrelation = autocorrelation / n

  ! Write autocorrelation data to file
  write(5, '(F20.7, F20.7, F20.7, F20.7)') tmin + it*dt*sigma*(m/ep)**(1/2.)*1e12, &
    autocorrelation(1), autocorrelation(2), autocorrelation(3)
  close(5)
end subroutine
subroutine write_energy(it, pos)
  implicit none
  integer, intent(in) :: it
  real :: potential_energy
  real, intent(in) :: pos(n,3,2)
  real :: kinetic_energy, total_energy
  integer :: i

  ! Calculate kinetic energy
  kinetic_energy = sum(pos(:,:,2)**2.) * 0.5 * ep ! Convert to erg
  ! Calculate potential energy
  potential_energy = 0.0
  call calpotent(potential_energy, pos(:,:,1)) !potential energy in erg
  ! Calculate total energy
  total_energy = kinetic_energy + potential_energy
  ! Write energy data to file
  if (it == 0) then
    open(unit = 3, file = "Data/Energy.txt", status = "replace", action = "write")
    write(3,*) " Time (ps)", " Kinetic Energy (erg) ", " Potential Energy (1e-12erg) ", " Total Energy (erg) "
  else
    open(unit = 3, file = "Data/Energy.txt", action = "write", position = "append")
  end if
  write(3, '(F20.7, F20.7, F20.7, F20.7)') tmin + it*dt*sigma*(m/ep)**(1/2.)*1e12, &
  kinetic_energy*1e12, potential_energy*1e12, total_energy*1e12
  close(3)
end subroutine
 real function LJpotent(rrs) ! rrs is in sigma
  implicit none
  real :: rrs
  real*8 :: rrinverse, rrinverse2, rrinverse6, V
  rrinverse = 1./rrs
  rrinverse2 = rrinverse**2.
  rrinverse6 = rrinverse2**3.
  LJpotent =  ep * (4.0 * rrinverse6 * (rrinverse6 - 1.0)) ! Energy in erg
 end function
end module
program liquid_argon
 use liquid
 implicit none
 integer :: it, is
 real :: Tp
 call cal()
 call Initial()
 do it = 0, ntrelax
  call verlexp(it, pos)
 end do
 Tp = sum(pos(:,:,2)**2.)*m/(3.*kb*n) * (ep/m) ! Temperature in Kelvin
 !pos(:,:,2) = pos(:,:,2) * sqrt(T / Tp)
 print*, "Relaxation: Done"
 open(unit = 1, file = "Data/Inivelocity.txt", status = "replace", action = "write")
    do is = 1, n
     write(1,*) pos(is,1,2), pos(is,2,2), pos(is,3,2), &
                sqrt(pos(is,1,2)**2. + pos(is,2,2)**2. + pos(is,3,2)**2.)
    end do
 close(1)
 call write_distant()
 print*, "Done write distant"
 do it = 0 , nt
  if (mod(it,100) == 0) then
   print '(A,F10.3,x,A)', "time: ", tmin + it*dt*sigma*(m/ep)**(1/2.)*1e12, "ps"
  end if
  call writedata(it)
  call verlexp(it, pos)
  call writetraject(it)
  call write_energy(it, pos)
  call write_velocity_autocorrelation(it)
 end do
end program