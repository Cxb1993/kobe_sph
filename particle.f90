module particle
  !integer, parameter :: n_particle = 125
  !integer, parameter :: n_particle = 225
  !integer, parameter :: n_particle = 1000
  integer, parameter :: n_particle = 3375
  !integer, parameter :: n_particle = 40000
  !integer, parameter :: nloop = 300
  integer, parameter :: nloop = 600
  integer, parameter :: ndim = 3
  !real(8), parameter :: h = 0.01d0 ! m
  !real(8), parameter :: h = 0.1d0 ! m
  !real(8), parameter :: h = 0.2d0 ! m
  !real(8), parameter :: h = 0.3d0 ! m
  !real(8), parameter :: h = 0.5d0 ! m
  !real(8), parameter :: h = 0.6d0 ! m
  !real(8), parameter :: h = 0.7d0 ! m
  !real(8), parameter :: h = 0.8d0 ! m
  !real(8), parameter :: h = 1.0d0 ! m
  !real(8), parameter :: h = 2.0d0 ! m
  !real(8), parameter :: h = 3.0d0 ! m
  !real(8), parameter :: h = 4.0d0 ! m
  !real(8), parameter :: h = 5.0d0 ! m
  !real(8), parameter :: h = 8.0d0 ! m
  real(8), parameter :: h = 10.0d0 ! m
  !real(8), parameter :: h = 15.0d0 ! m
  !real(8), parameter :: h = 20.0d0 ! m
  !real(8), parameter :: h = 30.0d0 ! m
  !real(8), parameter :: area = 0.125d0 ! m
  !real(8), parameter :: area = 0.0333d0 ! m
  !real(8), parameter :: area = 0.15d0 ! m
  !real(8), parameter :: area = 0.1d0 ! m
  real(8), parameter :: area = 1.0d0 ! m
  !real(8), parameter :: area = 10.0d0 ! m
  !real(8), parameter :: Lx = 0.5d0 ! m
  real(8), parameter :: Lx = 10.0d0 ! m
  !real(8), parameter :: Ly = 1.0d0 ! m
  real(8), parameter :: Ly = 10.0d0 ! m
  !real(8), parameter :: Lz = 1.0d0 ! m
  real(8), parameter :: Lz = 15.0d0 ! m
  real(8), parameter :: rho0 = 1000.0d0 ! kg / m^3
  !real(8), parameter :: rho0 = 600.0d0 ! kg / m^3
  !real(8), parameter :: rho0 = 100.0d0 ! kg / m^3
  !real(8), parameter :: rho0 = 0.01d0 ! kg / m^3
  !real(8), parameter :: rho0 = 0.05d0 ! kg / m^3
  !real(8), parameter :: rho0 = 0.1d0 ! kg / m^3
  !real(8), parameter :: rho0 = 0.3d0 ! kg / m^3
  !real(8), parameter :: rho0 = 0.5d0 ! kg / m^3
  !real(8), parameter :: rho0 = 1.0d0 ! kg / m^3
  !real(8), parameter :: rho0 = 0.001d0 ! kg / m^3
  !real(8), parameter :: rho0 = 0.002d0 ! kg / m^3
  !real(8), parameter :: rho0 = 0.005d0 ! kg / m^3
  !real(8), parameter :: rho0 = 0.0001d0 ! kg / m^3
  !real(8), parameter :: rho0 = 0.0d0 ! kg / m^3
  !real(8), parameter :: mass = 0.00020543d0 ! kg
  real(8), parameter :: mass = rho0 * Lx * Ly * Lz / n_particle ! kg
  real(8), parameter :: k = 1.0d0 ! gas constant
  real(8), parameter :: mu = 0.2 ! pascal-second (Pa.s) = 1 kg m^-1 s^-1
  !real(8), parameter :: dt = 0.004
  !real(8), parameter :: dt = 0.001
  real(8), parameter :: dt = 0.01
  !real(8), parameter :: dt = 0.0001
  real(8), parameter :: radius = 0.004 ! m
  real(8), parameter :: ext_stiff = 10000.0
  !real(8), parameter :: ext_stiff = 30000.0
  !real(8), parameter :: ext_stiff = 100000.0
  !real(8), parameter :: ext_stiff = 1000000.0
  !real(8), parameter :: ext_stiff = 1500000.0
  !real(8), parameter :: ext_stiff = 1800000.0
  !real(8), parameter :: ext_stiff = 2000000.0
  !real(8), parameter :: ext_stiff = 2200000.0
  !real(8), parameter :: ext_stiff = 0.0
  !real(8), parameter :: ext_damp = 256.0
  !real(8), parameter :: ext_damp = 750.0
  !real(8), parameter :: ext_damp = 1300.0
  !real(8), parameter :: ext_damp = 1500.0
  !real(8), parameter :: ext_damp = 2560.0
  !real(8), parameter :: ext_damp = 25600.0
  !real(8), parameter :: ext_damp = 0.0
  real(8), parameter :: ext_damp = 50.0
  real(8), parameter :: eps = 0.00001
  real(8) :: x_min(ndim) ! boundary condition
  real(8) :: x_max(ndim) ! boundary condition
  real(8) :: x(ndim, n_particle) ! position
  real(8) :: v(ndim, n_particle) ! velosity
  real(8) :: rho(n_particle) ! den0sity
  real(8) :: p(n_particle) ! pressure
  real(8) :: f(ndim, n_particle) ! force 
end module particle
