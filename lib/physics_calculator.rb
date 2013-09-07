require "physics_calculator/version"

module PhysicsCalculator

  include Math
  # methods for matrix and vector manipulation
  require 'matrix'

  # extend standard library vector class

  class ::Vector

    # dot_product as alias for inner product
    def dot_product(v_)
      self.inner_product(v_)
    end

    # define cross product
    def cross_product_(v_)
      raise ArgumentError, 'Vectors must be 3D' if (self.size != 3 || v_.size != 3)
      Vector.[](self[1]*v_[2] - self[2]*v_[1], self[2]*v_[0] - self[0]*v_[2], self[0]*v_[1] - self[1]*v_[0])
    end

    # define distance between two vectors
    def distance(v_)
      (self - v_).r
    end

    # define unit vector pointing from v to self
    def unit_(v_)
      (self - v_).normalize
    end

  end

  # Physical Constants (SI Units)


  # speed of light in a vacuum (m/s)
  SPEED_LIGHT = 299792458.0

  # planck's constant (J*s)
  PLANCK = 6.62606957e-34

  # reduced planck's constant
  PLANCK_REDUCED = 1.054571726e-34

  # vacuum permittivity (F/m)
  PERMITTIVITY = 8.854187817e-12

  # vacuum permeability ((V*s)/(A*m))
  PERMEABILITY = 1.2566370614e-6

  # elementary charge (C)
  CHARGE_ELEMENTARY = 1.602176565e-19

  # fine structure constant (dimensionless)
  FINE_STRUCTURE = 7.2973525698e-3

  # electron rest mass (kg)
  MASS_ELECTRON = 9.10938215e-31

  # proton rest mass (kg)
  MASS_PROTON = 1.672621777e-27

  # neutron rest mass (kg)
  MASS_NEUTRON = 1.674927351e-27

  # mass of the Earth (kg)
  MASS_EARTH = 5.97217e24

  # mass of the Sun (kg)
  MASS_SUN = 1.98855e30

  # mass of the Moon (kg)
  MASS_MOON = 7.3477e22

  # gravitational constant ((N*m^2)/kg^2))
  GRAVITATIONAL = 6.67384e-11

  # acceleration due to gravity on Earth's surface (m/s^2)
  ACCELERATION_GRAVITY = 9.80665

  # Boltzmann constant
  BOLTZMANN = 1.3806488e-23

  # energy of ground state of electron in hydrogen atom (Joules)
  ENERGY_BINDING = -2.1799e-19

  # Rydberg constant (1/m)
  RYDBERG = 1.09739e7

  # Aliases for convenience


  C = SPEED_LIGHT
  G = GRAVITATIONAL
  H = PLANCK_REDUCED
  I = Complex::I
  X_ = Vector.[](1,0,0)
  Y_ = Vector.[](0,1,0)
  Z_ = Vector.[](0,0,1)
  ACCELERATION_GRAVITY_ = -ACCELERATION_GRAVITY*Z_


  # arguments with suffix _ expect 3 dimensional vectors
  # arguments expect SI units


  # classical mechanics


  # gravitational force on mass1 at position1_ due to mass2 at position2_ (vector in Newtons)
  def force_gravity_(mass1, position1_, mass2, position2_)
    ((-G * mass1 * mass2) / (position1_.distance(position2_)**2)) * position1_.unit_(position2_)
  end

  # weight of mass at surface of Earth with -Z_ pointing at the center of Earth (vector in Newtons)
  def force_weight_(mass)
    ACCELERATION_GRAVITY_ * mass
  end

  # position of object under constant acceleration_ at time given initial conditions position0_ and velocity0_ (vector in meters)
  def position_constant_acceleration_(position0_, velocity0_, acceleration_, time)
    0.5*acceleration_*time**2 + velocity0_ * time + position0_
  end

  # velocity of object under constant acceleration_ at time given initial condition velocity0_ (vector in m / s)
  def velocity_constant_acceleration_(velocity0_, acceleration_, time)
    velocity0_ + acceleration_ * time
  end

  # position of object under acceleration due to gravity at time given initial conditions position0_ and velocity0_ (vector in meters)
  def position_gravity_acceleration_(position0_, velocity0_, time)
    position_constant_acceleration_(position0_, velocity0_, ACCELERATION_GRAVITY_, time)
  end

  # velocity of object under acceleration due to gravity at time given initial condition velocity0_ (vector in m / s)
  def velocity_gravity_acceleration_(velocity0_, time)
    velocity_constant_acceleration_(velocity0_, ACCELERATION_GRAVITY_, time)
  end

  # linear momentum of mass at velocity_ (vector in kg * m^2)
  def momentum_(mass, velocity_)
    mass * velocity_
  end

  # kinetic energy of mass at velocity_ (scalar in Joules)
  def energy_kinetic(mass, velocity_)
    0.5*mass*(velocity_.r)**2
  end

  # torque from applying f_ at position_ around axis_  (vector in N * m)
  def torque_(axis_, position_, force_)
    (position_ - axis_).cross_product_(force_)
  end

  # angular momentum from momentum mass*velocity_ around axis_ (vector in J * s)
  def angular_momentum_(axis_, position_, mass, velocity_)
    (position_ - axis_).cross_product_(momentum_(mass,velocity_))
  end

  # angular momentum of rotating rigid body with moment_inertia rotating at angular_velocity_ (vector in J * s)
  def angular_momentum_rigid_(moment_inertia, angular_velocity_)
    moment_inertia * angular_velocity_
  end

  # reduced mass of mass1 and mass2 (scalar in kilograms)
  def mass_reduced(mass1, mass2)
    (Float(mass1*mass2))/(mass1 + mass2)
  end

  # moment of inertia using parallel axis theorem using mass with moment_inertia_cm with axis_distance between the center of mass and axis of rotation (scalar in kg * m^2)
  def moment_inertia_parallel_axis(moment_inertia_cm, mass, axis_distance)
    moment_inertia_cm + mass*axis_distance**2
  end

  # moment of inertia of mass at distance away from axis of rotation (scalar in kg * m^2)
  def moment_inertia_point_mass(mass, distance)
    mass * distance**2
  end

  # moment of inertia of mass1 and mass2 at distance away from axis of rotation (scalar in kg * m^2)
  def moment_inertia_2_point_mass(mass1, mass2, distance)
    mass_reduced(mass1, mass2)*distance**2
  end

  # moment of inertia of rod with length and mass rotating about center of mass (scalar in kg * m^2)
  def moment_inertia_rod(mass, length)
    (mass*length**2) / 12.0
  end

  # moment of inertia of hollow sphere with radius and mass rotating about center of mass (scalar in kg * m^2)
  def moment_inertia_sphere(mass, radius)
    (2*mass*radius**2) / 3.0
  end

  # moment of inertia of solid ball with radius and mass rotating about center of mass (scalar in kg * m^2)
  def moment_inertia_ball(mass, radius)
    (2*mass*radius**2) / 5.0
  end


  # quantum mechanics


  # de Broglie wavelength for a particle with momentum (scalar in meters)
  def wavelength_de_broglie(momentum)
    PLANCK / momentum
  end

  # Hermite polynomial of degree n evaluated at x (helper for QHO eigenstates)
  def hermite_polynomial(n, x)
    if n == 0
      1
    elsif n == 1
      2*x
    else
      (2*x*hermite_polynomial(n-1, x) - 2*(n-1)*hermite_polynomial(n-2, x))
    end
  end

  # one dimensional quantum harmonic oscillator nth eigenstate for particle with mass, frequency, and position (scalar dimensionless)
  def eigenstate_qho(n, mass, frequency, position)
    xi = Math.sqrt((mass*frequency) / H)*position
    (((mass*frequency)/(PI*H))**0.25)*(((2**n)*Math.gamma(n+1))**0.5)*hermite_polynomial(n, xi)*Math.exp((-xi**2) / 2)
  end

  # one dimensional quantum harmonic oscillator nth eigenstate energy for particle with mass, and frequency (scalar in Joules)
  def energy_qho(n, frequency)
    (n + 0.5) * H * frequency
  end

  # eigenstate of particle in one dimensional infinite square well potential from 0 to width with mass and position (scalar dimensionless)
  def eigenstate_infinite_well(n, width, position)
    sqrt(2.0 / width) * sin((n * PI * position) / width)
  end

  # energy of particle in one dimensional infinite square well potential from 0 to width with mass (scalar in Joules)
  def energy_infinite_well(n, width, mass)
    (n**2 * PI**2 * H**2) / (2 * mass * width**2)
  end

  # energy of the nth (n = 1, 2, 3 ...) state of the electron in a hydrogen atom according to the Bohr model (scalar in Joules)
  def energy_bohr(n)
    ENERGY_BINDING / n**2
  end

  # wavelength of light emitted in transition of electron from n_initial state to n_final state according to the Rydberg formula (scalar in meters)
  def wavelength_rydberg(n_initial, n_final)
    1.0 / (RYDBERG*(1.0 / (n_final**2) - 1.0 / (n_initial**2))).abs
  end


  # electricity and magnetism


  # Lorentz force on particle with charge moving at velocity_ through electric_field_ and magnetic_field_
  def force_lorentz_(charge, velocity_, electric_field_, magnetic_field_)
    charge * (electric_field_ + velocity_.cross_product_(magnetic_field_))
  end

  # force felt by charge1 at position1_ due to charge2 at position2_
  def force_coulombs_law_(charge1, position1_, charge2, position2_)
    (((1/(4*PI*PERMITTIVITY))*charge1*charge2)/(position1_.distance(position2_)**2))*position1_.unit_(position2_)
  end

  # electric field at position_ due to point charge at charge_position_ (vector in N / C)
  def electric_field_point_charge_(charge, charge_position_, position_)
    (((1/(4*PI*PERMITTIVITY))*charge)/(position_.distance(charge_position_)**2))*position_.unit_(charge_position_)
  end

  # magnetic field at position_ due to point charge at charge_position_ moving at charge_velocity_ (vector in Teslas)
  def magnetic_field_point_charge_(charge, charge_position_, charge_velocity_, position_)
    (PERMEABILITY*charge)/(4*PI*position_.distance(charge_position_)**2)*charge_velocity_.cross_product_(position_.unit_(charge_position_))
  end

  # equivalent resistance of resistance1 and resistance2 in series (scalar in ohms)
  def resistance_series(resistance1, resistance2)
    resistance1+resistance2
  end

  # equivalent resistance of resistance1 and resistance2 in parallel (scalar in ohms)
  def resistance_parallel(resistance1, resistance2)
    Float(resistance1*resistance2)/(resistance1+resistance2)
  end

  # equivalent inductance of inductance1 and inductance2 in series (scalar in henry)
  def inductance_series(inductance1, inductance2)
    inductance1+inductance2
  end

  # equivalent inductance of inductance1 and inductance2 in parallel (scalar in henry)
  def inductance_parallel(inductance1, inductance2)
    Float(inductance1*inductance2)/(inductance1+inductance2)
  end

  # equivalent capacitance of capacitance1 and capacitance2 in series (scalar in Farads)
  def capacitance_series(capacitance1, capacitance2)
    Float(capacitance1*capacitance2)/(capacitance1+capacitance2)
  end

  # equivalent capacitance of capacitance1 and capacitance2 in parallel (scalar in Farads)
  def capacitance_parallel(capacitance1, capacitance2)
    capacitance1+capacitance2
  end

  # voltage due to current across resistance (scalar in Volts)
  def voltage_ohm_law(current, resistance)
    current*resistance
  end


  # special relativity


  # lorentz factor for object moving at velocity (scalar dimensionless)
  def lorentz_factor(speed)
    (1.0 / (Math.sqrt(1.0-(speed/C)**2)))
  end

  # time dilation for an object moving at velocity over time_interval (scalar in seconds)
  def time_dilation(time_interval, speed)
    lorentz_factor(speed)*time_interval
  end

  # length contraction for an object with length moving at speed (scalar in meters)
  def length_contraction(length, speed)
    length / lorentz_factor(speed)
  end

  # relativistic mass of an object moving at speed (scalar in kilograms)
  def mass_relativistic(mass, speed)
    lorentz_factor(speed)*mass
  end

  # relativistic momentum of an object with mass moving at speed (scalar in kg * m / s)
  def momentum_relativistic(mass, speed)
    lorentz_factor(speed) * mass * speed
  end

  # rest energy of mass (scalar in Joules)
  def energy_rest(mass)
    mass*C**2
  end

  # energy of a photon with frequency (scalar in Joules)
  def energy_photon(frequency)
    PLANCK*frequency
  end

  # momentum of a photon with frequency (scalar in kg * m / s)
  def momentum_photon(frequency)
    energy_photon(frequency)/C
  end

end
