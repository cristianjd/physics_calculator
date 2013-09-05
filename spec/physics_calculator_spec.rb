require 'rspec'
require 'physics_calculator'
include PhysicsCalculator

RSpec::Matchers.define :be_close_vector_to do |expected|
  match do |actual|
    n = 0
    (0..2).each do |x|
      a_value = actual[x]
      e_value = expected[x]
      if e_value == 0
        tolerance = 0.0001
      else
        tolerance = e_value.abs / 10000.0
      end
      if (a_value > e_value - tolerance) and (a_value < e_value + tolerance)
        n += 1
      end
    end
    n == 3
  end
end

RSpec::Matchers.define :be_close_scalar_to do |expected|
  match do |actual|
    if expected == 0
      tolerance = 0.0001
    else
      tolerance = expected.abs / 10000.0
    end
    (actual > expected - tolerance) and (actual < expected + tolerance)
  end
end

describe PhysicsCalculator do


  vec1_ = Vector.[](1, 2, 3)
  vec2_ = Vector.[](4, 5, 6)
  vec3_ = Vector.[](1, 0, -1)

  # added vector methods

  describe Vector do

    describe '#dot_product' do
      it 'returns dot product of two vectors' do
        expect(vec1_.dot_product(vec2_)).to eq 32
      end
    end

    describe '#cross_product_' do
      it 'returns cross product of two vectors' do
        expect(vec1_.cross_product_(vec2_)).to be_close_vector_to [-3,6,-3]
      end

      it 'raises error if vectors are not 3D' do
        expect {Vector.[](4,2).cross_product_(Vector.[](3,0))}.to raise_error ArgumentError
      end
    end

    describe '#distance' do
      it 'returns magnitude of difference of two vectors' do
        expect(vec2_.distance(vec1_)).to be_close_scalar_to 5.19615
      end
    end

    describe '#unit_' do
      it 'returns correct unit vector pointing from first vector to second' do
        expect(vec2_.unit_(vec1_)).to be_close_vector_to [0.57735, 0.57735, 0.57735]
      end
    end

  end

  # classical mechanics

  describe '#force_gravity_' do
    it 'returns correct force of gravity' do
      expect(force_gravity_(10, vec1_, 20, vec2_)).to be_close_vector_to [2.85418e-10, 2.85418e-10, 2.85418e-10]
    end
  end

  describe '#force_weight_' do
    it 'returns the correct weight force' do
      expect(force_weight_(10)).to be_close_vector_to [0, 0, -98.0665]
    end
  end

  describe '#position_constant_acceleration_' do
    it 'returns correct position under constant acceleration' do
      expect(position_constant_acceleration_(vec1_, vec2_, vec3_, 10)).to be_close_vector_to [91, 52, 13]
    end
  end

  describe '#velocity_constant_acceleration_' do
    it 'returns correct velocity under constant acceleration' do
      expect(velocity_constant_acceleration_(vec1_, vec2_, 10)).to be_close_vector_to [41, 52, 63]
    end
  end

  describe '#position_gravity_acceleration_' do
    it 'returns correct position under acceleration due to gravity' do
      expect(position_gravity_acceleration_(vec1_, vec2_, 10)).to be_close_vector_to [41, 52, -427.3275]
    end
  end

  describe '#velocity_gravity_acceleration_' do
    it 'returns correct velocity under acceleration due to gravity' do
      expect(velocity_gravity_acceleration_(vec1_, 10)).to be_close_vector_to [1, 2, -95.0665]
    end
  end

  describe '#momentum_' do
    it 'returns correct linear momentum' do
      expect(momentum_(10, vec1_)).to be_close_vector_to [10, 20, 30]
    end
  end

  describe '#energy_kinetic' do
    it 'returns correct kinetic energy' do
      expect(energy_kinetic(10, vec1_)).to be_close_scalar_to 70
    end
  end

  describe '#torque_' do
    it 'returns correct torque' do
      expect(torque_(vec1_, vec2_, vec3_)).to be_close_vector_to [-3, 6, -3]
    end
  end

  describe '#angular_momentum_' do
    it 'returns correct angular momentum' do
      expect(angular_momentum_(vec1_, vec2_, 10, vec3_)).to be_close_vector_to [-30, 60, -30]
    end
  end

  describe '#angular_momentum_rigid_' do
    it 'returns correct angular momentum' do
      expect(angular_momentum_rigid_(10, vec1_)).to be_close_vector_to [10, 20, 30]
    end
  end

  describe '#mass_reduced' do
    it 'returns correct reduced mass' do
      expect(mass_reduced(3, 4)).to be_close_scalar_to 1.714286
    end
  end

  describe '#moment_inertia_parallel_axis' do
    it 'returns correct moment of inertia' do
      expect(moment_inertia_parallel_axis(20,10,5)).to be_close_scalar_to 270
    end
  end

  describe '#moment_inertia_point_mass' do
    it 'returns correct moment of inertia' do
      expect(moment_inertia_point_mass(2,10)).to be_close_scalar_to 200
    end
  end

  describe '#moment_inertia_2_point_mass' do
    it 'returns correct moment of inertia' do
      expect(moment_inertia_2_point_mass(4,6,3)).to be_close_scalar_to 21.6
    end
  end

  describe '#moment_inertia_rod' do
    it 'returns correct moment of inertia' do
      expect(moment_inertia_rod(6, 10)).to be_close_scalar_to 50
    end
  end

  describe '#moment_inertia_sphere' do
    it 'returns correct moment of inertia' do
      expect(moment_inertia_sphere(3,5)).to be_close_scalar_to 50
    end
  end

  describe '#moment_inertia_ball' do
    it 'returns correct moment of inertia' do
      expect(moment_inertia_ball(3,5)).to be_close_scalar_to 30
    end
  end

  # quantum  mechanics

  describe '#wavelength_de_broglie' do
    it 'returns correct de Broglie wavelength' do
      expect(wavelength_de_broglie(10)).to be_close_scalar_to 6.62607e-35
    end
  end

  describe '#hermite_polynomial' do
    it 'returns correct hermite polynomial value' do
      expect(hermite_polynomial(3, 5)).to be_close_scalar_to 940
    end
  end

  describe '#eigenstate_qho' do
    it 'returns correct value of qho eigenstate' do
      expect(eigenstate_qho(0, 5, 10, 0)).to be_close_scalar_to  6.23284e8
    end
  end

  describe '#energy_qho' do
    it 'returns correct qho energy' do
      expect(energy_qho(4, 10)).to be_close_scalar_to 4.746e-33
    end
  end

  describe '#eigenstate_infinite_well' do
    it 'returns correct infinite well eigenstate' do
      expect(eigenstate_infinite_well(4,2,10.1)).to be_close_scalar_to 0.587785
    end
  end


  describe '#energy_infinite_well' do
    it 'returns correct infinite well energy' do
      expect(energy_infinite_well(3, 10, 5)).to be_close_scalar_to 9.878585e-70
    end
  end

  describe '#energy_bohr' do
    it 'returns correct Bohr energy' do
      expect(energy_bohr(6)).to be_close_scalar_to -6.055278e-21
    end
  end

  describe '#wavelength_rydberg' do
    it 'returns correct wavelength using Rydberg formula' do
      expect(wavelength_rydberg(2,5)).to be_close_scalar_to 4.339e-7
    end
  end

  # electricity and magnetism

  describe '#force_lorentz_' do
    it 'returns correct Lorentz force' do
      expect(force_lorentz_(10, vec1_, vec2_, vec3_)).to be_close_vector_to [20,90,40]
    end
  end

  describe '#force_coulombs_law_' do
    it 'returns correct force using Coulomb law' do
      expect(force_coulombs_law_(5,vec1_,-10,vec2_)).to be_close_vector_to [9.6092e9, 9.6092e9, 9.6092e9]
    end
  end

  describe '#electric_field_point_charge_' do
    it 'returns correct electric field due to point charge' do
      expect(electric_field_point_charge_(10,vec1_, vec2_)).to be_close_vector_to [1.92184e9, 1.92184e9, 1.92184e9]
    end
  end

  describe '#magnetic_field_point_charge_' do
    it 'returns correct magnetic field due to point charge' do
      expect(magnetic_field_point_charge_(5,vec1_,vec2_,vec3_)).to be_close_vector_to [-4.47214e-8,8.94427e-8,-4.47214e-8]
    end
  end

  describe '#resistance_series' do
    it 'returns correct resistance' do
      expect(resistance_series(4,5)).to be_close_scalar_to 9
    end
  end

  describe '#resistance_parallel' do
    it 'returns correct resistance' do
      expect(resistance_parallel(4,5)).to be_close_scalar_to 2.22222
    end
  end

  describe '#inductance_series' do
    it 'returns correct inductance' do
      expect(inductance_series(4,5)).to be_close_scalar_to 9
    end
  end

  describe '#inductance_parallel' do
    it 'returns correct inductance' do
      expect(inductance_parallel(4,5)).to be_close_scalar_to 2.22222
    end
  end

  describe '#capacitance_series' do
    it 'returns correct capacitance' do
      expect(capacitance_series(4,5)).to be_close_scalar_to 2.22222
    end
  end

  describe '#capacitance_parallel' do
    it 'returns correct capacitance' do
      expect(capacitance_parallel(4,5)).to be_close_scalar_to 9
    end
  end


  # special relativity

  speed = 259620269

  describe  '#lorentz_factor' do
    it 'returns correct lorentz factor value' do
      expect(lorentz_factor(speed)).to be_within(0.001).of 2.0
    end
  end

  describe '#time_dilation' do
    it 'returns correct time dilation' do
      expect(time_dilation(10, speed)).to be_within(0.01).percent_of 20.0
    end
  end

  describe '#length_contraction' do
    it 'returns correct length contraction' do
      expect(length_contraction(10, speed)).to be_within(0.01).percent_of 5.0
    end
  end

  describe '#mass_relativistic' do
    it 'returns correct relativistic mass' do
      expect(mass_relativistic(10, speed)).to be_within(0.01).percent_of 20.0
    end
  end

  describe '#momentum_relativistic' do
    it 'returns correct relativistic momentum' do
      expect(momentum_relativistic(10, speed)).to be_within(0.01).percent_of 5191948531
    end
  end

  describe '#energy_rest' do
    it 'returns correct rest energy' do
      expect(energy_rest(10)).to be_within(0.01).percent_of 8.9875518e17
    end
  end

  describe '#energy_photon' do
    it 'returns correct photon energy' do
      expect(energy_photon(1e10)).to be_within(0.01).percent_of 6.62606957e-24
    end
  end

  describe '#momentum_photon' do
    it 'returns correct photon momentum' do
      expect(momentum_photon(1e10)).to be_within(0.01).percent_of 2.2102189e-32
    end
  end

end