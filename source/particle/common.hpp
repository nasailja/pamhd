/*
Common particle functions of PAMHD.

Copyright 2014, 2015 Ilja Honkonen
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice, this
  list of conditions and the following disclaimer in the documentation and/or
  other materials provided with the distribution.

* Neither the names of the copyright holders nor the names of their contributors
  may be used to endorse or promote products derived from this software
  without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#ifndef PAMHD_PARTICLE_COMMON_HPP
#define PAMHD_PARTICLE_COMMON_HPP


#include "cmath"
#include "iostream"
#include "random"
#include "vector"


namespace pamhd {
namespace particle {


/*!
Returns kinetic energy of given particle.

Assumes Vector provides an API identical to Eigen vectors.
*/
template <class Vector> double get_kinetic_energy(
	const double mass,
	const Vector& velocity
) {
	return 0.5 * mass * velocity.squaredNorm();
}


/*!
Returns momentum of given particle.

Assumes Vector provides an API identical to Eigen vectors.
*/
template <class Vector> Vector get_momentum(
	const double mass,
	const Vector& velocity
) {
	return mass * velocity;
}


/*!
Returns radius and frequency of gyration of a particle in a magnetic field.

Assumes Vector provides an API identical to Eigen vectors.

Returns positive values, pair.first == radius, pair.second == frequency.
*/
template<class Vector> std::pair<double, double> get_gyro_info(
	const double charge_to_mass_ratio,
	const Vector& velocity,
	const Vector& magnetic_field
) {
	using std::fabs;
	using std::make_pair;

	const double
		B_magnitude = magnetic_field.norm(),
		c2m_abs = fabs(charge_to_mass_ratio);

	const Vector
		unit_B = magnetic_field / B_magnitude,
		perpendicular_velocity = velocity - velocity.dot(unit_B) * unit_B;

	return make_pair(
		perpendicular_velocity.norm() / c2m_abs / B_magnitude,
		c2m_abs * B_magnitude / 2 / M_PI
	);
}


/*!
Returns minimum and maximum allowed time step for a particle.

Maximum time step is a minimum value from:
	- flight time of particle through a cell
	- given fraction of particle's gyro period

Minimum time step is given fraction of maximum time step.

Assumes all simulation cells are cubes of given length.
*/
template<class Vector> std::pair<double, double> get_minmax_step(
	const double max_step_fraction,
	const double gyro_period_fraction,
	const double cell_length,
	const double charge_to_mass_ratio,
	const Vector& velocity,
	const Vector& magnetic_field
) {
	using std::make_pair;
	using std::min;

	const auto gyro_info
		= get_gyro_info(
			charge_to_mass_ratio,
			velocity,
			magnetic_field
		);

	const auto max_step
		= min(
			cell_length / velocity.norm(),
			gyro_period_fraction / gyro_info.second
		);

	return make_pair(max_step_fraction * max_step, max_step);
}


/*!
Returns particles with given average temperature and velocity.

Temperature components are used for standard deviation
of velocity in respective dimensions.
Total temperature is a sum of the components.

total_mass is the total mass assigned to created particles.
species_mass is used for particles' velocity distribution.

random_source is assumed to be std::mt19937 or similar.
*/
template <
	class Particle,
	class Mass_T,
	class Charge_Mass_Ratio_T,
	class Position_T,
	class Velocity_T,
	class Random_Source
> std::vector<Particle> create_particles(
	const typename Velocity_T::data_type bulk_velocity,
	const typename Position_T::data_type volume_min,
	const typename Position_T::data_type volume_max,
	const typename Velocity_T::data_type temperature,
	const size_t nr_of_particles,
	const double charge_mass_ratio,
	const double total_mass,
	const double species_mass,
	const double particle_temp_nrj_ratio,
	Random_Source& random_source
) {
	using std::sqrt;

	const Mass_T Mas{};
	const Position_T Pos{};
	const Velocity_T Vel{};
	const Charge_Mass_Ratio_T C2M{};

	const double particle_mass = total_mass / nr_of_particles;

	const auto
		std_dev_x = sqrt(particle_temp_nrj_ratio * temperature[0] / species_mass),
		std_dev_y = sqrt(particle_temp_nrj_ratio * temperature[1] / species_mass),
		std_dev_z = sqrt(particle_temp_nrj_ratio * temperature[2] / species_mass);

	std::normal_distribution<>
		velocity_generator_x(bulk_velocity[0], std_dev_x),
		velocity_generator_y(bulk_velocity[1], std_dev_y),
		velocity_generator_z(bulk_velocity[2], std_dev_z);

	std::uniform_real_distribution<>
		position_generator_x(volume_min[0], volume_max[0]),
		position_generator_y(volume_min[1], volume_max[1]),
		position_generator_z(volume_min[2], volume_max[2]);

	std::vector<Particle> particles(nr_of_particles);
	for (auto& particle: particles) {
		particle[Mas] = particle_mass;
		particle[C2M] = charge_mass_ratio;
		particle[Pos][0] = position_generator_x(random_source);
		particle[Pos][1] = position_generator_y(random_source);
		particle[Pos][2] = position_generator_z(random_source);
		particle[Vel][0] = velocity_generator_x(random_source);
		particle[Vel][1] = velocity_generator_y(random_source);
		particle[Vel][2] = velocity_generator_z(random_source);
	}

	return particles;
}


}} // namespaces


#endif // ifndef PAMHD_PARTICLE_COMMON_HPP
