/*
Tests bulk temparature calculation of particles in PAMHD.

Copyright 2015, 2016 Ilja Honkonen
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

#include "cmath"
#include "cstdlib"
#include "iostream"
#include "random"

#include "Eigen/Core"
#include "Eigen/Geometry"

#include "particle/common.hpp"
#include "particle/variables.hpp"


using namespace std;
using namespace Eigen;
using namespace pamhd::particle;


template<class Particle> Eigen::Vector3d get_bulk_momentum(const std::vector<Particle>& particles)
{
	Eigen::Vector3d bulk_momentum{0, 0, 0};
	for (const auto& particle: particles){
		bulk_momentum += particle[Mass()] * particle[Velocity()];
	}
	return bulk_momentum;
}


int main()
{
	using std::fabs;
	using std::min;
	using std::pow;
	using std::sqrt;

	const size_t nr_particles = 100000;
	const Vector3d
		volume_min(-5, 4, 9),
		volume_max(-1, 8, 33);
	const double particle_temp_nrj_ratio = 0.25;

	std::mt19937 random_source;
	random_source.seed(997);


	const Vector3d
		bulk_velocity1(0, 0, 0),
		temperature1(1, 1, 1);
	const double
		ref_temp1 = temperature1.sum() / 3,
		species_mass1 = 1,
		total_mass1 = nr_particles,
		species_particles1 = total_mass1 / species_mass1;
	const auto particles1
		= create_particles<
			Particle_Internal,
			Mass,
			Charge_Mass_Ratio,
			Position,
			Velocity,
			Particle_ID,
			Species_Mass
		>(
			bulk_velocity1,
			volume_min,
			volume_max,
			temperature1,
			nr_particles,
			10,
			total_mass1,
			species_mass1,
			particle_temp_nrj_ratio,
			random_source
		);
	const auto scalar_temp1
		= get_temperature<Mass, Velocity, Species_Mass>(particles1, particle_temp_nrj_ratio);
	if (abs(ref_temp1 - scalar_temp1) / ref_temp1 > 1e-2) {
		std::cerr <<  __FILE__ << "(" << __LINE__ << "): "
			<< ref_temp1 << " " << scalar_temp1
			<< std::endl;
		abort();
	}


	const Vector3d bulk_velocity2(-3, 0, 0);
	const auto particles2
		= create_particles<
			Particle_Internal,
			Mass,
			Charge_Mass_Ratio,
			Position,
			Velocity,
			Particle_ID,
			Species_Mass
		>(
			bulk_velocity2,
			volume_min,
			volume_max,
			temperature1,
			nr_particles,
			10,
			total_mass1,
			species_mass1,
			particle_temp_nrj_ratio,
			random_source
		);
	const auto scalar_temp2
		= get_temperature<Mass, Velocity, Species_Mass>(particles2, particle_temp_nrj_ratio);
	if (abs(ref_temp1 - scalar_temp2) / ref_temp1 > 1e-2) {
		std::cerr <<  __FILE__ << "(" << __LINE__ << "): "
			<< ref_temp1 << " " << scalar_temp1
			<< std::endl;
		abort();
	}

	auto particles3 = particles1;
	particles3.insert(particles3.end(), particles2.cbegin(), particles2.cend());
	const auto scalar_temp3
		= get_temperature<Mass, Velocity, Species_Mass>(particles3, particle_temp_nrj_ratio);
	const double ref_temp3
		= ref_temp1
		+ bulk_velocity2.squaredNorm() / 12 / particle_temp_nrj_ratio;
	if (abs(ref_temp3 - scalar_temp3) / ref_temp3 > 1e-2) {
		std::cerr <<  __FILE__ << "(" << __LINE__ << "): "
			<< ref_temp3 << " " << scalar_temp3
			<< std::endl;
		abort();
	}


	const Vector3d
		bulk_velocity4(30, 8, 5),
		temperature4(13, 17, 23);
	const auto particles4
		= create_particles<
			Particle_Internal,
			Mass,
			Charge_Mass_Ratio,
			Position,
			Velocity,
			Particle_ID,
			Species_Mass
		>(
			bulk_velocity4,
			volume_min,
			volume_max,
			temperature4,
			nr_particles,
			10,
			total_mass1,
			species_mass1,
			particle_temp_nrj_ratio,
			random_source
		);

	auto particles5 = particles2;
	particles5.insert(particles5.end(), particles4.cbegin(), particles4.cend());
	const auto scalar_temp5
		= get_temperature<Mass, Velocity, Species_Mass>(particles5, particle_temp_nrj_ratio);
	const double ref_temp5
		= (temperature1 + temperature4).sum() / 6
		+ (bulk_velocity2 - bulk_velocity4).squaredNorm()
			/ 12 / particle_temp_nrj_ratio;
	if (abs(ref_temp5 - scalar_temp5) / ref_temp5 > 1e-2) {
		std::cerr <<  __FILE__ << "(" << __LINE__ << "): "
			<< ref_temp5 << " " << scalar_temp5
			<< std::endl;
		abort();
	}


	const double
		total_mass6 = 1000 * nr_particles,
		species_particles6 = total_mass6 / species_mass1;
	const auto particles6
		= create_particles<
			Particle_Internal,
			Mass,
			Charge_Mass_Ratio,
			Position,
			Velocity,
			Particle_ID,
			Species_Mass
		>(
			bulk_velocity4,
			volume_min,
			volume_max,
			temperature4,
			nr_particles,
			10,
			total_mass6,
			species_mass1,
			particle_temp_nrj_ratio,
			random_source
		);

	auto particles7 = particles2;
	particles7.insert(particles7.end(), particles6.cbegin(), particles6.cend());
	const auto scalar_temp7
		= get_temperature<Mass, Velocity, Species_Mass>(particles7, particle_temp_nrj_ratio);
	const double ref_temp7
		= (temperature1 * species_particles1 + temperature4 * species_particles6).sum()
			/ (species_particles1 + species_particles6) / 3
		+ (bulk_velocity2 - bulk_velocity4).squaredNorm()
			/ 3 / particle_temp_nrj_ratio
			* species_particles1 * species_particles6
			/ pow(species_particles1 + species_particles6, 2);
	if (abs(ref_temp7 - scalar_temp7) / ref_temp7 > 1e-2) {
		std::cerr <<  __FILE__ << "(" << __LINE__ << "): "
			<< ref_temp7 << " " << scalar_temp7
			<< std::endl;
		abort();
	}


	return EXIT_SUCCESS;
}
