/*
Initializes particle solution of PAMHD.

Copyright 2015 Ilja Honkonen
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

#ifndef PAMHD_PARTICLE_INITIALIZE_HPP
#define PAMHD_PARTICLE_INITIALIZE_HPP


#include "iostream"
#include "random"

#include "particle/common.hpp"


namespace pamhd {
namespace particle {


/*!
Creates particles in given cells as defined by given initial conditions.

Returns the total number of particles created.

\param [Init_Cond] Compatible with pamhd::boundaries::Initial_Condition
\param [Grid] Compatible with DCCRG (github.com/fmihpc/dccrg)
\param [Number_Density_T] Used to access particle number density in init_cond
\param [Temperature_T] ...access bulk particle temperature in init_cond
\param [Charge_Mass_Ratio_T] ...charge to mass ratio of particles in init_cond
\param [Species_Mass_T] ...mass of particle species in init_cond to create
\param [Bulk_Velocity_T] ...bulk particle velocity in init_cond
\param [Nr_Particles_In_Cell_T] ...number of particle / cell to create in init_cond
\param [Particles_T] ...variable to assign the created particles to in given cells
\param [Particle] Type of particles stored in a container accessed by Particles_T
\param [Particle_Mass_T] Used to access the mass data of a particle
\param [Particle_Charge_Mass_Ratio_T] ...charge to mass ratio data of a particle
\param [Particle_Position_T] ...position data of a particle
\param [Particle_Velocity_T] ...velocity data of a particle
\param [Particle_ID_T] ...id data of a particle
\param [init_cond] Initial condition for particles obtained from user
\param [simulation_time] Time to give to init_cond when querying data
\param [cells] Cells into which to create particles
\param [random_source] Source to use when creating particle distributions
\param [particle_temp_nrj_ratio] https://en.wikipedia.org/wiki/Boltzmann_constant
\param [replace] Whether to replace/add particles in/to given cells
\param [verbose] Whether to print what is being done to stdout on MPI rank 0
*/
template<
	class Number_Density_T,
	class Temperature_T,
	class Charge_Mass_Ratio_T,
	class Bulk_Velocity_T,
	class Nr_Particles_In_Cell_T,
	class Particles_T,
	class Particle,
	class Particle_Mass_T,
	class Particle_Charge_Mass_Ratio_T,
	class Particle_Position_T,
	class Particle_Velocity_T,
	class Particle_ID_T,
	class Particle_Species_Mass_T,
	class Init_Cond,
	class Grid
> size_t initialize(
	Init_Cond& init_cond,
	const double simulation_time,
	const std::vector<uint64_t>& cells,
	Grid& grid,
	std::mt19937_64& random_source,
	const double particle_temp_nrj_ratio,
	const unsigned long long int first_particle_id,
	const unsigned long long int particle_id_increase,
	const bool replace,
	const bool verbose
) {
	if (verbose && grid.get_rank() == 0) {
		std::cout << "Setting default particle state... ";
		std::cout.flush();
	}

	size_t nr_particles_created = 0;
	auto current_id_start = first_particle_id;
	for (const auto cell_id: cells) {
		random_source.seed(cell_id);

		const auto
			cell_start = grid.geometry.get_min(cell_id),
			cell_end = grid.geometry.get_max(cell_id),
			cell_length = grid.geometry.get_length(cell_id),
			cell_center = grid.geometry.get_center(cell_id);

		// classify cells for setting non-default initial state
		init_cond.add_cell(
			cell_id,
			cell_start,
			cell_end
		);

		auto* const cell_data = grid[cell_id];
		if (cell_data == NULL) {
			std::cerr <<  __FILE__ << "(" << __LINE__ << ") No data for cell: "
				<< cell_id
				<< std::endl;
			abort();
		}

		// set default state
		const auto
			number_density
				= init_cond.default_data.get_data(
					Number_Density_T(),
					cell_center,
					simulation_time
				),
			temperature
				= init_cond.default_data.get_data(
					Temperature_T(),
					cell_center,
					simulation_time
				),
			charge_mass_ratio
				= init_cond.default_data.get_data(
					Charge_Mass_Ratio_T(),
					cell_center,
					simulation_time
				),
			species_mass
				= init_cond.default_data.get_data(
					Particle_Species_Mass_T(),
					cell_center,
					simulation_time
				);
		const auto bulk_velocity
			= init_cond.default_data.get_data(
				Bulk_Velocity_T(),
				cell_center,
				simulation_time
			);
		const auto nr_particles
			= init_cond.default_data.get_data(
				Nr_Particles_In_Cell_T(),
				cell_center,
				simulation_time
			);

		auto new_particles
			= create_particles<
				Particle,
				Particle_Mass_T,
				Particle_Charge_Mass_Ratio_T,
				Particle_Position_T,
				Particle_Velocity_T,
				Particle_ID_T,
				Particle_Species_Mass_T
			>(
				bulk_velocity,
				Eigen::Vector3d{cell_start[0], cell_start[1], cell_start[2]},
				Eigen::Vector3d{cell_end[0], cell_end[1], cell_end[2]},
				Eigen::Vector3d{temperature, temperature, temperature},
				nr_particles,
				charge_mass_ratio,
				species_mass * number_density * cell_length[0] * cell_length[1] * cell_length[2],
				species_mass,
				particle_temp_nrj_ratio,
				random_source,
				current_id_start,
				particle_id_increase
			);
		nr_particles_created += nr_particles;

		if (replace) {
			(*cell_data)[Particles_T()] = std::move(new_particles);
		} else {
			(*cell_data)[Particles_T()].insert(
				(*cell_data)[Particles_T()].end(),
				new_particles.begin(),
				new_particles.end()
			);
		}

		current_id_start += nr_particles * particle_id_increase;
	}

	// set non-default initial conditions
	if (verbose && grid.get_rank() == 0) {
		std::cout << "done\nSetting non-default initial particle state... ";
		std::cout.flush();
	}
	for (size_t bdy_id = 0; bdy_id < init_cond.get_number_of_boundaries(); bdy_id++) {
		for (const auto& cell_id: init_cond.get_cells(bdy_id)) {
			const auto
				cell_start = grid.geometry.get_min(cell_id),
				cell_end = grid.geometry.get_max(cell_id),
				cell_length = grid.geometry.get_length(cell_id),
				cell_center = grid.geometry.get_center(cell_id);

			const auto
				number_density
					= init_cond.get_data(
						Number_Density_T(),
						bdy_id,
						cell_center,
						simulation_time
					),
				temperature
					= init_cond.get_data(
						Temperature_T(),
						bdy_id,
						cell_center,
						simulation_time
					),
				charge_mass_ratio
					= init_cond.get_data(
						Charge_Mass_Ratio_T(),
						bdy_id,
						cell_center,
						simulation_time
					),
				species_mass
					= init_cond.get_data(
						Particle_Species_Mass_T(),
						bdy_id,
						cell_center,
						simulation_time
					);
			const auto bulk_velocity
				= init_cond.get_data(
					Bulk_Velocity_T(),
					bdy_id,
					cell_center,
					simulation_time
				);
			const auto nr_particles
				= init_cond.get_data(
					Nr_Particles_In_Cell_T(),
					bdy_id,
					cell_center,
					simulation_time
				);

			auto* const cell_data = grid[cell_id];
			if (cell_data == NULL) {
				std::cerr <<  __FILE__ << "(" << __LINE__ << ") No data for cell: "
					<< cell_id
					<< std::endl;
				abort();
			}

			auto new_particles
				= create_particles<
					Particle,
					Particle_Mass_T,
					Particle_Charge_Mass_Ratio_T,
					Particle_Position_T,
					Particle_Velocity_T,
					Particle_ID_T,
					Particle_Species_Mass_T
				>(
					bulk_velocity,
					Eigen::Vector3d{cell_start[0], cell_start[1], cell_start[2]},
					Eigen::Vector3d{cell_end[0], cell_end[1], cell_end[2]},
					Eigen::Vector3d{temperature, temperature, temperature},
					nr_particles,
					charge_mass_ratio,
					species_mass * number_density * cell_length[0] * cell_length[1] * cell_length[2],
					species_mass,
					particle_temp_nrj_ratio,
					random_source,
					current_id_start,
					particle_id_increase
				);
			nr_particles_created += nr_particles;

			if (replace) {
				(*cell_data)[Particles_T()] = std::move(new_particles);
			} else {
				(*cell_data)[Particles_T()].insert(
					(*cell_data)[Particles_T()].end(),
					new_particles.begin(),
					new_particles.end()
				);
			}

			current_id_start += nr_particles * particle_id_increase;
		}
	}
	if (verbose && grid.get_rank() == 0) {
		std::cout << "done" << std::endl;
	}

	return nr_particles_created;
}

}} // namespaces

#endif // ifndef PAMHD_PARTICLE_INITIALIZE_HPP
