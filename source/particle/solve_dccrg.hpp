/*
Particle propagator of PAMHD for DCCRG grid.

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

#ifndef PAMHD_PARTICLE_SOLVE_DCCRG_HPP
#define PAMHD_PARTICLE_SOLVE_DCCRG_HPP


#include "type_traits"
#include "utility"

#include "dccrg.hpp"
#include "Eigen/Core"
#include "Eigen/Geometry"

#include "interpolate.hpp"


namespace pamhd {
namespace particle {


/*!
State type used in boost::numeric::odeint particle propagation.
state[0] = particle position, state[1] = velocity
*/
using state_t = std::array<Eigen::Vector3d, 2>;


/*!
Object given to boost::numeric::odeint to propagate a particle.

Interpolates electric and magnetic fields from given values to
particle position at each call to operator().
*/
class Particle_Propagator
{
private:

	const double& charge_to_mass_ratio;

	const Eigen::Vector3d
		&data_start, &data_end;

	const std::array<Eigen::Vector3d, 27>
		&electric_field, &magnetic_field;


public:

	/*!
	Arguments except charge to mass ratio are passed to interpolate().
	*/
	Particle_Propagator(
		const double& given_charge_to_mass_ratio,
		const Eigen::Vector3d& given_data_start,
		const Eigen::Vector3d& given_data_end,
		const std::array<Eigen::Vector3d, 27>& given_electric_field,
		const std::array<Eigen::Vector3d, 27>& given_magnetic_field
	) :
		charge_to_mass_ratio(given_charge_to_mass_ratio),
		data_start(given_data_start),
		data_end(given_data_end),
		electric_field(given_electric_field),
		magnetic_field(given_magnetic_field)
	{}

	void operator()(
		const state_t& state,
		state_t& change,
		const double
	) const {
		const auto
			E_at_pos
				= interpolate(
					state[0],
					this->data_start,
					this->data_end,
					this->electric_field
				),
			B_at_pos
				= interpolate(
					state[0],
					this->data_start,
					this->data_end,
					this->magnetic_field
				);

		change[0] = state[1];
		change[1]
			= this->charge_to_mass_ratio
			* (E_at_pos + state[1].cross(B_at_pos));
	};
};


/*!
Propagates particles in given cells for a given amount of time.

Returns the longest allowed time step for given cells
and their neighbors. Particles which propagate outside of the
cell in which they are stored are moved to the External_Particles_T
list of their previous cell and added to Particle_Destinations_T
information.
*/
template<
	class Cell_T,
	class Electric_Field_T,
	class Magnetic_Field_T,
	class Nr_Particles_External_T,
	class Particles_Internal_T,
	class Particles_External_T,
	class Particle_Position_T,
	class Particle_Velocity_T,
	class Charge_To_Mass_Ratio_T,
	class Particle_Mass_T,
	class Particle_Destination_T,
	class Stepper
> double solve(
	const double dt,
	const std::vector<uint64_t>& cell_ids,
	dccrg::Dccrg<Cell_T, dccrg::Cartesian_Geometry>& grid
) {
	namespace odeint = boost::numeric::odeint;
	using std::isnan;
	using std::is_same;

	static_assert(
		is_same<Stepper, odeint::euler<state_t>>::value
		or is_same<Stepper, odeint::modified_midpoint<state_t>>::value
		or is_same<Stepper, odeint::runge_kutta4<state_t>>::value
		or is_same<Stepper, odeint::runge_kutta_cash_karp54<state_t>>::value
		or is_same<Stepper, odeint::runge_kutta_fehlberg78<state_t>>::value,
		"Only odeint steppers without internal state are supported."
	);

	constexpr Electric_Field_T Ele{};
	constexpr Magnetic_Field_T Mag{};
	constexpr Nr_Particles_External_T Nr_Ext{};
	constexpr Particles_Internal_T Part_Int{};
	constexpr Particles_External_T Part_Ext{};
	constexpr Particle_Position_T Pos{};
	constexpr Particle_Velocity_T Vel{};
	constexpr Charge_To_Mass_Ratio_T C2M{};
	constexpr Particle_Mass_T Mas{};
	constexpr Particle_Destination_T Des{};


	Stepper stepper;
	double max_time_step = std::numeric_limits<double>::max();
	for (const auto& cell_id: cell_ids) {

		// get field data from neighborhood for interpolation
		std::array<Eigen::Vector3d, 27> electric_fields, magnetic_fields;

		const auto* const neighbor_ids = grid.get_neighbors_of(cell_id);
		if (neighbor_ids == nullptr) {
			std::cerr << __FILE__ << "(" << __LINE__ << "): " << std::endl;
			abort();
		}

		if (neighbor_ids->size() != 26) {
			std::cerr << __FILE__ << "(" << __LINE__ << "): "
				<< "Unsupported neighborhood: " << neighbor_ids->size()
				<< std::endl;
			abort();
		}

		auto* const cell_data = grid[cell_id];
		if (cell_data == NULL) {
			std::cerr << __FILE__ << "(" << __LINE__ << "): " << std::endl;
			abort();
		}

		// put current cell's data into middle
		electric_fields[14] = (*cell_data)[Ele];
		magnetic_fields[14] = (*cell_data)[Mag];

		for (size_t i = 0; i < neighbor_ids->size(); i++) {
			const auto neighbor_id = (*neighbor_ids)[i];

			// use same values in missing neighbors
			if (neighbor_id == dccrg::error_cell) {
				if (i <= 13) {
					electric_fields[i] = electric_fields[14];
					magnetic_fields[i] = magnetic_fields[14];
				} else {
					electric_fields[i + 1] = electric_fields[14];
					magnetic_fields[i + 1] = magnetic_fields[14];
				}

				continue;
			}

			const auto* const neighbor_data = grid[neighbor_id];
			if (neighbor_data == nullptr) {
				std::cerr << __FILE__ << "(" << __LINE__ << "): " << std::endl;
				abort();
			}

			if (i <= 13) {
				electric_fields[i] = (*neighbor_data)[Ele];
				magnetic_fields[i] = (*neighbor_data)[Mag];
			} else {
				electric_fields[i + 1] = (*neighbor_data)[Ele];
				magnetic_fields[i + 1] = (*neighbor_data)[Mag];
			}
		}

		const auto
			cell_min = grid.geometry.get_min(cell_id),
			cell_max = grid.geometry.get_max(cell_id),
			cell_center = grid.geometry.get_center(cell_id),
			cell_length = grid.geometry.get_length(cell_id);

		const Eigen::Vector3d
			interpolation_start{
				cell_center[0] - cell_length[0],
				cell_center[1] - cell_length[1],
				cell_center[2] - cell_length[2]
			},
			interpolation_end{
				cell_center[0] + cell_length[0],
				cell_center[1] + cell_length[1],
				cell_center[2] + cell_length[2]
			};

		for (size_t i = 0; i < (*cell_data)[Part_Int].size(); i++) {
			const auto& particle = (*cell_data)[Part_Int][i];

			// check max length of time step for next step
			max_time_step =
				std::min(max_time_step,
				std::min(fabs(cell_length[0] / particle[Vel][0]),
				std::min(fabs(cell_length[1] / particle[Vel][1]),
				         fabs(cell_length[2] / particle[Vel][2])
				)));

			const Particle_Propagator propagator(
				particle[C2M],
				interpolation_start,
				interpolation_end,
				electric_fields,
				magnetic_fields
			);

			state_t state{{particle[Pos], particle[Vel]}};
			stepper.do_step(propagator, state, 0.0, dt);
			(*cell_data)[Part_Int][i][Vel] = state[1];


			// take into account periodic grid
			const std::array<double, 3> real_pos
				= grid.geometry.get_real_coordinate({{
					state[0][0],
					state[0][1],
					state[0][2]
				}});

			// remove from simulation if particle not inside of grid
			if (
				isnan(real_pos[0])
				or isnan(real_pos[1])
				or isnan(real_pos[2])
			) {

				(*cell_data)[Part_Int].erase((*cell_data)[Part_Int].begin() + i);
				i--;

			// move to ext list if particle outside of current cell
			} else if (
				real_pos[0] < cell_min[0]
				or real_pos[0] > cell_max[0]
				or real_pos[1] < cell_min[1]
				or real_pos[1] > cell_max[1]
				or real_pos[2] < cell_min[2]
				or real_pos[2] > cell_max[2]
			) {
				uint64_t destination = dccrg::error_cell;

				for (const auto& neighbor_id: *neighbor_ids) {
					if (neighbor_id == dccrg::error_cell) {
						continue;
					}

					const auto
						neighbor_min = grid.geometry.get_min(neighbor_id),
						neighbor_max = grid.geometry.get_max(neighbor_id);

					if (
						real_pos[0] >= neighbor_min[0]
						and real_pos[0] <= neighbor_max[0]
						and real_pos[1] >= neighbor_min[1]
						and real_pos[1] <= neighbor_max[1]
						and real_pos[2] >= neighbor_min[2]
						and real_pos[2] <= neighbor_max[2]
					) {
						destination = neighbor_id;
						break;
					}
				}

				if (destination != dccrg::error_cell) {
					const auto index = (*cell_data)[Part_Ext].size();

					(*cell_data)[Part_Ext].resize(index + 1);
					assign(
						(*cell_data)[Part_Ext][index],
						(*cell_data)[Part_Int][i]
					);
					(*cell_data)[Part_Ext][index][Des] = destination;
					(*cell_data)[Part_Ext][index][Pos] = {
						real_pos[0],
						real_pos[1],
						real_pos[2]
					};

					(*cell_data)[Part_Int].erase((*cell_data)[Part_Int].begin() + i);
					i--;
				} else {
					std::cerr << __FILE__ << "(" << __LINE__ << ")" << std::endl;
					abort();
				}

			// update particle position in this cell's internal list
			} else {
				(*cell_data)[Part_Int][i][Pos] = {
					real_pos[0],
					real_pos[1],
					real_pos[2]
				};
			}
		}

		(*cell_data)[Nr_Ext] = (*cell_data)[Part_Ext].size();
	}

	// only allow particles to propagate through half a
	// cell in order to not break field interpolation
	return max_time_step / 2.0;
}


template<
	class Cell_T,
	class Nr_Particles_External_T,
	class Particles_External_T
> void resize_receiving_containers(
	const std::vector<uint64_t>& cell_ids,
	dccrg::Dccrg<Cell_T, dccrg::Cartesian_Geometry>& grid
) {
	for (const auto& cell_id: cell_ids) {

		auto* const cell_ptr = grid[cell_id];
		if (cell_ptr == NULL) {
			std::cerr << __FILE__ << "(" << __LINE__ << ")" << std::endl;
			abort();
		}
		auto& cell_data = *cell_ptr;

		cell_data[Particles_External_T()]
			.resize(cell_data[Nr_Particles_External_T()]);
	}
}


template<
	class Cell_T,
	class Nr_Particles_External_T,
	class Particles_Internal_T,
	class Particles_External_T,
	class Particle_Destination_T
> void incorporate_external_particles(
	const std::vector<uint64_t>& cell_ids,
	dccrg::Dccrg<Cell_T, dccrg::Cartesian_Geometry>& grid
) {
	constexpr Nr_Particles_External_T Nr_Ext{};
	constexpr Particles_Internal_T Part_Int{};
	constexpr Particles_External_T Part_Ext{};
	constexpr Particle_Destination_T Dest{};

	for (const auto& cell_id: cell_ids) {

		auto* const cell_ptr = grid[cell_id];
		if (cell_ptr == nullptr) {
			std::cerr << __FILE__ << "(" << __LINE__ << ")" << std::endl;
			abort();
		}
		auto& cell_data = *cell_ptr;

		// assign particles from neighbors' external list to this cell
		const auto* const neighbor_ids = grid.get_neighbors_of(cell_id);
		if (neighbor_ids == nullptr) {
			std::cerr << __FILE__ << "(" << __LINE__ << ")" << std::endl;
			abort();
		}

		for (const auto& neighbor_id: *neighbor_ids) {
			if (neighbor_id == dccrg::error_cell) {
				continue;
			}

			auto* const neighbor_ptr = grid[neighbor_id];
			if (neighbor_ptr == nullptr) {
				std::cerr << __FILE__ << "(" << __LINE__ << ")" << std::endl;
				abort();
			}
			auto& neighbor_data = *neighbor_ptr;

			for (auto& particle: neighbor_data[Part_Ext]) {
				if (particle[Dest] == dccrg::error_cell) {
					continue;
				}

				if (particle[Dest] == cell_id) {
					particle[Dest] = dccrg::error_cell;

					const auto index = cell_data[Part_Int].size();

					cell_data[Part_Int].resize(index + 1);

					assign(cell_data[Part_Int][index], particle);
				}
			}
		}
	}
}


template<
	class Cell_T,
	class Nr_Particles_External_T,
	class Particles_External_T
> void remove_external_particles(
	const std::vector<uint64_t>& cell_ids,
	dccrg::Dccrg<Cell_T, dccrg::Cartesian_Geometry>& grid
) {
	for (auto cell_id: cell_ids) {

		auto* const cell_data = grid[cell_id];
		if (cell_data == nullptr) {
			std::cerr << __FILE__ << "(" << __LINE__ << ")" << std::endl;
			abort();
		}

		(*cell_data)[Nr_Particles_External_T()] = 0;
		(*cell_data)[Particles_External_T()].clear();
	}
}


}} // namespaces

#endif // ifndef PAMHD_PARTICLE_SOLVE_DCCRG_HPP
