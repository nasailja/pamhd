/*
Particle data accumulator of PAMHD built on top of DCCRG.

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

#ifndef PAMHD_PARTICLE_ACCUMULATE_DCCRG_HPP
#define PAMHD_PARTICLE_ACCUMULATE_DCCRG_HPP


#include "iostream"
#include "vector"

#include "dccrg.hpp"
#include "dccrg_cartesian_geometry.hpp"
#include "Eigen/Core"

#include "particle/accumulate.hpp"
#include "particle/common.hpp"
#include "particle/variables.hpp"


namespace pamhd {
namespace particle {


/*!
Accumulates particle data in given cells to those cells and their neighbors.

Accumulated data from particles in the same cell is not zeroed before accumulating.

Accumulated data is written to local cells directly, accumulated data to
remote neighbors is stored in local cells' accumulation list.

Updates the number of remote accumulated values.
*/
template<
	class Particles_Getter,
	class Particle_Position_Getter,
	class Particle_Value_Getter,
	class Bulk_Value_Getter,
	class Bulk_Value_In_List_Getter,
	class Target_In_List_Getter,
	class Accumulation_List_Length_Getter,
	class Accumulation_List_Getter,
	class Cell
> void accumulate(
	const std::vector<uint64_t>& cell_ids,
	dccrg::Dccrg<Cell, dccrg::Cartesian_Geometry>& grid,
	Particles_Getter Part,
	Particle_Position_Getter Part_Pos,
	Particle_Value_Getter Part_Val,
	Bulk_Value_Getter Bulk_Val,
	Bulk_Value_In_List_Getter List_Bulk_Val,
	Target_In_List_Getter List_Target,
	Accumulation_List_Length_Getter List_Len,
	Accumulation_List_Getter Accu_List,
	const bool clear_at_start = true
) {
	for (const auto& cell_id: cell_ids) {
		const auto
			cell_min_tmp = grid.geometry.get_min(cell_id),
			cell_max_tmp = grid.geometry.get_max(cell_id),
			cell_length_tmp = grid.geometry.get_length(cell_id),
			cell_center_tmp = grid.geometry.get_center(cell_id);
		const Eigen::Vector3d
			cell_min{cell_min_tmp[0], cell_min_tmp[1], cell_min_tmp[2]},
			cell_max{cell_max_tmp[0], cell_max_tmp[1], cell_max_tmp[2]},
			cell_length{cell_length_tmp[0], cell_length_tmp[1], cell_length_tmp[2]},
			cell_center{cell_center_tmp[0], cell_center_tmp[1], cell_center_tmp[2]};

		auto* const cell_data = grid[cell_id];
		if (cell_data == nullptr) {
			std::cerr << __FILE__ << "(" << __LINE__ << ")" << std::endl;
			abort();
		}

		if (clear_at_start) {
			Accu_List(*cell_data).clear();
		}

		// cache required neighbor data
		std::vector<bool> is_locals;
		std::vector<uint64_t> neighbor_ids;
		std::vector<decltype(grid[cell_id])> neighbor_datas;
		std::vector<Eigen::Vector3d> neighbor_mins, neighbor_maxs;

		for (const auto& x_offset: {-1, 0, 1}) {
		for (const auto& y_offset: {-1, 0, 1}) {
		for (const auto& z_offset: {-1, 0, 1}) {
			const auto neighbors
				= grid.get_neighbors_of_at_offset(cell_id, x_offset, y_offset, z_offset);
			if (
				neighbors.size() == 0
				and (x_offset != 0 or y_offset != 0 or z_offset != 0)
			) {
				std::cerr << __FILE__ << "(" << __LINE__ << ")" << std::endl;
				abort();
			}

			for (const auto& neighbor_id: neighbors) {
				if (neighbor_id == dccrg::error_cell) {
					continue;
				}

				if (
					grid.get_refinement_level(neighbor_id) != grid.get_refinement_level(cell_id)
				) {
					std::cerr << __FILE__ << "(" << __LINE__ << ") "
						<< "Different cell refinement levels not supported"
						<< std::endl;
					abort();
				}

				neighbor_ids.emplace_back(neighbor_id);

				auto* const neighbor_data = grid[neighbor_id];
				if (neighbor_data == nullptr) {
					std::cerr << __FILE__ << "(" << __LINE__ << "): "
						<< "No data for neighbor " << neighbor_id
						<< " of cell " << cell_id
						<< std::endl;
					abort();
				}
				neighbor_datas.emplace_back(neighbor_data);

				if (grid.is_local(neighbor_id)) {
					is_locals.push_back(true);
				} else {
					is_locals.push_back(false);
				}

				const auto
					neighbor_center = grid.geometry.get_center(neighbor_id),
					neighbor_length = grid.geometry.get_length(neighbor_id);

				Eigen::Vector3d
					neighbor_min{
						neighbor_center[0] - neighbor_length[0] / 2,
						neighbor_center[1] - neighbor_length[1] / 2,
						neighbor_center[2] - neighbor_length[2] / 2
					},
					neighbor_max{
						neighbor_center[0] + neighbor_length[0] / 2,
						neighbor_center[1] + neighbor_length[1] / 2,
						neighbor_center[2] + neighbor_length[2] / 2
					};

				// handle periodic grid geometry
				if (x_offset < 0 and neighbor_center[0] > cell_min[0]) {
					neighbor_max[0] = cell_min[0];
					neighbor_min[0] = cell_min[0] - neighbor_length[0];
				} else if (x_offset > 0 and neighbor_center[0] < cell_max[0]) {
					neighbor_min[0] = cell_max[0];
					neighbor_max[0] = cell_max[0] + neighbor_length[0];
				}

				if (y_offset < 0 and neighbor_center[1] > cell_min[1]) {
					neighbor_max[1] = cell_min[1];
					neighbor_min[1] = cell_min[1] - neighbor_length[1];
				} else if (y_offset > 0 and neighbor_center[1] < cell_max[1]) {
					neighbor_min[1] = cell_max[1];
					neighbor_max[1] = cell_max[1] + neighbor_length[1];
				}

				if (z_offset < 0 and neighbor_center[2] > cell_min[2]) {
					neighbor_max[2] = cell_min[2];
					neighbor_min[2] = cell_min[2] - neighbor_length[2];
				} else if (z_offset > 0 and neighbor_center[2] < cell_max[2]) {
					neighbor_min[2] = cell_max[2];
					neighbor_max[2] = cell_max[2] + neighbor_length[2];
				}

				neighbor_mins.push_back(neighbor_min);
				neighbor_maxs.push_back(neighbor_max);
			}
		}}}

		for (auto& particle: Part(*cell_data)) {
			auto& position = Part_Pos(particle);
			const Eigen::Vector3d
				value_box_min{
					position[0] - cell_length[0] / 2,
					position[1] - cell_length[1] / 2,
					position[2] - cell_length[2] / 2
				},
				value_box_max{
					position[0] + cell_length[0] / 2,
					position[1] + cell_length[1] / 2,
					position[2] + cell_length[2] / 2
				};

			// accumulate to current cell
			Bulk_Val(*cell_data)
				+= get_accumulated_value(
					Part_Val(*cell_data, particle),
					value_box_min,
					value_box_max,
					cell_min,
					cell_max
				);

			// accumulate to neighbors
			for (size_t i = 0; i < is_locals.size(); i++) {
				const auto accumulated_value
					= get_accumulated_value(
						Part_Val(*(neighbor_datas[i]), particle),
						value_box_min,
						value_box_max,
						neighbor_mins[i],
						neighbor_maxs[i]
					);

				// same as for current cell above
				if (is_locals[i]) {
					Bulk_Val(*(neighbor_datas[i])) += accumulated_value;
				// accumulate values to a list in current cell
				} else {
					// use this index in the list
					size_t accumulation_index = 0;

					// find the index of target neighbor
					const auto neighbor_id = neighbor_ids[i];
					auto iter
						= std::find_if(
							Accu_List(*cell_data).begin(),
							Accu_List(*cell_data).end(),
							[&neighbor_id, &List_Target](
								const decltype(*Accu_List(*cell_data).begin()) candidate_item
							) {
								if (List_Target(candidate_item) == neighbor_id) {
									return true;
								} else {
									return false;
								}
							}
						);

					// found
					if (iter != Accu_List(*cell_data).end()) {
						List_Bulk_Val(*iter) += accumulated_value;
					// create the item
					} else {
						const auto old_size = Accu_List(*cell_data).size();
						Accu_List(*cell_data).resize(old_size + 1);
						auto& new_item = Accu_List(*cell_data)[old_size];
						List_Target(new_item) = neighbor_id;
						List_Bulk_Val(new_item) = accumulated_value;
					}
				}
			}
		}

		List_Len(*cell_data) = Accu_List(*cell_data).size();
	}
}


/*!
Allocates memory for lists of accumulated particle data in copies of remote neighbors of local cells.
*/
template<
	class Grid,
	class Accumulation_List_Getter,
	class Accumulation_List_Length_Getter
> void allocate_accumulation_lists(
	Grid& grid,
	Accumulation_List_Getter Accu_List,
	Accumulation_List_Length_Getter List_Len
) {
	for (const auto& remote_cell_id: grid.get_remote_cells_on_process_boundary()) {
		auto* const cell_data = grid[remote_cell_id];
		if (cell_data == nullptr) {
			std::cerr << __FILE__ << "(" << __LINE__ << ")" << std::endl;
			abort();
		}

		Accu_List(*cell_data).resize(List_Len(*cell_data));
	}
}


/*!
Adds accumulated particle data from remote neighbors' to local cells.
*/
template<
	class Grid,
	class Bulk_Value_Getter,
	class Value_In_List_Getter,
	class Target_In_List_Getter,
	class Accumulation_List_Getter
> void accumulate_from_remote_neighbors(
	Grid& grid,
	Bulk_Value_Getter Bulk_Val,
	Value_In_List_Getter Value_In_List,
	Target_In_List_Getter Target_In_List,
	Accumulation_List_Getter Accu_List
) {
	for (const auto& remote_cell_id: grid.get_remote_cells_on_process_boundary()) {
		auto* const source_data = grid[remote_cell_id];
		if (source_data == nullptr) {
			std::cerr << __FILE__ << "(" << __LINE__ << ")" << std::endl;
			abort();
		}

		for (auto& item: Accu_List(*source_data)) {
			if (not grid.is_local(Target_In_List(item))) {
				continue;
			}

			auto* const target_data = grid[Target_In_List(item)];
			if (target_data == nullptr) {
				std::cerr << __FILE__ << "(" << __LINE__ << ")" << std::endl;
				abort();
			}

			Bulk_Val(*target_data) += Value_In_List(item);
		}
	}
}


template<
	class Cell,
	class Particles_Getter,
	class Particle_Position_Getter,
	class Particle_Mass_Getter,
	class Particle_Species_Mass_Getter,
	class Particle_Momentum_Getter,
	class Particle_Relative_Velocity2_Getter,
	class Number_Of_Particles_Getter,
	class Bulk_Mass_Getter,
	class Bulk_Momentum_Getter,
	class Bulk_Relative_Velocity2_Getter,
	class Number_Of_Particles_In_List_Getter,
	class Bulk_Mass_In_List_Getter,
	class Bulk_Momentum_In_List_Getter,
	class Bulk_Relative_Velocity2_In_List_Getter,
	class Bulk_Velocity_Getter,
	class Target_In_List_Getter,
	class Accumulation_List_Length_Getter,
	class Accumulation_List_Getter,
	class Accumulation_List_Length_Variable,
	class Accumulation_List_Variable,
	class Bulk_Velocity_Variable
> void accumulate_mhd_data(
	const std::vector<uint64_t>& inner_cell_ids,
	const std::vector<uint64_t>& outer_cell_ids,
	dccrg::Dccrg<Cell, dccrg::Cartesian_Geometry>& grid,
	Particles_Getter Particles,
	Particle_Position_Getter Particle_Position,
	Particle_Mass_Getter Particle_Mass,
	Particle_Species_Mass_Getter Particle_Species_Mass,
	Particle_Momentum_Getter Particle_Momentum,
	Particle_Relative_Velocity2_Getter Particle_Relative_Velocity2,
	Number_Of_Particles_Getter Number_Of_Particles,
	Bulk_Mass_Getter Bulk_Mass,
	Bulk_Momentum_Getter Bulk_Momentum,
	Bulk_Relative_Velocity2_Getter Bulk_Relative_Velocity2,
	Bulk_Velocity_Getter Bulk_Velocity,
	Number_Of_Particles_In_List_Getter Accu_List_Number_Of_Particles,
	Bulk_Mass_In_List_Getter Accu_List_Bulk_Mass,
	Bulk_Momentum_In_List_Getter Accu_List_Bulk_Momentum,
	Bulk_Relative_Velocity2_In_List_Getter Accu_List_Bulk_Relative_Velocity2,
	Target_In_List_Getter Accu_List_Target,
	Accumulation_List_Length_Getter Accu_List_Length,
	Accumulation_List_Getter Accu_List,
	Accumulation_List_Length_Variable accu_list_len_var,
	Accumulation_List_Variable accu_list_var,
	Bulk_Velocity_Variable bulk_vel_var
) {
	auto cell_ids = inner_cell_ids;
	cell_ids.insert(cell_ids.end(), outer_cell_ids.cbegin(), outer_cell_ids.cend());

	for (const auto& cell: cell_ids) {
		auto* const cell_data = grid[cell];
		if (cell_data == nullptr) {
			std::cerr <<  __FILE__ << "(" << __LINE__ << ")" << std::endl;
			abort();
		}
		Bulk_Mass(*cell_data) = 0;
		Bulk_Momentum(*cell_data) = {0, 0, 0};
	}

	accumulate(
		outer_cell_ids,
		grid,
		Particles,
		Particle_Position,
		Particle_Mass,
		Bulk_Mass,
		Accu_List_Bulk_Mass,
		Accu_List_Target,
		Accu_List_Length,
		Accu_List
	);
	accumulate(
		outer_cell_ids,
		grid,
		Particles,
		Particle_Position,
		Particle_Momentum,
		Bulk_Momentum,
		Accu_List_Bulk_Momentum,
		Accu_List_Target,
		Accu_List_Length,
		Accu_List,
		false // keep previous data in accumulation lists
	);

	Cell::set_transfer_all(true, accu_list_len_var);
	grid.start_remote_neighbor_copy_updates();

	accumulate(
		inner_cell_ids,
		grid,
		Particles,
		Particle_Position,
		Particle_Mass,
		Bulk_Mass,
		Accu_List_Bulk_Mass,
		Accu_List_Target,
		Accu_List_Length,
		Accu_List
	);
	accumulate(
		inner_cell_ids,
		grid,
		Particles,
		Particle_Position,
		Particle_Momentum,
		Bulk_Momentum,
		Accu_List_Bulk_Momentum,
		Accu_List_Target,
		Accu_List_Length,
		Accu_List
	);

	grid.wait_remote_neighbor_copy_update_receives();

	allocate_accumulation_lists(
		grid,
		Accu_List,
		Accu_List_Length
	);

	grid.wait_remote_neighbor_copy_update_sends();
	Cell::set_transfer_all(false, accu_list_len_var);

	Cell::set_transfer_all(true, accu_list_var);
	grid.start_remote_neighbor_copy_updates();
	grid.wait_remote_neighbor_copy_update_receives();

	accumulate_from_remote_neighbors(
		grid,
		Bulk_Mass,
		Accu_List_Bulk_Mass,
		Accu_List_Target,
		Accu_List
	);
	accumulate_from_remote_neighbors(
		grid,
		Bulk_Momentum,
		Accu_List_Bulk_Momentum,
		Accu_List_Target,
		Accu_List
	);

	grid.wait_remote_neighbor_copy_update_sends();
	Cell::set_transfer_all(false, accu_list_var);

	/*
	Accumulate particle data required for pressure
	*/

	// needs remote neighbors' bulk velocity
	for (const auto& cell: outer_cell_ids) {
		auto* const cell_data = grid[cell];
		if (cell_data == nullptr) {
			std::cerr <<  __FILE__ << "(" << __LINE__ << ")" << std::endl;
			abort();
		}
		Number_Of_Particles(*cell_data)     =
		Bulk_Relative_Velocity2(*cell_data) = 0;
		Bulk_Velocity(*cell_data) = Bulk_Momentum(*cell_data) / Bulk_Mass(*cell_data);
	}

	Cell::set_transfer_all(true, bulk_vel_var);
	grid.start_remote_neighbor_copy_updates();

	for (const auto& cell: inner_cell_ids) {
		auto* const cell_data = grid[cell];
		if (cell_data == nullptr) {
			std::cerr <<  __FILE__ << "(" << __LINE__ << ")" << std::endl;
			abort();
		}
		Number_Of_Particles(*cell_data)     =
		Bulk_Relative_Velocity2(*cell_data) = 0;
		Bulk_Velocity(*cell_data) = Bulk_Momentum(*cell_data) / Bulk_Mass(*cell_data);
	}

	grid.wait_remote_neighbor_copy_update_receives();

	// returns number of particles of mass species mass
	const auto particle_counter
		= [
			&Particle_Mass,
			&Particle_Species_Mass
		](
			Cell& cell,
			decltype(*Particles(*(grid[0])).begin())& particle
		) ->double {
			return Particle_Mass(cell, particle) / Particle_Species_Mass(cell, particle);
		};
	accumulate(
		outer_cell_ids,
		grid,
		Particles,
		Particle_Position,
		particle_counter,
		Number_Of_Particles,
		Accu_List_Number_Of_Particles,
		Accu_List_Target,
		Accu_List_Length,
		Accu_List
	);
	accumulate(
		outer_cell_ids,
		grid,
		Particles,
		Particle_Position,
		Particle_Relative_Velocity2,
		Bulk_Relative_Velocity2,
		Accu_List_Bulk_Relative_Velocity2,
		Accu_List_Target,
		Accu_List_Length,
		Accu_List,
		false
	);

	grid.wait_remote_neighbor_copy_update_sends();
	Cell::set_transfer_all(false, bulk_vel_var);


	Cell::set_transfer_all(true, accu_list_len_var);
	grid.start_remote_neighbor_copy_updates();

	accumulate(
		inner_cell_ids,
		grid,
		Particles,
		Particle_Position,
		particle_counter,
		Number_Of_Particles,
		Accu_List_Number_Of_Particles,
		Accu_List_Target,
		Accu_List_Length,
		Accu_List
	);
	accumulate(
		inner_cell_ids,
		grid,
		Particles,
		Particle_Position,
		Particle_Relative_Velocity2,
		Bulk_Relative_Velocity2,
		Accu_List_Bulk_Relative_Velocity2,
		Accu_List_Target,
		Accu_List_Length,
		Accu_List
	);

	grid.wait_remote_neighbor_copy_update_receives();

	allocate_accumulation_lists(
		grid,
		Accu_List,
		Accu_List_Length
	);

	grid.wait_remote_neighbor_copy_update_sends();
	Cell::set_transfer_all(false, accu_list_len_var);

	Cell::set_transfer_all(true, accu_list_var);
	grid.start_remote_neighbor_copy_updates();
	grid.wait_remote_neighbor_copy_update_receives();

	accumulate_from_remote_neighbors(
		grid,
		Number_Of_Particles,
		Accu_List_Number_Of_Particles,
		Accu_List_Target,
		Accu_List
	);
	accumulate_from_remote_neighbors(
		grid,
		Bulk_Relative_Velocity2,
		Accu_List_Bulk_Relative_Velocity2,
		Accu_List_Target,
		Accu_List
	);

	grid.wait_remote_neighbor_copy_update_sends();
	Cell::set_transfer_all(false, accu_list_var);
}


/*!
Skips a cell if it doesn't have particles and no_particles_allowed == true.
*/
template<
	class Cell,
	class Geometry,
	class Number_Of_Particles_Getter,
	class Particle_Bulk_Mass_Getter,
	class Particle_Bulk_Momentum_Getter,
	class Particle_Bulk_Relative_Velocity2_Getter,
	class Particle_List_Getter,
	class Total_Mass_Getter,
	class MHD_Mass_Getter,
	class MHD_Momentum_Getter,
	class MHD_Energy_Getter,
	class MHD_Magnetic_Field_Getter
> void fill_mhd_fluid_values(
	const std::vector<uint64_t>& cell_ids,
	dccrg::Dccrg<Cell, Geometry>& grid,
	const double adiabatic_index,
	const double vacuum_permeability,
	const double particle_temp_nrj_ratio,
	Number_Of_Particles_Getter Number_Of_Particles,
	Particle_Bulk_Mass_Getter Particle_Bulk_Mass,
	Particle_Bulk_Momentum_Getter Particle_Bulk_Momentum,
	Particle_Bulk_Relative_Velocity2_Getter Particle_Bulk_Relative_Velocity2,
	Particle_List_Getter Particle_List,
	Total_Mass_Getter Total_Mass,
	MHD_Mass_Getter MHD_Mass,
	MHD_Momentum_Getter MHD_Momentum,
	MHD_Energy_Getter MHD_Energy,
	MHD_Magnetic_Field_Getter MHD_Magnetic_Field
) {
	for (const auto& cell: cell_ids) {
		auto* const cell_data = grid[cell];
		if (cell_data == nullptr) {
			std::cerr <<  __FILE__ << "(" << __LINE__ << ")" << std::endl;
			abort();
		}

		if (Particle_Bulk_Mass(*cell_data) <= 0) {
			MHD_Mass(*cell_data) = 0;
			MHD_Momentum(*cell_data) = {0, 0, 0};
			MHD_Energy(*cell_data) = 0;
			continue;
		}

		const auto length = grid.geometry.get_length(cell);
		const auto volume = length[0] * length[1] * length[2];

		MHD_Mass(*cell_data) = Particle_Bulk_Mass(*cell_data) / volume;
		MHD_Momentum(*cell_data) = Particle_Bulk_Momentum(*cell_data) / volume;

		const double
			pressure
				= Number_Of_Particles(*cell_data)
				* Particle_Bulk_Relative_Velocity2(*cell_data)
				/ Particle_Bulk_Mass(*cell_data)
				/ 3 / volume,
			mass_frac = MHD_Mass(*cell_data) / Total_Mass(*cell_data);

		MHD_Energy(*cell_data)
			= pressure / (adiabatic_index - 1)
			+ MHD_Momentum(*cell_data).squaredNorm() / MHD_Mass(*cell_data) / 2
			+ mass_frac * MHD_Magnetic_Field(*cell_data).squaredNorm() / vacuum_permeability / 2;
	}
}


}} // namespaces

#endif // ifndef PAMHD_PARTICLE_ACCUMULATE_DCCRG_HPP
