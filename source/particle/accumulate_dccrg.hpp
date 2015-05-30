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


#include "vector"

#include "dccrg.hpp"
#include "dccrg_cartesian_geometry.hpp"
#include "Eigen/Core"

#include "particle/accumulate.hpp"


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
	Accumulation_List_Getter Accu_List
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

		Accu_List(*cell_data).clear();

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

				if (grid.is_local(neighbor_id)) {
					auto* const neighbor_data = grid[neighbor_id];
					if (neighbor_data == nullptr) {
						std::cerr << __FILE__ << "(" << __LINE__ << ")" << std::endl;
						abort();
					}
					neighbor_datas.emplace_back(neighbor_data);
					is_locals.push_back(true);
				} else {
					neighbor_datas.emplace_back(nullptr);
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
					Part_Val(particle),
					value_box_min,
					value_box_max,
					cell_min,
					cell_max
				);

			// accumulate to neighbors
			for (size_t i = 0; i < is_locals.size(); i++) {
				const auto accumulated_value
					= get_accumulated_value(
						Part_Val(particle),
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
	class Accumulation_Target_T,
	class Bulk_Value_T,
	class Grid,
	class Bulk_Value_Getter,
	class Accumulation_List_Getter
> void accumulate_from_remote_neighbors(
	Grid& grid,
	Bulk_Value_Getter Bulk_Val,
	Accumulation_List_Getter Accu_List
) {
	for (const auto& remote_cell_id: grid.get_remote_cells_on_process_boundary()) {
		auto* const source_data = grid[remote_cell_id];
		if (source_data == nullptr) {
			std::cerr << __FILE__ << "(" << __LINE__ << ")" << std::endl;
			abort();
		}

		for (const auto& item: Accu_List(*source_data)) {
			const auto target = item[Accumulation_Target_T()];
			if (not grid.is_local(target)) {
				continue;
			}

			auto* const target_data = grid[target];
			if (target_data == nullptr) {
				std::cerr << __FILE__ << "(" << __LINE__ << ")" << std::endl;
				abort();
			}

			Bulk_Val(*target_data) += item[Bulk_Value_T()];
		}
	}
}


}} // namespaces

#endif // ifndef PAMHD_PARTICLE_ACCUMULATE_DCCRG_HPP
