/*
Handles boundary cell classification logic of MHD part of PAMHD.

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

#ifndef PAMHD_MHD_BOUNDARIES_HPP
#define PAMHD_MHD_BOUNDARIES_HPP


#include "cmath"
#include "limits"
#include "utility"
#include "vector"

#include "boost/optional.hpp"
#include "dccrg.hpp"

#include "mhd/N_solve.hpp"
#include "mhd/variables.hpp"


namespace pamhd {
namespace mhd {


/*!
Handles classification of boundary cells.

Boundary cells without non-boundary neighbors
are instead classified as dont_solve cells.

Assumes grid library provides a DCCRG compatible API.
*/
template<class Cell_Type_Variable> class Boundary_Classifier
{
public:

	//! cell types
	static constexpr typename Cell_Type_Variable::data_type
		normal_cell = 0,
		dont_solve_cell = 1,
		value_boundary_cell = 2,
		copy_boundary_cell = 3;

	/*!
	Source_Of_Copy_Cell_T variable stores the id of cell that
	is source of data for copy boundary cell.
	*/
	template<
		class Value_Boundary_Id_T,
		class Source_Of_Copy_Cell_T,
		class Cell_Data,
		class Geometry,
		class Value_Boundaries,
		class Copy_Boundaries
	> void classify(
		const double simulation_time,
		dccrg::Dccrg<Cell_Data, Geometry>& grid,
		Value_Boundaries& value_boundaries,
		Copy_Boundaries& copy_boundaries
	) {
		using std::get;

		constexpr Cell_Type_Variable cell_type{};
		constexpr Value_Boundary_Id_T value_bdy_id{};
		constexpr Source_Of_Copy_Cell_T source{};

		const auto cell_data_pointers = grid.get_cell_data_pointers();

		for (const auto& cell_item: cell_data_pointers) {
			boost::optional<size_t> classification_result;

			const auto& cell_id = get<0>(cell_item);

			// process inner and outer cells
			if (cell_id == dccrg::error_cell) {
				continue;
			}

			const auto& offset = get<2>(cell_item);

			// skip neighbors
			if (offset[0] != 0 or offset[1] != 0 or offset[2] != 0) {
				continue;
			}

			const auto
				cell_start = grid.geometry.get_min(cell_id),
				cell_end = grid.geometry.get_max(cell_id);

			classification_result = value_boundaries.add_cell(
				simulation_time,
				cell_id,
				cell_start,
				cell_end
			);
			if (not classification_result) {
				std::cerr <<  __FILE__ << "(" << __LINE__ << "): "
					"Couldn't add cell " << cell_id << " to value boundaries."
					<< std::endl;
				abort();
			}

			auto* const cell_data = get<1>(cell_item);

			if (*classification_result > 0) {
				(*cell_data)[cell_type] = value_boundary_cell;
				continue;
			}

			classification_result = copy_boundaries.add_cell(
				cell_id,
				cell_start,
				cell_end
			);
			if (not classification_result) {
				std::cerr <<  __FILE__ << "(" << __LINE__ << "): "
					"Couldn't add cell " << cell_id << " to copy boundaries."
					<< std::endl;
				abort();
			}
			if (*classification_result > 0) {
				(*cell_data)[cell_type] = copy_boundary_cell;
				continue;
			}

			(*cell_data)[cell_type] = normal_cell;
		}
		Cell_Data::set_transfer_all(true, cell_type);
		grid.update_copies_of_remote_neighbors();
		Cell_Data::set_transfer_all(false, cell_type);

		// figure out dont solve cells
		for (size_t i = 0; i < cell_data_pointers.size(); i++) {
			const auto& cell_id = get<0>(cell_data_pointers[i]);

			if (cell_id == dccrg::error_cell) {
				continue;
			}

			const auto& offset = get<2>(cell_data_pointers[i]);
			if (offset[0] != 0 or offset[1] != 0 or offset[2] != 0) {
				continue;
			}

			auto* const cell_data = get<1>(cell_data_pointers[i]);
			if (
				(*cell_data)[cell_type] == normal_cell
				or (*cell_data)[cell_type] == dont_solve_cell
			) {
				continue;
			}

			bool has_normal_neighbor = false;
			i++;
			while (i < cell_data_pointers.size()) {
				const auto& neighbor_id = get<0>(cell_data_pointers[i]);

				if (neighbor_id == dccrg::error_cell) {
					i--;
					break;
				}

				const auto& neigh_offset = get<2>(cell_data_pointers[i]);
				if (neigh_offset[0] == 0 and neigh_offset[1] == 0 and neigh_offset[2] == 0) {
					i--;
					break;
				}

				auto* const neighbor_data = get<1>(cell_data_pointers[i]);
				if ((*neighbor_data)[cell_type] == normal_cell) {
					has_normal_neighbor = true;
					break;
				}

				i++;
			}

			if (not has_normal_neighbor) {
				(*cell_data)[cell_type] = dont_solve_cell;
			}
		}

		// tell other processes about dont_solve_cells
		Cell_Data::set_transfer_all(true, cell_type);
		grid.update_copies_of_remote_neighbors();
		Cell_Data::set_transfer_all(false, cell_type);

		// set boundary ids of value boundary cells
		for (size_t bdy_id = 0; bdy_id < value_boundaries.get_number_of_boundaries(); bdy_id++) {
			for (const auto& cell_id: value_boundaries.get_cells(bdy_id)) {
				auto* const cell_data = grid[cell_id];
				if (cell_data == nullptr) {
					std::cerr <<  __FILE__ << "(" << __LINE__ << ")" << std::endl;
					abort();
				}
				(*cell_data)[value_bdy_id] = bdy_id;
			}
		}

		// set sources of copy boundary cells
		for (size_t i = 0; i < cell_data_pointers.size(); i++) {
			const auto& cell_id = get<0>(cell_data_pointers[i]);

			if (cell_id == dccrg::error_cell) {
				continue;
			}

			const auto& offset = get<2>(cell_data_pointers[i]);
			if (offset[0] != 0 or offset[1] != 0 or offset[2] != 0) {
				continue;
			}

			auto* const cell_data = get<1>(cell_data_pointers[i]);
			if ((*cell_data)[cell_type] != copy_boundary_cell) {
				continue;
			}

			bool source_found = false;
			i++;
			while (i < cell_data_pointers.size()) {
				const auto& neighbor_id = get<0>(cell_data_pointers[i]);

				if (neighbor_id == dccrg::error_cell) {
					i--;
					break;
				}

				const auto& neigh_offset = get<2>(cell_data_pointers[i]);
				if (neigh_offset[0] == 0 and neigh_offset[1] == 0 and neigh_offset[2] == 0) {
					i--;
					break;
				}

				auto* const neigh_data = get<1>(cell_data_pointers[i]);
				if ((*neigh_data)[cell_type] == normal_cell) {
					source_found = true;
					(*cell_data)[source] = neighbor_id;
					break;
				}

				i++;
			}

			// make copy cell without sources into dont solve
			if (not source_found) {
				(*cell_data)[cell_type] = dont_solve_cell;
			}
		}

		// tell other processes about dont_solve_cells
		Cell_Data::set_transfer_all(true, cell_type);
		grid.update_copies_of_remote_neighbors();
		Cell_Data::set_transfer_all(false, cell_type);
	}
};


}} // namespaces


#endif // ifndef PAMHD_MHD_BOUNDARIES_HPP
