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
template<class Cell> class Boundary_Classifier
{
private:

	std::vector<Cell>
		// cells given to classify but without dont_solves
		inner_solve_cells_, outer_solve_cells_,
		// classified cells but without boundaries and dont_solves
		inner_normal_cells_, outer_normal_cells_;

	std::vector<
		std::pair<
			Cell, // cell id
			size_t // id of 1st boundary of any type that cell belongs to
		>
	> value_boundary_cells_;

	std::vector<
		std::pair<
			Cell, // cell id
			Cell // source cell id for data
		>
	> copy_boundary_cells_;

	std::vector<Cell> dont_solve_cells_;


public:

	Boundary_Classifier() :
		inner_solve_cells(inner_solve_cells_),
		outer_solve_cells(outer_solve_cells_),
		inner_normal_cells(inner_normal_cells_),
		outer_normal_cells(outer_normal_cells_),
		value_boundary_cells(value_boundary_cells_),
		copy_boundary_cells(copy_boundary_cells_),
		dont_solve_cells(dont_solve_cells_)
	{}

	const std::vector<Cell>
		//! inner cells given to classify() but without dont_solve_cells
		&inner_solve_cells,
		//! outer cells given to classify() but without dont_solve_cells
		&outer_solve_cells,
		//! inner cells given to classify() but without dont_solve and boundary cells
		&inner_normal_cells,
		//! outer cells given to classify() but without dont_solve and boundary cells
		&outer_normal_cells;

	//! cells given to classify classified as value boundary cells
	const std::vector<
		std::pair<
			Cell,
			size_t
		>
	>& value_boundary_cells;

	//! cells given to classify classified as copy boundary cells
	const std::vector<
		std::pair<
			Cell,
			Cell
		>
	>& copy_boundary_cells;

	//! cells given to classify that don't have non-boundary face neighbors
	const std::vector<Cell>& dont_solve_cells;


	template<
		class Cell_Type_V,
		class Cell_Data,
		class Geometry,
		class Value_Boundaries,
		class Copy_Boundaries
	> void classify(
		const double simulation_time,
		const std::vector<Cell>& given_inner_cells,
		const std::vector<Cell>& given_outer_cells,
		dccrg::Dccrg<Cell_Data, Geometry>& grid,
		Value_Boundaries& value_boundaries,
		Copy_Boundaries& copy_boundaries
	) {
		constexpr auto
			normal_value = 0,
			value_boundary_value = 1,
			copy_boundary_value = 2,
			dont_solve_value
				= std::numeric_limits<typename Cell_Type_V::data_type>::max();

		static_assert(
			copy_boundary_value < dont_solve_value,
			"Given type representing cell boundary type is too small."
		);

		std::vector<Cell> cells(
			given_inner_cells.cbegin(),
			given_inner_cells.cend()
		);
		cells.insert(
			cells.end(),
			given_outer_cells.cbegin(),
			given_outer_cells.cend()
		);

		value_boundaries.clear_cells();
		copy_boundaries.clear_cells();
		for (const auto& cell_id: cells) {
			boost::optional<size_t> classification_result;

			auto* const cell_data = grid[cell_id];
			if (cell_data == nullptr) {
				std::cerr <<  __FILE__ << "(" << __LINE__ << "): " << cell_id << std::endl;
				abort();
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
			if (*classification_result > 0) {
				(*cell_data)[Cell_Type_V()] = value_boundary_value;
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
				(*cell_data)[Cell_Type_V()] = copy_boundary_value;
				continue;
			}

			(*cell_data)[Cell_Type_V()] = normal_value;
		}
		Cell_Data::set_transfer_all(true, Cell_Type_V());
		grid.update_copies_of_remote_neighbors();
		Cell_Data::set_transfer_all(false, Cell_Type_V());

		// figure out dont_solve_cells_
		const auto has_normal_neighbor
			= [&grid](const Cell& cell){
				bool ret_val = false;
				const auto face_neighbors_of = grid.get_face_neighbors_of(cell);
				for (const auto& item: face_neighbors_of) {
					const auto& neighbor = item.first;
					auto* const neighbor_data = grid[neighbor];
					if (neighbor_data == nullptr) {
						std::cerr <<  __FILE__ << "(" << __LINE__ << "): " << neighbor << std::endl;
						abort();
					}
					if ((*neighbor_data)[Cell_Type_V()] == normal_value) {
						ret_val = true;
						break;
					}
				}
				return ret_val;
			};
		for (const auto& cell: cells) {
			auto* const cell_data = grid[cell];
			if (cell_data == nullptr) {
				std::cerr <<  __FILE__ << "(" << __LINE__ << ")" << std::endl;
				abort();
			}
			if (
				(*cell_data)[Cell_Type_V()] != normal_value
				and not has_normal_neighbor(cell)
			) {
				(*cell_data)[Cell_Type_V()] = dont_solve_value;
			}
		}

		// populate copy_boundary_cells_ and set their source
		this->copy_boundary_cells_.clear();
		for (const auto& cell: cells) {
			auto* const cell_data = grid[cell];
			if (cell_data == nullptr) {
				std::cerr <<  __FILE__ << "(" << __LINE__ << ")" << std::endl;
				abort();
			}
			if ((*cell_data)[Cell_Type_V()] != copy_boundary_value) {
				continue;
			}

			bool source_found = false;
			const auto face_neighbors_of = grid.get_face_neighbors_of(cell);
			for (const auto& item: face_neighbors_of) {
				const auto& neighbor = item.first;

				auto* const neighbor_data = grid[neighbor];
				if (neighbor_data == nullptr) {
					std::cerr <<  __FILE__ << "(" << __LINE__ << ")" << std::endl;
					abort();
				}

				if ((*neighbor_data)[Cell_Type_V()] == normal_value) {
					source_found = true;
					this->copy_boundary_cells_.emplace_back(cell, neighbor);
					break;
				}
			}
			if (not source_found) {
				// make copy cell without sources into dont solves
				(*cell_data)[Cell_Type_V()] = dont_solve_value;
			}
		}

		// populate value_boundary_cells_
		this->value_boundary_cells_.clear();
		for (size_t bdy_id = 0; bdy_id < value_boundaries.get_number_of_boundaries(); bdy_id++) {
			for (const auto& cell_id: value_boundaries.get_cells(bdy_id)) {
				auto* const cell_data = grid[cell_id];
				if (cell_data == nullptr) {
					std::cerr <<  __FILE__ << "(" << __LINE__ << ")" << std::endl;
					abort();
				}
				if (
					(*cell_data)[Cell_Type_V()] != value_boundary_value
					and (*cell_data)[Cell_Type_V()] != dont_solve_value
				) {
					std::cerr <<  __FILE__ << "(" << __LINE__ << ")" << std::endl;
					abort();
				}
				this->value_boundary_cells_.emplace_back(cell_id, bdy_id);
			}
		}

		this->dont_solve_cells_.clear();
		this->inner_solve_cells_.clear();
		this->inner_normal_cells_.clear();
		for (const auto& cell: given_inner_cells) {
			auto* const cell_data = grid[cell];
			if (cell_data == nullptr) {
				std::cerr <<  __FILE__ << "(" << __LINE__ << ")" << std::endl;
				abort();
			}
			if ((*cell_data)[Cell_Type_V()] == dont_solve_value) {
				this->dont_solve_cells_.push_back(cell);
			} else {
				this->inner_solve_cells_.push_back(cell);
				if ((*cell_data)[Cell_Type_V()] == normal_value) {
					this->inner_normal_cells_.push_back(cell);
				}
			}
		}
		this->outer_solve_cells_.clear();
		this->outer_normal_cells_.clear();
		for (const auto& cell: given_outer_cells) {
			auto* const cell_data = grid[cell];
			if (cell_data == nullptr) {
				std::cerr <<  __FILE__ << "(" << __LINE__ << ")" << std::endl;
				abort();
			}
			if ((*cell_data)[Cell_Type_V()] == dont_solve_value) {
				this->dont_solve_cells_.push_back(cell);
			} else {
				this->outer_solve_cells_.push_back(cell);
				if ((*cell_data)[Cell_Type_V()] == normal_value) {
					this->outer_normal_cells_.push_back(cell);
				}
			}
		}
	}
};


}} // namespaces


#endif // ifndef PAMHD_MHD_BOUNDARIES_HPP
