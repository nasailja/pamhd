/*
Tests vector field divergence removal of PAMHD in 3d.

Copyright 2014, 2015, 2016 Ilja Honkonen
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


#include "array"
#include "cstdlib"
#include "iostream"
#include "iomanip"
#include "limits"
#include "tuple"
#include "vector"

#include "dccrg.hpp"
#include "dccrg_cartesian_geometry.hpp"
#include "prettyprint.hpp"

#include "divergence/remove.hpp"


int Poisson_Cell::transfer_switch = Poisson_Cell::INIT;


std::array<double, 3> function(const std::array<double, 3>& r)
{
	return {{0.1 * r[0], 1 + 0.1 * std::sin(r[1]), 1 - 0.1 * std::cos(r[2])}};
}

double div_of_function(const std::array<double, 3>& r)
{
	return 0.1 * (1 + std::cos(r[1]) + std::sin(r[2]));
}


struct Vector_Field {
	using data_type = std::array<double, 3>;
};

struct Divergence_Before {
	using data_type = double;
};

struct Divergence_After {
	using data_type = double;
};

struct Gradient {
	using data_type = std::array<double, 3>;
};

using Cell = gensimcell::Cell<
	gensimcell::Always_Transfer,
	Vector_Field,
	Divergence_Before,
	Divergence_After,
	Gradient
>;


int main(int argc, char* argv[])
{
	if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
		std::cerr << "Couldn't initialize MPI." << std::endl;
		abort();
	}

	MPI_Comm comm = MPI_COMM_WORLD;

	int rank = 0, comm_size = 0;
	if (MPI_Comm_rank(comm, &rank) != MPI_SUCCESS) {
		std::cerr << "Couldn't obtain MPI rank." << std::endl;
		abort();
	}
	if (MPI_Comm_size(comm, &comm_size) != MPI_SUCCESS) {
		std::cerr << "Couldn't obtain size of MPI communicator." << std::endl;
		abort();
	}


	// intialize Zoltan
	float zoltan_version;
	if (Zoltan_Initialize(argc, argv, &zoltan_version) != ZOLTAN_OK) {
		std::cerr << "Zoltan_Initialize failed." << std::endl;
		abort();
	}

	const unsigned int neighborhood_size = 0;
	const int max_refinement_level = 0;

	double old_div = std::numeric_limits<double>::max();
	size_t old_nr_of_cells = 0;
	for (size_t nr_of_cells = 8; nr_of_cells <= 64; nr_of_cells *= 2) {

		dccrg::Dccrg<Cell, dccrg::Cartesian_Geometry> grid;

		const std::array<uint64_t, 3> grid_size{{
			nr_of_cells + 2,
			nr_of_cells + 2,
			nr_of_cells + 2
		}};

		if (
			not grid.initialize(
				grid_size,
				comm,
				"RANDOM",
				neighborhood_size,
				max_refinement_level,
				false,
				false,
				false
			)
		) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< ": Couldn't initialize grid."
				<< std::endl;
			abort();
		}

		const std::array<double, 3>
			cell_length{{
				double(1) / (grid_size[0] - 2),
				2 * M_PI / (grid_size[1] - 2),
				2 * M_PI / (grid_size[2] - 2)
			}},
			grid_start{{
				-cell_length[0], -cell_length[1], -cell_length[2]
			}};

		dccrg::Cartesian_Geometry::Parameters geom_params;
		geom_params.start = grid_start;
		geom_params.level_0_cell_length = cell_length;

		if (not grid.set_geometry(geom_params)) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< ": Couldn't set grid geometry."
				<< std::endl;
			abort();
		}

		grid.balance_load();

		const auto all_cells = grid.get_cells();
		for (const auto& cell: all_cells) {
			auto* const cell_data = grid[cell];
			if (cell_data == NULL) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< ": No data for cell " << cell
					<< std::endl;
				abort();
			}

			const auto center = grid.geometry.get_center(cell);
			(*cell_data)[Vector_Field()] = function(center);
		}
		grid.update_copies_of_remote_neighbors();

		// classify cells
		std::vector<uint64_t>
			solve_cells,
			boundary_cells;
		for (const auto& cell: all_cells) {
			const auto index = grid.mapping.get_indices(cell);
			if (
				index[0] > 0
				and index[0] < grid_size[0] - 1
				and index[1] > 0
				and index[1] < grid_size[1] - 1
				and index[2] > 0
				and index[2] < grid_size[2] - 1
			) {
				solve_cells.push_back(cell);
			} else {
				boundary_cells.push_back(cell);
			}
		}

		// apply copy boundaries
		for (const auto& cell: boundary_cells) {
			const auto index = grid.mapping.get_indices(cell);
			auto neighbor_index = index;

			if (index[0] == 0) {
				neighbor_index[0] = index[0] + 1;
			} else if (index[0] == grid_size[0] - 1) {
				neighbor_index[0] = index[0] - 1;
			} else if (index[1] == 0) {
				neighbor_index[1] = index[1] + 1;
			} else if (index[1] == grid_size[1] - 1) {
				neighbor_index[1] = index[1] - 1;
			} else if (index[2] == 0) {
				neighbor_index[2] = index[2] + 1;
			} else if (index[2] == grid_size[2] - 1) {
				neighbor_index[2] = index[2] - 1;
			}
			const auto neighbor = grid.mapping.get_cell_from_indices(neighbor_index, 0);

			const auto* const neighbor_data = grid[neighbor];
			auto* const cell_data = grid[cell];
			if (cell_data == NULL or neighbor_data == NULL) {
				std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
				abort();
			}

			(*cell_data)[Vector_Field()] = (*neighbor_data)[Vector_Field()];
		}
		grid.update_copies_of_remote_neighbors();

		auto Vector_Getter = [](Cell& cell_data) -> Vector_Field::data_type& {
			return cell_data[Vector_Field()];
		};
		auto Div_After_Getter = [](Cell& cell_data) -> Divergence_After::data_type& {
			return cell_data[Divergence_After()];
		};
		const double div_before = pamhd::divergence::get_divergence(
			solve_cells,
			grid,
			Vector_Getter,
			[](Cell& cell_data) -> Divergence_Before::data_type& {
				return cell_data[Divergence_Before()];
			}
		);

		pamhd::divergence::remove(
			solve_cells,
			boundary_cells,
			{},
			grid,
			Vector_Getter,
			Div_After_Getter,
			[](Cell& cell_data) -> Gradient::data_type& {
				return cell_data[Gradient()];
			},
			1000, 0, 1e-2, 2, 10, false
		);
		grid.update_copies_of_remote_neighbors();

		// update copy boundaries to correspond to removed divergence
		for (const auto& cell: boundary_cells) {
			const auto index = grid.mapping.get_indices(cell);
			auto neighbor_index = index;

			if (index[0] == 0) {
				neighbor_index[0] = index[0] + 1;
			} else if (index[0] == grid_size[0] - 1) {
				neighbor_index[0] = index[0] - 1;
			} else if (index[1] == 0) {
				neighbor_index[1] = index[1] + 1;
			} else if (index[1] == grid_size[1] - 1) {
				neighbor_index[1] = index[1] - 1;
			} else if (index[2] == 0) {
				neighbor_index[2] = index[2] + 1;
			} else if (index[2] == grid_size[2] - 1) {
				neighbor_index[2] = index[2] - 1;
			}
			const auto neighbor = grid.mapping.get_cell_from_indices(neighbor_index, 0);

			const auto* const neighbor_data = grid[neighbor];
			auto* const cell_data = grid[cell];
			if (cell_data == NULL or neighbor_data == NULL) {
				std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
				abort();
			}

			(*cell_data)[Vector_Field()] = (*neighbor_data)[Vector_Field()];
		}
		grid.update_copies_of_remote_neighbors();

		const double div_after = pamhd::divergence::get_divergence(
			solve_cells,
			grid,
			Vector_Getter,
			Div_After_Getter
		);

		if (div_after > div_before) {
			if (grid.get_rank() == 0) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< ": Divergence after removal " << div_after
					<< " is larger than before " << div_before
					<< " with " << nr_of_cells << " cells."
					<< std::endl;
			}
			abort();
		}

		if (old_nr_of_cells > 0) {
			const double
				order_of_accuracy
					= -log(div_after / old_div)
					/ log(double(nr_of_cells) / old_nr_of_cells);

			if (
				(nr_of_cells == 16 and order_of_accuracy < 0.4)
				or (nr_of_cells == 32 and order_of_accuracy < 0.5)
				or (nr_of_cells == 64 and order_of_accuracy < 0.7)
			) {
				if (grid.get_rank() == 0) {
					std::cerr << __FILE__ << ":" << __LINE__
						<< ": Order of accuracy from "
						<< old_nr_of_cells << " to " << nr_of_cells
						<< " is too low: " << order_of_accuracy
						<< std::endl;
				}
				abort();
			}
		}

		old_nr_of_cells = nr_of_cells;
		old_div = div_after;
	}

	MPI_Finalize();

	return EXIT_SUCCESS;
}
