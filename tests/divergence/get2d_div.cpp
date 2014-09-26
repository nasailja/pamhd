/*
Tests vector field divergence calculation of PAMHD in 2d.

Copyright 2014 Ilja Honkonen
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
#include "limits"
#include "tuple"
#include "vector"

#include "dccrg.hpp"
#include "dccrg_cartesian_geometry.hpp"

#include "divergence/remove.hpp"


// assigned to y component of vector field
double function(const double x, const double y)
{
	return x*x/y;
}

double div_of_function(const double x, const double y)
{
	return -std::pow(x/y, 2);
}


struct Vector_Field {
	using data_type = std::array<double, 3>;
};

struct Divergence {
	using data_type = double;
};

using Cell = gensimcell::Cell<
	gensimcell::Always_Transfer,
	Vector_Field,
	Divergence
>;


/*!
Returns maximum norm if p == 0
*/
template<class Grid_T> double get_diff_lp_norm(
	const std::vector<uint64_t>& cells,
	const Grid_T& grid,
	const double p,
	const double cell_volume
) {
	double local_norm = 0, global_norm = 0;
	for (const auto& cell: cells) {
		const auto* const cell_data = grid[cell];
		if (cell_data == NULL) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< ": No data for cell " << cell
				<< std::endl;
			abort();
		}

		const auto center = grid.geometry.get_center(cell);

		const auto div_of = div_of_function(center[0], center[1]);

		if (p == 0) {
			local_norm = std::max(
				local_norm,
				std::fabs((*cell_data)[Divergence()] - div_of)
			);
		} else {
			local_norm += std::pow(
				std::fabs((*cell_data)[Divergence()] - div_of),
				p
			);
		}
	}
	local_norm *= cell_volume;

	if (p == 0) {
		MPI_Comm comm = grid.get_communicator();
		MPI_Allreduce(&local_norm, &global_norm, 1, MPI_DOUBLE, MPI_MAX, comm);
		MPI_Comm_free(&comm);
		return global_norm;
	} else {
		MPI_Comm comm = grid.get_communicator();
		MPI_Allreduce(&local_norm, &global_norm, 1, MPI_DOUBLE, MPI_SUM, comm);
		MPI_Comm_free(&comm);
		return std::pow(global_norm, 1.0 / p);
	}
}


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

	const std::array<double, 3> grid_length{{5, 5, 1}};

	double old_norm = std::numeric_limits<double>::max();
	size_t old_nr_of_cells = 0;
	for (size_t nr_of_cells = 8; nr_of_cells <= 256; nr_of_cells *= 2) {

		dccrg::Dccrg<Cell, dccrg::Cartesian_Geometry> grid;

		const std::array<uint64_t, 3> grid_size{{nr_of_cells + 2, nr_of_cells + 2, 1}};

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
				grid_length[0] / grid_size[0],
				grid_length[1] / grid_size[1],
				grid_length[2] / grid_size[2]
			}},
			grid_start{{
				1, 1, 0
			}};

		const double cell_volume
			= cell_length[0] * cell_length[1] * cell_length[2];

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

			auto& vec = (*cell_data)[Vector_Field()];
			vec[0] = 0;
			vec[1] = function(center[0], center[1]);
			vec[2] = 0;
		}
		grid.update_copies_of_remote_neighbors();

		std::vector<uint64_t> solve_cells;
		for (const auto& cell: all_cells) {
			const auto index = grid.mapping.get_indices(cell);
			if (
				index[0] > 0
				and index[0] < grid_size[0] - 1
				and index[1] > 0
				and index[1] < grid_size[1] - 1
			) {
				solve_cells.push_back(cell);
			}
		}

		pamhd::divergence::get_divergence(
			solve_cells,
			grid,
			std::tuple<Vector_Field>(),
			std::tuple<Divergence>()
		);

		const double
			p_of_norm = 2,
			norm = get_diff_lp_norm(solve_cells, grid, p_of_norm, cell_volume);

		if (norm > old_norm) {
			if (grid.get_rank() == 0) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< ": Norm with " << nr_of_cells
					<< " cells " << norm
					<< " is larger than with " << nr_of_cells / 2
					<< " cells " << old_norm
					<< std::endl;
			}
			abort();
		}

		if (old_nr_of_cells >= 64) {
			const double order_of_accuracy
				= -log(norm / old_norm)
				/ log(double(nr_of_cells) / old_nr_of_cells);

			if (order_of_accuracy < 1.7) {
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
		old_norm = norm;
	}

	MPI_Finalize();

	return EXIT_SUCCESS;
}
