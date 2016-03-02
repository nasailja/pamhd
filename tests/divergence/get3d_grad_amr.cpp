/*
Tests scalar field gradient calculation of PAMHD in 3d.

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
#include "limits"
#include "vector"

#include "dccrg.hpp"
#include "dccrg_cartesian_geometry.hpp"
#include "prettyprint.hpp"

#include "divergence/remove.hpp"


double function(const std::array<double, 3>& r)
{
	const double x = r[0], y = r[1], z = r[2];

	return
		std::exp(2*x) / 2 + 5 * x*x / 2 + 6 * x
		- std::pow(std::cos(y), 2) / 2 + 2 * y*y + 5 * y
		+ 3 * z*z*z*z / 4 - 2 * z*z*z / 3 + z*z / 2 - 2 * z;
}

std::array<double, 3> grad_of_function(const std::array<double, 3>& r)
{
	const double x = r[0], y = r[1], z = r[2];

	return {{
		std::exp(2*x) + 5*x + 6,
		std::sin(y) * std::cos(y) + 4*y + 5,
		3*z*z*z - 2*z*z + z - 2
	}};
}


struct Scalar {
	using data_type = double;
};

struct Gradient {
	using data_type = std::array<double, 3>;
};

using Cell = gensimcell::Cell<
	gensimcell::Always_Transfer,
	Scalar,
	Gradient
>;


/*!
Returns maximum norm
*/
template<class Grid_T> std::array<double, 3> get_max_norm(
	const std::vector<uint64_t>& cells,
	const Grid_T& grid
) {
	std::array<double, 3>
		local_norm{{0, 0, 0}},
		global_norm{{0, 0, 0}};

	for (const auto& cell: cells) {
		const auto* const cell_data = grid[cell];
		if (cell_data == NULL) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< ": No data for cell " << cell
				<< std::endl;
			abort();
		}

		const auto center = grid.geometry.get_center(cell);

		const auto grad_of = grad_of_function(center);

		for (size_t i = 0; i < 3; i++) {
			local_norm[i] = std::max(
				local_norm[i],
				std::fabs((*cell_data)[Gradient()][i] - grad_of[i])
			);
		}
	}

	MPI_Comm comm = grid.get_communicator();
	MPI_Allreduce(
		local_norm.data(),
		global_norm.data(),
		3,
		MPI_DOUBLE,
		MPI_MAX,
		comm
	);
	MPI_Comm_free(&comm);

	return global_norm;
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
	const int max_refinement_level = 1;

	std::array<double, 3> old_norm{{
		std::numeric_limits<double>::max(),
		std::numeric_limits<double>::max(),
		std::numeric_limits<double>::max()
	}};
	size_t old_nr_of_cells = 0;
	for (size_t nr_of_cells = 4; nr_of_cells <= 16; nr_of_cells *= 2) {

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
				false, false, false
			)
		) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< ": Couldn't initialize grid."
				<< std::endl;
			abort();
		}

		const std::array<double, 3>
			cell_length{{
				double(5) / (grid_size[0] - 2),
				double(5) / (grid_size[1] - 2),
				double(1) / (grid_size[2] - 2),
			}},
			grid_start{{1, 1, 0}};

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

		for (int i = 0; i < max_refinement_level; i++) {
			for (const auto& cell: grid.get_cells()) {
				const auto center = grid.geometry.get_center(cell);
				if (
					    center[0] > grid_start[0] + cell_length[0] + 5.0*2/4
					and center[0] < grid_start[0] + cell_length[0] + 5.0*4/4
					and center[1] > grid_start[1] + cell_length[1] + 5.0*1/4
					and center[1] < grid_start[1] + cell_length[1] + 5.0*3/4
					and center[2] > grid_start[2] + cell_length[2] + 2.0/4
					and center[2] < grid_start[2] + cell_length[2] + 4.0/4
				) {
					grid.refine_completely(cell);
				}
			}
			grid.stop_refining();
		}

		std::vector<uint64_t> solve_cells;
		for (const auto& cell: grid.get_cells()) {
			auto* const cell_data = grid[cell];
			if (cell_data == nullptr) {
				std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
				abort();
			}

			const auto center = grid.geometry.get_center(cell);
			(*cell_data)[Scalar()] = function(center);

			if (
				    center[0] > grid_start[0] + cell_length[0]
				and center[0] < grid_start[0] + cell_length[0] + 5
				and center[1] > grid_start[1] + cell_length[1]
				and center[1] < grid_start[1] + cell_length[1] + 5
				and center[2] > grid_start[2] + cell_length[2]
				and center[2] < grid_start[2] + cell_length[2] + 1
			) {
				solve_cells.push_back(cell);
			}
		}
		grid.update_copies_of_remote_neighbors();

		pamhd::divergence::get_gradient(
			solve_cells,
			grid,
			[](Cell& cell_data) -> Scalar::data_type& {
				return cell_data[Scalar()];
			},
			[](Cell& cell_data) -> Gradient::data_type& {
				return cell_data[Gradient()];
			}
		);

		const auto norm = get_max_norm(solve_cells, grid);

		for (size_t dim = 0; dim < 3; dim++) {
			if (norm[dim] > old_norm[dim]) {
				if (grid.get_rank() == 0) {
					std::cerr << __FILE__ << ":" << __LINE__
						<< " dim " << dim
						<< ": Norm with " << solve_cells.size()
						<< " cells " << norm[dim]
						<< " is larger than with " << old_nr_of_cells
						<< " cells " << old_norm[dim]
						<< std::endl;
				}
				abort();
			}

			if (old_nr_of_cells > 0) {
				const double order_of_accuracy
					= -log(norm[dim] / old_norm[dim])
					/ log(double(solve_cells.size()) / old_nr_of_cells);

				if (order_of_accuracy < 0.15) {
					if (grid.get_rank() == 0) {
						std::cerr << __FILE__ << ":" << __LINE__
							<< " dim " << dim
							<< ": Order of accuracy from "
							<< old_nr_of_cells << " to " << solve_cells.size()
							<< " is too low: " << order_of_accuracy
							<< std::endl;
					}
					abort();
				}
			}
		}

		old_nr_of_cells = solve_cells.size();
		old_norm = norm;
	}

	MPI_Finalize();

	return EXIT_SUCCESS;
}
