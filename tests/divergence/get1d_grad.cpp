/*
Tests scalar field gradient calculation of PAMHD in 1d.

Copyright 2014, 2015 Ilja Honkonen
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
#include "functional"
#include "iostream"
#include "limits"
#include "tuple"
#include "vector"

#include "dccrg.hpp"
#include "dccrg_cartesian_geometry.hpp"

#include "divergence/remove.hpp"


double function(const double x)
{
	return std::pow(std::sin(x / 2), 2);
}

// gradient only in x dimension
double grad_of_function(const double x)
{
	return std::cos(x / 2) * std::sin(x / 2);
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
Returns maximum norm if p == 0
*/
template<class Grid_T> double get_diff_lp_norm(
	const std::vector<uint64_t>& cells,
	const Grid_T& grid,
	const double p,
	const double cell_volume,
	const size_t dimension
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

		if (p == 0) {
			local_norm = std::max(
				local_norm,
				std::fabs(
					(*cell_data)[Gradient()][dimension]
					- grad_of_function(center[dimension])
				)
			);
		} else {
			local_norm += std::pow(
				std::fabs(
					(*cell_data)[Gradient()][dimension]
					- grad_of_function(center[dimension])
				),
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

	double
		old_norm_x = std::numeric_limits<double>::max(),
		old_norm_y = std::numeric_limits<double>::max(),
		old_norm_z = std::numeric_limits<double>::max();

	size_t old_nr_of_cells = 0;
	for (size_t nr_of_cells = 8; nr_of_cells <= 4096; nr_of_cells *= 2) {

		dccrg::Dccrg<Cell, dccrg::Cartesian_Geometry> grid_x, grid_y, grid_z;

		const std::array<uint64_t, 3>
			grid_size_x{{nr_of_cells + 2, 1, 1}},
			grid_size_y{{1, nr_of_cells + 2, 1}},
			grid_size_z{{1, 1, nr_of_cells + 2}};

		if (not grid_x.initialize(grid_size_x,comm,"RANDOM",0,0,false,false,false)) {
			std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
			abort();
		}
		if (not grid_y.initialize(grid_size_y,comm,"RANDOM",0,0,false,false,false)) {
			std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
			abort();
		}
		if (not grid_z.initialize(grid_size_z,comm,"RANDOM",0,0,false,false,false)) {
			std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
			abort();
		}

		const std::array<double, 3>
			cell_length_x{{2 * M_PI / (grid_size_x[0] - 2), 1, 1}},
			cell_length_y{{1, 2 * M_PI / (grid_size_y[1] - 2), 1}},
			cell_length_z{{1, 1, 2 * M_PI / (grid_size_z[2] - 2)}},
			grid_start_x{{-cell_length_x[0], 0, 0}},
			grid_start_y{{0, -cell_length_y[1], 0}},
			grid_start_z{{0, 0, -cell_length_z[2]}};

		const double cell_volume
			= cell_length_x[0] * cell_length_y[1] * cell_length_z[2];

		dccrg::Cartesian_Geometry::Parameters geom_params_x,geom_params_y,geom_params_z;
		geom_params_x.start = grid_start_x;
		geom_params_x.level_0_cell_length = cell_length_x;
		geom_params_y.start = grid_start_y;
		geom_params_y.level_0_cell_length = cell_length_y;
		geom_params_z.start = grid_start_z;
		geom_params_z.level_0_cell_length = cell_length_z;


		if (not grid_x.set_geometry(geom_params_x)) {
			std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
			abort();
		}
		if (not grid_y.set_geometry(geom_params_y)) {
			std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
			abort();
		}
		if (not grid_z.set_geometry(geom_params_z)) {
			std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
			abort();
		}

		const auto all_cells = grid_x.get_cells();
		for (const auto& cell: all_cells) {
			auto
				*const cell_data_x = grid_x[cell],
				*const cell_data_y = grid_y[cell],
				*const cell_data_z = grid_z[cell];

			if (cell_data_x == NULL or cell_data_y == NULL or cell_data_z == NULL) {
				std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
				abort();
			}

			const auto center = grid_x.geometry.get_center(cell);

			(*cell_data_x)[Scalar()] = function(center[0]);
			(*cell_data_y)[Scalar()] = function(center[0]);
			(*cell_data_z)[Scalar()] = function(center[0]);
		}
		grid_x.update_copies_of_remote_neighbors();
		grid_y.update_copies_of_remote_neighbors();
		grid_z.update_copies_of_remote_neighbors();

		// exclude one layer of boundary cells
		std::vector<uint64_t> solve_cells;
		for (const auto& cell: all_cells) {
			const auto index = grid_x.mapping.get_indices(cell);
			if (index[0] > 0 and index[0] < grid_size_x[0] - 1) {
				solve_cells.push_back(cell);
			}
		}

		auto Scalar_Getter = [](Cell& cell_data) -> Scalar::data_type& {
			return cell_data[Scalar()];
		};
		auto Gradient_Getter = [](Cell& cell_data) -> Gradient::data_type& {
			return cell_data[Gradient()];
		};
		pamhd::divergence::get_gradient(
			solve_cells, grid_x, Scalar_Getter, Gradient_Getter
		);
		pamhd::divergence::get_gradient(
			solve_cells, grid_y, Scalar_Getter, Gradient_Getter
		);
		pamhd::divergence::get_gradient(
			solve_cells, grid_z, Scalar_Getter, Gradient_Getter
		);

		const double
			p_of_norm = 2,
			norm_x = get_diff_lp_norm(solve_cells, grid_x, p_of_norm, cell_volume, 0),
			norm_y = get_diff_lp_norm(solve_cells, grid_y, p_of_norm, cell_volume, 1),
			norm_z = get_diff_lp_norm(solve_cells, grid_z, p_of_norm, cell_volume, 2);

		if (norm_x > old_norm_x) {
			if (grid_x.get_rank() == 0) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< ": X norm with " << nr_of_cells
					<< " cells " << norm_x
					<< " is larger than with " << nr_of_cells / 2
					<< " cells " << old_norm_x
					<< std::endl;
			}
			abort();
		}
		if (norm_y > old_norm_y) {
			if (grid_y.get_rank() == 0) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< ": X norm with " << nr_of_cells
					<< " cells " << norm_y
					<< " is larger than with " << nr_of_cells / 2
					<< " cells " << old_norm_y
					<< std::endl;
			}
			abort();
		}
		if (norm_z > old_norm_z) {
			if (grid_z.get_rank() == 0) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< ": X norm with " << nr_of_cells
					<< " cells " << norm_z
					<< " is larger than with " << nr_of_cells / 2
					<< " cells " << old_norm_z
					<< std::endl;
			}
			abort();
		}

		if (old_nr_of_cells > 0) {
			const double
				order_of_accuracy_x
					= -log(norm_x / old_norm_x)
					/ log(double(nr_of_cells) / old_nr_of_cells),
				order_of_accuracy_y
					= -log(norm_y / old_norm_y)
					/ log(double(nr_of_cells) / old_nr_of_cells),
				order_of_accuracy_z
					= -log(norm_z / old_norm_z)
					/ log(double(nr_of_cells) / old_nr_of_cells);

			if (order_of_accuracy_x < 2.95) {
				if (grid_x.get_rank() == 0) {
					std::cerr << __FILE__ << ":" << __LINE__
						<< ": Order of accuracy from "
						<< old_nr_of_cells << " to " << nr_of_cells
						<< " is too low: " << order_of_accuracy_x
						<< std::endl;
				}
				abort();
			}
			if (order_of_accuracy_y < 2.95) {
				if (grid_y.get_rank() == 0) {
					std::cerr << __FILE__ << ":" << __LINE__
						<< ": Order of accuracy from "
						<< old_nr_of_cells << " to " << nr_of_cells
						<< " is too low: " << order_of_accuracy_y
						<< std::endl;
				}
				abort();
			}
			if (order_of_accuracy_z < 2.95) {
				if (grid_z.get_rank() == 0) {
					std::cerr << __FILE__ << ":" << __LINE__
						<< ": Order of accuracy from "
						<< old_nr_of_cells << " to " << nr_of_cells
						<< " is too low: " << order_of_accuracy_z
						<< std::endl;
				}
				abort();
			}
		}

		old_nr_of_cells = nr_of_cells;
		old_norm_x = norm_x;
		old_norm_y = norm_y;
		old_norm_z = norm_z;
	}

	MPI_Finalize();

	return EXIT_SUCCESS;
}
