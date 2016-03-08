/*
Tests vector field divergence removal of PAMHD in 1d.

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
#include "gensimcell.hpp"

#include "divergence/remove.hpp"


int Poisson_Cell::transfer_switch = Poisson_Cell::INIT;


double function(const double x)
{
	return 1 + 0.1 * std::sin(x);
}

double div_removed_function()
{
	return 1;
}


struct Vector_Field {
	using data_type = std::array<double, 3>;
};

struct Divergence {
	using data_type = double;
};

struct Gradient {
	using data_type = std::array<double, 3>;
};

using Cell = gensimcell::Cell<
	gensimcell::Always_Transfer,
	Vector_Field,
	Divergence,
	Gradient
>;


/*!
Returns maximum norm.
*/
template<class Grid_T> double get_max_norm(
	const std::vector<uint64_t>& cells,
	const Grid_T& grid,
	const size_t dimension
) {
	double local_norm = 0, global_norm = 0;
	for (const auto& cell: cells) {
		const auto* const cell_data = grid[cell];
		if (cell_data == nullptr) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< ": No data for cell " << cell
				<< std::endl;
			abort();
		}

		local_norm = std::max(
			local_norm,
			std::fabs((*cell_data)[Vector_Field()][dimension] - div_removed_function())
		);
	}

	MPI_Comm comm = grid.get_communicator();
	MPI_Allreduce(&local_norm, &global_norm, 1, MPI_DOUBLE, MPI_MAX, comm);
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

	std::array<double, 3> old_norms{{
		std::numeric_limits<double>::max(),
		std::numeric_limits<double>::max(),
		std::numeric_limits<double>::max()
	}};

	size_t old_nr_of_cells = 0;
	for (size_t nr_of_cells = 64; nr_of_cells <= 512; nr_of_cells *= 2) {

		std::array<dccrg::Dccrg<Cell, dccrg::Cartesian_Geometry>, 3> grids;

		const std::array<std::array<uint64_t, 3>, 3> grid_sizes{{
			{{nr_of_cells + 2, 1, 1}},
			{{1, nr_of_cells + 2, 1}},
			{{1, 1, nr_of_cells + 2}}
		}};

		const int max_ref_lvl = 1;
		for (size_t dim = 0; dim < 3; dim++) {
			if (
				not grids[dim].initialize(
					grid_sizes[dim],
					comm,
					"RANDOM",
					0,
					max_ref_lvl,
					false, false, false
				)
			) {
				std::cerr << __FILE__ << ":" << __LINE__ << " " << dim << std::endl;
				abort();
			}
		}

		const std::array<std::array<double, 3>, 3>
			cell_lengths{{
				{{2 * M_PI / (grid_sizes[0][0] - 2), 1, 1}},
				{{1, 2 * M_PI / (grid_sizes[1][1] - 2), 1}},
				{{1, 1, 2 * M_PI / (grid_sizes[2][2] - 2)}}
			}},
			grid_starts{{
				{{-cell_lengths[0][0], 0, 0}},
				{{0, -cell_lengths[1][1], 0}},
				{{0, 0, -cell_lengths[2][2]}}
			}};

		std::array<dccrg::Cartesian_Geometry::Parameters, 3> geom_params;
		for (size_t dim = 0; dim < 3; dim++) {
			geom_params[dim].start = grid_starts[dim];
			geom_params[dim].level_0_cell_length = cell_lengths[dim];

			if (not grids[dim].set_geometry(geom_params[dim])) {
				std::cerr << __FILE__ << ":" << __LINE__ << " " << dim << std::endl;
				abort();
			}
		}

		for (size_t dim = 0; dim < 3; dim++) {
			for (int i = 0; i < max_ref_lvl; i++) {
				for (const auto& cell: grids[dim].get_cells()) {
					const auto center = grids[dim].geometry.get_center(cell);
					if (
						(center[dim] > 2.0*M_PI*1/8 and center[dim] < 2.0*M_PI*3/8)
						or (center[dim] > 2.0*M_PI*5/8 and center[dim] < 2.0*M_PI*7/8)
					) {
						grids[dim].refine_completely(cell);
					}
				}
				grids[dim].stop_refining();
			}
		}

		std::array<std::vector<uint64_t>, 3>
			solve_cells,
			boundary_cells;
		for (size_t dim = 0; dim < 3; dim++) {
			for (const auto& cell: grids[dim].get_cells()) {
				const auto center = grids[dim].geometry.get_center(cell);
				if (center[dim] > 0 and center[dim] < 2 * M_PI) {
					solve_cells[dim].push_back(cell);
				} else {
					boundary_cells[dim].push_back(cell);
				}

				auto* const cell_data = grids[dim][cell];
				if (cell_data == nullptr) {
					std::cerr << __FILE__ << ":" << __LINE__ << " " << dim << std::endl;
					abort();
				}
				auto& vec = (*cell_data)[Vector_Field()];
				vec[0] = vec[1] = vec[2] = 0;
				vec[dim] = function(center[dim]);
			}
			grids[dim].update_copies_of_remote_neighbors();
		}

		auto Vector_Getter = [](Cell& cell_data) -> Vector_Field::data_type& {
			return cell_data[Vector_Field()];
		};
		auto Divergence_Getter = [](Cell& cell_data) -> Divergence::data_type& {
			return cell_data[Divergence()];
		};
		auto Gradient_Getter = [](Cell& cell_data) -> Gradient::data_type& {
			return cell_data[Gradient()];
		};

		for (size_t dim = 0; dim < 3; dim++) {
			pamhd::divergence::remove(
				solve_cells[dim],
				boundary_cells[dim],
				{},
				grids[dim],
				Vector_Getter,
				Divergence_Getter,
				Gradient_Getter,
				1000, 0, 1e-15, 2, 100, 5, false, false
			);
		}

		const std::array<double, 3> norms{{
			get_max_norm(solve_cells[0], grids[0], 0),
			get_max_norm(solve_cells[1], grids[1], 1),
			get_max_norm(solve_cells[2], grids[2], 2)
		}};

		for (size_t dim = 0; dim < 3; dim++) {
			if (old_nr_of_cells > 1000 and norms[dim] > old_norms[dim]) {
				if (grids[dim].get_rank() == 0) {
					std::cerr << __FILE__ << ":" << __LINE__
						<< " dim " << dim
						<< ": max norm with " << solve_cells[dim].size()
						<< " cells " << norms[dim]
						<< " is larger than with " << old_nr_of_cells
						<< " cells " << old_norms[dim]
						<< std::endl;
				}
				abort();
			}
		}

		for (size_t dim = 0; dim < 3; dim++) {
			if (old_nr_of_cells > 0) {
				const double order_of_accuracy
					= -log(norms[dim] / old_norms[dim])
					/ log(double(solve_cells[dim].size()) / old_nr_of_cells);

				if (old_nr_of_cells > 1000 and order_of_accuracy < 1) {
					if (grids[dim].get_rank() == 0) {
						std::cerr << __FILE__ << ":" << __LINE__
							<< " dim " << dim
							<< ": order of accuracy from "
							<< old_nr_of_cells << " to " << solve_cells[dim].size()
							<< " is too low: " << order_of_accuracy
							<< std::endl;
					}
					abort();
				}
			}
		}

		old_nr_of_cells = solve_cells[0].size();
		old_norms = norms;
	}

	MPI_Finalize();

	return EXIT_SUCCESS;
}
