/*
Tests vector field double curl calculation of PAMHD in 3d.

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


#include "array"
#include "cmath"
#include "cstdlib"
#include "iomanip"
#include "iostream"

#include "dccrg.hpp"
#include "dccrg_cartesian_geometry.hpp"
#include "divergence/remove.hpp"
#include "prettyprint.hpp"

using namespace std;


struct Vector {
	using data_type = std::array<double, 3>;
};

struct Result {
	using data_type = std::array<double, 3>;
};

using Cell = gensimcell::Cell<
	gensimcell::Always_Transfer,
	Vector,
	Result
>;
using Grid = dccrg::Dccrg<Cell, dccrg::Cartesian_Geometry>;

const auto Vector_Getter
	= [](Cell& cell_data) -> Vector::data_type& {
		return cell_data[Vector()];
	};
const auto Result_Getter
	= [](Cell& cell_data) -> Result::data_type& {
		return cell_data[Result()];
	};

// returns true on success
bool check(
	Grid& grid,
	const std::array<std::array<double, 3>, 7>& init_cond,
	const std::array<double, 3>& reference_result
) {
	auto cells = grid.get_cells();
	for (const auto& cell: cells) {
		auto* const cell_data = grid[cell];
		if (cell_data == nullptr) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< ": No data for cell " << cell
				<< std::endl;
			abort();
		}
		(*cell_data)[Result()] = {0, 0, 0};

		switch(cell) {
		case 5:
			(*cell_data)[Vector()] = init_cond[0];
			break;
		case 11:
			(*cell_data)[Vector()] = init_cond[1];
			break;
		case 13:
			(*cell_data)[Vector()] = init_cond[2];
			break;
		case 14:
			(*cell_data)[Vector()] = init_cond[3];
			break;
		case 15:
			(*cell_data)[Vector()] = init_cond[4];
			break;
		case 17:
			(*cell_data)[Vector()] = init_cond[5];
			break;
		case 23:
			(*cell_data)[Vector()] = init_cond[6];
			break;
		default:
			(*cell_data)[Vector()] = {0, 0, 0};
			break;
		}
	}
	grid.update_copies_of_remote_neighbors();

	if (grid.is_local(14)) {
		pamhd::divergence::get_curl_curl(
			{14},
			grid,
			Vector_Getter,
			Result_Getter
		);
	} else {
		pamhd::divergence::get_curl_curl(
			{},
			grid,
			Vector_Getter,
			Result_Getter
		);
	}

	int failure = 0;
	if (grid.is_local(14)) {
		const auto* const cell_data = grid[14];
		if (cell_data == nullptr) {
			std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
			abort();
		}
		if ((*cell_data)[Result()][0] != reference_result[0]) {
			std::cerr << __FILE__ << ":" << __LINE__ << std::setprecision(20)
				<< ": X, got " << (*cell_data)[Result()][0]
				<< " but should be " << reference_result[0]
				<< std::endl;
			failure++;
		}
		if ((*cell_data)[Result()][1] != reference_result[1]) {
			std::cerr << __FILE__ << ":" << __LINE__ << std::setprecision(20)
				<< ": Y, got " << (*cell_data)[Result()][1]
				<< " but should be " << reference_result[1]
				<< std::endl;
			failure++;
		}
		if ((*cell_data)[Result()][2] != reference_result[2]) {
			std::cerr << __FILE__ << ":" << __LINE__ << std::setprecision(20)
				<< ": Z, got " << (*cell_data)[Result()][2]
				<< " but should be " << reference_result[2]
				<< std::endl;
			failure++;
		}
	}

	int failure_global = 0;
	MPI_Allreduce(&failure, &failure_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	if (failure_global == 0) {
		return true;
	} else {
		return false;
	}
}


int main(int argc, char* argv[])
{
	using std::sqrt;

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

	std::array<uint64_t, 3> grid_size{3, 3, 3};
	const auto init_grid
		= [&](Grid& grid){
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
		};

	const std::array<double, 3> cell_length{1, 2, 3}, grid_start{3, 2, 1};
	const auto init_geom
		= [&](Grid& grid){
			dccrg::Cartesian_Geometry::Parameters geom_params;
			geom_params.start = grid_start;
			geom_params.level_0_cell_length = cell_length;
			if (not grid.set_geometry(geom_params)) {
				std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
				abort();
			}
		};

	// only Vx components
	{
		Grid grid;
		init_grid(grid);
		init_geom(grid); // dx = 1, dy = 2, dz = 3
		if (
			not check(
				grid,
				{{ // Vx, Vy, Vz
					{{7, 0, 0}}, // -z: cell 5
					{{6, 0, 0}}, // -y: cell 11
					{{5, 0, 0}}, // -x: cell 13
					{{4, 0, 0}}, // cell 14
					{{3, 0, 0}}, // +x: cell 15
					{{2, 0, 0}}, // +y: cell 17
					{{1, 0, 0}}  // +z: cell 23
				}},
				{
					// dydxVy - dydyVx - dzdzVx + dzdxVz
					-2*(-4/2.0/2 + 2/2.0/2/2 + 6/2.0/2/2) - 2*(-4/3.0/3 + 1/2.0/3/3 + 7/2.0/3/3),
					// dzdyVz - dzdzVy - dxdxVy + dxdyVx
					(3+2-6-5)/2.0/sqrt(1*1+2*2),
					// dxdzVx - dxdxVz - dydyVz + dydzVy
					(1+3-5-7)/2.0/sqrt(1*1+3*3)
				}
			)
		) {
			if (rank == 0) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " Vx failed"
					<< std::endl;
			}
			MPI_Finalize();
			return EXIT_FAILURE;
		}
	}
	// only Vy components
	{
		Grid grid;
		init_grid(grid);
		init_geom(grid); // dx = 1, dy = 2, dz = 3
		if (
			not check(
				grid,
				{{ // Vx, Vy, Vz
					{{0, 1, 0}}, // -z: cell 5
					{{0, 7, 0}}, // -y: cell 11
					{{0, 2, 0}}, // -x: cell 13
					{{0, 6, 0}}, // cell 14
					{{0, 3, 0}}, // +x: cell 15
					{{0, 5, 0}}, // +y: cell 17
					{{0, 4, 0}}  // +z: cell 23
				}},
				{
					// dydxVy - dydyVx - dzdzVx + dzdxVz
					(3+5-7-2)/2.0/sqrt(1*1+2*2),
					// dzdyVz - dzdzVy - dxdxVy + dxdyVx
					-2*(-6/3.0/3 + 4/2.0/3/3 + 1/2.0/3/3) - 2*(-6/1.0/1 + 3/2.0/1/1 + 2/2.0/1/1),
					// dxdzVx - dxdxVz - dydyVz + dydzVy
					(5+4-1-7)/2.0/sqrt(2*2+3*3)
				}
			)
		) {
			if (rank == 0) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " Vy failed"
					<< std::endl;
			}
			MPI_Finalize();
			return EXIT_FAILURE;
		}
	}
	// only Vz components
	{
		Grid grid;
		init_grid(grid);
		init_geom(grid); // dx = 1, dy = 2, dz = 3
		if (
			not check(
				grid,
				{{ // Vx, Vy, Vz
					{{0, 0, 1}}, // -z: cell 5
					{{0, 0, 2}}, // -y: cell 11
					{{0, 0, 3}}, // -x: cell 13
					{{0, 0, 4}}, // cell 14
					{{0, 0, 5}}, // +x: cell 15
					{{0, 0, 6}}, // +y: cell 17
					{{0, 0, 7}}  // +z: cell 23
				}},
				{
					// dydxVy - dydyVx - dzdzVx + dzdxVz
					(5+7-1-3)/2.0/sqrt(1*1+3*3),
					// dzdyVz - dzdzVy - dxdxVy + dxdyVx
					(7+6-1-2)/2.0/sqrt(2*2+3*3),
					// dxdzVx - dxdxVz - dydyVz + dydzVy
					-2*(-4/1.0/1 + 7/2.0/1/1 + 1/2.0/1/1) - 2*(-4/2.0/2 + 6/2.0/2/2 + 2/2.0/2/2)
				}
			)
		) {
			if (rank == 0) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " Vz failed"
					<< std::endl;
			}
			MPI_Finalize();
			return EXIT_FAILURE;
		}
	}
	// rotation in z plane 1
	{
		Grid grid;
		grid_size[2] = 1;
		init_grid(grid);
		init_geom(grid); // dx = 1, dy = 2, dz = 3
		if (
			not check(
				grid,
				{{ // Vx, Vy, Vz
					{{ 0,  0, 0}}, // -z: cell 5
					{{-1,  2, 0}}, // -y: cell 11
					{{ 1,  2, 0}}, // -x: cell 13
					{{ 0,  0, 0}}, // cell 14
					{{-1, -2, 0}}, // +x: cell 15
					{{ 1, -2, 0}}, // +y: cell 17
					{{ 0,  0, 0}}  // +z: cell 23
				}},
				{
					// dydxVy - dydyVx - dzdzVx + dzdxVz
					(-2-2-2-2)/2.0/sqrt(1*1+2*2) - 2*(0+1/2.0/2/2-1/2.0/2/2),
					// dzdyVz - dzdzVy - dxdxVy + dxdyVx
					-2*(0-2/2.0/2/2+2/2.0/2/2) + (1-1-1+1)/2.0/sqrt(1*1+2*2),
					// dxdzVx - dxdxVz - dydyVz + dydzVy
					0
				}
			)
		) {
			if (rank == 0) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " rotation 1 failed"
					<< std::endl;
			}
			MPI_Finalize();
			return EXIT_FAILURE;
		}
	}

	// rotation in z plane 2
	{
		Grid grid;
		grid_size[2] = 1;
		init_grid(grid);
		init_geom(grid); // dx = 1, dy = 2, dz = 3
		if (
			not check(
				grid,
				{{ // Vx, Vy, Vz
					{{ 0,  0, 0}}, // -z: cell 5
					{{ 1,  2, 0}}, // -y: cell 11
					{{ 1, -2, 0}}, // -x: cell 13
					{{ 0,  0, 0}}, // cell 14
					{{-1,  2, 0}}, // +x: cell 15
					{{-1, -2, 0}}, // +y: cell 17
					{{ 0,  0, 0}}  // +z: cell 23
				}},
				{
					// dydxVy - dydyVx - dzdzVx + dzdxVz
					-(-2-2-2-2)/2.0/sqrt(1*1+2*2) + 2*(0+1/2.0/2/2-1/2.0/2/2),
					// dzdyVz - dzdzVy - dxdxVy + dxdyVx
					2*(0-2/2.0/2/2+2/2.0/2/2) - (1-1-1+1)/2.0/sqrt(1*1+2*2),
					// dxdzVx - dxdxVz - dydyVz + dydzVy
					0
				}
			)
		) {
			if (rank == 0) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " rotation 1 failed"
					<< std::endl;
			}
			MPI_Finalize();
			return EXIT_FAILURE;
		}
	}

	MPI_Finalize();

	return EXIT_SUCCESS;
}
