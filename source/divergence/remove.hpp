/*
Functions for working with divergence of vector field.

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

#ifndef PAMHD_DIVERGENCE_REMOVE_HPP
#define PAMHD_DIVERGENCE_REMOVE_HPP

#include "cmath"
#include "limits"
#include "tuple"
#include "vector"

#include "dccrg.hpp"
#include "dccrg_cartesian_geometry.hpp"
#include "gensimcell.hpp"
#include "mpi.h"
#include "tests/poisson/poisson_solve.hpp" // part of dccrg


namespace pamhd {
namespace divergence {


/*!
Calculates divergence of vector variable in given cells.

Returns the total divergence i.e. sum of absolute
divergence in given cells of all processes divided
by number of cells on all processes in which divergence
was calculated.

Vector variable must have a defined value in face
neighbors of given cells. In given cells contribution
to divergence is 0 from dimensions with at least one
missing neighbor.

Vector must be an object that, when given a reference
to the data of one grid cell, returns a reference to
the data from which divergnce should be calculated.
Similarly Divergence should return a reference to data
in which calculated divergence should be stored.

Assumes Grid_T API is compatible with dccrg.

Example:

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

using Grid_T = dccrg::Dccrg<Cell, Cartesian_Geometry>;

pamhd::divergence::get_divergence(
	std::vector<uint64_t>(),
	grid,
	[](Cell& cell_data)->Vector_Field::data_type& {
		return cell_data[Vector_Field()];
	},
	[](Cell& cell_data)->Divergence::data_type& {
		return cell_data[Divergence()];
	}
);
*/
template <
	class Grid_T,
	class Vector_Getter,
	class Divergence_Getter
> double get_divergence(
	const std::vector<uint64_t>& cells,
	Grid_T& grid,
	Vector_Getter Vector,
	Divergence_Getter Divergence
) {
	double local_divergence = 0, global_divergence = 0;
	uint64_t local_calculated_cells = 0, global_calculated_cells = 0;
	for (const auto& cell: cells) {
		// get distance between neighbors in same dimension
		std::array<double, 3>
			// distance from current cell on neg and pos side (> 0)
			neigh_neg_dist{{0, 0, 0}},
			neigh_pos_dist{{0, 0, 0}};

		// number of neighbors in each dimension
		std::array<size_t, 3> nr_neighbors{{0, 0, 0}};

		const auto cell_length = grid.geometry.get_length(cell);
		const auto face_neighbors_of = grid.get_face_neighbors_of(cell);
		for (const auto& item: face_neighbors_of) {
			const auto neighbor = item.first;
			const auto direction = item.second;
			if (direction == 0 or std::abs(direction) > 3) {
				std::cerr << __FILE__ << "(" << __LINE__<< ")" << std::endl;
				abort();
			}
			const size_t dim = std::abs(direction) - 1;

			nr_neighbors[dim]++;

			const auto neighbor_length
				= grid.geometry.get_length(neighbor);

			const double distance
				= (cell_length[dim] + neighbor_length[dim]) / 2.0;

			if (direction < 0) {
				neigh_neg_dist[dim] = distance;
			} else {
				neigh_pos_dist[dim] = distance;
			}
		}

		bool have_enough_neighbors = false;
		for (auto dim = 0; dim < 3; dim++) {
			if (
				nr_neighbors[dim] == 2
				or nr_neighbors[dim] == 5
				or nr_neighbors[dim] == 8
			) {
				have_enough_neighbors = true;
				break;
			}
		}

		if (not have_enough_neighbors) {
			continue;
		}
		local_calculated_cells++;

		auto* const cell_data = grid[cell];
		if (cell_data == nullptr) {
			std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
				<< "No data for simulation cell " << cell
				<< std::endl;
			abort();
		}

		auto& div = Divergence(*cell_data);
		div = 0;

		const auto cell_ref_lvl = grid.get_refinement_level(cell);

		for (const auto& item: face_neighbors_of) {
			const auto neighbor = item.first;
			const auto direction = item.second;
			const size_t dim = std::abs(direction) - 1;

			if (
				nr_neighbors[dim] != 2
				and nr_neighbors[dim] != 5
				and nr_neighbors[dim] != 8
			) {
				continue;
			}

			double multiplier = 1 / (neigh_pos_dist[dim] + neigh_neg_dist[dim]);
			if (direction < 0) {
				multiplier *= -1;
			}

			const auto neigh_ref_lvl = grid.get_refinement_level(neighbor);
			if (neigh_ref_lvl > cell_ref_lvl) {
				multiplier /= 4;
			}

			auto* const neighbor_data = grid[neighbor];
			if (neighbor_data == nullptr) {
				std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
					<< "No data for neighbor " << neighbor
					<< " of cell " << cell
					<< std::endl;
				abort();
			}

			div += multiplier * Vector(*neighbor_data)[dim];
		}

		local_divergence += std::fabs(div);
	}

	MPI_Comm comm = grid.get_communicator();
	MPI_Allreduce(
		&local_divergence,
		&global_divergence,
		1,
		MPI_DOUBLE,
		MPI_SUM,
		comm
	);
	MPI_Allreduce(
		&local_calculated_cells,
		&global_calculated_cells,
		1,
		MPI_UINT64_T,
		MPI_SUM,
		comm
	);
	MPI_Comm_free(&comm);

	return global_divergence / global_calculated_cells;
}


/*!
Calculates gradient of scalar variable in given cells.

See get_divergence() for more info on arguments, etc.
*/
template <
	class Grid_T,
	class Scalar_Getter,
	class Gradient_Getter
> void get_gradient(
	const std::vector<uint64_t>& cells,
	Grid_T& grid,
	Scalar_Getter Scalar,
	Gradient_Getter Gradient
) {
	for (const auto& cell: cells) {
		// get distance between neighbors in same dimension
		std::array<double, 3>
			// distance from current cell on neg and pos side (> 0)
			neigh_neg_dist{{0, 0, 0}},
			neigh_pos_dist{{0, 0, 0}};

		// number of neighbors in each dimension
		std::array<size_t, 3> nr_neighbors{{0, 0, 0}};

		const auto cell_length = grid.geometry.get_length(cell);
		const auto face_neighbors_of = grid.get_face_neighbors_of(cell);
		for (const auto& item: face_neighbors_of) {
			const auto neighbor = item.first;
			const auto direction = item.second;
			if (direction == 0 or std::abs(direction) > 3) {
				std::cerr << __FILE__ << "(" << __LINE__<< ")" << std::endl;
				abort();
			}
			const size_t dim = std::abs(direction) - 1;

			nr_neighbors[dim]++;

			const auto neighbor_length
				= grid.geometry.get_length(neighbor);

			const double distance
				= (cell_length[dim] + neighbor_length[dim]) / 2.0;

			if (direction < 0) {
				neigh_neg_dist[dim] = distance;
			} else {
				neigh_pos_dist[dim] = distance;
			}
		}

		bool have_enough_neighbors = false;
		for (auto dim = 0; dim < 3; dim++) {
			if (
				nr_neighbors[dim] == 2
				or nr_neighbors[dim] == 5
				or nr_neighbors[dim] == 8
			) {
				have_enough_neighbors = true;
				break;
			}
		}

		auto* const cell_data = grid[cell];
		if (cell_data == nullptr) {
			std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
				<< "No data for simulation cell " << cell
				<< std::endl;
			abort();
		}

		auto& gradient = Gradient(*cell_data);
		gradient[0] =
		gradient[1] =
		gradient[2] = 0;

		if (not have_enough_neighbors) {
			continue;
		}

		const auto cell_ref_lvl = grid.get_refinement_level(cell);

		for (const auto& item: face_neighbors_of) {
			const auto neighbor = item.first;
			const auto direction = item.second;
			const size_t dim = std::abs(direction) - 1;

			if (
				nr_neighbors[dim] != 2
				and nr_neighbors[dim] != 5
				and nr_neighbors[dim] != 8
			) {
				continue;
			}

			double multiplier = 1 / (neigh_pos_dist[dim] + neigh_neg_dist[dim]);
			if (direction < 0) {
				multiplier *= -1;
			}

			const auto neigh_ref_lvl = grid.get_refinement_level(neighbor);
			if (neigh_ref_lvl > cell_ref_lvl) {
				multiplier /= 4;
			}

			auto* const neighbor_data = grid[neighbor];
			if (neighbor_data == nullptr) {
				std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
					<< "No data for neighbor " << neighbor
					<< " of cell " << cell
					<< std::endl;
				abort();
			}

			gradient[dim] += multiplier * Scalar(*neighbor_data);
		}
	}
}


/*!
Calculates curl of vector variable in given cells.

See get_divergence() for more info on arguments, etc.
*/
template <
	class Grid_T,
	class Vector_Getter,
	class Curl_Getter
> void get_curl(
	const std::vector<uint64_t>& cells,
	Grid_T& grid,
	Vector_Getter Vector,
	Curl_Getter Curl
) {
	for (const auto& cell: cells) {
		if (grid.get_refinement_level(cell) != 0) {
			std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
				<< "Adaptive mesh refinement not supported"
				<< std::endl;
			abort();
		}

		const auto cell_length = grid.geometry.get_length(cell);

		const auto face_neighbors_of = grid.get_face_neighbors_of(cell);

		// get distance between neighbors in same dimension
		std::array<double, 3>
			// distance from current cell on neg and pos side (> 0)
			neigh_neg_dist{{0, 0, 0}},
			neigh_pos_dist{{0, 0, 0}};

		// number of neighbors in each dimension
		std::array<size_t, 3> nr_neighbors{{0, 0, 0}};

		for (const auto& item: face_neighbors_of) {
			const auto neighbor = item.first;
			const auto direction = item.second;
			if (direction == 0 or std::abs(direction) > 3) {
				std::cerr << __FILE__ << "(" << __LINE__<< ")" << std::endl;
				abort();
			}
			const size_t dim = std::abs(direction) - 1;

			nr_neighbors[dim]++;

			const auto neighbor_length
				= grid.geometry.get_length(neighbor);

			const double distance
				= (cell_length[dim] + neighbor_length[dim]) / 2.0;

			if (direction < 0) {
				neigh_neg_dist[dim] = distance;
			} else {
				neigh_pos_dist[dim] = distance;
			}
		}

		bool have_enough_neighbors = false;
		for (auto dim = 0; dim < 3; dim++) {
			if (nr_neighbors[dim] == 2) {
				have_enough_neighbors = true;
			}
		}

		auto* const cell_data = grid[cell];
		if (cell_data == nullptr) {
			std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
				<< "No data for simulation cell " << cell
				<< std::endl;
			abort();
		}

		auto
			&vec = Vector(*cell_data),
			&curl = Curl(*cell_data);

		curl[0] =
		curl[1] =
		curl[2] = 0;

		if (not have_enough_neighbors) {
			continue;
		}

		/*
		curl_0 = diff_vec_2 / diff_pos_1 - diff_vec_1 / diff_pos_2
		curl_1 = diff_vec_0 / diff_pos_2 - diff_vec_2 / diff_pos_0
		curl_2 = diff_vec_1 / diff_pos_0 - diff_vec_0 / diff_pos_1
		*/
		for (auto dim0 = 0; dim0 < 3; dim0++) {

			const auto
				dim1 = (dim0 + 1) % 3,
				dim2 = (dim0 + 2) % 3;

			// zero in dimensions with missing neighbor(s)
			if (nr_neighbors[dim1] == 2) {
				curl[dim0] += vec[dim2] * (neigh_pos_dist[dim1] - neigh_neg_dist[dim1]);
			}
			if (nr_neighbors[dim2] == 2) {
				curl[dim0] -= vec[dim1] * (neigh_pos_dist[dim2] - neigh_neg_dist[dim2]);
			}
		}

		for (const auto& item: face_neighbors_of) {
			const auto neighbor = item.first;
			const auto direction = item.second;
			const size_t dim0 = std::abs(direction) - 1;

			if (nr_neighbors[dim0] != 2) {
				continue;
			}

			const auto
				dim1 = (dim0 + 1) % 3,
				dim2 = (dim0 + 2) % 3;

			double multiplier = 0;
			if (direction < 0) {
				multiplier = -neigh_pos_dist[dim0] / neigh_neg_dist[dim0];
			} else {
				multiplier = neigh_neg_dist[dim0] / neigh_pos_dist[dim0];
			}
			multiplier /= (neigh_pos_dist[dim0] + neigh_neg_dist[dim0]);

			auto* const neighbor_data = grid[neighbor];
			if (neighbor_data == nullptr) {
				std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
					<< "No data for neighbor " << neighbor
					<< " of cell " << cell
					<< std::endl;
				abort();
			}
			auto& neigh_vec = Vector(*neighbor_data);

			curl[dim2] += multiplier * neigh_vec[dim1];
			curl[dim1] -= multiplier * neigh_vec[dim2];
		}
	}
}


/*!
Removes divergence of a vector variable from given cells.

Returns total divergence of vector variable
(from get_divergence()) before removing divergence.

Solves phi from div(grad(phi)) = div(vector) and assigns
vector = vector - grad(phi) after which div(vector) -> 0.

Vector variable must have been updated between processes
before calling this function.

Vector, Divergence and Gradient should return a reference
to the variable's data.

See get_divergence() for more info on arguments, etc.

The transfer of Vector and Divergence variables' data
must have been switched on in cells of simulation grid
before calling this function.

The arguments max_iterations to verbose are given
directly to the constructor of dccrg Poisson equation
solver class:
https://gitorious.org/dccrg/dccrg/source/master:tests/poisson/poisson_solve.hpp

Poisson solver is applied to previous results retries
number of times in a row, i.e. solver is used once
if retries == 0.
*/
template <
	class Cell_T,
	class Geometry_T,
	class Vector_Getter,
	class Divergence_Getter,
	class Gradient_Getter
	// TODO: add possibility to reuse solution from previous call
> double remove(
	const std::vector<uint64_t>& cells,
	const std::vector<uint64_t>& boundary_cells,
	const std::vector<uint64_t>& skip_cells,
	dccrg::Dccrg<Cell_T, Geometry_T>& simulation_grid,
	Vector_Getter Vector,
	Divergence_Getter Divergence,
	Gradient_Getter Gradient,
	const unsigned int max_iterations = 1000,
	const unsigned int min_iterations = 0,
	const double stop_residual = 1e-15,
	const double p_of_norm = 2,
	const double stop_after_residual_increase = 10,
	const unsigned int retries = 0,
	const bool verbose = false
) {
	std::vector<uint64_t> solve_cells;
	solve_cells.reserve(cells.size() + boundary_cells.size());

	solve_cells.insert(
		solve_cells.end(),
		cells.cbegin(),
		cells.cend()
	);
	solve_cells.insert(
		solve_cells.end(),
		boundary_cells.cbegin(),
		boundary_cells.cend()
	);

	/*
	Prepare solution grid and source term
	*/

	simulation_grid.update_copies_of_remote_neighbors();

	// rhs for poisson solver = div(vec)
	const double ret_val = get_divergence(
		cells,
		simulation_grid,
		Vector,
		Divergence
	);
	// zero divergence in boundary cells
	for (const auto& cell: boundary_cells) {
		auto* const cell_data = simulation_grid[cell];
		if (cell_data == nullptr) { abort(); }
		Divergence(*cell_data) = 0;
	}

	dccrg::Dccrg<Poisson_Cell, Geometry_T> poisson_grid(simulation_grid);

	// transfer rhs to poisson grid
	for (const auto& cell: solve_cells) {
		auto* const poisson_data = poisson_grid[cell];
		if (poisson_data == nullptr) {
			std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
				<< "No data for poisson cell " << cell
				<< std::endl;
			abort();
		}

		auto* const simulation_data = simulation_grid[cell];
		if (simulation_data == nullptr) {
			std::cerr << __FILE__ << "(" << __LINE__<< "): "
				<< "No data for simulation cell " << cell
				<< std::endl;
			abort();
		}

		poisson_data->solution = 0;
		poisson_data->rhs = Divergence(*simulation_data);
	}

	Poisson_Solve solver(
		max_iterations,
		min_iterations,
		stop_residual,
		p_of_norm,
		stop_after_residual_increase,
		verbose
	);

	// solve phi in div(grad(phi)) = rhs
	size_t iters = 0;
	while (iters <= retries) {
		solver.solve(
			cells,
			poisson_grid,
			skip_cells,
			[iters](){
				if (iters == 0) {
					return false;
				} else {
					return true;
				}
			}()
		);
		iters++;
	}

	/*
	Remove divergence with Vec = Vec - grad(phi)
	*/

	// store phi (solution) in divergence variable
	for (const auto& cell: cells) {
		auto* const poisson_data = poisson_grid[cell];
		if (poisson_data == nullptr) {
			std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
				<< "No data for poisson cell " << cell
				<< std::endl;
			abort();
		}

		auto* const simulation_data = simulation_grid[cell];
		if (simulation_data == nullptr) {
			std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
				<< "No data for simulation cell " << cell
				<< std::endl;
			abort();
		}

		Divergence(*simulation_data) = poisson_data->solution;
	}

	simulation_grid.update_copies_of_remote_neighbors();

	get_gradient(
		cells,
		simulation_grid,
		Divergence,
		Gradient
	);

	for (const auto& cell: cells) {
		auto* const cell_data = simulation_grid[cell];
		if (cell_data == nullptr) {
			std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
				<< "No data for simulation cell " << cell
				<< std::endl;
			abort();
		}

		for (size_t dim = 0; dim < 3; dim++) {
			Vector(*cell_data)[dim] -= Gradient(*cell_data)[dim];
		}
	}

	// clean up divergence variable
	for (const auto& cell: boundary_cells) {
		auto* const simulation_data = simulation_grid[cell];
		if (simulation_data == nullptr) {
			std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
				<< "No data for simulation cell " << cell
				<< std::endl;
			abort();
		}

		Divergence(*simulation_data) = 0;
	}

	return ret_val;
}


/*!
Calculates curl(curl(vector)) of vector variable in given cells.

See get_divergence() for more info on arguments, etc.
*/
template <
	class Cell,
	class Vector_Getter,
	class Result_Getter
> void get_curl_curl(
	const std::vector<uint64_t>& cells,
	dccrg::Dccrg<Cell, dccrg::Cartesian_Geometry>& grid,
	Vector_Getter Vector,
	Result_Getter Result
) {
	for (const auto& cell: cells) {
		if (grid.get_refinement_level(cell) != 0) {
			std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
				<< "Adaptive mesh refinement not supported"
				<< std::endl;
			abort();
		}

		const auto cell_length = grid.geometry.get_length(cell);

		const auto face_neighbors_of = grid.get_face_neighbors_of(cell);

		// get distance between neighbors in same dimension
		std::array<double, 3>
			// distance from current cell on neg and pos side (> 0)
			neigh_neg_dist{{0, 0, 0}},
			neigh_pos_dist{{0, 0, 0}};

		// number of neighbors in each dimension
		std::array<size_t, 3> nr_neighbors{{0, 0, 0}};
		// pointers to neighbor data, -x, +x, -y, +y, -z, +z
		std::array<Cell*, 6> neighbor_datas{nullptr};

		for (const auto& item: face_neighbors_of) {
			const auto neighbor = item.first;
			const auto direction = item.second;
			if (direction == 0 or std::abs(direction) > 3) {
				std::cerr << __FILE__ << "(" << __LINE__<< ")" << std::endl;
				abort();
			}
			const size_t dim = std::abs(direction) - 1;

			nr_neighbors[dim]++;

			const auto neighbor_length
				= grid.geometry.get_length(neighbor);

			const double distance
				= (cell_length[dim] + neighbor_length[dim]) / 2.0;

			if (direction < 0) {
				neigh_neg_dist[dim] = distance;
			} else {
				neigh_pos_dist[dim] = distance;
			}

			auto* const neighbor_data = grid[neighbor];
			if (neighbor_data == nullptr) {
				std::cerr << __FILE__ << "(" << __LINE__<< ")" << std::endl;
				abort();
			}
			size_t index = dim * 2;
			if (direction > 0) {
				index++;
			}
			neighbor_datas[index] = neighbor_data;
		}

		bool have_enough_neighbors = false;
		for (auto dim = 0; dim < 3; dim++) {
			if (nr_neighbors[dim] == 2) {
				have_enough_neighbors = true;
			}
		}

		auto* const cell_data = grid[cell];
		if (cell_data == nullptr) {
			std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
				<< "No data for simulation cell " << cell
				<< std::endl;
			abort();
		}

		auto
			&vec = Vector(*cell_data),
			&result = Result(*cell_data);

		result[0] =
		result[1] =
		result[2] = 0;

		if (not have_enough_neighbors) {
			continue;
		}

		/*
		result_x = dydxVy - dydyVx - dzdzVx + dzdxVz
		result_y = dzdyVz - dzdzVy - dxdxVy + dxdyVx
		result_z = dxdzVx - dxdxVz - dydyVz + dydzVy

		dxdyA = dydxA = [A(0, 1) + A(1, 0) - A(-1, 0) - A(0, -1)] / 2 / sqrt(dx*dx+dy*dy)
		dxdxA = 2*[-A(0, 0) + A(1, 0)/2 + A(-1, 0)/2]/dx/dx
		...
		*/

		// result_x
		if (nr_neighbors[1] == 2) {
			// -dydyVx
			result[0] -= 2 * (
				-vec[0] / neigh_neg_dist[1] / neigh_pos_dist[1]
				+ Vector(*neighbor_datas[3])[0] / neigh_pos_dist[1] / (neigh_pos_dist[1] + neigh_neg_dist[1])
				+ Vector(*neighbor_datas[2])[0] / neigh_neg_dist[1] / (neigh_pos_dist[1] + neigh_neg_dist[1])
			);
			// dydxVy
			if (nr_neighbors[0] == 2) {
				result[0] += 0.5 * (
					Vector(*neighbor_datas[3])[1]
					+ Vector(*neighbor_datas[1])[1]
					- Vector(*neighbor_datas[0])[1]
					- Vector(*neighbor_datas[2])[1]
				) / std::sqrt(
					neigh_pos_dist[0]*neigh_pos_dist[0]
					+ neigh_pos_dist[1]*neigh_pos_dist[1]
				);
			}
		}
		if (nr_neighbors[2] == 2) {
			// -dzdzVx
			result[0] -= 2 * (
				-vec[0] / neigh_neg_dist[2] / neigh_pos_dist[2]
				+ Vector(*neighbor_datas[5])[0] / neigh_pos_dist[2] / (neigh_pos_dist[2] + neigh_neg_dist[2])
				+ Vector(*neighbor_datas[4])[0] / neigh_neg_dist[2] / (neigh_pos_dist[2] + neigh_neg_dist[2])
			);
			// dzdxVz
			if (nr_neighbors[0] == 2) {
				result[0] += 0.5 * (
					Vector(*neighbor_datas[5])[2]
					+ Vector(*neighbor_datas[1])[2]
					- Vector(*neighbor_datas[0])[2]
					- Vector(*neighbor_datas[4])[2]
				) / std::sqrt(
					neigh_pos_dist[0]*neigh_pos_dist[0]
					+ neigh_pos_dist[2]*neigh_pos_dist[2]
				);
			}
		}
		// result_y
		if (nr_neighbors[2] == 2) {
			// -dzdzVy
			result[1] -= 2 * (
				-vec[1] / neigh_neg_dist[2] / neigh_pos_dist[2]
				+ Vector(*neighbor_datas[5])[1] / neigh_pos_dist[2] / (neigh_pos_dist[2] + neigh_neg_dist[2])
				+ Vector(*neighbor_datas[4])[1] / neigh_neg_dist[2] / (neigh_pos_dist[2] + neigh_neg_dist[2])
			);
			// dzdyVz
			if (nr_neighbors[1] == 2) {
				result[1] += 0.5 * (
					Vector(*neighbor_datas[5])[2]
					+ Vector(*neighbor_datas[3])[2]
					- Vector(*neighbor_datas[2])[2]
					- Vector(*neighbor_datas[4])[2]
				) / std::sqrt(
					neigh_pos_dist[1]*neigh_pos_dist[1]
					+ neigh_pos_dist[2]*neigh_pos_dist[2]
				);
			}
		}
		if (nr_neighbors[0] == 2) {
			// -dxdxVy
			result[1] -= 2 * (
				-vec[1] / neigh_neg_dist[0] / neigh_pos_dist[0]
				+ Vector(*neighbor_datas[1])[1] / neigh_pos_dist[0] / (neigh_pos_dist[0] + neigh_neg_dist[0])
				+ Vector(*neighbor_datas[0])[1] / neigh_neg_dist[0] / (neigh_pos_dist[0] + neigh_neg_dist[0])
			);
			// dydxVx
			if (nr_neighbors[1] == 2) {
				result[1] += 0.5 * (
					Vector(*neighbor_datas[3])[0]
					+ Vector(*neighbor_datas[1])[0]
					- Vector(*neighbor_datas[0])[0]
					- Vector(*neighbor_datas[2])[0]
				) / std::sqrt(
					neigh_pos_dist[0]*neigh_pos_dist[0]
					+ neigh_pos_dist[1]*neigh_pos_dist[1]
				);
			}
		}
		// result_z
		if (nr_neighbors[0] == 2) {
			// -dxdxVz
			result[2] -= 2 * (
				-vec[2] / neigh_neg_dist[0] / neigh_pos_dist[0]
				+ Vector(*neighbor_datas[1])[2] / neigh_pos_dist[0] / (neigh_pos_dist[0] + neigh_neg_dist[0])
				+ Vector(*neighbor_datas[0])[2] / neigh_neg_dist[0] / (neigh_pos_dist[0] + neigh_neg_dist[0])
			);
			// dxdzVx
			if (nr_neighbors[2] == 2) {
				result[2] += 0.5 * (
					Vector(*neighbor_datas[5])[0]
					+ Vector(*neighbor_datas[1])[0]
					- Vector(*neighbor_datas[0])[0]
					- Vector(*neighbor_datas[4])[0]
				) / std::sqrt(
					neigh_pos_dist[0]*neigh_pos_dist[0]
					+ neigh_pos_dist[2]*neigh_pos_dist[2]
				);
			}
		}
		if (nr_neighbors[1] == 2) {
			// -dydyVz
			result[2] -= 2 * (
				-vec[2] / neigh_neg_dist[1] / neigh_pos_dist[1]
				+ Vector(*neighbor_datas[3])[2] / neigh_pos_dist[1] / (neigh_pos_dist[1] + neigh_neg_dist[1])
				+ Vector(*neighbor_datas[2])[2] / neigh_neg_dist[1] / (neigh_pos_dist[1] + neigh_neg_dist[1])
			);
			// dydzVy
			if (nr_neighbors[2] == 2) {
				result[2] += 0.5 * (
					Vector(*neighbor_datas[5])[1]
					+ Vector(*neighbor_datas[3])[1]
					- Vector(*neighbor_datas[2])[1]
					- Vector(*neighbor_datas[4])[1]
				) / std::sqrt(
					neigh_pos_dist[1]*neigh_pos_dist[1]
					+ neigh_pos_dist[2]*neigh_pos_dist[2]
				);
			}
		}
	}
}


}} // namespaces

#endif // ifndef PAMHD_DIVERGENCE_REMOVE_HPP
