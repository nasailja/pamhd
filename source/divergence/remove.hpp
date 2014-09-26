/*
Functions for working with divergence of vector field.

Copyright 2014 Ilja Honkonen
All rights reserved.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef PAMHD_DIVERGENCE_REMOVE_HPP
#define PAMHD_DIVERGENCE_REMOVE_HPP


#include "cmath"
#include "limits"
#include "tuple"
#include "vector"

#include "dccrg.hpp"
#include "gensimcell.hpp"
#include "mpi.h"
#include "tests/poisson/poisson_solve.hpp" // part of dccrg


namespace pamhd {
namespace divergence {


/*!
Calculates divergence of vector variable in given cells.

Returns the total divergence i.e. sum of absolute
divergence in given cells.

Uses gensimcell::get() for accessing vector variable from
which divergence is calculated and scalar variable where
divergence is stored.

Vector variable must have a defined value in face
neighbors of given cells. In given cells contribution
to divergence is 0 from dimensions with at least one
missing neighbor.

Assumes Grid_T API is compatible with dccrg.

Vec_ and Div_Containers are dummy classes, store the
Nested_Vars_To_ in e.g. std::tuple.

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
	std::tuple<Vector_Field>(),
	std::tuple<Divergence>()
);
*/
template <
	class Grid_T,
	template<class... Nested_Vars_To_Vec> class Vec_Container,
	class... Nested_Vars_To_Vec,
	template<class... Nested_Vars_To_Div> class Div_Container,
	class... Nested_Vars_To_Div
> double get_divergence(
	const std::vector<uint64_t>& cells,
	Grid_T& grid,
	const Vec_Container<Nested_Vars_To_Vec...>&,
	const Div_Container<Nested_Vars_To_Div...>&
) {
	double local_divergence = 0, global_divergence = 0;
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
				std::cerr <<  __FILE__ << "(" << __LINE__<< ")" << std::endl;
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

		// calculate divergence
		auto* const cell_data = grid[cell];
		if (cell_data == NULL) {
			std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
				<< "No data for simulation cell " << cell
				<< std::endl;
			abort();
		}

		const auto& vec = gensimcell::get(*cell_data, Nested_Vars_To_Vec()...);
		auto& div = gensimcell::get(*cell_data, Nested_Vars_To_Div()...);

		div = 0;
		for (auto dim = 0; dim < 3; dim++) {
			// zero contribution to gradient from dims with missing neighbors
			if (nr_neighbors[dim] == 2) {
				div += vec[dim] * (neigh_pos_dist[dim] - neigh_neg_dist[dim]);
			}
		}

		for (const auto& item: face_neighbors_of) {
			const auto neighbor = item.first;
			const auto direction = item.second;
			const size_t dim = std::abs(direction) - 1;

			if (nr_neighbors[dim] != 2) {
				continue;
			}

			const auto* const neighbor_data = grid[neighbor];
			if (neighbor_data == NULL) {
				std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
					<< "No data for neighbor " << neighbor
					<< " of cell " << cell
					<< std::endl;
				abort();
			}
			const auto& neigh_vec
				= gensimcell::get(*neighbor_data, Nested_Vars_To_Vec()...);

			double multiplier = 0;
			if (direction < 0) {
				multiplier = -neigh_pos_dist[dim] / neigh_neg_dist[dim];
			} else {
				multiplier = neigh_neg_dist[dim] / neigh_pos_dist[dim];
			}
			multiplier /= (neigh_pos_dist[dim] + neigh_neg_dist[dim]);

			div += multiplier * neigh_vec[dim];
		}

		local_divergence += std::fabs(div);
	}

	MPI_Comm comm = grid.get_communicator();
	MPI_Allreduce(&local_divergence, &global_divergence, 1, MPI_DOUBLE, MPI_SUM, comm);
	MPI_Comm_free(&comm);

	return global_divergence;
}


/*!
Calculates gradient of scalar variable in given cells.

See get_divergence() for more info.
*/
template <
	class Grid_T,
	template<class... Nested_Vars_To_Sca> class Sca_Container,
	class... Nested_Vars_To_Sca,
	template<class... Nested_Vars_To_Grad> class Grad_Container,
	class... Nested_Vars_To_Grad
> void get_gradient(
	const std::vector<uint64_t>& cells,
	Grid_T& grid,
	const Sca_Container<Nested_Vars_To_Sca...>&,
	const Grad_Container<Nested_Vars_To_Grad...>&
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
				std::cerr <<  __FILE__ << "(" << __LINE__<< ")" << std::endl;
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

		// calculate gradient
		auto* const cell_data = grid[cell];
		if (cell_data == NULL) {
			std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
				<< "No data for simulation cell " << cell
				<< std::endl;
			abort();
		}

		const auto& scalar = gensimcell::get(*cell_data, Nested_Vars_To_Sca()...);
		auto& gradient = gensimcell::get(*cell_data, Nested_Vars_To_Grad()...);

		gradient[0] =
		gradient[1] =
		gradient[2] = 0;
		for (auto dim = 0; dim < 3; dim++) {
			// gradient is zero in dimensions with missing neighbor(s)
			if (nr_neighbors[dim] == 2) {
				gradient[dim] += scalar * (neigh_pos_dist[dim] - neigh_neg_dist[dim]);
			}
		}

		for (const auto& item: face_neighbors_of) {
			const auto neighbor = item.first;
			const auto direction = item.second;
			const size_t dim = std::abs(direction) - 1;

			if (nr_neighbors[dim] != 2) {
				continue;
			}

			const auto* const neighbor_data = grid[neighbor];
			if (neighbor_data == NULL) {
				std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
					<< "No data for neighbor " << neighbor
					<< " of cell " << cell
					<< std::endl;
				abort();
			}
			const auto& neigh_scalar
				= gensimcell::get(*neighbor_data, Nested_Vars_To_Sca()...);

			double multiplier = 0;
			if (direction < 0) {
				multiplier = -neigh_pos_dist[dim] / neigh_neg_dist[dim];
			} else {
				multiplier = neigh_neg_dist[dim] / neigh_pos_dist[dim];
			}
			multiplier /= (neigh_pos_dist[dim] + neigh_neg_dist[dim]);

			gradient[dim] += multiplier * neigh_scalar;
		}
	}
}


/*!
Removes divergence of a vector variable from given cells.

Solves phi from div(grad(phi)) = div(vector) and assigns
vector = vector - grad(phi) after which div(vector) == 0.

Vector variable must have been updated between processes
before calling this function.

Vec_Container represents the path to the vector variable
in cells of simulation_grid whose divergence is removed.

Div_Container represents the path to the temporary
scalar variable in cells of simulation_grid to use
by this function.

Grad_Container represents the path to the temporary
vector variable in cells of simulation_grid to use
by this function.

The transfer of Vec and Div variables must have been
switched on in cells of simulation grid before
calling this function.
*/
template <
	class Cell_T,
	class Geometry_T,
	template<class... Nested_Vars_To_Vec> class Vec_Container,
	class... Nested_Vars_To_Vec,
	template<class... Nested_Vars_To_Div> class Div_Container,
	class... Nested_Vars_To_Div,
	template<class... Nested_Vars_To_Grad> class Grad_Container,
	class... Nested_Vars_To_Grad
> void remove(
	const std::vector<uint64_t>& cells,
	dccrg::Dccrg<Cell_T, Geometry_T>& simulation_grid,
	const Vec_Container<Nested_Vars_To_Vec...>&,
	const Div_Container<Nested_Vars_To_Div...>&,
	const Grad_Container<Nested_Vars_To_Grad...>&
) {
	/*
	Prepare solution grid and source term
	*/

	simulation_grid.update_copies_of_remote_neighbors();

	// rhs for poisson solver = div(vec)
	get_divergence(
		cells,
		simulation_grid,
		std::tuple<Nested_Vars_To_Vec...>(),
		std::tuple<Nested_Vars_To_Div...>()
	);

	dccrg::Dccrg<Poisson_Cell, Geometry_T> poisson_grid(simulation_grid);

	// transfer rhs to poisson grid
	for (const auto& cell: cells) {
		auto* const poisson_data = poisson_grid[cell];
		if (poisson_data == NULL) {
			std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
				<< "No data for poisson cell " << cell
				<< std::endl;
			abort();
		}

		const auto* const simulation_data = simulation_grid[cell];
		if (simulation_data == NULL) {
			std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
				<< "No data for simulation cell " << cell
				<< std::endl;
			abort();
		}


		poisson_data->solution = 0;
		poisson_data->rhs = gensimcell::get(*simulation_data, Nested_Vars_To_Div()...);
	}

	// solve phi in div(grad(phi)) = rhs
	Poisson_Solve solver;
	solver.solve(cells, poisson_grid);

	/*
	Remove divergence with B = B - grad(phi)
	*/

	// store phi in divergence variable
	for (const auto& cell: cells) {
		const auto* const poisson_data = poisson_grid[cell];
		if (poisson_data == NULL) {
			std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
				<< "No data for poisson cell " << cell
				<< std::endl;
			abort();
		}

		auto* const simulation_data = simulation_grid[cell];
		if (simulation_data == NULL) {
			std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
				<< "No data for simulation cell " << cell
				<< std::endl;
			abort();
		}

		gensimcell::get(*simulation_data, Nested_Vars_To_Div()...)
			= poisson_data->solution;
	}

	simulation_grid.update_copies_of_remote_neighbors();

	get_gradient(
		cells,
		simulation_grid,
		std::tuple<Nested_Vars_To_Div...>(),
		std::tuple<Nested_Vars_To_Grad...>()
	);

	for (const auto& cell: cells) {
		auto* const data = simulation_grid[cell];
		if (data == NULL) {
			std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
				<< "No data for simulation cell " << cell
				<< std::endl;
			abort();
		}

		const auto& grad = gensimcell::get(*data, Nested_Vars_To_Grad()...);
		auto& vec = gensimcell::get(*data, Nested_Vars_To_Vec()...);
		for (size_t dim = 0; dim < 3; dim++) {
			vec[dim] -= grad[dim];
		}
	}
}


}} // namespaces

#endif // ifndef PAMHD_DIVERGENCE_REMOVE_HPP
