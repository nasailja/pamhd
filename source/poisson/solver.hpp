/*
Serial PAMHD solvers of Poisson's equation.

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

#ifndef PAMHD_POISSON_SOLVER_HPP
#define PAMHD_POISSON_SOLVER_HPP


#include "cmath"
#include "iostream"
#include "stdexcept"

#include "Eigen/Core"
#include "Eigen/IterativeLinearSolvers"


namespace pamhd {
namespace poisson {

/*!
Wraps given index + offset to index within grid_size [0, grid_size[.
*/
size_t wrap(const size_t index, const int offset, const size_t grid_size)
{
	if (offset < 0) {
		return (index + grid_size - (std::abs(offset) % grid_size)) % grid_size;
	} else {
		return (index + offset) % grid_size;
	}
}

size_t map_3_to_1(
	const size_t x_i,
	const size_t y_i,
	const size_t z_i,
	const size_t grid_size_x,
	const size_t grid_size_y
) {
	return
		x_i
		+ y_i * grid_size_x
		+ z_i * grid_size_x * grid_size_y;
}


/*!
Failsafe Poisson solver used as a last resort.

From http://www.rsmas.miami.edu/personal/miskandarani/Courses/MSC321/Projects/prjpoisson.pdf
with corrected geometric factors in equations 15-18.
*/
template<class Scalar_T> std::vector<Scalar_T> solve_failsafe(
	const size_t grid_size_x,
	const size_t grid_size_y,
	const size_t grid_size_z,
	const Scalar_T cell_length_x,
	const Scalar_T cell_length_y,
	const Scalar_T cell_length_z,
	const std::vector<Scalar_T>& rhs
) {
	const size_t nr_of_cells = grid_size_x * grid_size_y * grid_size_z;

	if (rhs.size() == 0 or rhs.size() != nr_of_cells) {
		return {};
	}

	const Scalar_T
		cell_len_x2 = cell_length_x * cell_length_x,
		cell_len_y2 = cell_length_y * cell_length_y,
		cell_len_z2 = cell_length_z * cell_length_z,
		common_factor // incorrect in eqs 15-18
			= 0.5 / (
				cell_len_x2 * cell_len_y2
				+ cell_len_x2 * cell_len_z2
				+ cell_len_y2 * cell_len_z2
			);

	std::vector<Scalar_T>
		solution(nr_of_cells, 0),
		new_solution(nr_of_cells, 0);

	size_t iterations = 0;
	Scalar_T solution_norm = std::numeric_limits<Scalar_T>::max();
	while (iterations++ < 100000 and solution_norm > 1e-10) {
		solution_norm = 0;

		for (size_t z_i = 0; z_i < grid_size_z; z_i++) {
		for (size_t y_i = 0; y_i < grid_size_y; y_i++) {
		for (size_t x_i = 0; x_i < grid_size_x; x_i++) {

			const auto cell_index
				= map_3_to_1(x_i, y_i, z_i, grid_size_x, grid_size_y);

			new_solution[cell_index]
				= -rhs[cell_index]
				* cell_len_x2
				* cell_len_y2
				* cell_len_z2
				* common_factor;

			const size_t
				neigh_x_neg_index
					= map_3_to_1(
						wrap(x_i, -1, grid_size_x), y_i, z_i,
						grid_size_x, grid_size_y
					),
				neigh_x_pos_index
					= map_3_to_1(
						wrap(x_i, +1, grid_size_x), y_i, z_i,
						grid_size_x, grid_size_y
					),
				neigh_y_neg_index
					= map_3_to_1(
						x_i, wrap(y_i, -1, grid_size_y), z_i,
						grid_size_x, grid_size_y
					),
				neigh_y_pos_index
					= map_3_to_1(
						x_i, wrap(y_i, +1, grid_size_y), z_i,
						grid_size_x, grid_size_y
					),
				neigh_z_neg_index
					= map_3_to_1(
						x_i, y_i, wrap(z_i, -1, grid_size_z),
						grid_size_x, grid_size_y
					),
				neigh_z_pos_index
					= map_3_to_1(
						x_i, y_i, wrap(z_i, +1, grid_size_z),
						grid_size_x, grid_size_y
					);

			const Scalar_T
				new_solution_x
					= (solution[neigh_x_neg_index] + solution[neigh_x_pos_index])
					* cell_len_y2 * cell_len_z2,
				new_solution_y
					= (solution[neigh_y_neg_index] + solution[neigh_y_pos_index])
					* cell_len_x2 * cell_len_z2,
				new_solution_z
					= (solution[neigh_z_neg_index] + solution[neigh_z_pos_index])
					* cell_len_x2 * cell_len_y2;

			new_solution[cell_index]
				+= common_factor
				* (new_solution_x + new_solution_y + new_solution_z);

			solution_norm
				= std::max(
					solution_norm,
					std::fabs(new_solution[cell_index] - solution[cell_index])
				);
		}}}

		for (size_t i = 0; i < solution.size(); i++) {
			solution[i] = new_solution[i];
		}
	}

	return solution;
}


/*!
Uses Eigen::BiCGSTAB.
*/
template<class Scalar_T> std::vector<Scalar_T> solve_bicgstab(
	const size_t grid_size_x,
	const size_t grid_size_y,
	const size_t grid_size_z,
	const Scalar_T cell_length_x,
	const Scalar_T cell_length_y,
	const Scalar_T cell_length_z,
	const std::vector<Scalar_T>& rhs
) {
	const size_t nr_of_cells = grid_size_x * grid_size_y * grid_size_z;

	if (rhs.size() == 0 or rhs.size() != nr_of_cells) {
		throw std::invalid_argument("Number of cells and rhs size do not match");
	}

	const Scalar_T
		cell_len_x2 = cell_length_x * cell_length_x,
		cell_len_y2 = cell_length_y * cell_length_y,
		cell_len_z2 = cell_length_z * cell_length_z;

	// solve Ax = b
	Eigen::VectorXd x(nr_of_cells), b(nr_of_cells);
	Eigen::SparseMatrix<Scalar_T> A(nr_of_cells, nr_of_cells);
	std::vector<Eigen::Triplet<Scalar_T>> A_non_zeros;

	// fill A and b
	for (size_t z_i = 0; z_i < grid_size_z; z_i++) {
	for (size_t y_i = 0; y_i < grid_size_y; y_i++) {
	for (size_t x_i = 0; x_i < grid_size_x; x_i++) {

		const size_t cell_index
			= map_3_to_1(x_i, y_i, z_i, grid_size_x, grid_size_y);

		b[cell_index] = rhs[cell_index];

		A_non_zeros.push_back(
			Eigen::Triplet<Scalar_T>(
				cell_index,
				cell_index,
				-2.0 / cell_len_x2 - 2.0 / cell_len_y2 - 2.0 / cell_len_z2
			)
		);

		const size_t
			neigh_x_neg_index
				= map_3_to_1(
					wrap(x_i, -1, grid_size_x), y_i, z_i,
					grid_size_x, grid_size_y
				),
			neigh_x_pos_index
				= map_3_to_1(
					wrap(x_i, +1, grid_size_x), y_i, z_i,
					grid_size_x, grid_size_y
				),
			neigh_y_neg_index
				= map_3_to_1(
					x_i, wrap(y_i, -1, grid_size_y), z_i,
					grid_size_x, grid_size_y
				),
			neigh_y_pos_index
				= map_3_to_1(
					x_i, wrap(y_i, +1, grid_size_y), z_i,
					grid_size_x, grid_size_y
				),
			neigh_z_neg_index
				= map_3_to_1(
					x_i, y_i, wrap(z_i, -1, grid_size_z),
					grid_size_x, grid_size_y
				),
			neigh_z_pos_index
				= map_3_to_1(
					x_i, y_i, wrap(z_i, +1, grid_size_z),
					grid_size_x, grid_size_y
				);

		const Eigen::Triplet<Scalar_T>
			neigh_x_neg_triplet(cell_index, neigh_x_neg_index, 1.0 / cell_len_x2),
			neigh_x_pos_triplet(cell_index, neigh_x_pos_index, 1.0 / cell_len_x2),
			neigh_y_neg_triplet(cell_index, neigh_y_neg_index, 1.0 / cell_len_y2),
			neigh_y_pos_triplet(cell_index, neigh_y_pos_index, 1.0 / cell_len_y2),
			neigh_z_neg_triplet(cell_index, neigh_z_neg_index, 1.0 / cell_len_z2),
			neigh_z_pos_triplet(cell_index, neigh_z_pos_index, 1.0 / cell_len_z2);

		A_non_zeros.push_back(neigh_x_neg_triplet);
		A_non_zeros.push_back(neigh_x_pos_triplet);
		A_non_zeros.push_back(neigh_y_neg_triplet);
		A_non_zeros.push_back(neigh_y_pos_triplet);
		A_non_zeros.push_back(neigh_z_neg_triplet);
		A_non_zeros.push_back(neigh_z_pos_triplet);
	}}}
	A.setFromTriplets(A_non_zeros.begin(), A_non_zeros.end());

	Eigen::BiCGSTAB<Eigen::SparseMatrix<Scalar_T>> solver;
	solver.compute(A);
	if(solver.info() != Eigen::Success) {
		throw std::runtime_error("Couldn't compute eigendecomposition");
	}

	x = solver.solve(b);

	if(solver.info() != Eigen::Success) {
		std::cout << A << std::endl;
		throw std::runtime_error("Couldn't solve equations");
	}

	std::vector<Scalar_T> solution(nr_of_cells);
	for (size_t i = 0; i < nr_of_cells; i++) {
		solution[i] = x[i];
	}

	return solution;
}


/*!
Solves 2-dimensional Poisson's equation
(https://en.wikipedia.org/wiki/Poisson%27s_equation)
with rhs == f and returns a pointer to the result
allocated with malloc which is as long as rhs.

Solution is calculated on a periodic 2d grid with
given length, the length of rhs must be
grid_size_x * grid_size_y.

Caller is responsible for allocating and
deallocating both rhs and the returned result.

Interprets rhs as row-major, i.e. cell at indices
x_i and y_i is stored at address rhs + x_i + y_i * grid_size_x
where x_i is column index and y_i row index of cell.

Returns NULL on failure.
*/
template<class Scalar_T> Scalar_T* solve(
	const size_t grid_size_x,
	const size_t grid_size_y,
	const size_t grid_size_z,
	const double cell_length_x,
	const double cell_length_y,
	const double cell_length_z,
	const Scalar_T* const rhs
) {
	if (
		rhs == NULL
		or grid_size_x == 0
		or grid_size_y == 0
		or grid_size_z == 0
	) {
		std::cerr << "Invalid input: "
			<< rhs << ", "
			<< grid_size_x << ", "
			<< grid_size_y << ", "
			<< grid_size_z << "."
			<< std::endl;
		return NULL;
	}

	const size_t nr_of_cells = grid_size_x * grid_size_y * grid_size_z;

	std::vector<Scalar_T> rhs_impl(nr_of_cells);
	for (size_t i = 0; i < nr_of_cells; i++) {
		rhs_impl[i] = *(rhs + i);
	}

	std::vector<Scalar_T> solution
		= solve_bicgstab(
			grid_size_x,
			grid_size_y,
			grid_size_z,
			cell_length_x,
			cell_length_y,
			cell_length_z,
			rhs_impl
		);

	if (solution.size() == 0) {
		solution = solve_failsafe(
			grid_size_x,
			grid_size_y,
			grid_size_z,
			cell_length_x,
			cell_length_y,
			cell_length_z,
			rhs_impl
		);
	}

	Scalar_T* ret_val = (Scalar_T*)malloc(sizeof(Scalar_T) * nr_of_cells);
	if (ret_val == NULL) {
		std::cerr << "Couldn't allocate memory for output." << std::endl;
		return NULL;
	}

	for (size_t i = 0; i < nr_of_cells; i++) {
		*(ret_val + i) = solution[i];
	}

	return ret_val;
}

}} // namespaces

#endif // ifndef PAMHD_POISSON_SOLVER_HPP
