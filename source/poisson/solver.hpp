/*


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

#ifndef PAMHD_POISSON_SOLVER_HPP
#define PAMHD_POISSON_SOLVER_HPP


#include "iostream"

#include "Eigen/Core"
#include "Eigen/IterativeLinearSolvers"


namespace pamhd {
namespace poisson {


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
	const double cell_length_x,
	const double cell_length_y,
	const Scalar_T* const rhs
) {
	if (
		rhs == NULL
		or grid_size_x == 0
		or grid_size_y == 0
	) {
		std::cerr << "Invalid input: "
			<< rhs << ", " << grid_size_x << ", " << grid_size_y << "."
			<< std::endl;
		return NULL;
	}

	const size_t n = grid_size_x * grid_size_y;

	// solve Ax = b
	Eigen::VectorXd x(n), b(n);
	Eigen::SparseMatrix<Scalar_T> A(n, n);
	std::vector<Eigen::Triplet<Scalar_T>> A_non_zeros;

	// fill A and b
	for (size_t y_i = 0; y_i < grid_size_y; y_i++) {
	for (size_t x_i = 0; x_i < grid_size_x; x_i++) {

		const size_t cell_index = x_i + y_i * grid_size_x;

		b[cell_index] = *(rhs + cell_index);
		A_non_zeros.push_back(
			Eigen::Triplet<Scalar_T>(
				cell_index,
				cell_index,
				-4 / (cell_length_x * cell_length_y)
			)
		);

		for (long long int y_offset: {1ull, 0ull, grid_size_y - 1ull}) {
		for (long long int x_offset: {1ull, 0ull, grid_size_x - 1ull}) {
			if (
				(x_offset == 0 and y_offset == 0)
				or (x_offset != 0 and y_offset != 0)
			) {
				continue;
			}

			const size_t
				x_final = (x_i + x_offset) % grid_size_x,
				y_final = (y_i + y_offset) % grid_size_y,
				neighbor_index = x_final + y_final * grid_size_x;

			const Scalar_T geometry_factor
				= [&](){
					if (x_offset != 0) {
						return 1 / std::pow(cell_length_x, 2);
					} else {
						return 1 / std::pow(cell_length_y, 2);
					}
				}();

			A_non_zeros.push_back(
				Eigen::Triplet<Scalar_T>(
					cell_index,
					neighbor_index,
					geometry_factor
				)
			);
		}}
	}}
	A.setFromTriplets(A_non_zeros.begin(), A_non_zeros.end());

	Eigen::BiCGSTAB<Eigen::SparseMatrix<Scalar_T>> solver;
	//Eigen::ConjugateGradient<Eigen::SparseMatrix<Scalar_T>> solver;
	solver.compute(A);
	if(solver.info() != Eigen::Success) {
		std::cerr << "Decomposition failed." << std::endl;
		return NULL;
	}

	x = solver.solve(b);
	/*std::cout << "error: " << solver.error() << std::endl;
	std::cout << "iterations: " << solver.iterations() << std::endl;*/

	if(solver.info() != Eigen::Success) {
		std::cerr << "Could not solve." << std::endl;
		return NULL;
	}

	Scalar_T* ret_val = (Scalar_T*)malloc(sizeof(Scalar_T) * n);
	if (ret_val == NULL) {
		std::cerr << "Couldn't allocate memory for output." << std::endl;
		return NULL;
	}

	for (size_t j = 0; j < grid_size_y; j++) {
	for (size_t i = 0; i < grid_size_x; i++) {
		const size_t cell_index = i + j * grid_size_x;
		*(ret_val + cell_index) = x[cell_index];
	}}

	return ret_val;
}

}} // namespaces

#endif // ifndef PAMHD_POISSON_SOLVER_HPP
