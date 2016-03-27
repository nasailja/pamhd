/*
Tests serial failsafe Poisson solver of PAMHD in 1d.

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


#include "cstdlib"
#include "iostream"
#include "vector"

#include "poisson/solver.hpp"

//! size == number of cells, length == physical length
template <class Scalar_T> std::vector<Scalar_T> get_analytic_solution(
	const Scalar_T grid_length_x,
	const size_t grid_size_x,
	const size_t grid_size_y,
	const size_t grid_size_z
) {
	std::vector<Scalar_T> solution(grid_size_x * grid_size_y * grid_size_z, 0);

	for (size_t z_i = 0; z_i < grid_size_z; z_i++) {
	for (size_t y_i = 0; y_i < grid_size_y; y_i++) {
	for (size_t x_i = 0; x_i < grid_size_x; x_i++) {
		const Scalar_T x = grid_length_x * (x_i + 0.5) / grid_size_x;

		const size_t cell_index
			= pamhd::poisson::map_3_to_1(x_i, y_i, z_i, grid_size_x, grid_size_y);
		solution[cell_index] = std::sin(x / 2 + 3);
	}}}

	return solution;
}


template <class Scalar_T> std::vector<Scalar_T> get_rhs(
	const Scalar_T grid_length_x,
	const size_t grid_size_x,
	const size_t grid_size_y,
	const size_t grid_size_z
) {
	std::vector<Scalar_T> rhs(grid_size_x * grid_size_y * grid_size_z, 0);

	for (size_t z_i = 0; z_i < grid_size_z; z_i++) {
	for (size_t y_i = 0; y_i < grid_size_y; y_i++) {
	for (size_t x_i = 0; x_i < grid_size_x; x_i++) {
		const Scalar_T x = grid_length_x * (x_i + 0.5) / grid_size_x;

		const size_t cell_index
			= pamhd::poisson::map_3_to_1(x_i, y_i, z_i, grid_size_x, grid_size_y);
		rhs[cell_index] = -std::sin(x / 2 + 3) / 4;
	}}}

	return rhs;
}


template <class Scalar_T> std::vector<Scalar_T> normalize_solution(
	const Scalar_T grid_length_x,
	const size_t grid_size_x,
	const size_t grid_size_y,
	const size_t grid_size_z,
	const std::vector<Scalar_T>& solution
) {
	Scalar_T avg_analytic = 0;
	for (
		const auto i:
		get_analytic_solution(
			grid_length_x,
			grid_size_x,
			grid_size_y,
			grid_size_z
		)
	) {
		avg_analytic += i;
	}
	avg_analytic /= solution.size();

	Scalar_T avg_solution = 0;
	for (const auto i: solution) {
		avg_solution += i;
	}
	avg_solution /= solution.size();

	std::vector<Scalar_T> ret_val;
	ret_val.reserve(solution.size());
	for (const auto i: solution) {
		ret_val.push_back(i - avg_solution + avg_analytic);
	}

	return ret_val;
}


/*!
Returns maximum norm if p == 0
*/
template <class Scalar_T> Scalar_T get_diff_lp_norm(
	const std::vector<Scalar_T>& data1,
	const std::vector<Scalar_T>& data2,
	const Scalar_T& p
) {
	if (data1.size() != data2.size()) {
		abort();
	}

	if (p == 0) {
		Scalar_T maximum = std::numeric_limits<Scalar_T>::lowest();

		for (size_t i = 0; i < data1.size(); i++) {
			maximum = std::max(maximum, std::fabs(data1[i] - data2[i]));
		}

		return maximum;
	}

	Scalar_T norm = 0;
	for (size_t i = 0; i < data1.size(); i++) {
		norm += std::pow(std::fabs(data1[i] - data2[i]), p);
	}

	return std::pow(norm, Scalar_T(1) / p);
}


int main()
{
	using scalar_type = double;

	constexpr double
		grid_length_x = 4 * M_PI,
		grid_length_y = 1,
		grid_length_z = 1;

	constexpr size_t
		grid_size_x = 32,
		grid_size_y = 1,
		grid_size_z = 1;

	const auto rhs
		= get_rhs<scalar_type>(
			grid_length_x,
			grid_size_x,
			grid_size_y,
			grid_size_z
		);

	auto solution
		= pamhd::poisson::solve_failsafe(
			grid_size_x,
			grid_size_y,
			grid_size_z,
			grid_length_x / grid_size_x,
			grid_length_y / grid_size_y,
			grid_length_z / grid_size_z,
			rhs
		);

	solution = normalize_solution(
			grid_length_x,
			grid_size_x,
			grid_size_y,
			grid_size_z,
			solution
	);

	const auto analytic
		= get_analytic_solution(
			grid_length_x,
			grid_size_x,
			grid_size_y,
			grid_size_z
		);

	const double
		diff_l1_norm = get_diff_lp_norm(solution, analytic, scalar_type(1)),
		diff_l2_norm = get_diff_lp_norm(solution, analytic, scalar_type(2)),
		diff_linf_norm = get_diff_lp_norm(solution, analytic, scalar_type(0));

	if (diff_l1_norm > 0.1) {
		std::cerr <<  __FILE__ << "(" << __LINE__ << "): "
			<< "L1 norm too large: " << diff_l1_norm
			<< std::endl;
		abort();
	}
	if (diff_l2_norm > 0.02) {
		std::cerr <<  __FILE__ << "(" << __LINE__ << "): "
			<< "L2 norm too large: " << diff_l2_norm
			<< std::endl;
		abort();
	}
	if (diff_linf_norm > 0.01) {
		std::cerr <<  __FILE__ << "(" << __LINE__ << "): "
			<< "Infinite L norm too large: " << diff_linf_norm
			<< std::endl;
		abort();
	}

	return EXIT_SUCCESS;
}
