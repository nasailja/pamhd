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


#include "cstdlib"
#include "iostream"
#include "vector"

#include "poisson/solver.hpp"

template <class Scalar_T> std::vector<Scalar_T> get_rhs(
	const size_t grid_size_x,
	const size_t grid_size_y
) {
	std::vector<Scalar_T> rhs(grid_size_x * grid_size_y, 0);

	for (size_t y_i = 0; y_i < grid_size_y; y_i++) {
		for (size_t x_i = 0; x_i < grid_size_x; x_i++) {
			const size_t cell_index = x_i + y_i * grid_size_x;

			rhs[cell_index]
				= std::sin(2 * M_PI * (x_i + 0.5) / grid_size_x)
				* std::cos(4 * M_PI * (y_i + 0.5) / grid_size_y);
		}
	}

	return rhs;
}


template <class Scalar_T> std::vector<Scalar_T> get_analytic_solution(
	const size_t grid_size_x,
	const size_t grid_size_y
) {
	std::vector<double> solution(grid_size_x * grid_size_y, 0);

	for (size_t y_i = 0; y_i < grid_size_y; y_i++) {
		for (size_t x_i = 0; x_i < grid_size_x; x_i++) {
			const size_t cell_index = x_i + y_i * grid_size_x;

			solution[cell_index]
				= (std::sin(2 * M_PI * (x_i + 0.5) / grid_size_x)
					* std::cos(4 * M_PI * (y_i + 0.5) / grid_size_y))
				/ (-20 * M_PI * M_PI);
		}
	}

	return solution;
}


template <class Scalar_T> void print_grid(
	const size_t grid_size_x,
	const size_t grid_size_y,
	const std::vector<Scalar_T>& data
) {
	for (size_t y_i = 0; y_i < grid_size_y; y_i++) {
		for (size_t x_i = 0; x_i < grid_size_x; x_i++) {
			const size_t cell_index = x_i + y_i * grid_size_x;

			std::cout << data[cell_index] << " ";
		}
		std::cout << "\n";
	}
	std::cout << std::endl;
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

	constexpr size_t
		grid_size_x = 10,
		grid_size_y = 10,
		n = grid_size_x * grid_size_y;

	const auto rhs = get_rhs<scalar_type>(grid_size_x, grid_size_y);

	scalar_type* solution_tmp
		= pamhd::poisson::solve(
			grid_size_x,
			grid_size_y,
			1.0 / grid_size_x,
			1.0 / grid_size_y,
			rhs.data()
		);
	if (solution_tmp == NULL) {
		std::cerr << "Solution failed." << std::endl;
		return EXIT_FAILURE;
	}
	std::vector<scalar_type> solution(n, 0);
	for (size_t i = 0; i < n; i++) {
		solution[i] = *(solution_tmp + i);
	}
	free(solution_tmp);

	// offset solution to average == 0
	scalar_type avg = 0;
	for (const auto i: solution) {
		avg += i;
	}
	avg /= n;
	for (auto i: solution) {
		i -= avg;
	}

	/*print_grid(grid_size_x, grid_size_y, rhs);
	std::cout << "Solution:\n";
	print_grid(grid_size_x, grid_size_y, solution);*/

	const auto analytic = get_analytic_solution<scalar_type>(grid_size_x, grid_size_y);
	/*std::cout << "Analytic:\n";
	print_grid(
		grid_size_x,
		grid_size_y,
		analytic
	);*/

	std::cout << "Norms:\n"
		<< "  1: " << get_diff_lp_norm(solution, analytic, scalar_type(1)) << "\n"
		<< "  2: " << get_diff_lp_norm(solution, analytic, scalar_type(2)) << "\n"
		<< "Inf: " << get_diff_lp_norm(solution, analytic, scalar_type(0)) << "\n"
		<< std::endl;

	return EXIT_SUCCESS;
}
