//! Wrappers for instances of pamhd::poisson::solve
#include "poisson/solver.hpp"

double* pamhd_poisson_solve_double(
	const size_t* grid_size_x,
	const size_t* grid_size_y,
	const double* cell_length_x,
	const double* cell_length_y,
	const double** const rhs
) {
	return pamhd::poisson::solve(
		*grid_size_x,
		*grid_size_y,
		*cell_length_x,
		*cell_length_y,
		*rhs
	);
}
