program pamhd_poisson_test
use, intrinsic :: iso_c_binding
implicit none

	interface
		function pamhd_poisson_solve_double(&
			grid_size_x,&
			grid_size_y,&
			cell_size_x,&
			cell_size_y,&
			rhs&
		) result(ret_val) bind(c)
			use, intrinsic :: iso_c_binding
			integer(c_size_t) :: grid_size_x, grid_size_y
			real(c_double) :: cell_length_x, cell_length_y
			real(c_double), dimension(*) :: rhs
		end function
	end interface

	integer(c_size_t), parameter :: grid_size_x = 10, grid_size_y = 10
	real(c_double), parameter ::&
		cell_length_x = 1.0 / grid_size_x,&
		cell_length_y = 1.0 / grid_size_y,&
		pi = 3.1415926535897932384626

	real(c_double), dimension(grid_size_x * grid_size_y) :: rhs, solution
	real :: x, y
	integer x_i, y_i

	do y_i = 0, grid_size_y
	do x_i = 0, grid_size_x

		x = (x_i + 0.5) / grid_size_x
		y = (y_i + 0.5) / grid_size_y

		rhs(x_i + y_i * grid_size_x) = sin(2*pi*x) * cos(4*pi*y)

	end do
	end do

	solution&
		= pamhd_poisson_solve_double(&
			grid_size_x,&
			grid_size_y,&
			1.0 / grid_size_x,&
			1.0 / grid_size_y,&
			rhs&
		)

end program
