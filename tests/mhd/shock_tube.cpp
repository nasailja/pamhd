/*
Shock tube test program for MHD mhd_solvers of PAMHD.

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


#include "array"
#include "boost/lexical_cast.hpp"
#include "boost/program_options.hpp"
#include "cmath"
#include "cstdlib"
#include "iostream"
#include "string"
#include "fstream"
#include "iomanip"
#include "iostream"
#include "limits"
#include "mpi.h"
#include "sstream"
#include "string"
#include "type_traits"

#include "dccrg.hpp"
#include "dccrg_cartesian_geometry.hpp"

/*
Eigen amd MPI must be included before generic cell
to get MPI support for Eigen types in generic cell
*/
#include "Eigen/Core"
#include "gensimcell.hpp"

#include "mhd/common.hpp"
#include "mhd/hll_athena.hpp"
#include "mhd/initialize.hpp"
#include "mhd/save.hpp"
#include "mhd/solve.hpp"
#include "mhd/variables.hpp"


using namespace std;


int main(int argc, char* argv[])
{
	using Cell_T = gensimcell::Cell<
		pamhd::mhd::MHD_State_Conservative,
		pamhd::mhd::MHD_Flux_Conservative
	>;
	using Grid_T = dccrg::Dccrg<Cell_T, dccrg::Cartesian_Geometry>;

	// always transfer all variables of MHD state
	pamhd::mhd::MHD_State_Conservative::data_type::set_transfer_all(
		true,
		pamhd::mhd::Mass_Density(),
		pamhd::mhd::Momentum_Density(),
		pamhd::mhd::Total_Energy_Density(),
		pamhd::mhd::Magnetic_Field()
	);

	constexpr double
		adiabatic_index = 5.0 / 3.0,    
		vacuum_permeability = 4e-7 * M_PI,
		proton_mass = 1.672621777e-27;

	/*
	Initialize MPI
	*/

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


	/*
	Get program options
	*/

	bool verbose = false;
	double save_mhd_n = -1;
	std::string mhd_solver;
	size_t cells_n = 0;

	boost::program_options::options_description options(
		"Usage: program_name [options], where options are:"
	);
	options.add_options()
		("help", "Print this help message")
		("cells",
			boost::program_options::value<size_t>(&cells_n)
				->default_value(1000),
			"Number of grid cells along the tube excluding boundary cells")
		("mhd-solver",
			boost::program_options::value<std::string>(&mhd_solver)
				->default_value("hll_athena"),
			"MHD solver to use, one of: hll_athena")
		("save-mhd-n",
			boost::program_options::value<double>(&save_mhd_n)->default_value(-1),
			"Save results every arg seconds, 0 saves "
			"initial and final states, -1 doesn't save")
		("verbose", "Print run time information");

	boost::program_options::variables_map option_variables;
	boost::program_options::store(
		boost::program_options::parse_command_line(argc, argv, options),
		option_variables
	);
	boost::program_options::notify(option_variables);

	if (option_variables.count("help") > 0) {
		if (rank == 0) {
			cout << options << endl;
		}
		MPI_Finalize();
		return EXIT_SUCCESS;
	}

	if (option_variables.count("verbose") > 0) {
		verbose = true;
	}

	/*
	Transform parameters to internal units
	*/

	// boundary cells
	cells_n += 2;


	/*
	Initialize simulation grid
	*/

	Grid_T grid;

	const std::array<uint64_t, 3> grid_length = {cells_n, 1, 1};
	const unsigned int neighborhood_size = 1;
	if (not grid.initialize(
		grid_length,
		comm,
		"RANDOM",
		neighborhood_size,
		0,
		false, false, false
	)) {
		std::cerr << __FILE__ << ":" << __LINE__
			<< ": Couldn't initialize grid."
			<< std::endl;
		abort();
	}

	// set the grid's geometry
	const double cell_length = 1e5 / grid_length[0];
	dccrg::Cartesian_Geometry::Parameters geom_params;
	geom_params.start[0] =
	geom_params.start[1] =
	geom_params.start[2] = -cell_length;
	geom_params.level_0_cell_length[0] =
	geom_params.level_0_cell_length[1] =
	geom_params.level_0_cell_length[2] = cell_length;
	if (not grid.set_geometry(geom_params)) {
		std::cerr << __FILE__ << ":" << __LINE__
			<< ": Couldn't set grid geometry."
			<< std::endl;
		abort();
	}

	grid.balance_load();


	/*
	Initialize MHD
	*/

	std::vector<uint64_t>
		cells = grid.get_cells(),
		inner_cells = grid.get_local_cells_not_on_process_boundary(),
		outer_cells = grid.get_local_cells_on_process_boundary(),
		boundary_cells;

	if (grid.is_local(1)) {
		boundary_cells.push_back(1);
	}
	if (grid.is_local(cells_n)) {
		boundary_cells.push_back(cells_n);
	}

	if (verbose and rank == 0) {
		cout << "Initializing MHD" << endl;
	}
	pamhd::mhd::initialize<
		Grid_T,
		pamhd::mhd::MHD_State_Conservative,
		pamhd::mhd::MHD_Flux_Conservative,
		pamhd::mhd::Mass_Density,
		pamhd::mhd::Momentum_Density,
		pamhd::mhd::Total_Energy_Density,
		pamhd::mhd::Magnetic_Field
	>(
		grid,
		cells,
		adiabatic_index,
		vacuum_permeability,
		proton_mass
	);

	/*
	Simulate
	*/

	const double simulation_duration = 1;
	double
		max_dt = 0,
		simulation_time = 0,
		next_mhd_save = save_mhd_n;
	while (simulation_time < simulation_duration) {

		/*
		Save simulation to disk
		*/

		Cell_T::set_transfer_all(
			false,
			pamhd::mhd::MHD_State_Conservative()
		);

		if (
			(save_mhd_n >= 0 and simulation_time == 0)
			or (save_mhd_n > 0 and simulation_time >= next_mhd_save)
		) {
			next_mhd_save += save_mhd_n;

			if (verbose and rank == 0) {
				cout << "Saving MHD at time " << simulation_time << endl;
			}

			pamhd::mhd::save<
				Grid_T,
				Cell_T,
				pamhd::mhd::MHD_State_Conservative
			>(grid, simulation_time);
		}

		if (simulation_time >= simulation_duration) {
			break;
		}


		/*
		Get maximum allowed time step
		*/

		double
			time_step_factor = 0.5,
			// don't step over the final simulation time
			until_end = simulation_duration - simulation_time,
			local_time_step = std::min(time_step_factor * max_dt, until_end),
			time_step = -1;

		if (
			MPI_Allreduce(
				&local_time_step,
				&time_step,
				1,
				MPI_DOUBLE,
				MPI_MIN,
				comm
			) != MPI_SUCCESS
		) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< ": Couldn't set reduce time step."
				<< std::endl;
			abort();
		}


		/*
		Solve
		*/

		max_dt = std::numeric_limits<double>::max();

		if (verbose and rank == 0) {
			cout << "Solving MHD at time " << simulation_time
				<< " s with time step " << time_step << " s" << endl;
		}

		Cell_T::set_transfer_all(
			true,
			pamhd::mhd::MHD_State_Conservative()
		);

		grid.start_remote_neighbor_copy_updates();

		max_dt = std::min(
			max_dt,
			pamhd::mhd::solve<
				Grid_T,
				pamhd::mhd::MHD_State_Conservative,
				pamhd::mhd::MHD_Flux_Conservative,
				pamhd::mhd::Mass_Density,
				pamhd::mhd::Momentum_Density,
				pamhd::mhd::Total_Energy_Density,
				pamhd::mhd::Magnetic_Field
			>(
				"hll_athena",
				grid,
				inner_cells,
				time_step,
				adiabatic_index,
				vacuum_permeability
			)
		);

		grid.wait_remote_neighbor_copy_update_receives();

		max_dt = std::min(
			max_dt,
			pamhd::mhd::solve<
				Grid_T,
				pamhd::mhd::MHD_State_Conservative,
				pamhd::mhd::MHD_Flux_Conservative,
				pamhd::mhd::Mass_Density,
				pamhd::mhd::Momentum_Density,
				pamhd::mhd::Total_Energy_Density,
				pamhd::mhd::Magnetic_Field
			>(
				"hll_athena",
				grid,
				outer_cells,
				time_step,
				adiabatic_index,
				vacuum_permeability
			)
		);

		pamhd::mhd::apply_solution<
			Grid_T,
			Cell_T,
			pamhd::mhd::MHD_State_Conservative,
			pamhd::mhd::MHD_Flux_Conservative,
			pamhd::mhd::Mass_Density,
			pamhd::mhd::Momentum_Density,
			pamhd::mhd::Total_Energy_Density,
			pamhd::mhd::Magnetic_Field
		>(grid, inner_cells);

		grid.wait_remote_neighbor_copy_update_sends();

		pamhd::mhd::apply_solution<
			Grid_T,
			Cell_T,
			pamhd::mhd::MHD_State_Conservative,
			pamhd::mhd::MHD_Flux_Conservative,
			pamhd::mhd::Mass_Density,
			pamhd::mhd::Momentum_Density,
			pamhd::mhd::Total_Energy_Density,
			pamhd::mhd::Magnetic_Field
		>(grid, outer_cells);

		// set boundary conditions
		pamhd::mhd::initialize<
			Grid_T,
			pamhd::mhd::MHD_State_Conservative,
			pamhd::mhd::MHD_Flux_Conservative,
			pamhd::mhd::Mass_Density,
			pamhd::mhd::Momentum_Density,
			pamhd::mhd::Total_Energy_Density,
			pamhd::mhd::Magnetic_Field
		>(
			grid,
			boundary_cells,
			adiabatic_index,
			vacuum_permeability,
			proton_mass
		);

		simulation_time += time_step;
	}

	Cell_T::set_transfer_all(
		false,
		pamhd::mhd::MHD_State_Conservative()
	);

	if (save_mhd_n >= 0) {
		if (verbose and rank == 0) {
			cout << "Saving MHD at time " << simulation_time << endl;
		}
		pamhd::mhd::save<
			Grid_T,
			Cell_T,
			pamhd::mhd::MHD_State_Conservative
		>(grid, simulation_time);
	}

	if (verbose and rank == 0) {
		cout << "Simulation finished at time " << simulation_time << endl;
	}

	MPI_Finalize();

	return EXIT_SUCCESS;
}
