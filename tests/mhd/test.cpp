/*
MHD test program of PAMHD.

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
Eigen and MPI must be included before generic cell
to get MPI support for Eigen types in generic cell
*/
#include "Eigen/Core"
#include "gensimcell.hpp"

#include "mhd/common.hpp"
#include "mhd/initialize.hpp"
#include "mhd/save.hpp"
#include "mhd/solve.hpp"
#include "mhd/variables.hpp"
#include "boundaries/initial_condition.hpp"
#include "boundaries/time_dependent_boundary.hpp"


using namespace std;


int main(int argc, char* argv[])
{
	using Cell_T = pamhd::mhd::Cell;
	using Grid_T = dccrg::Dccrg<Cell_T, dccrg::Cartesian_Geometry>;

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


	using Init_Cond_T
		= pamhd::boundaries::Initial_Condition<
			uint64_t,
			std::array<double, 3>,
			pamhd::mhd::Number_Density,
			pamhd::mhd::Velocity,
			pamhd::mhd::Pressure,
			pamhd::mhd::Magnetic_Field
		>;
	Init_Cond_T initial_condition;

	using Boundary_T
		= pamhd::boundaries::Time_Dependent_Boundary<
			uint64_t,
			std::array<double, 3>,
			double,
			pamhd::mhd::Number_Density,
			pamhd::mhd::Velocity,
			pamhd::mhd::Pressure,
			pamhd::mhd::Magnetic_Field
		>;
	std::array<Boundary_T, 6> boundaries; // order: -x, +x, -y, +y, -z, +z


	/*
	Program options
	*/

	bool verbose = false, save_fluxes = false;
	double
		save_mhd_n = -1,
		start_time = 0,
		end_time = 1,
		adiabatic_index = 5.0 / 3.0,    
		vacuum_permeability = 4e-7 * M_PI,
		proton_mass = 1.672621777e-27;
	std::string
		mhd_solver("hll_athena"),
		config_file_name("");

	std::array<uint64_t, 3> number_of_cells{1, 1, 1};
	std::array<double, 3>
		simulation_volume{1, 1, 1},
		simulation_volume_start{0, 0, 0};
	std::array<bool, 3> is_periodic{false, false, false};

	boost::program_options::options_description
		options(
			"All supported options"
		),
		// grouped options for printing help
		basic_options(
			"Usage: program_name [options], where options are:"
		),
		initial_condition_options(
			"Options for initial condition: default plasma state "
			"and state on the negative side of the tube"
		),
		boundary_options(
			"Options for plasma state at negative and positive ends of shock tube"
		);

	basic_options.add_options()
		("help", "Print this help message")
		("verbose", "Print run time information")
		("initial-help", "Print help for initial condition options")
		("boundary-help", "Print help for boundary condition options")
		("config-file",
			boost::program_options::value<std::string>(&config_file_name)
				->default_value(config_file_name),
			"Read options from file arg (added to command line, not read if name empty)")
		("cells",
			boost::program_options::value<std::array<uint64_t, 3>>(&number_of_cells)
				->default_value(number_of_cells),
			"Number of grid cells, including any boundary cells")
		("periodic",
			boost::program_options::value<std::array<bool, 3>>(&is_periodic)
				->default_value(is_periodic),
			"Whether grid is periodic (1 == yes, 0 == no")
		("time-start",
			boost::program_options::value<double>(&start_time)
				->default_value(start_time),
			"Start time of simulation (s)")
		("time-length",
			boost::program_options::value<double>(&end_time)
				->default_value(end_time),
			"Length of simulation (s)")
		("volume-start",
			boost::program_options::value<std::array<double, 3>>(&simulation_volume_start)
				->default_value(simulation_volume_start),
			"Start coordinate of simulated volume, including any boundaries (m)")
		("volume-length",
			boost::program_options::value<std::array<double, 3>>(&simulation_volume)
				->default_value(simulation_volume),
			"Length of simulated volume, including any boundaries (m)")
		("solver-mhd",
			boost::program_options::value<std::string>(&mhd_solver)
				->default_value(mhd_solver),
			"MHD solver to use, one of: hll_athena")
		("save-mhd-n",
			boost::program_options::value<double>(&save_mhd_n)
				->default_value(save_mhd_n),
			"Save results every arg seconds, 0 saves "
			"initial and final states, -1 doesn't save")
		("save-mhd-fluxes", "Save fluxes of MHD variables");

	options.add(basic_options);
	initial_condition.add_options(options, "initial.");
	boundaries[0].add_options(options, "x-neg.");
	boundaries[1].add_options(options, "x-pos.");
	boundaries[2].add_options(options, "y-neg.");
	boundaries[3].add_options(options, "y-pos.");
	boundaries[4].add_options(options, "z-neg.");
	boundaries[5].add_options(options, "z-pos.");

	// group options for printing help
	initial_condition.add_options(initial_condition_options, "initial.");
	boundaries[0].add_options(boundary_options, "x-neg.");
	boundaries[1].add_options(boundary_options, "x-pos.");
	boundaries[2].add_options(boundary_options, "y-neg.");
	boundaries[3].add_options(boundary_options, "y-pos.");
	boundaries[4].add_options(boundary_options, "z-neg.");
	boundaries[5].add_options(boundary_options, "z-pos.");

	boost::program_options::variables_map option_variables;
	boost::program_options::store(
		boost::program_options::parse_command_line(argc, argv, options),
		option_variables
	);
	boost::program_options::notify(option_variables);

	if (option_variables.count("help") > 0) {
		if (rank == 0) {
			cout << basic_options << endl;
		}
		MPI_Finalize();
		return EXIT_SUCCESS;
	}

	if (option_variables.count("initial-help") > 0) {
		if (rank == 0) {
			cout << initial_condition_options << endl;
		}
		MPI_Finalize();
		return EXIT_SUCCESS;
	}

	if (option_variables.count("boundary-help") > 0) {
		if (rank == 0) {
			cout << boundary_options << endl;
		}
		MPI_Finalize();
		return EXIT_SUCCESS;
	}

	if (option_variables.count("verbose") > 0) {
		verbose = true;
	}

	if (option_variables.count("save-mhd-fluxes") > 0) {
		save_fluxes = true;
	}

	if (config_file_name != "") {
		boost::program_options::store(
			boost::program_options::parse_config_file<char>(
				config_file_name.c_str(),
				options
			),
			option_variables
		);
		boost::program_options::notify(option_variables);
	}

	if (option_variables.count("verbose") > 0) {
		verbose = true;
	}

	if (option_variables.count("save-mhd-fluxes") > 0) {
		save_fluxes = true;
	}

	// update derived grid parameters based on given options
	std::array<double, 3> cell_volume{
		simulation_volume[0] / number_of_cells[0],
		simulation_volume[1] / number_of_cells[1],
		simulation_volume[2] / number_of_cells[2]
	};


	/*
	Initialize simulation grid
	*/

	Grid_T grid;

	const unsigned int neighborhood_size = 1;
	if (not grid.initialize(
		number_of_cells,
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

	// set grid geometry
	dccrg::Cartesian_Geometry::Parameters geom_params;
	geom_params.start = simulation_volume_start;
	geom_params.level_0_cell_length = cell_volume;

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
		outer_cells = grid.get_local_cells_on_process_boundary();

	if (verbose and rank == 0) {
		cout << "Initializing MHD... " << endl;
	}
	pamhd::mhd::initialize<
		Grid_T,
		Init_Cond_T,
		pamhd::mhd::MHD_State_Conservative,
		pamhd::mhd::MHD_Flux_Conservative,
		pamhd::mhd::Mass_Density,
		pamhd::mhd::Momentum_Density,
		pamhd::mhd::Total_Energy_Density,
		pamhd::mhd::Magnetic_Field
	>(
		grid,
		initial_condition,
		cells,
		adiabatic_index,
		vacuum_permeability,
		proton_mass
	);
	if (verbose and rank == 0) {
		cout << "Done initializing MHD" << endl;
	}

	/*
	Simulate
	*/

	double
		max_dt = 0,
		simulation_time = start_time,
		next_mhd_save = save_mhd_n;
	while (simulation_time < end_time) {

		Cell_T::set_transfer_all(
			true,
			pamhd::mhd::MHD_State_Conservative()
		);
		/*grid.balance_load();

		cells = grid.get_cells();
		inner_cells = grid.get_local_cells_not_on_process_boundary();
		outer_cells = grid.get_local_cells_on_process_boundary();*/

		/*
		Apply boundaries
		*/

		// remove classifications from previous steps
		for (size_t i = 0; i < boundaries.size(); i++) {
			boundaries[i].clear_cells();
		}

		for (const auto cell_id: cells) {
			const auto
				cell_min = grid.geometry.get_min(cell_id),
				cell_max = grid.geometry.get_max(cell_id);

			for (size_t i = 0; i < boundaries.size(); i++) {
				if (not boundaries[i].exists()) {
					continue;
				}

				boundaries[i].add_cell(
					cell_id,
					cell_min,
					cell_max,
					simulation_time
				);
			}
		}

		const pamhd::mhd::Mass_Density Rho{};
		const pamhd::mhd::Number_Density N{};
		const pamhd::mhd::Velocity V{};
		const pamhd::mhd::Pressure P{};
		const pamhd::mhd::Magnetic_Field B{};

		for (size_t i = 0; i < boundaries.size(); i++) {
			if (not boundaries[i].exists()) {
				continue;
			}

			for (const auto cell_id: boundaries[i].get_boundary_cells(simulation_time)) {
				auto* cell_data = grid[cell_id];
				if (cell_data == NULL) {
					std::cerr <<  __FILE__ << "(" << __LINE__ << "): "
						"No data for cell: " << cell_id
						<< std::endl;
					abort();
				}

				pamhd::mhd::MHD_Primitive temp;
				temp[Rho]
					= boundaries[i].get_boundary_data(N, simulation_time)
					* proton_mass;
				temp[V] = boundaries[i].get_boundary_data(V, simulation_time);
				temp[P] = boundaries[i].get_boundary_data(P, simulation_time);
				temp[B] = boundaries[i].get_boundary_data(B, simulation_time);

				auto& state = (*cell_data)[pamhd::mhd::MHD_State_Conservative()];
				state = pamhd::mhd::get_conservative<
					pamhd::mhd::MHD_Conservative,
					pamhd::mhd::MHD_Primitive,
					pamhd::mhd::Momentum_Density,
					pamhd::mhd::Total_Energy_Density,
					pamhd::mhd::Mass_Density,
					pamhd::mhd::Velocity,
					pamhd::mhd::Pressure,
					pamhd::mhd::Magnetic_Field
				>(temp, adiabatic_index, vacuum_permeability);
			}
		}


		/*
		Get maximum allowed time step
		*/

		double
			time_step_factor = 0.5,
			// don't step over the final simulation time
			until_end = end_time - simulation_time,
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

		pamhd::mhd::zero_fluxes<
			Grid_T,
			Cell_T,
			pamhd::mhd::MHD_Flux_Conservative,
			pamhd::mhd::Mass_Density,
			pamhd::mhd::Momentum_Density,
			pamhd::mhd::Total_Energy_Density,
			pamhd::mhd::Magnetic_Field
		>(grid, cells);

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
				mhd_solver,
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
				mhd_solver,
				grid,
				outer_cells,
				time_step,
				adiabatic_index,
				vacuum_permeability
			)
		);

		pamhd::mhd::apply_fluxes<
			Grid_T,
			Cell_T,
			pamhd::mhd::MHD_State_Conservative,
			pamhd::mhd::MHD_Flux_Conservative
		>(grid, inner_cells);

		grid.wait_remote_neighbor_copy_update_sends();

		pamhd::mhd::apply_fluxes<
			Grid_T,
			Cell_T,
			pamhd::mhd::MHD_State_Conservative,
			pamhd::mhd::MHD_Flux_Conservative
		>(grid, outer_cells);

		simulation_time += time_step;

		/*
		Save simulation to disk
		*/

		Cell_T::set_transfer_all(
			false,
			pamhd::mhd::MHD_State_Conservative()
		);

		if (
			(save_mhd_n >= 0 and (simulation_time == 0 or simulation_time >= end_time))
			or (save_mhd_n > 0 and simulation_time >= next_mhd_save)
		) {
			if (next_mhd_save <= simulation_time) {
				next_mhd_save += save_mhd_n;
			}

			if (verbose and rank == 0) {
				cout << "Saving MHD at time " << simulation_time << endl;
			}

			pamhd::mhd::Save::save<
				Grid_T,
				Cell_T,
				pamhd::mhd::MHD_State_Conservative,
				pamhd::mhd::MHD_Flux_Conservative
			>(grid, simulation_time, save_fluxes);
		}
	}

	if (verbose and rank == 0) {
		cout << "Simulation finished at time " << simulation_time << endl;
	}

	MPI_Finalize();

	return EXIT_SUCCESS;
}
