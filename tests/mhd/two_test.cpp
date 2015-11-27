/*
Two-fluid MHD test program of PAMHD.

Copyright 2014, 2015 Ilja Honkonen
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
#include "cmath"
#include "cstdlib"
#include "exception"
#include "iostream"
#include "string"
#include "fstream"
#include "iomanip"
#include "iostream"
#include "limits"
#include "sstream"
#include "string"
#include "type_traits"

#include "boost/filesystem.hpp"
#include "boost/lexical_cast.hpp"
#include "boost/program_options.hpp"
#include "dccrg.hpp"
#include "dccrg_cartesian_geometry.hpp"
#include "Eigen/Core" // must be included before gensimcell.hpp
#include "mpi.h" // must be included before gensimcell.hpp
#include "gensimcell.hpp"

#include "boundaries/copy_boundary.hpp"
#include "boundaries/initial_condition.hpp"
#include "boundaries/value_boundaries.hpp"
#include "divergence/remove.hpp"
#include "grid_options.hpp"
#include "mhd/boundaries.hpp"
#include "mhd/common.hpp"
#include "mhd/save.hpp"
#include "mhd/solve.hpp"
#include "mhd/hll_athena.hpp"
#include "mhd/hlld_athena.hpp"
#include "mhd/roe_athena.hpp"
#include "mhd/variables.hpp"
#include "pamhd/initialize.hpp"


using namespace std;

/*
Controls transfer of variables in poisson solver
which doesn't use generic cell
*/
int Poisson_Cell::transfer_switch = Poisson_Cell::INIT;

using Cell = gensimcell::Cell<
	gensimcell::Optional_Transfer,
	pamhd::mhd::HD1_State,
	pamhd::mhd::HD2_State,
	pamhd::mhd::MHD_State_Conservative,
	pamhd::mhd::Electric_Current_Density,
	pamhd::mhd::Cell_Type,
	pamhd::mhd::Copy_Source,
	pamhd::mhd::Value_Boundary_Id,
	pamhd::mhd::Cell_Type_Fluid2,
	pamhd::mhd::Copy_Source_Fluid2,
	pamhd::mhd::Value_Boundary_Id_Fluid2,
	pamhd::mhd::Cell_Type_Field,
	pamhd::mhd::Value_Boundary_Id_Field,
	pamhd::mhd::Copy_Source_Field,
	pamhd::mhd::MPI_Rank,
	pamhd::mhd::Resistivity,
	pamhd::mhd::Magnetic_Field_Resistive,
	pamhd::mhd::Magnetic_Field_Temp,
	pamhd::mhd::Magnetic_Field_Divergence,
	pamhd::mhd::Scalar_Potential_Gradient,
	pamhd::mhd::HD1_Flux,
	pamhd::mhd::HD2_Flux,
	pamhd::mhd::MHD_Flux_Conservative
>;
using Grid = dccrg::Dccrg<Cell, dccrg::Cartesian_Geometry>;

// reference to total mass density of all fluids in given cell
const auto Mas
	= [](Cell& cell_data)->typename pamhd::mhd::Mass_Density::data_type&{
		return cell_data[pamhd::mhd::MHD_State_Conservative()][pamhd::mhd::Mass_Density()];
	};
const auto Mom
	= [](Cell& cell_data)->typename pamhd::mhd::Momentum_Density::data_type&{
		return cell_data[pamhd::mhd::MHD_State_Conservative()][pamhd::mhd::Momentum_Density()];
	};
const auto Nrj
	= [](Cell& cell_data)->typename pamhd::mhd::Total_Energy_Density::data_type&{
		return cell_data[pamhd::mhd::MHD_State_Conservative()][pamhd::mhd::Total_Energy_Density()];
	};
const auto Mag
	= [](Cell& cell_data)->typename pamhd::mhd::Magnetic_Field::data_type&{
		return cell_data[pamhd::mhd::MHD_State_Conservative()][pamhd::mhd::Magnetic_Field()];
	};

// field before divergence removal in case removal fails
const auto Mag_tmp
	= [](Cell& cell_data)->typename pamhd::mhd::Magnetic_Field_Temp::data_type&{
		return cell_data[pamhd::mhd::Magnetic_Field_Temp()];
	};
// divergence of magnetic field
const auto Mag_div
	= [](Cell& cell_data)->typename pamhd::mhd::Magnetic_Field_Divergence::data_type&{
		return cell_data[pamhd::mhd::Magnetic_Field_Divergence()];
	};
// adjustment to magnetic field due to resistivity
const auto Mag_res
	= [](Cell& cell_data)->typename pamhd::mhd::Magnetic_Field_Resistive::data_type&{
		return cell_data[pamhd::mhd::Magnetic_Field_Resistive()];
	};
// curl of magnetic field
const auto Cur
	= [](Cell& cell_data)->typename pamhd::mhd::Electric_Current_Density::data_type&{
		return cell_data[pamhd::mhd::Electric_Current_Density()];
	};

// flux / total change of magnetic field over one time step
const auto Mag_f
	= [](Cell& cell_data)->typename pamhd::mhd::Magnetic_Field::data_type&{
		return cell_data[pamhd::mhd::MHD_Flux_Conservative()][pamhd::mhd::Magnetic_Field()];
	};

// reference to mass density of fluid 1 in given cell
const auto Mas1
	= [](Cell& cell_data)->typename pamhd::mhd::Mass_Density::data_type&{
		return cell_data[pamhd::mhd::HD1_State()][pamhd::mhd::Mass_Density()];
	};
const auto Mom1
	= [](Cell& cell_data)->typename pamhd::mhd::Momentum_Density::data_type&{
		return cell_data[pamhd::mhd::HD1_State()][pamhd::mhd::Momentum_Density()];
	};
const auto Nrj1
	= [](Cell& cell_data)->typename pamhd::mhd::Total_Energy_Density::data_type&{
		return cell_data[pamhd::mhd::HD1_State()][pamhd::mhd::Total_Energy_Density()];
	};
// reference to mass density of fluid 2 in given cell
const auto Mas2
	= [](Cell& cell_data)->typename pamhd::mhd::Mass_Density::data_type&{
		return cell_data[pamhd::mhd::HD2_State()][pamhd::mhd::Mass_Density()];
	};
const auto Mom2
	= [](Cell& cell_data)->typename pamhd::mhd::Momentum_Density::data_type&{
		return cell_data[pamhd::mhd::HD2_State()][pamhd::mhd::Momentum_Density()];
	};
const auto Nrj2
	= [](Cell& cell_data)->typename pamhd::mhd::Total_Energy_Density::data_type&{
		return cell_data[pamhd::mhd::HD2_State()][pamhd::mhd::Total_Energy_Density()];
	};

// flux of mass density of fluid 1 over one time step
const auto Mas1_f
	= [](Cell& cell_data)->typename pamhd::mhd::Mass_Density::data_type&{
		return cell_data[pamhd::mhd::HD1_Flux()][pamhd::mhd::Mass_Density()];
	};
const auto Mom1_f
	= [](Cell& cell_data)->typename pamhd::mhd::Momentum_Density::data_type&{
		return cell_data[pamhd::mhd::HD1_Flux()][pamhd::mhd::Momentum_Density()];
	};
const auto Nrj1_f
	= [](Cell& cell_data)->typename pamhd::mhd::Total_Energy_Density::data_type&{
		return cell_data[pamhd::mhd::HD1_Flux()][pamhd::mhd::Total_Energy_Density()];
	};
const auto Mas2_f
	= [](Cell& cell_data)->typename pamhd::mhd::Mass_Density::data_type&{
		return cell_data[pamhd::mhd::HD2_Flux()][pamhd::mhd::Mass_Density()];
	};
const auto Mom2_f
	= [](Cell& cell_data)->typename pamhd::mhd::Momentum_Density::data_type&{
		return cell_data[pamhd::mhd::HD2_Flux()][pamhd::mhd::Momentum_Density()];
	};
const auto Nrj2_f
	= [](Cell& cell_data)->typename pamhd::mhd::Total_Energy_Density::data_type&{
		return cell_data[pamhd::mhd::HD2_Flux()][pamhd::mhd::Total_Energy_Density()];
	};
const auto Cell_t
	= [](Cell& cell_data)->typename pamhd::mhd::Cell_Type::data_type&{
		return cell_data[pamhd::mhd::Cell_Type()];
	};


int main(int argc, char* argv[])
{
	using std::min;
	using std::pow;

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


	pamhd::boundaries::Initial_Condition<
		uint64_t,
		double,
		std::array<double, 3>,
		pamhd::mhd::Magnetic_Field
	> initial_field;

	// one for each fluid
	pamhd::boundaries::Initial_Condition<
		uint64_t,
		double,
		std::array<double, 3>,
		pamhd::mhd::Number_Density,
		pamhd::mhd::Velocity,
		pamhd::mhd::Pressure
	> init_cond_fluid1, init_cond_fluid2;

	pamhd::boundaries::Value_Boundaries<
		uint64_t,
		double,
		double,
		std::array<double, 3>,
		pamhd::mhd::Magnetic_Field
	> value_bdy_field;

	pamhd::boundaries::Value_Boundaries<
		uint64_t,
		double,
		double,
		std::array<double, 3>,
		pamhd::mhd::Number_Density,
		pamhd::mhd::Velocity,
		pamhd::mhd::Pressure
	> value_bdy_fluid1, value_bdy_fluid2;

	pamhd::boundaries::Copy_Boundary<
		uint64_t,
		double,
		std::array<double, 3>
	> copy_bdy_fluid1, copy_bdy_fluid2, copy_bdy_fields;

	pamhd::grid::Options grid_options;
	grid_options.data.set_expression(pamhd::grid::Number_Of_Cells(), "{1, 1, 1}");
	grid_options.data.set_expression(pamhd::grid::Volume(), "{1, 1, 1}");
	grid_options.data.set_expression(pamhd::grid::Start(), "{0, 0, 0}");
	grid_options.data.set_expression(pamhd::grid::Periodic(), "{false, false, false}");

	/*
	Program options
	*/

	bool verbose = false, save_fluxes = false;
	size_t
		poisson_iterations_max = 1000,
		poisson_iterations_min = 0;
	double
		save_mhd_n = -1,
		start_time = 0,
		end_time = 1,
		time_step_factor = 0.5,
		resistivity = 0,
		remove_div_B_n = 0.1,
		poisson_norm_stop = 1e-10,
		poisson_norm_increase_max = 10,
		adiabatic_index = 5.0 / 3.0,    
		vacuum_permeability = 4e-7 * M_PI,
		proton_mass = 1.672621777e-27;
	std::string
		mhd_solver_str("roe_athena"),
		config_file_name(""),
		boundary_file_name(""),
		lb_name("RCB"),
		output_directory("./");

	boost::program_options::options_description
		options(
			"Usage: program_name [options], where options are"
		),
		// grouped options for printing help
		initial_condition_help(
			"Options for initial conditions which set cell data to specified "
			"values before the simulation is started of the simulation."
		),
		value_boundary_help(
			"Options for value boundary conditions which set cell data "
			"to specified values after each simulation time step. "
			"Boundaries of different fluids have to be in identical positions."
		),
		copy_boundary_help(
			"Options for copy boundaries which set cell data to average "
			"value of neighboring non-boundary cells. "
			"Boundaries of different fluids have to be in identical positions."
		);

	// handle general options
	options.add_options()
		("help", "Print this help message")
		("verbose", "Print run time information")
		("initial-help", "Print help for initial condition options")
		("value-boundary-help", "Print help for value boundary condition options")
		("copy-boundary-help", "Print help for copy boundary condition options")
		("config-file",
			boost::program_options::value<std::string>(&config_file_name)
				->default_value(config_file_name),
			"Read options also from file arg (command line has priority, "
			"not read if empty string)")
		("output-directory",
			boost::program_options::value<std::string>(&output_directory)
				->default_value(output_directory),
			"Output simulation results into directory arg (relative to "
			"current working directory")
		("time-start",
			boost::program_options::value<double>(&start_time)
				->default_value(start_time),
			"Start time of simulation (s)")
		("time-length",
			boost::program_options::value<double>(&end_time)
				->default_value(end_time),
			"Length of simulation (s)")
		("time-step-factor",
			boost::program_options::value<double>(&time_step_factor)
				->default_value(time_step_factor),
			"Multiply maximum allowed time step (CFL condition) with factor arg")
		("vacuum-permeability",
			boost::program_options::value<double>(&vacuum_permeability)
				->default_value(vacuum_permeability),
			"https://en.wikipedia.org/wiki/Vacuum_permeability in V*s/(A*m)")
		("resistivity",
			boost::program_options::value<double>(&resistivity)
				->default_value(resistivity),
			"Use arg as plasma resistivity (V*m/A, E += arg * J")
		("adiabatic-index",
			boost::program_options::value<double>(&adiabatic_index)
				->default_value(adiabatic_index),
			"https://en.wikipedia.org/wiki/Heat_capacity_ratio")
		("proton-mass",
			boost::program_options::value<double>(&proton_mass)
				->default_value(proton_mass),
			"Mass of a proton in kg")
		("load-balancer",
			boost::program_options::value<std::string>(&lb_name)
				->default_value(lb_name),
			"Load balancing algorithm to use (for example RANDOM, "
			"RCB, HSFC, HYPERGRAPH; for a list of available algorithms see "
			"http://www.cs.sandia.gov/zoltan/ug_html/ug_alg.html)")
		("solver-mhd",
			boost::program_options::value<std::string>(&mhd_solver_str)
				->default_value(mhd_solver_str),
			"MHD solver to use, one of: hll_athena")
		("remove-div-B-n",
			boost::program_options::value<double>(&remove_div_B_n)
				->default_value(remove_div_B_n),
			"Remove divergence of magnetic field every arg seconds (<= 0 doesn't remove)")
		("save-mhd-n",
			boost::program_options::value<double>(&save_mhd_n)
				->default_value(save_mhd_n),
			"Save results every arg seconds, 0 saves "
			"initial and final states, -1 doesn't save")
		("save-mhd-fluxes", "Save fluxes of MHD variables");

	grid_options.add_options("grid.", options);
	initial_field.add_initialization_options("initial-field.", options);
	init_cond_fluid1.add_initialization_options("initial-fluid1.", options);
	init_cond_fluid2.add_initialization_options("initial-fluid2.", options);
	value_bdy_field.add_initialization_options("value-boundary-field.", options);
	value_bdy_fluid1.add_initialization_options("value-boundary-fluid1.", options);
	value_bdy_fluid2.add_initialization_options("value-boundary-fluid2.", options);
	copy_bdy_fluid1.add_initialization_options("copy-boundary-fluid1.", options);
	copy_bdy_fluid2.add_initialization_options("copy-boundary-fluid2.", options);
	copy_bdy_fields.add_initialization_options("copy-boundary-fields.", options);

	boost::program_options::variables_map option_variables;
	try {
		boost::program_options::store(
			boost::program_options::command_line_parser(argc, argv)
				.options(options)
				.allow_unregistered()
				.run(),
			option_variables
		);
	} catch (std::exception& e) {
		if (rank == 0) {
			std::cerr <<  __FILE__ << "(" << __LINE__ << "): "
				<< "Couldn't parse command line options: " << e.what()
				<< std::endl;
		}
		abort();
	}
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

	if (option_variables.count("save-mhd-fluxes") > 0) {
		save_fluxes = true;
	}

	if (config_file_name != "") {
		try {
			boost::program_options::store(
				boost::program_options::parse_config_file<char>(
					config_file_name.c_str(),
					options,
					true
				),
				option_variables
			);
		} catch (std::exception& e) {
			if (rank == 0) {
				std::cerr <<  __FILE__ << "(" << __LINE__ << "): "
					<< "Couldn't parse general options from file "
					<< config_file_name << ": " << e.what()
					<< std::endl;
			}
			abort();
		}
		boost::program_options::notify(option_variables);
	}

	initial_field.add_options("initial-field.", options);
	init_cond_fluid1.add_options("initial-fluid1.", options);
	init_cond_fluid1.add_options("initial-fluid1.", initial_condition_help);
	init_cond_fluid2.add_options("initial-fluid2.", options);
	init_cond_fluid2.add_options("initial-fluid2.", initial_condition_help);
	value_bdy_field.add_options("value-boundary-field.", options);
	value_bdy_field.add_options("value-boundary-field.", value_boundary_help);
	value_bdy_fluid1.add_options("value-boundary-fluid1.", options);
	value_bdy_fluid1.add_options("value-boundary-fluid1.", value_boundary_help);
	value_bdy_fluid2.add_options("value-boundary-fluid2.", options);
	value_bdy_fluid2.add_options("value-boundary-fluid2.", value_boundary_help);
	copy_bdy_fluid1.add_options("copy-boundary-fluid1.", options);
	copy_bdy_fluid1.add_options("copy-boundary-fluid1.", copy_boundary_help);
	copy_bdy_fluid2.add_options("copy-boundary-fluid2.", options);
	copy_bdy_fluid2.add_options("copy-boundary-fluid2.", copy_boundary_help);
	copy_bdy_fields.add_options("copy-boundary-fields.", options);
	copy_bdy_fields.add_options("copy-boundary-fields.", copy_boundary_help);

	// parse again to get boundaries' options
	try {
		boost::program_options::store(
			boost::program_options::command_line_parser(argc, argv)
				.options(options)
				.allow_unregistered()
				.run(),
			option_variables
		);
	} catch (std::exception& e) {
		if (rank == 0) {
			std::cerr <<  __FILE__ << "(" << __LINE__ << "): "
				<< "Couldn't parse command line options: " << e.what()
				<< std::endl;
		}
		abort();
	}

	if (config_file_name != "") {
		try {
			boost::program_options::store(
				boost::program_options::parse_config_file<char>(
					config_file_name.c_str(),
					options,
					true
				),
				option_variables
			);
		} catch (std::exception& e) {
			if (rank == 0) {
				std::cerr <<  __FILE__ << "(" << __LINE__ << "): "
					<< "Couldn't parse boundary options from file "
					<< config_file_name << ": " << e.what()
					<< std::endl;
			}
			abort();
		}
	}
	boost::program_options::notify(option_variables);

	if (option_variables.count("initial-help") > 0) {
		if (rank == 0) {
			cout << initial_condition_help << endl;
		}
		MPI_Finalize();
		return EXIT_SUCCESS;
	}

	if (option_variables.count("value-boundary-help") > 0) {
		if (rank == 0) {
			cout << value_boundary_help << endl;
		}
		MPI_Finalize();
		return EXIT_SUCCESS;
	}

	if (option_variables.count("copy-boundary-help") > 0) {
		if (rank == 0) {
			cout << copy_boundary_help << endl;
		}
		MPI_Finalize();
		return EXIT_SUCCESS;
	}


	if (not isnormal(adiabatic_index) or adiabatic_index < 0) {
		if (rank == 0) {
			std::cerr <<  __FILE__ << "(" << __LINE__ << "): "
				<< "Invalid adiabatic index: " << adiabatic_index
				<< std::endl;
		}
		abort();
	}
	if (not isnormal(vacuum_permeability) or vacuum_permeability < 0) {
		if (rank == 0) {
			std::cerr <<  __FILE__ << "(" << __LINE__ << "): "
				<< "Invalid vacuum permeability: " << vacuum_permeability
				<< std::endl;
		}
		abort();
	}


	if (option_variables.count("verbose") > 0) {
		verbose = true;
	}

	if (option_variables.count("save-mhd-fluxes") > 0) {
		save_fluxes = true;
	}

	const auto mhd_solver
		= [&mhd_solver_str](){
			if (mhd_solver_str == "hll_athena") {

				return pamhd::mhd::athena::get_flux_hll<
					pamhd::mhd::MHD_Conservative,
					pamhd::mhd::Mass_Density,
					pamhd::mhd::Momentum_Density,
					pamhd::mhd::Total_Energy_Density,
					pamhd::mhd::Magnetic_Field
				>;

			} else if (mhd_solver_str == "hlld_athena") {

				return pamhd::mhd::athena::get_flux_hlld<
					pamhd::mhd::MHD_Conservative,
					pamhd::mhd::Mass_Density,
					pamhd::mhd::Momentum_Density,
					pamhd::mhd::Total_Energy_Density,
					pamhd::mhd::Magnetic_Field
				>;

			} else if (mhd_solver_str == "roe_athena") {

				return pamhd::mhd::athena::get_flux_roe<
					pamhd::mhd::MHD_Conservative,
					pamhd::mhd::Mass_Density,
					pamhd::mhd::Momentum_Density,
					pamhd::mhd::Total_Energy_Density,
					pamhd::mhd::Magnetic_Field
				>;

			} else {

				std::cerr <<  __FILE__ << "(" << __LINE__ << ") Invalid solver: "
					<< mhd_solver_str << ", use --help to list available solvers"
					<< std::endl;
				abort();
			}
		}();


	if (rank == 0 and output_directory != "") {
		try {
			boost::filesystem::create_directories(output_directory);
		} catch (const boost::filesystem::filesystem_error& e) {
			std::cerr <<  __FILE__ << "(" << __LINE__ << ") "
				"Couldn't create output directory " << output_directory << ": "
				<< e.what()
				<< std::endl;
			abort();
		}
	}


	/*
	Set/update derived grid parameters based on given options
	*/

	// muparserx doesn't support uint64_t so convert from int
	const auto nr_cells_tmp
		= grid_options.data.get_data(pamhd::grid::Number_Of_Cells());

	const std::array<uint64_t, 3> number_of_cells{{
		uint64_t(nr_cells_tmp[0]),
		uint64_t(nr_cells_tmp[1]),
		uint64_t(nr_cells_tmp[2])
	}};

	const std::array<double, 3>
		simulation_volume
			= grid_options.data.get_data(pamhd::grid::Volume()),
		cell_volume{{
			simulation_volume[0] / number_of_cells[0],
			simulation_volume[1] / number_of_cells[1],
			simulation_volume[2] / number_of_cells[2]
		}};


	/*
	Initialize simulation grid
	*/
	Grid grid;

	const unsigned int neighborhood_size = 0;
	const auto periodic = grid_options.data.get_data(pamhd::grid::Periodic());
	if (not grid.initialize(
		number_of_cells,
		comm,
		lb_name.c_str(),
		neighborhood_size,
		0,
		periodic[0],
		periodic[1],
		periodic[2]
	)) {
		std::cerr << __FILE__ << ":" << __LINE__
			<< ": Couldn't initialize grid."
			<< std::endl;
		abort();
	}

	// set grid geometry
	dccrg::Cartesian_Geometry::Parameters geom_params;
	geom_params.start = grid_options.data.get_data(pamhd::grid::Start());
	geom_params.level_0_cell_length = cell_volume;

	if (not grid.set_geometry(geom_params)) {
		std::cerr << __FILE__ << ":" << __LINE__
			<< ": Couldn't set grid geometry."
			<< std::endl;
		abort();
	}

	grid.balance_load();


	/*
	Simulate
	*/

	double
		max_dt = 0,
		simulation_time = start_time,
		next_mhd_save = save_mhd_n,
		next_rem_div_B = remove_div_B_n;

	std::vector<uint64_t>
		cells = grid.get_cells(),
		inner_cells = grid.get_local_cells_not_on_process_boundary(),
		outer_cells = grid.get_local_cells_on_process_boundary();

	// initialize MHD
	if (verbose and rank == 0) {
		cout << "Initializing MHD... " << endl;
	}

	// zero fluxes
	for (const auto& cell_id: cells) {
		auto* const cell_data = grid[cell_id];
		if (cell_data == nullptr) {
			std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
			abort();
		}
		Mas1_f(*cell_data)    =
		Mas2_f(*cell_data)    =
		Mom1_f(*cell_data)[0] =
		Mom1_f(*cell_data)[1] =
		Mom1_f(*cell_data)[2] =
		Mom2_f(*cell_data)[0] =
		Mom2_f(*cell_data)[1] =
		Mom2_f(*cell_data)[2] =
		Nrj1_f(*cell_data)    =
		Nrj2_f(*cell_data)    =
		Mag_f(*cell_data)[0]  =
		Mag_f(*cell_data)[1]  =
		Mag_f(*cell_data)[2]  = 0;
	}

	pamhd::initialize_field(initial_field, grid, cells, 0, Mag);
	pamhd::initialize_fluid(
		init_cond_fluid1,
		grid, cells, 0,
		adiabatic_index, vacuum_permeability, proton_mass,
		Mas1, Mom1, Nrj1
	);
	pamhd::initialize_fluid(
		init_cond_fluid2,
		grid, cells, 0,
		adiabatic_index, vacuum_permeability, proton_mass,
		Mas2, Mom2, Nrj2
	);
	// add magnetic field contribution to total energy densities
	for (const auto& cell_id: cells) {
		auto* const cell_data = grid[cell_id];
		if (cell_data == nullptr) {
			std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
			abort();
		}
		const auto
			total_mass = Mas1(*cell_data) + Mas2(*cell_data),
			mass_frac1 = Mas1(*cell_data) / total_mass,
			mass_frac2 = Mas2(*cell_data) / total_mass;
		Nrj1(*cell_data) += mass_frac1 * 0.5 * Mag(*cell_data).squaredNorm() / vacuum_permeability;
		Nrj2(*cell_data) += mass_frac2 * 0.5 * Mag(*cell_data).squaredNorm() / vacuum_permeability;
	}
	if (verbose and rank == 0) {
		cout << "Done initializing MHD" << endl;
	}

	pamhd::mhd::Boundary_Classifier<pamhd::mhd::Cell_Type> bdy_classifier_fluid1;
	pamhd::mhd::Boundary_Classifier<pamhd::mhd::Cell_Type_Fluid2> bdy_classifier_fluid2;
	pamhd::mhd::Boundary_Classifier<pamhd::mhd::Cell_Type_Field> bdy_classifier_fields;

	bdy_classifier_fluid1.classify<
		pamhd::mhd::Value_Boundary_Id,
		pamhd::mhd::Copy_Source
	>(
		simulation_time,
		grid,
		value_bdy_fluid1,
		copy_bdy_fluid1
	);
	bdy_classifier_fluid2.classify<
		pamhd::mhd::Value_Boundary_Id_Fluid2,
		pamhd::mhd::Copy_Source_Fluid2
	>(
		simulation_time,
		grid,
		value_bdy_fluid2,
		copy_bdy_fluid2
	);
	bdy_classifier_fields.classify<
		pamhd::mhd::Value_Boundary_Id_Field,
		pamhd::mhd::Copy_Source_Field
	>(
		simulation_time,
		grid,
		value_bdy_field,
		copy_bdy_fields
	);

	size_t simulated_steps = 0;
	while (simulation_time < end_time) {
		simulated_steps++;

		/*
		Get maximum allowed time step
		*/
		double
			// don't step over the final simulation time
			until_end = end_time - simulation_time,
			local_time_step = min(time_step_factor * max_dt, until_end),
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

		// calculate total values
		for (const auto& cell_id: outer_cells) {
			auto* const cell_data = grid[cell_id];
			if (cell_data == nullptr) {
				std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
				abort();
			}
			Mas(*cell_data) = Mas1(*cell_data) + Mas2(*cell_data);
			Mom(*cell_data) = Mom1(*cell_data) + Mom2(*cell_data);
			Nrj(*cell_data) = Nrj1(*cell_data) + Nrj2(*cell_data);
		}

		Cell::set_transfer_all(
			true,
			pamhd::mhd::MHD_State_Conservative(),
			pamhd::mhd::HD1_State(),
			pamhd::mhd::HD2_State(),
			pamhd::mhd::Cell_Type()
		);
		grid.start_remote_neighbor_copy_updates();

		for (const auto& cell_id: inner_cells) {
			auto* const cell_data = grid[cell_id];
			if (cell_data == nullptr) {
				std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
				abort();
			}
			Mas(*cell_data) = Mas1(*cell_data) + Mas2(*cell_data);
			Mom(*cell_data) = Mom1(*cell_data) + Mom2(*cell_data);
			Nrj(*cell_data) = Nrj1(*cell_data) + Nrj2(*cell_data);
		}

		pamhd::divergence::get_curl(
			inner_cells,
			grid,
			Mag,
			Cur
		);
		for (const auto& cell: inner_cells) {
			auto* const cell_data = grid[cell];
			if (cell_data == nullptr) {
				std::cerr <<  __FILE__ << "(" << __LINE__ << ")" << std::endl;
				abort();
			}
			Cur(*cell_data) /= vacuum_permeability;
		}

		double solve_max_dt = -1;
		size_t solve_index = 0;
		std::tie(
			solve_max_dt,
			solve_index
		) = pamhd::mhd::N_solve(
			mhd_solver,
			0,
			grid,
			time_step,
			adiabatic_index,
			vacuum_permeability,
			Mas1, Mas, Mom, Nrj, Mag,
			Mas1_f, Mom1_f, Nrj1_f, Mag_f,
			Cell_t,
			bdy_classifier_fluid1.normal_cell,
			bdy_classifier_fluid1.dont_solve_cell
		);
		max_dt = min(
			max_dt,
			solve_max_dt
		);

		std::tie(
			solve_max_dt,
			std::ignore
		) = pamhd::mhd::N_solve(
			mhd_solver,
			0,
			grid,
			time_step,
			adiabatic_index,
			vacuum_permeability,
			Mas2, Mas, Mom, Nrj, Mag,
			Mas2_f, Mom2_f, Nrj2_f, Mag_f,
			Cell_t,
			bdy_classifier_fluid1.normal_cell,
			bdy_classifier_fluid1.dont_solve_cell
		);
		max_dt = min(
			max_dt,
			solve_max_dt
		);

		grid.wait_remote_neighbor_copy_update_receives();

		std::tie(
			solve_max_dt,
			std::ignore
		) = pamhd::mhd::N_solve(
			mhd_solver,
			solve_index + 1,
			grid,
			time_step,
			adiabatic_index,
			vacuum_permeability,
			Mas1, Mas, Mom, Nrj, Mag,
			Mas1_f, Mom1_f, Nrj1_f, Mag_f,
			Cell_t,
			bdy_classifier_fluid1.normal_cell,
			bdy_classifier_fluid1.dont_solve_cell
		);
		max_dt = min(
			max_dt,
			solve_max_dt
		);

		std::tie(
			solve_max_dt,
			solve_index
		) = pamhd::mhd::N_solve(
			mhd_solver,
			solve_index + 1,
			grid,
			time_step,
			adiabatic_index,
			vacuum_permeability,
			Mas2, Mas, Mom, Nrj, Mag,
			Mas2_f, Mom2_f, Nrj2_f, Mag_f,
			Cell_t,
			bdy_classifier_fluid1.normal_cell,
			bdy_classifier_fluid1.dont_solve_cell
		);
		max_dt = min(
			max_dt,
			solve_max_dt
		);

		pamhd::divergence::get_curl(
			outer_cells,
			grid,
			Mag,
			Cur
		);
		for (const auto& cell: outer_cells) {
			auto* const cell_data = grid[cell];
			if (cell_data == nullptr) {
				std::cerr <<  __FILE__ << "(" << __LINE__ << ")" << std::endl;
				abort();
			}
			Cur(*cell_data) /= vacuum_permeability;
		}

		grid.wait_remote_neighbor_copy_update_sends();
		Cell::set_transfer_all(
			false,
			pamhd::mhd::MHD_State_Conservative(),
			pamhd::mhd::HD1_State(),
			pamhd::mhd::HD2_State(),
			pamhd::mhd::Cell_Type()
		);

		// transfer J for calculating additional contributions to B
		Cell::set_transfer_all(true, pamhd::mhd::Electric_Current_Density());
		grid.start_remote_neighbor_copy_updates();

		// add contribution to change of B from resistivity
		pamhd::divergence::get_curl(
			inner_cells,
			grid,
			Cur,
			Mag_res
		);
		for (const auto& cell: inner_cells) {
			auto* const cell_data = grid[cell];
			if (cell_data == nullptr) {
				std::cerr <<  __FILE__ << "(" << __LINE__ << ")" << std::endl;
				abort();
			}
			Mag_res(*cell_data) *= -resistivity;
			Mag_f(*cell_data) += Mag_res(*cell_data);
		}

		grid.wait_remote_neighbor_copy_update_receives();

		pamhd::divergence::get_curl(
			outer_cells,
			grid,
			Cur,
			Mag_res
		);
		for (const auto& cell: outer_cells) {
			auto* const cell_data = grid[cell];
			if (cell_data == nullptr) {
				std::cerr <<  __FILE__ << "(" << __LINE__ << ")" << std::endl;
				abort();
			}
			Mag_res(*cell_data) *= -resistivity;
			Mag_f(*cell_data) += Mag_res(*cell_data);
		}

		pamhd::mhd::apply_fluxes(
			grid,
			adiabatic_index,
			vacuum_permeability,
			Mas1, Mom1, Nrj1, Mag,
			Mas1_f, Mom1_f, Nrj1_f, Mag_f,
			Cell_t,
			bdy_classifier_fluid1.normal_cell,
			false
		);
		pamhd::mhd::apply_fluxes(
			grid,
			adiabatic_index,
			vacuum_permeability,
			Mas2, Mom2, Nrj2, Mag,
			Mas2_f, Mom2_f, Nrj2_f, Mag_f,
			Cell_t,
			bdy_classifier_fluid1.normal_cell,
			false
		);

		grid.wait_remote_neighbor_copy_update_sends();
		Cell::set_transfer_all(false, pamhd::mhd::Electric_Current_Density());

		simulation_time += time_step;


		/*
		Remove divergence of magnetic field
		*/

		if (remove_div_B_n > 0 and simulation_time >= next_rem_div_B) {
			next_rem_div_B += remove_div_B_n;

			if (verbose and rank == 0) {
				cout << "Removing divergence of B at time "
					<< simulation_time << "...  ";
				cout.flush();
			}

			// save old B in case div removal fails
			for (const auto& cell: cells) {
				auto* const cell_data = grid[cell];
				if (cell_data == nullptr) {
					std::cerr <<  __FILE__ << "(" << __LINE__ << "): "
						"No data for cell " << cell
						<< std::endl;
					abort();
				}

				Mag_tmp(*cell_data) = Mag(*cell_data);
			}

			Cell::set_transfer_all(
				true,
				pamhd::mhd::MHD_State_Conservative(),
				pamhd::mhd::Magnetic_Field_Divergence()
			);
			const auto div_before
				= pamhd::divergence::remove(
					cells,
					{},
					{},
					grid,
					Mag,
					Mag_div,
					[](Cell& cell_data)
						-> pamhd::mhd::Scalar_Potential_Gradient::data_type&
					{
						return cell_data[pamhd::mhd::Scalar_Potential_Gradient()];
					},
					poisson_iterations_max,
					poisson_iterations_min,
					poisson_norm_stop,
					2,
					poisson_norm_increase_max,
					false
				);
			Cell::set_transfer_all(false, pamhd::mhd::Magnetic_Field_Divergence());

			grid.update_copies_of_remote_neighbors();
			Cell::set_transfer_all(false, pamhd::mhd::MHD_State_Conservative());
			const double div_after
				= pamhd::divergence::get_divergence(
					cells,
					grid,
					Mag,
					Mag_div
				);

			// restore old B
			if (div_after > div_before) {
				if (verbose and rank == 0) {
					cout << "failed (" << div_after
						<< "), restoring previous value (" << div_before << ")."
						<< endl;
				}
				for (const auto& cell: cells) {
					auto* const cell_data = grid[cell];
					if (cell_data == nullptr) {
						std::cerr <<  __FILE__ << "(" << __LINE__ << "): "
							"No data for cell " << cell
							<< std::endl;
						abort();
					}

					Mag(*cell_data) = Mag_tmp(*cell_data);
				}
			} else {

				if (verbose and rank == 0) {
					cout << div_before << " -> " << div_after << endl;
				}

				// keep pressure/temperature constant over div removal
				for (auto& cell: cells) {
					auto* const cell_data = grid[cell];
					if (cell_data == nullptr) {
						std::cerr <<  __FILE__ << "(" << __LINE__ << "): "
							"No data for cell " << cell
							<< std::endl;
						abort();
					}

					const auto mag_nrj_diff
						= (
							Mag(*cell_data).squaredNorm()
							- Mag_tmp(*cell_data).squaredNorm()
						) / (2 * vacuum_permeability);

					const auto
						total_mass = Mas1(*cell_data) + Mas2(*cell_data),
						mass_frac1 = Mas1(*cell_data) / total_mass,
						mass_frac2 = Mas2(*cell_data) / total_mass;
					Nrj1(*cell_data) += mass_frac1 * mag_nrj_diff;
					Nrj2(*cell_data) += mass_frac2 * mag_nrj_diff;
				}
			}
		}


		// set boundary data
		const pamhd::mhd::Number_Density N{};
		const pamhd::mhd::Velocity V{};
		const pamhd::mhd::Pressure P{};
		const pamhd::mhd::Magnetic_Field B{};

		const auto& cell_data_pointers = grid.get_cell_data_pointers();

		// value boundaries
		for (const auto& cell_item: cell_data_pointers) {
			const auto& cell_id = get<0>(cell_item);
			if (cell_id == dccrg::error_cell) {
				continue;
			}

			// skip neighbor cells
			const auto& offset = get<2>(cell_item);
			if (offset[0] != 0 or offset[1] != 0 or offset[2] != 0) {
				continue;
			}

			auto* const cell_data = get<1>(cell_item);

			// fluid 1
			if ((*cell_data)[pamhd::mhd::Cell_Type()] == bdy_classifier_fluid1.value_boundary_cell) {

				const auto& bdy_id = (*cell_data)[pamhd::mhd::Value_Boundary_Id()];
				const auto cell_center = grid.geometry.get_center(cell_id);

				const auto mass_density
					= value_bdy_fluid1.get_data(N, bdy_id, cell_center, simulation_time)
					* proton_mass;
				const auto velocity
					= value_bdy_fluid1.get_data(V, bdy_id, cell_center, simulation_time);
				const auto pressure
					= value_bdy_fluid1.get_data(P, bdy_id, cell_center, simulation_time);

				Mas1(*cell_data) = mass_density;
				Mas(*cell_data) = Mas1(*cell_data) + Mas2(*cell_data);
				Mom1(*cell_data) = mass_density * velocity;
				if (mass_density > 0 and pressure > 0) {
					Nrj1(*cell_data) = pamhd::mhd::get_total_energy_density(
						mass_density,
						velocity,
						pressure,
						std::array<double, 3>{0, 0, 0},
						adiabatic_index,
						vacuum_permeability
					);
				} else {
					Nrj1(*cell_data) = 0;
				}
			}

			// fluid 2
			if ((*cell_data)[pamhd::mhd::Cell_Type_Fluid2()] == bdy_classifier_fluid2.value_boundary_cell) {

				const auto& bdy_id = (*cell_data)[pamhd::mhd::Value_Boundary_Id_Fluid2()];
				const auto cell_center = grid.geometry.get_center(cell_id);

				const auto mass_density
					= value_bdy_fluid2.get_data(N, bdy_id, cell_center, simulation_time)
					* proton_mass;
				const auto velocity
					= value_bdy_fluid2.get_data(V, bdy_id, cell_center, simulation_time);
				const auto pressure
					= value_bdy_fluid2.get_data(P, bdy_id, cell_center, simulation_time);

				Mas2(*cell_data) = mass_density;
				Mas(*cell_data) = Mas1(*cell_data) + Mas2(*cell_data);
				Mom2(*cell_data) = mass_density * velocity;
				if (mass_density > 0 and pressure > 0) {
					Nrj2(*cell_data) = pamhd::mhd::get_total_energy_density(
						mass_density,
						velocity,
						pressure,
						std::array<double, 3>{0, 0, 0},
						adiabatic_index,
						vacuum_permeability
					);
				} else {
					Nrj2(*cell_data) = 0;
				}
			}

			// field(s)
			if ((*cell_data)[pamhd::mhd::Cell_Type_Field()] == bdy_classifier_fields.value_boundary_cell) {
				const auto& bdy_id = (*cell_data)[pamhd::mhd::Value_Boundary_Id_Field()];
				const auto cell_center = grid.geometry.get_center(cell_id);
				const auto magnetic_field
					= value_bdy_field.get_data(B, bdy_id, cell_center, simulation_time);
				Mag(*cell_data) = magnetic_field;
			}
		}

		/*
		Add magnetic nrj to fluid boundary cells
		*/
		for (const auto& cell_item: cell_data_pointers) {
			const auto& cell_id = get<0>(cell_item);
			if (cell_id == dccrg::error_cell) {
				continue;
			}

			// skip neighbor cells
			const auto& offset = get<2>(cell_item);
			if (offset[0] != 0 or offset[1] != 0 or offset[2] != 0) {
				continue;
			}

			auto* const cell_data = get<1>(cell_item);

			if ((*cell_data)[pamhd::mhd::Cell_Type()] == bdy_classifier_fluid1.value_boundary_cell) {
				const auto mass_frac = Mas1(*cell_data) / Mas(*cell_data);
				Nrj1(*cell_data) += mass_frac * 0.5 * Mag(*cell_data).squaredNorm() / vacuum_permeability;
			}

			if ((*cell_data)[pamhd::mhd::Cell_Type_Fluid2()] == bdy_classifier_fluid2.value_boundary_cell) {
				const auto mass_frac = Mas2(*cell_data) / Mas(*cell_data);
				Nrj2(*cell_data) += mass_frac * 0.5 * Mag(*cell_data).squaredNorm() / vacuum_permeability;
			}
		}

		// copy up-to-date data
		Cell::set_transfer_all(
			true,
			pamhd::mhd::MHD_State_Conservative(),
			pamhd::mhd::HD1_State(),
			pamhd::mhd::HD2_State()
		);
		grid.update_copies_of_remote_neighbors();
		Cell::set_transfer_all(
			false,
			pamhd::mhd::MHD_State_Conservative(),
			pamhd::mhd::HD1_State(),
			pamhd::mhd::HD2_State()
		);

		// copy boundaries
		for (const auto& cell_item: cell_data_pointers) {
			const auto& cell_id = get<0>(cell_item);
			if (cell_id == dccrg::error_cell) {
				continue;
			}

			const auto& offset = get<2>(cell_item);
			if (offset[0] != 0 or offset[1] != 0 or offset[2] != 0) {
				continue;
			}

			auto* const target_data = get<1>(cell_item);

			if ((*target_data)[pamhd::mhd::Cell_Type()] == bdy_classifier_fluid1.copy_boundary_cell) {
				const auto& source_id = (*target_data)[pamhd::mhd::Copy_Source()];
				auto* const source_data = grid[source_id];
				if (source_data == nullptr) {
					std::cerr <<  __FILE__ << "(" << __LINE__ << "): "
						<< "No data for source cell " << source_id
						<< " of cell " << cell_id
						<< std::endl;
					abort();
				}

				Mas1(*target_data) = Mas1(*source_data);
				Mas(*target_data) = Mas1(*target_data) + Mas2(*target_data);
				Mom1(*target_data) = Mom1(*source_data);
				Mom(*target_data) = Mom1(*target_data) + Mom2(*target_data);
				Nrj1(*target_data) = Nrj1(*source_data);
				Nrj(*target_data) = Nrj1(*target_data) + Nrj2(*target_data);
			}

			if ((*target_data)[pamhd::mhd::Cell_Type_Fluid2()] == bdy_classifier_fluid2.copy_boundary_cell) {
				const auto& source_id = (*target_data)[pamhd::mhd::Copy_Source_Fluid2()];
				auto* const source_data = grid[source_id];
				if (source_data == nullptr) {
					std::cerr <<  __FILE__ << "(" << __LINE__ << "): "
						<< "No data for source cell: " << source_id
						<< " of cell " << cell_id
						<< std::endl;
					abort();
				}

				Mas2(*target_data) = Mas2(*source_data);
				Mas(*target_data) = Mas1(*target_data) + Mas2(*target_data);
				Mom2(*target_data) = Mom2(*source_data);
				Mom(*target_data) = Mom1(*target_data) + Mom2(*target_data);
				Nrj2(*target_data) = Nrj2(*source_data);
				Nrj(*target_data) = Nrj1(*target_data) + Nrj2(*target_data);
			}

			if ((*target_data)[pamhd::mhd::Cell_Type_Field()] == bdy_classifier_fields.copy_boundary_cell) {
				const auto& source_id = (*target_data)[pamhd::mhd::Cell_Type_Field()];
				auto* const source_data = grid[source_id];
				if (source_data == nullptr) {
					std::cerr <<  __FILE__ << "(" << __LINE__ << "): "
						<< "No data for source cell: " << source_id
						<< " of cell " << cell_id
						<< std::endl;
					abort();
				}
				Mag(*target_data) = Mag(*source_data);
			}
		}


		/*
		Save simulation to disk
		*/

		if (
			(save_mhd_n >= 0 and (simulation_time == start_time or simulation_time >= end_time))
			or (save_mhd_n > 0 and simulation_time >= next_mhd_save)
		) {
			if (next_mhd_save <= simulation_time) {
				next_mhd_save += save_mhd_n;
			}

			if (verbose and rank == 0) {
				cout << "Saving (M)HD at time " << simulation_time << endl;
			}

			if (
				not pamhd::mhd::save(
					boost::filesystem::canonical(
						boost::filesystem::path(output_directory)
					).append("mhd_").generic_string(),
					grid,
					1,
					simulation_time,
					adiabatic_index,
					proton_mass,
					vacuum_permeability,
					pamhd::mhd::MHD_State_Conservative(),
					pamhd::mhd::Electric_Current_Density(),
					pamhd::mhd::Cell_Type(),
					pamhd::mhd::MPI_Rank(),
					pamhd::mhd::Resistivity()
				)
			) {
				std::cerr <<  __FILE__ << "(" << __LINE__ << "): "
					"Couldn't save mhd result."
					<< std::endl;
				abort();
			}

			if (
				not pamhd::mhd::save(
					boost::filesystem::canonical(
						boost::filesystem::path(output_directory)
					).append("hd1_").generic_string(),
					grid,
					1,
					simulation_time,
					adiabatic_index,
					proton_mass,
					vacuum_permeability,
					pamhd::mhd::HD1_State()
				)
			) {
				std::cerr <<  __FILE__ << "(" << __LINE__ << "): "
					"Couldn't save fluid1 result."
					<< std::endl;
				abort();
			}

			if (
				not pamhd::mhd::save(
					boost::filesystem::canonical(
						boost::filesystem::path(output_directory)
					).append("hd2_").generic_string(),
					grid,
					1,
					simulation_time,
					adiabatic_index,
					proton_mass,
					vacuum_permeability,
					pamhd::mhd::HD2_State()
				)
			) {
				std::cerr <<  __FILE__ << "(" << __LINE__ << "): "
					"Couldn't save fluid2 result."
					<< std::endl;
				abort();
			}
		}

	}

	if (verbose and rank == 0) {
		cout << "Simulation finished at time " << simulation_time << endl;
	}

	MPI_Finalize();

	return EXIT_SUCCESS;
}
