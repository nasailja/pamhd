/*
Propagates test particles (0 mass) in prescribed electric and magnetic fields.

Copyright 2015 Ilja Honkonen
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

#include "array"
#include "cmath"
#include "cstdlib"
#include "fstream"
#include "iostream"
#include "random"
#include "string"

#include "boost/numeric/odeint.hpp"
#include "dccrg.hpp"
#include "dccrg_cartesian_geometry.hpp"
#include "Eigen/Core" // must be included before gensimcell
#include "Eigen/Geometry"
#include "mpi.h" // must be included before gensimcell
#include "gensimcell.hpp"


#include "boundaries/initial_condition.hpp"
#include "grid_options.hpp"
#include "particle/save.hpp"
#include "particle/solve_dccrg.hpp"
#include "particle/variables.hpp"

using namespace std;
using namespace pamhd::particle;


/*
Variables for bulk particle parameters of initial and value boundaries.
*/
struct Mass_Density {
	using data_type = double;
	static std::string get_name() { return std::string("Bulk mass density"); }
	static std::string get_option_name() { return std::string("bulk-mass-density"); }
	static std::string get_option_help() { return std::string(""); }
};
struct Temperature {
	using data_type = double;
	static std::string get_name() { return std::string("Bulk temperature"); }
	static std::string get_option_name() { return std::string("bulk-temperature"); }
	static std::string get_option_help() { return std::string(""); }
};
struct Bulk_Velocity {
	using data_type = Eigen::Vector3d;
	static std::string get_name() { return std::string("Bulk velocity"); }
	static std::string get_option_name() { return std::string("bulk-velocity"); }
	static std::string get_option_help() { return std::string(""); }
};
struct Species_Mass {
	using data_type = double;
	static std::string get_name() { return std::string("Species mass"); }
	static std::string get_option_name() { return std::string("species-mass"); }
	static std::string get_option_help() { return std::string("Mass of individual constituents of particle species"); }
};
// number of macroparticles / simulation cell
struct Nr_Particles_In_Cell {
	using data_type = double;
	static std::string get_name() { return std::string("Nr particles in cell"); }
	static std::string get_option_name() { return std::string("nr-particles-in-cell"); }
	static std::string get_option_help() { return std::string(""); }
};


/*
Initializes electric and magnetic fields in given grid.

Fields are initialized either from given file or from
given options of given file == "".

File format (ascii data):
x, y and z coordinate of 1st E & B location in meters,
x, y and z components of 1st electric field in volt/meter,
x, y and z components of 1st magnetic field in tesla,
x, y and z coordinates of 2nd E & B location in meters,
...
Each local cell in given grid is initialized from the closest value
in given file.
*/
template<class Init_Cond_T> void initialize_fields(
	Init_Cond_T& init_cond,
	const std::vector<uint64_t>& cell_ids,
	dccrg::Dccrg<pamhd::particle::Cell, dccrg::Cartesian_Geometry>& grid,
	const std::string& fields_file
) {
	using std::pow;

	// set initial fields from file
	if (fields_file != "") {
		std::ifstream infile(fields_file);

		/*
		Read in file data
		*/

		std::vector<
			// r, E, B
			std::array<
				Eigen::Vector3d,
				3
			>
		> file_data;
		file_data.reserve(cell_ids.size()); // guess

		while (infile) {
			Eigen::Vector3d r, E, B;

			infile >> r[0];

			if (infile.eof()) {
				break;
			} else if (not infile.good()) {
				std::cerr << "Couldn't extract first coordinate of location from "
					<< fields_file << " after " << file_data.size() << " items"
					<< std::endl;
				abort();
			}

			infile >> r[1] >> r[2];
			if (not infile.good()) {
				std::cerr << "Couldn't extract other coordinates from "
					<< fields_file << " after " << file_data.size() << " items"
					<< std::endl;
				abort();
			}

			infile >> E[0] >> E[1] >> E[2];
			if (not infile.good()) {
				std::cerr << "Couldn't extract electric field from "
					<< fields_file << " after " << file_data.size() << " items"
					<< std::endl;
				abort();
			}

			infile >> B[0] >> B[1] >> B[2];
			if (not infile.good()) {
				std::cerr << "Couldn't extract magnetic field from "
					<< fields_file << " after " << file_data.size() << " items"
					<< std::endl;
				abort();
			}

			file_data.push_back({{r, E, B}});
		}

		// set field values in local cells
		for (const auto& cell_id: cell_ids) {
			const auto cell_center = grid.geometry.get_center(cell_id);

			const auto closest_item = std::min_element(
				file_data.cbegin(),
				file_data.cend(),
				[&cell_center](
					const std::array<Eigen::Vector3d, 3>& a,
					const std::array<Eigen::Vector3d, 3>& b
				) {
					const auto
						a_distance2
							= pow(a[0][0] - cell_center[0], 2)
							+ pow(a[0][1] - cell_center[1], 2)
							+ pow(a[0][2] - cell_center[2], 2),
						b_distance2
							= pow(b[0][0] - cell_center[0], 2)
							+ pow(b[0][1] - cell_center[1], 2)
							+ pow(b[0][2] - cell_center[2], 2);

					return a_distance2 < b_distance2;
				}
			);

			if (closest_item == file_data.cend()) {
				std::cerr << "Couldn't find field data item closest to center of cell "
					<< cell_id
					<< std::endl;
				abort();
			}

			auto* const cell_data = grid[cell_id];
			if (cell_data == nullptr) {
				std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
				abort();
			}

			(*cell_data)[Electric_Field()] = (*closest_item)[1];
			(*cell_data)[Magnetic_Field()] = (*closest_item)[2];
		}

	// set initial fields from command line
	} else {

		for (const auto cell_id: cell_ids) {
			const auto
				cell_start = grid.geometry.get_min(cell_id),
				cell_end = grid.geometry.get_max(cell_id),
				cell_center = grid.geometry.get_center(cell_id);

			init_cond.add_cell(
				cell_id,
				cell_start,
				cell_end
			);

			auto* const cell_data = grid[cell_id];
			if (cell_data == NULL) {
				std::cerr <<  __FILE__ << "(" << __LINE__ << ") No data for cell: "
					<< cell_id
					<< std::endl;
				abort();
			}

			// set default state
			(*cell_data)[Electric_Field()]
				= init_cond.default_data.get_data(
					Electric_Field(),
					cell_center,
					0
				);
			(*cell_data)[Magnetic_Field()]
				= init_cond.default_data.get_data(
					Magnetic_Field(),
					cell_center,
					0
				);
		}
	}
}


/*
Creates particles in given cells as defined initial conditions obtained from user.
*/
template<class Init_Cond_T> void initialize_particles(
	Init_Cond_T& init_cond,
	const std::vector<uint64_t>& cells,
	dccrg::Dccrg<pamhd::particle::Cell, dccrg::Cartesian_Geometry>& grid,
	std::mt19937& random_source,
	const double particle_temp_nrj_ratio
) {
	if (grid.get_rank() == 0) {
		std::cout << "Setting default particle state... ";
		std::cout.flush();
	}
	for (const auto cell_id: cells) {
		const auto
			cell_start = grid.geometry.get_min(cell_id),
			cell_end = grid.geometry.get_max(cell_id),
			cell_length = grid.geometry.get_length(cell_id);

		init_cond.add_cell(
			cell_id,
			cell_start,
			cell_end
		);

		auto* const cell_data = grid[cell_id];
		if (cell_data == NULL) {
			std::cerr <<  __FILE__ << "(" << __LINE__ << ") No data for cell: "
				<< cell_id
				<< std::endl;
			abort();
		}

		// set default state
		const auto
			mass_density = init_cond.default_data.get_data(Mass_Density(), {{0, 0, 0}}, 0),
			temperature = init_cond.default_data.get_data(Temperature(), {{0, 0, 0}}, 0),
			charge_mass_ratio = init_cond.default_data.get_data(Charge_Mass_Ratio(), {{0, 0, 0}}, 0),
			species_mass = init_cond.default_data.get_data(Species_Mass(), {{0, 0, 0}}, 0);
		const auto bulk_velocity = init_cond.default_data.get_data(Bulk_Velocity(), {{0, 0, 0}}, 0);
		const auto nr_particles = init_cond.default_data.get_data(Nr_Particles_In_Cell(), {{0, 0, 0}}, 0);
		
		(*cell_data)[Particles_Internal()]
			= create_particles<
				Particle_T<>,
				Mass,
				Charge_Mass_Ratio,
				Position,
				Velocity,
				std::mt19937
			>(
				bulk_velocity,
				Eigen::Vector3d{cell_start[0], cell_start[1], cell_start[2]},
				Eigen::Vector3d{cell_end[0], cell_end[1], cell_end[2]},
				Eigen::Vector3d{temperature / 3, temperature / 3, temperature / 3},
				nr_particles,
				charge_mass_ratio,
				mass_density * cell_length[0] * cell_length[1] * cell_length[2],
				species_mass,
				particle_temp_nrj_ratio,
				random_source
			);
	}

	// set non-default initial conditions
	if (grid.get_rank() == 0) {
		std::cout << "done\nSetting non-default initial particle state... ";
		std::cout.flush();
	}
	for (size_t bdy_i = 0; bdy_i < init_cond.get_number_of_boundaries(); bdy_i++) {
		const auto
			mass_density = init_cond.get_data(Mass_Density(), bdy_i, {{0, 0, 0}}, 0),
			temperature = init_cond.get_data(Temperature(), bdy_i, {{0, 0, 0}}, 0),
			charge_mass_ratio = init_cond.get_data(Charge_Mass_Ratio(), bdy_i, {{0, 0, 0}}, 0),
			species_mass = init_cond.get_data(Species_Mass(), bdy_i, {{0, 0, 0}}, 0);
		const auto bulk_velocity = init_cond.get_data(Bulk_Velocity(), bdy_i, {{0, 0, 0}}, 0);
		const auto nr_particles = init_cond.get_data(Nr_Particles_In_Cell(), bdy_i, {{0, 0, 0}}, 0);

		for (const auto& cell_id: init_cond.get_cells(bdy_i)) {
			const auto
				cell_start = grid.geometry.get_min(cell_id),
				cell_end = grid.geometry.get_max(cell_id),
				cell_length = grid.geometry.get_length(cell_id);

			auto* const cell_data = grid[cell_id];
			if (cell_data == NULL) {
				std::cerr <<  __FILE__ << "(" << __LINE__ << ") No data for cell: "
					<< cell_id
					<< std::endl;
				abort();
			}

			(*cell_data)[Particles_Internal()]
				= create_particles<
					Particle_T<>,
					Mass,
					Charge_Mass_Ratio,
					Position,
					Velocity,
					std::mt19937
				>(
					bulk_velocity,
					Eigen::Vector3d{cell_start[0], cell_start[1], cell_start[2]},
					Eigen::Vector3d{cell_end[0], cell_end[1], cell_end[2]},
					Eigen::Vector3d{temperature / 3, temperature / 3, temperature / 3},
					nr_particles,
					charge_mass_ratio,
					mass_density * cell_length[0] * cell_length[1] * cell_length[2],
					species_mass,
					particle_temp_nrj_ratio,
					random_source
				);
		}
	}
	if (grid.get_rank() == 0) {
		std::cout << "done" << std::endl;
	}
}


int main(int argc, char* argv[])
{
	using Cell_T = pamhd::particle::Cell;
	using Grid_T = dccrg::Dccrg<Cell_T, dccrg::Cartesian_Geometry>;

	constexpr double Re = 6.371e6; // radius of earth

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
		Mass_Density,
		Temperature,
		Bulk_Velocity,
		Nr_Particles_In_Cell,
		Charge_Mass_Ratio,
		Species_Mass
	> initial_particles;

	pamhd::boundaries::Initial_Condition<
		uint64_t,
		double,
		std::array<double, 3>,
		Electric_Field,
		Magnetic_Field
	> initial_fields;


	pamhd::grid::Options grid_options;
	grid_options.data.set_expression(pamhd::grid::Number_Of_Cells(), "{1, 1, 1}");
	grid_options.data.set_expression(pamhd::grid::Volume(), "{1, 1, 1}");
	grid_options.data.set_expression(pamhd::grid::Start(), "{0, 0, 0}");
	grid_options.data.set_expression(pamhd::grid::Periodic(), "{false, false, false}");

	bool verbose = false;
	double
		save_particle_n = -1,
		start_time = 0,
		end_time = 1,
		time_step_factor = 0.5,
		adiabatic_index = 5.0 / 3.0,    
		vacuum_permeability = 4e-7 * M_PI,
		particle_temp_nrj_ratio = 1.3806488e-23;

	std::string
		config_file_name(""),
		lb_name("RCB"),
		output_directory(""),
		fields_file_name("");

	boost::program_options::options_description
		options(
			"Usage: program_name [options], where options are"
		),
		// grouped options for printing help
		initial_condition_help(
			"Options for initial condition which sets cell data to values given "
			"in --config-file at start of simulation (not read if empty strging)"
		);

	// handle general options
	options.add_options()
		("help", "Print this help message")
		("verbose", "Print run time information")
		("initial-help", "Print help for initial condition options")
		("config-file",
			boost::program_options::value<std::string>(&config_file_name)
				->default_value(config_file_name),
			"Read options also from file arg (command line has priority, "
			"not read if empty string)")
		("fields-file",
			boost::program_options::value<std::string>(&fields_file_name)
				->default_value(fields_file_name),
			"Initialize electric and magnetic fields from "
			"file arg instead of config-file")
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
		("save-particle-n",
			boost::program_options::value<double>(&save_particle_n)
				->default_value(save_particle_n),
			"Save results every arg seconds, 0 saves "
			"initial and final states, -1 doesn't save")
		("boltzmann",
			boost::program_options::value<double>(&particle_temp_nrj_ratio)
				->default_value(particle_temp_nrj_ratio),
			"Ratio of energy and temperature when calculating particle velocities (J/K, Boltzmann constant)");


	grid_options.add_options("grid.", options);
	initial_particles.add_initialization_options("initial-particles.", options);
	initial_fields.add_initialization_options("initial-fields.", options);

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

	initial_particles.add_options("initial-particles.", options);
	initial_particles.add_options("initial-particles.", initial_condition_help);
	initial_fields.add_options("initial-fields.", options);
	initial_fields.add_options("initial-fields.", initial_condition_help);

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

	if (option_variables.count("verbose") > 0) {
		verbose = true;
	}


	Grid_T grid;

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

	const unsigned int neighborhood_size = 1;
	const auto periodic = grid_options.data.get_data(pamhd::grid::Periodic());
	if (
		not grid.initialize(
			number_of_cells,
			comm,
			lb_name.c_str(),
			neighborhood_size,
			0,
			periodic[0],
			periodic[1],
			periodic[2]
		)
	) {
		std::cerr << __FILE__ << ":" << __LINE__
			<< ": Couldn't initialize grid."
			<< std::endl;
		abort();
	}

	dccrg::Cartesian_Geometry::Parameters geom_params;
	geom_params.start = grid_options.data.get_data(pamhd::grid::Start());
	geom_params.level_0_cell_length = cell_volume;

	if (not grid.set_geometry(geom_params)) {
		std::cerr << __FILE__ << "(" << __LINE__ << "): "
			<< "Couldn't set geometry of one or more grids."
			<< std::endl;
		abort();
	}

	grid.balance_load();


	const auto
		inner_cell_ids = grid.get_local_cells_not_on_process_boundary(),
		outer_cell_ids = grid.get_local_cells_on_process_boundary(),
		remote_cell_ids = grid.get_remote_cells_on_process_boundary(),
		cell_ids = grid.get_cells();

	// set initial condition
	std::mt19937 random_source;
	random_source.seed(grid.get_rank());

	initialize_fields(
		initial_fields,
		cell_ids,
		grid,
		fields_file_name
	);

	initialize_particles(
		initial_particles,
		cell_ids,
		grid,
		random_source,
		particle_temp_nrj_ratio
	);

	// allocate copies of remote neighbor cells
	grid.update_copies_of_remote_neighbors();

	double
		max_dt = 0,
		next_particle_save = save_particle_n,
		simulation_time = start_time;

	size_t simulated_steps = 0;
	while (simulation_time < end_time) {
		simulated_steps++;

		double
			// don't step over the final simulation time
			until_end = end_time - simulation_time,
			local_time_step = std::min(0.5 * max_dt, until_end),
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

		if (grid.get_rank() == 0) {
			cout << "Solving particles at time " << simulation_time
				<< " s with time step " << time_step << " s" << endl;
		}

		max_dt = std::numeric_limits<double>::max();

		max_dt = std::min(
			max_dt,
			pamhd::particle::solve<
				Cell_T,
				pamhd::particle::Electric_Field,
				pamhd::particle::Magnetic_Field,
				pamhd::particle::Nr_Particles_External,
				pamhd::particle::Particles_Internal,
				pamhd::particle::Particles_External,
				pamhd::particle::Position,
				pamhd::particle::Velocity,
				pamhd::particle::Charge_Mass_Ratio,
				pamhd::particle::Mass,
				pamhd::particle::Destination_Cell,
				boost::numeric::odeint::runge_kutta_fehlberg78<state_t>
			>(time_step, outer_cell_ids, grid)
		);

		Cell_T::set_transfer_all(
			true,
			pamhd::particle::Electric_Field(),
			pamhd::particle::Magnetic_Field(),
			pamhd::particle::Nr_Particles_External()
		);
		grid.start_remote_neighbor_copy_updates();

		max_dt = std::min(
			max_dt,
			pamhd::particle::solve<
				Cell_T,
				pamhd::particle::Electric_Field,
				pamhd::particle::Magnetic_Field,
				pamhd::particle::Nr_Particles_External,
				pamhd::particle::Particles_Internal,
				pamhd::particle::Particles_External,
				pamhd::particle::Position,
				pamhd::particle::Velocity,
				pamhd::particle::Charge_Mass_Ratio,
				pamhd::particle::Mass,
				pamhd::particle::Destination_Cell,
				boost::numeric::odeint::runge_kutta_fehlberg78<state_t>
			>(time_step, inner_cell_ids, grid)
		);

		simulation_time += time_step;

		grid.wait_remote_neighbor_copy_update_receives();
		pamhd::particle::resize_receiving_containers<
			Cell_T,
			pamhd::particle::Nr_Particles_External,
			pamhd::particle::Particles_External
		>(remote_cell_ids, grid);

		grid.wait_remote_neighbor_copy_update_sends();

		Cell_T::set_transfer_all(
			false,
			pamhd::particle::Electric_Field(),
			pamhd::particle::Magnetic_Field(),
			pamhd::particle::Nr_Particles_External()
		);
		Cell_T::set_transfer_all(
			true,
			pamhd::particle::Particles_External()
		);

		grid.start_remote_neighbor_copy_updates();

		pamhd::particle::incorporate_external_particles<
			Cell_T,
			pamhd::particle::Nr_Particles_External,
			pamhd::particle::Particles_Internal,
			pamhd::particle::Particles_External,
			pamhd::particle::Destination_Cell
		>(inner_cell_ids, grid);

		grid.wait_remote_neighbor_copy_update_receives();

		pamhd::particle::incorporate_external_particles<
			Cell_T,
			pamhd::particle::Nr_Particles_External,
			pamhd::particle::Particles_Internal,
			pamhd::particle::Particles_External,
			pamhd::particle::Destination_Cell
		>(outer_cell_ids, grid);

		pamhd::particle::remove_external_particles<
			Cell_T,
			pamhd::particle::Nr_Particles_External,
			pamhd::particle::Particles_External
		>(inner_cell_ids, grid);

		grid.wait_remote_neighbor_copy_update_sends();
		Cell_T::set_transfer_all(
			false,
			pamhd::particle::Particles_External()
		);

		pamhd::particle::remove_external_particles<
			Cell_T,
			pamhd::particle::Nr_Particles_External,
			pamhd::particle::Particles_External
		>(outer_cell_ids, grid);


		if (
			(save_particle_n >= 0 and (simulation_time == 0 or simulation_time >= end_time))
			or (save_particle_n > 0 and simulation_time >= next_particle_save)
		) {
			if (next_particle_save <= simulation_time) {
				next_particle_save += save_particle_n;
			}

			if (rank == 0) {
				cout << "Saving particles at time " << simulation_time << endl;
			}

			if (
				not Save::save<
					Grid_T,
					Cell_T
				>(
					"",
					grid,
					simulation_time,
					vacuum_permeability
				)
			) {
				std::cerr <<  __FILE__ << "(" << __LINE__ << "): "
					"Couldn't save particle result."
					<< std::endl;
				abort();
			}
		}
	}

	MPI_Finalize();

	return EXIT_SUCCESS;
}
