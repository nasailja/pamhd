/*
MHD test program of PAMHD.

Copyright 2014, 2015, 2016 Ilja Honkonen
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
#include "fstream"
#include "iostream"
#include "string"
#include "fstream"
#include "iomanip"
#include "iostream"
#include "limits"
#include "streambuf"
#include "string"
#include "type_traits"

#include "boost/filesystem.hpp"
#include "boost/lexical_cast.hpp"
#include "dccrg.hpp"
#include "dccrg_cartesian_geometry.hpp"
#include "Eigen/Core" // must be included before gensimcell.hpp
#include "mpi.h" // must be included before gensimcell.hpp
#include "gensimcell.hpp"
#include "rapidjson/document.h"
#include "rapidjson/error/en.h"

#include "boundaries/geometries.hpp"
#include "boundaries/multivariable_boundaries.hpp"
#include "boundaries/multivariable_initial_conditions.hpp"
#include "divergence/remove.hpp"
#include "grid_options.hpp"
#include "mhd/background_magnetic_field.hpp"
#include "mhd/boundaries.hpp"
#include "mhd/common.hpp"
#include "mhd/hll_athena.hpp"
#include "mhd/hlld_athena.hpp"
#include "mhd/initialize.hpp"
#include "mhd/options.hpp"
#include "mhd/roe_athena.hpp"
#include "mhd/rusanov.hpp"
#include "mhd/save.hpp"
#include "mhd/solve.hpp"
#include "mhd/variables.hpp"


using namespace std;

/*
Controls transfer of variables in poisson solver
which doesn't use generic cell
*/
int Poisson_Cell::transfer_switch = Poisson_Cell::INIT;

using Cell = pamhd::mhd::Cell;
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

// references to background magnetic fields
const auto Bg_B_Pos_X
	= [](Cell& cell_data)->typename pamhd::mhd::Bg_Magnetic_Field_Pos_X::data_type&{
		return cell_data[pamhd::mhd::Bg_Magnetic_Field_Pos_X()];
	};
const auto Bg_B_Pos_Y
	= [](Cell& cell_data)->typename pamhd::mhd::Bg_Magnetic_Field_Pos_Y::data_type&{
		return cell_data[pamhd::mhd::Bg_Magnetic_Field_Pos_Y()];
	};
const auto Bg_B_Pos_Z
	= [](Cell& cell_data)->typename pamhd::mhd::Bg_Magnetic_Field_Pos_Z::data_type&{
		return cell_data[pamhd::mhd::Bg_Magnetic_Field_Pos_Z()];
	};

// solver info variable for boundary logic
const auto Sol_Info
	= [](Cell& cell_data)->typename pamhd::mhd::Solver_Info::data_type&{
		return cell_data[pamhd::mhd::Solver_Info()];
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
// electrical resistivity
const auto Res
	= [](Cell& cell_data)->typename pamhd::mhd::Resistivity::data_type&{
		return cell_data[pamhd::mhd::Resistivity()];
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

// flux / total change of mass density over one time step
const auto Mas_f
	= [](Cell& cell_data)->typename pamhd::mhd::Mass_Density::data_type&{
		return cell_data[pamhd::mhd::MHD_Flux_Conservative()][pamhd::mhd::Mass_Density()];
	};
const auto Mom_f
	= [](Cell& cell_data)->typename pamhd::mhd::Momentum_Density::data_type&{
		return cell_data[pamhd::mhd::MHD_Flux_Conservative()][pamhd::mhd::Momentum_Density()];
	};
const auto Nrj_f
	= [](Cell& cell_data)->typename pamhd::mhd::Total_Energy_Density::data_type&{
		return cell_data[pamhd::mhd::MHD_Flux_Conservative()][pamhd::mhd::Total_Energy_Density()];
	};
const auto Mag_f
	= [](Cell& cell_data)->typename pamhd::mhd::Magnetic_Field::data_type&{
		return cell_data[pamhd::mhd::MHD_Flux_Conservative()][pamhd::mhd::Magnetic_Field()];
	};


int main(int argc, char* argv[])
{
	using std::asin;
	using std::atan2;
	using std::get;
	using std::min;
	using std::sqrt;


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


	// read and parse json data from configuration file
	if (argc != 2) {
		if (argc < 2 and rank == 0) {
			std::cerr
				<< "Name of configuration file required."
				<< std::endl;
		}
		if (argc > 2 and rank == 0) {
			std::cerr
				<< "Too many arguments given to " << argv[0]
				<< ": " << argc - 1 << ", should be 1"
				<< std::endl;
		}
		MPI_Finalize();
		return EXIT_FAILURE;
	}

	std::ifstream json_file(argv[1]);
	if (not json_file.good()) {
		if (rank == 0) {
			std::cerr << "Couldn't open configuration file " << argv[1] << std::endl;
		}
		MPI_Finalize();
		return EXIT_FAILURE;
	}

	std::string json{
		std::istreambuf_iterator<char>(json_file),
		std::istreambuf_iterator<char>()
	};

	rapidjson::Document document;
	document.Parse(json.c_str());
	if (document.HasParseError()) {
		std::cerr << "Couldn't parse json data in file " << argv[1]
			<< " at character position " << document.GetErrorOffset()
			<< ": " << rapidjson::GetParseError_En(document.GetParseError())
			<< std::endl;
		MPI_Finalize();
		return EXIT_FAILURE;
	}

	// mhd options
	pamhd::mhd::Options options_mhd{document};

	if (rank == 0 and options_mhd.output_directory != "") {
		try {
			boost::filesystem::create_directories(options_mhd.output_directory);
		} catch (const boost::filesystem::filesystem_error& e) {
			std::cerr <<  __FILE__ << "(" << __LINE__ << ") "
				"Couldn't create output directory "
				<< options_mhd.output_directory << ": "
				<< e.what()
				<< std::endl;
			abort();
		}
	}

	using geometry_id_t = unsigned int;

	pamhd::boundaries::Geometries<
		geometry_id_t,
		std::array<double, 3>,
		double,
		uint64_t
	> geometries;
	geometries.set(document);

	pamhd::boundaries::Multivariable_Initial_Conditions<
		geometry_id_t,
		pamhd::mhd::Number_Density,
		pamhd::mhd::Velocity,
		pamhd::mhd::Pressure,
		pamhd::mhd::Magnetic_Field
	> initial_conditions;
	initial_conditions.set(document);

	pamhd::boundaries::Multivariable_Boundaries<
		uint64_t,
		geometry_id_t,
		pamhd::mhd::Number_Density,
		pamhd::mhd::Velocity,
		pamhd::mhd::Pressure,
		pamhd::mhd::Magnetic_Field
	> boundaries;
	boundaries.set(document);

	pamhd::mhd::Background_Magnetic_Field<
		pamhd::mhd::Magnetic_Field::data_type
	> background_B;
	background_B.set(document);


	const auto mhd_solver
		= [&options_mhd, &background_B, &rank](){
			if (options_mhd.solver == "rusanov") {

				return pamhd::mhd::get_flux_rusanov<
					pamhd::mhd::MHD_Conservative,
					pamhd::mhd::Magnetic_Field::data_type,
					pamhd::mhd::Mass_Density,
					pamhd::mhd::Momentum_Density,
					pamhd::mhd::Total_Energy_Density,
					pamhd::mhd::Magnetic_Field
				>;

			} else if (options_mhd.solver == "hll-athena") {

				return pamhd::mhd::athena::get_flux_hll<
					pamhd::mhd::MHD_Conservative,
					pamhd::mhd::Magnetic_Field::data_type,
					pamhd::mhd::Mass_Density,
					pamhd::mhd::Momentum_Density,
					pamhd::mhd::Total_Energy_Density,
					pamhd::mhd::Magnetic_Field
				>;

			} else if (options_mhd.solver == "hlld-athena") {

				if (background_B.exists() and rank == 0) {
					std::cout << "NOTE: background magnetic field ignored by hlld-athena solver." << std::endl;
				}

				return pamhd::mhd::athena::get_flux_hlld<
					pamhd::mhd::MHD_Conservative,
					pamhd::mhd::Magnetic_Field::data_type,
					pamhd::mhd::Mass_Density,
					pamhd::mhd::Momentum_Density,
					pamhd::mhd::Total_Energy_Density,
					pamhd::mhd::Magnetic_Field
				>;

			} else if (options_mhd.solver == "roe-athena") {

				if (background_B.exists() and rank == 0) {
					std::cout << "NOTE: background magnetic field ignored by roe-athena solver." << std::endl;
				}

				return pamhd::mhd::athena::get_flux_roe<
					pamhd::mhd::MHD_Conservative,
					pamhd::mhd::Magnetic_Field::data_type,
					pamhd::mhd::Mass_Density,
					pamhd::mhd::Momentum_Density,
					pamhd::mhd::Total_Energy_Density,
					pamhd::mhd::Magnetic_Field
				>;
			} else {
				std::cerr <<  __FILE__ << "(" << __LINE__ << "): "
					<< "Unsupported solver: " << options_mhd.solver
					<< std::endl;
				abort();
			}
		}();


	/*
	Prepare resistivity
	*/

	pamhd::boundaries::Math_Expression<pamhd::mhd::Resistivity> resistivity;
	mup::Value J_val;
	mup::Variable J_var(&J_val);
	resistivity.add_expression_variable("J", J_var);

	const auto resistivity_name = pamhd::mhd::Resistivity::get_option_name();
	if (not document.HasMember(resistivity_name.c_str())) {
		if (rank == 0) {
			std::cerr << __FILE__ "(" << __LINE__
				<< "): Configuration file doesn't have a "
				<< resistivity_name << " key."
				<< std::endl;
		};
		MPI_Finalize();
		return EXIT_FAILURE;
	}
	const auto& json_resistivity = document[resistivity_name.c_str()];
	if (not json_resistivity.IsString()) {
		if (rank == 0) {
			std::cerr << __FILE__ "(" << __LINE__
				<< "): Resistivity option is not of type string."
				<< std::endl;
		};
		MPI_Finalize();
		return EXIT_FAILURE;
	}

	resistivity.set_expression(json_resistivity.GetString());


	/*
	Initialize simulation grid
	*/
	Grid grid;

	pamhd::grid::Options grid_options;
	grid_options.set(document);

	const unsigned int neighborhood_size = 0;
	const auto& number_of_cells = grid_options.get_number_of_cells();
	const auto& periodic = grid_options.get_periodic();
	if (not grid.initialize(
		number_of_cells,
		comm,
		options_mhd.lb_name.c_str(),
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
	const std::array<double, 3>
		simulation_volume
			= grid_options.get_volume(),
		cell_volume{
			simulation_volume[0] / number_of_cells[0],
			simulation_volume[1] / number_of_cells[1],
			simulation_volume[2] / number_of_cells[2]
		};

	dccrg::Cartesian_Geometry::Parameters geom_params;
	geom_params.start = grid_options.get_start();
	geom_params.level_0_cell_length = cell_volume;

	if (not grid.set_geometry(geom_params)) {
		std::cerr << __FILE__ << ":" << __LINE__
			<< ": Couldn't set grid geometry."
			<< std::endl;
		abort();
	}

	grid.balance_load();

	// update owner process of cells for saving into file
	for (auto& cell: grid.cells) {
		(*cell.data)[pamhd::mhd::MPI_Rank()] = rank;
	}

	// assign cells into boundary geometries
	for (const auto& cell: grid.cells) {
		const auto
			start = grid.geometry.get_min(cell.id),
			end = grid.geometry.get_max(cell.id);
		geometries.overlaps(start, end, cell.id);
	}

	// pointer to data of every local cell and its neighbor(s)
	const auto& cell_data_pointers = grid.get_cell_data_pointers();

	// index of first outer cell in dccrg's cell data pointer cache
	size_t outer_cell_start_i = 0;
	for (const auto& item: cell_data_pointers) {
		outer_cell_start_i++;
		if (get<0>(item) == dccrg::error_cell) {
			break;
		}
	}


	/*
	Simulate
	*/

	const double time_end = options_mhd.time_start + options_mhd.time_length;
	double
		max_dt = 0,
		simulation_time = options_mhd.time_start,
		next_mhd_save = options_mhd.save_mhd_n,
		next_rem_div_B = options_mhd.remove_div_B_n;

	std::vector<uint64_t>
		cells = grid.get_cells(),
		inner_cells = grid.get_local_cells_not_on_process_boundary(),
		outer_cells = grid.get_local_cells_on_process_boundary();

	// initialize MHD
	const bool verbose = true;
	if (verbose and rank == 0) {
		cout << "Initializing MHD... " << endl;
	}
	pamhd::mhd::initialize(
		geometries,
		initial_conditions,
		background_B,
		grid,
		cells,
		simulation_time,
		options_mhd.adiabatic_index,
		options_mhd.vacuum_permeability,
		options_mhd.proton_mass,
		verbose,
		Mas, Mom, Nrj, Mag,
		Bg_B_Pos_X, Bg_B_Pos_Y, Bg_B_Pos_Z,
		Mas_f, Mom_f, Nrj_f, Mag_f
	);

	// update background field between processes
	Cell::set_transfer_all(
		true,
		pamhd::mhd::Bg_Magnetic_Field_Pos_X(),
		pamhd::mhd::Bg_Magnetic_Field_Pos_Y(),
		pamhd::mhd::Bg_Magnetic_Field_Pos_Z()
	);
	grid.update_copies_of_remote_neighbors();
	Cell::set_transfer_all(
		false,
		pamhd::mhd::Bg_Magnetic_Field_Pos_X(),
		pamhd::mhd::Bg_Magnetic_Field_Pos_Y(),
		pamhd::mhd::Bg_Magnetic_Field_Pos_Z()
	);

	// initialize resistivity
	for (auto& cell: grid.cells) {
		Res(*cell.data) = 0;
	}

	pamhd::mhd::apply_boundaries(
		grid,
		boundaries,
		geometries,
		simulation_time,
		Mas,
		Mom,
		Nrj,
		Mag,
		options_mhd.proton_mass,
		options_mhd.adiabatic_index,
		options_mhd.vacuum_permeability
	);
	if (verbose and rank == 0) {
		cout << "Done initializing MHD" << endl;
	}

	/*
	Classify cells into normal, boundary and dont_solve
	*/

	pamhd::mhd::set_solver_info<pamhd::mhd::Solver_Info>(
		grid, boundaries, geometries, Sol_Info
	);
	// make lists from above for divergence removal functions
	std::vector<uint64_t> solve_cells, bdy_cells, skip_cells;
	for (const auto& cell: grid.cells) {
		if ((Sol_Info(*cell.data) & pamhd::mhd::Solver_Info::dont_solve) > 0) {
			skip_cells.push_back(cell.id);
		} else if (Sol_Info(*cell.data) > 0) {
			bdy_cells.push_back(cell.id);
		} else {
			solve_cells.push_back(cell.id);
		}
	}

	size_t simulated_steps = 0;
	while (simulation_time < time_end) {
		simulated_steps++;

		/*
		Get maximum allowed time step
		*/
		double
			// don't step over the final simulation time
			until_end = time_end - simulation_time,
			local_time_step = min(options_mhd.time_step_factor * max_dt, until_end),
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

		Cell::set_transfer_all(true, pamhd::mhd::MHD_State_Conservative());
		grid.start_remote_neighbor_copy_updates();

		pamhd::divergence::get_curl(
			// TODO: consider boundaries
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
			Cur(*cell_data) /= options_mhd.vacuum_permeability;
		}

		double solve_max_dt = -1;
		size_t solve_index = 0;
		std::tie(
			solve_max_dt,
			solve_index
		) = pamhd::mhd::solve<pamhd::mhd::Solver_Info>(
			mhd_solver,
			0,
			grid,
			time_step,
			options_mhd.adiabatic_index,
			options_mhd.vacuum_permeability,
			Mas, Mom, Nrj, Mag,
			Bg_B_Pos_X, Bg_B_Pos_Y, Bg_B_Pos_Z,
			Mas_f, Mom_f, Nrj_f, Mag_f,
			Sol_Info
		);
		max_dt = min(
			max_dt,
			solve_max_dt
		);

		grid.wait_remote_neighbor_copy_update_receives();

		std::tie(
			solve_max_dt,
			solve_index
		) = pamhd::mhd::solve<pamhd::mhd::Solver_Info>(
			mhd_solver,
			solve_index + 1,
			grid,
			time_step,
			options_mhd.adiabatic_index,
			options_mhd.vacuum_permeability,
			Mas, Mom, Nrj, Mag,
			Bg_B_Pos_X, Bg_B_Pos_Y, Bg_B_Pos_Z,
			Mas_f, Mom_f, Nrj_f, Mag_f,
			Sol_Info
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
			Cur(*cell_data) /= options_mhd.vacuum_permeability;
		}

		grid.wait_remote_neighbor_copy_update_sends();
		Cell::set_transfer_all(false, pamhd::mhd::MHD_State_Conservative());


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

			const auto c = grid.geometry.get_center(cell);
			const auto r = sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2]);

			J_val = Cur(*cell_data).norm();
			Res(*cell_data) = resistivity.evaluate(
				simulation_time,
				c[0], c[1], c[2],
				r, asin(c[2] / r), atan2(c[1], c[0])
			);

			//TODO keep pressure/temperature constant despite electric resistivity
			Mag_res(*cell_data) *= -Res(*cell_data);
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

			const auto c = grid.geometry.get_center(cell);
			const auto r = sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2]);

			J_val = Cur(*cell_data).norm();
			Res(*cell_data) = resistivity.evaluate(
				simulation_time,
				c[0], c[1], c[2],
				r, asin(c[2] / r), atan2(c[1], c[0])
			);

			Mag_res(*cell_data) *= -Res(*cell_data);
			Mag_f(*cell_data) += Mag_res(*cell_data);
		}

		grid.wait_remote_neighbor_copy_update_sends();
		Cell::set_transfer_all(false, pamhd::mhd::Electric_Current_Density());


		pamhd::mhd::apply_fluxes<pamhd::mhd::Solver_Info>(
			grid,
			options_mhd.min_pressure,
			options_mhd.adiabatic_index,
			options_mhd.vacuum_permeability,
			Mas, Mom, Nrj, Mag,
			Mas_f, Mom_f, Nrj_f, Mag_f,
			Sol_Info
		);

		simulation_time += time_step;


		/*
		Remove divergence of magnetic field
		*/

		if (options_mhd.remove_div_B_n > 0 and simulation_time >= next_rem_div_B) {
			next_rem_div_B += options_mhd.remove_div_B_n;

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
					solve_cells,
					bdy_cells,
					skip_cells,
					grid,
					Mag,
					Mag_div,
					[](Cell& cell_data)
						-> pamhd::mhd::Scalar_Potential_Gradient::data_type&
					{
						return cell_data[pamhd::mhd::Scalar_Potential_Gradient()];
					},
					options_mhd.poisson_iterations_max,
					options_mhd.poisson_iterations_min,
					options_mhd.poisson_norm_stop,
					2,
					options_mhd.poisson_norm_increase_max,
					0,
					false
				);
			Cell::set_transfer_all(false, pamhd::mhd::Magnetic_Field_Divergence());

			grid.update_copies_of_remote_neighbors();
			Cell::set_transfer_all(false, pamhd::mhd::MHD_State_Conservative());
			const double div_after
				= pamhd::divergence::get_divergence(
					solve_cells,
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
						) / (2 * options_mhd.vacuum_permeability);

					Nrj(*cell_data) += mag_nrj_diff;
				}
			}
		}


		pamhd::mhd::apply_boundaries(
			grid,
			boundaries,
			geometries,
			simulation_time,
			Mas,
			Mom,
			Nrj,
			Mag,
			options_mhd.proton_mass,
			options_mhd.adiabatic_index,
			options_mhd.vacuum_permeability
		);


		/*
		Save simulation to disk
		*/

		if (
			(
				options_mhd.save_mhd_n >= 0
				and (
					simulation_time == options_mhd.time_start
					or simulation_time >= time_end
				)
			) or (options_mhd.save_mhd_n > 0 and simulation_time >= next_mhd_save)
		) {
			if (next_mhd_save <= simulation_time) {
				next_mhd_save
					+= options_mhd.save_mhd_n
					* ceil((simulation_time - next_mhd_save) / options_mhd.save_mhd_n);
			}

			if (verbose and rank == 0) {
				cout << "Saving MHD at time " << simulation_time << endl;
			}

			if (
				not pamhd::mhd::save(
					boost::filesystem::canonical(
						boost::filesystem::path(options_mhd.output_directory)
					).append("mhd_").generic_string(),
					grid,
					2,
					simulation_time,
					options_mhd.adiabatic_index,
					options_mhd.proton_mass,
					options_mhd.vacuum_permeability,
					pamhd::mhd::MHD_State_Conservative(),
					pamhd::mhd::Electric_Current_Density(),
					pamhd::mhd::Solver_Info(),
					pamhd::mhd::MPI_Rank(),
					pamhd::mhd::Resistivity(),
					pamhd::mhd::Bg_Magnetic_Field_Pos_X(),
					pamhd::mhd::Bg_Magnetic_Field_Pos_Y(),
					pamhd::mhd::Bg_Magnetic_Field_Pos_Z()
				)
			) {
				std::cerr <<  __FILE__ << "(" << __LINE__ << "): "
					"Couldn't save mhd result."
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
