/*
Reference test program for MHD solvers of PAMHD.

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
#include "sstream"
#include "string"
#include "type_traits"

#include "gensimcell.hpp"
#include "mhd/common.hpp"
#include "mhd/hll_athena.hpp"
#include "mhd/variables.hpp"


using namespace std;


constexpr double
	adiabatic_index = 5.0 / 3.0,    
	vacuum_permeability = 4e-7 * M_PI,
	proton_mass = 1.672621777e-27;


// 1 boundary cell in each direction
constexpr size_t grid_length = 1000 + 2;

using Cell_T = gensimcell::Cell<
	pamhd::mhd::MHD_State_Conservative,
	pamhd::mhd::MHD_Flux_Conservative
>;
using Grid_T = std::array<Cell_T, grid_length>;


/*!
Returns start coordinate of the grid.
*/
constexpr double get_grid_start()
{
	return 0.0;
}

/*!
Returns end coordinate of the grid.
*/
constexpr double get_grid_end()
{
	return 1e5;
}


/*!
Returns the size of a grid cell.
*/
constexpr double get_cell_size()
{
	static_assert(
		std::tuple_size<Grid_T>::value > 0,
		"Grid must have at least one cell"
	);
	return (get_grid_end() - get_grid_start()) / std::tuple_size<Grid_T>::value;
}


/*!
Returns the center of a cell located at given index in given grid.

Index starts from 0.
Returns a quiet NaN in case of error.
*/
double get_cell_center(const size_t index)
{
	static_assert(
		std::tuple_size<Grid_T>::value > 0,
		"Grid must have at least one cell"
	);

	if (index >= std::tuple_size<Grid_T>::value) {
		return std::numeric_limits<double>::quiet_NaN();
	}

	return (0.5 + index) * get_cell_size();
}


//! Sets the initial state of simulation, zeroes fluxes
template <
	class MHD_T,
	class MHD_Flux_T,
	class Mass_Density_T,
	class Momentum_Density_T,
	class Total_Energy_Density_T,
	class Magnetic_Field_T
> void initialize_mhd(
	Grid_T& grid,
	const double adiabatic_index,
	const double vacuum_permeability
) {
	static_assert(
		std::tuple_size<Grid_T>::value > 0,
		"Grid must have at least one cell"
	);

	static_assert(
		std::is_same<typename MHD_T::data_type, typename MHD_Flux_T::data_type>::value,
		"The data types of variables MHD_T and MHD_Flux_T must be equal"
	);

	// shorthand for referring to variables
	const Mass_Density_T Rho{};
	const Momentum_Density_T M{};
	const Total_Energy_Density_T E{};
	const Magnetic_Field_T B{};
	const pamhd::mhd::Velocity V{};
	const pamhd::mhd::Pressure P{};

	for (size_t cell_i = 0; cell_i < grid.size(); cell_i++) {

		auto
			&state = grid[cell_i][MHD_T()],
			&flux = grid[cell_i][MHD_Flux_T()];

		flux[Rho]  =
		flux[E]    =
		flux[M][0] =
		flux[M][1] =
		flux[M][2] =
		flux[B][0] =
		flux[B][1] =
		flux[B][2] = 0;

		// specify initial state in primitive variables
		pamhd::mhd::MHD_Primitive temp;

		state[M][0] =
		state[M][1] =
		state[M][2] =
		temp[V][0]  =
		temp[V][1]  =
		temp[V][2]  = 0;

		state[B][0] =
		temp[B][0]  = 1.5e-9;

		state[B][2] =
		temp[B][2]  = 0;

		const double center = get_cell_center(cell_i);
		if (center < (get_grid_end() - get_grid_start()) / 2) {
			state[Rho] =
			temp[Rho]  = proton_mass * 3e6;

			state[B][1] =
			temp[B][1]  = 1e-9;

			temp[P] = 3e-12;

		} else {

			state[Rho] =
			temp[Rho]  = proton_mass * 1e6;

			state[B][1] =
			temp[B][1]  = -1e-9;

			temp[P] = 1e-12;
		}

		state[E]
			= pamhd::mhd::get_total_energy_density<
				pamhd::mhd::MHD_Primitive,
				Mass_Density_T,
				pamhd::mhd::Velocity,
				pamhd::mhd::Pressure,
				Magnetic_Field_T
			>(
				temp,
				adiabatic_index,
				vacuum_permeability
			);
	}
}


/*!
Advances the MHD solution for one time step of length dt with given solver.

Returns the maximum allowed length of time step for the
next step.
*/
template <
	class MHD_T,
	class MHD_Flux_T,
	class Mass_Density_T,
	class Momentum_Density_T,
	class Total_Energy_Density_T,
	class Magnetic_Field_T
> double solve_mhd(
	const std::string& solver,
	Grid_T& grid,
	const double dt,
	const double adiabatic_index,
	const double vacuum_permeability
) {
	static_assert(
		std::tuple_size<Grid_T>::value >= 3,
		"Grid must have at least three cells"
	);

	static_assert(
		std::is_same<typename MHD_T::data_type, typename MHD_Flux_T::data_type>::value,
		"The data types of variables MHD_T and MHD_Flux_T must be equal"
	);

	const Mass_Density_T Rho{};
	const Momentum_Density_T M{};
	const Total_Energy_Density_T E{};
	const Magnetic_Field_T B{};

	const double
		cell_size = get_cell_size(),
		face_area = 1.0;

	double max_dt = std::numeric_limits<double>::max();

	// calculate fluxes
	for (size_t cell_i = 0; cell_i < std::tuple_size<Grid_T>::value - 1; cell_i++) {

		Cell_T
			&cell = grid[cell_i],
			&neighbor = grid[cell_i + 1];

		double max_vel = -1;
		typename MHD_Flux_T::data_type flux;

		if (solver == "hll_athena") {

			std::tie(
				flux,
				max_vel
			) = pamhd::mhd::athena::get_flux_hll<
				typename MHD_T::data_type,
				Mass_Density_T,
				Momentum_Density_T,
				Total_Energy_Density_T,
				Magnetic_Field_T
			>(
				cell[MHD_T()],
				neighbor[MHD_T()],
				face_area,
				dt,
				adiabatic_index,
				vacuum_permeability
			);

		} else {

			std::cerr <<  __FILE__ << "(" << __LINE__ << ") Invalid solver: "
				<< solver << ", use --help to list available solvers"
				<< std::endl;
			abort();
		}

		max_dt = std::min(max_dt, cell_size / max_vel);

		if (cell_i > 0) {
			// positive flux flows neg->pos, i.e. out of current cell
			cell[MHD_Flux_T()][Rho] -= flux[Rho];
			cell[MHD_Flux_T()][M] -= flux[M];
			cell[MHD_Flux_T()][E] -= flux[E];
			cell[MHD_Flux_T()][B] -= flux[B];
		}
		if (cell_i < std::tuple_size<Grid_T>::value - 2) {
			neighbor[MHD_Flux_T()][Rho] += flux[Rho];
			neighbor[MHD_Flux_T()][M] += flux[M];
			neighbor[MHD_Flux_T()][E] += flux[E];
			neighbor[MHD_Flux_T()][B] += flux[B];
		}
	}

	// apply fluxes
	const double inverse_volume = 1.0 / (face_area * cell_size);

	for (auto& cell: grid) {
		pamhd::mhd::apply_flux<
			Grid_T::value_type,
			MHD_T,
			MHD_Flux_T,
			Mass_Density_T,
			Momentum_Density_T,
			Total_Energy_Density_T,
			Magnetic_Field_T
		>(
			cell,
			inverse_volume
		);
	}

	return max_dt;
}

/*!
Writes given advection simulation to a file plottable with gnuplot.

Returns the name of the file written which
is derived from given simulation time.
*/
template <
	class MHD_T,
	class Mass_Density_T,
	class Momentum_Density_T,
	class Total_Energy_Density_T,
	class Magnetic_Field_T
> std::string plot_mhd(
	const Grid_T& grid,
	const double simulation_time,
	const double adiabatic_index,
	const double vacuum_permeability
) {
	std::ostringstream time_string, normalization_string;
	time_string
		<< std::setw(4)
		<< std::setfill('0')
		<< size_t(simulation_time * 1000);

	const std::string
		gnuplot_file_name(
			"mhd_"
			+ time_string.str()
			+ "_ms.dat"
		),
		plot_file_name(
			"mhd_"
			+ time_string.str()
			+ "_ms.png"
		),
		B_plot_file_name(
			"mhd_B_"
			+ time_string.str()
			+ "_ms.png"
		),
		v_plot_file_name(
			"mhd_v_"
			+ time_string.str()
			+ "_ms.png"
		);

	std::ofstream gnuplot_file(gnuplot_file_name);

	// scale outputs to -1 < value < 1
	std::array<double, 4> ranges;
	for (auto& range: ranges) {
		range = 0;
	}
	for (size_t cell_i = 0; cell_i < std::tuple_size<Grid_T>::value; cell_i++) {
		const auto& cell_data = grid[cell_i];
		const auto& state = cell_data[MHD_T()];
		const auto& mass = state[Mass_Density_T()];
		const auto& mom = state[Momentum_Density_T()];
		const auto& mag = state[Magnetic_Field_T()];
		const auto vel(mom / mass);
		const double p
			= pamhd::mhd::get_pressure<
				typename MHD_T::data_type,
				Mass_Density_T,
				Momentum_Density_T,
				Total_Energy_Density_T,
				Magnetic_Field_T
			>(state, adiabatic_index, vacuum_permeability);

		if (ranges[0] < mass) ranges[0] = mass;
		if (ranges[1] < std::fabs(vel[0])) ranges[1] = std::fabs(vel[0]);
		if (ranges[1] < std::fabs(vel[1])) ranges[1] = std::fabs(vel[1]);
		if (ranges[1] < std::fabs(vel[2])) ranges[1] = std::fabs(vel[2]);
		if (ranges[2] < std::fabs(mag[0])) ranges[2] = std::fabs(mag[0]);
		if (ranges[2] < std::fabs(mag[1])) ranges[2] = std::fabs(mag[1]);
		if (ranges[2] < std::fabs(mag[2])) ranges[2] = std::fabs(mag[2]);
		if (ranges[3] < p) ranges[3] = p;
	}
	// don't divide by 0
	for (auto& range: ranges) {
		if (range == 0) {
			range = 1;
		}
	}
	gnuplot_file
		<< "set term png enhanced\nset output '"
		<< plot_file_name
		<< "'\nset xlabel 'X / X_0 (X_0 = "
		<< boost::lexical_cast<std::string>(get_grid_end() - get_grid_start())
		<< " m)'\nset ylabel 'F / F_{max}'\n"
		   "set key horizontal outside bottom\n"
		   "set yrange [0 : 1.05]\n"
		   "plot "
		     "'-' using 1:2 with line linewidth 2 title 'density', "
		     "'-' u 1:2 w l lw 2 t 'pressure'\n"
		     ;

	// mass density
	for (size_t cell_i = 0; cell_i < std::tuple_size<Grid_T>::value; cell_i++) {
		const auto& cell = grid[cell_i];
		const double x
			= get_cell_center(cell_i)
			/ (get_grid_end() - get_grid_start());
		gnuplot_file << x << " "
			<< cell[MHD_T()][Mass_Density_T()] / ranges[0] << "\n";
	}
	gnuplot_file << "end\n";
	// pressure
	for (size_t cell_i = 0; cell_i < std::tuple_size<Grid_T>::value; cell_i++) {
		const auto& cell_data = grid[cell_i];
		const auto& state = cell_data[MHD_T()];
		const double x = get_cell_center(cell_i) / (get_grid_end() - get_grid_start());
		gnuplot_file
			<< x << " "
			<< pamhd::mhd::get_pressure<
				typename MHD_T::data_type,
				Mass_Density_T,
				Momentum_Density_T,
				Total_Energy_Density_T,
				Magnetic_Field_T
			>(state, adiabatic_index, vacuum_permeability)
			/ ranges[3]
			<< "\n";
	}
	gnuplot_file << "end\n";

	normalization_string << std::scientific << std::setprecision(1);
	normalization_string << ranges[1];
	gnuplot_file
		<< "set term png enhanced\nset output '"
		<< v_plot_file_name
		<< "'\nset xlabel 'X / X_{max} (X_{max} = "
		<< boost::lexical_cast<std::string>(get_grid_end() - get_grid_start())
		<< " m)'\nset ylabel 'V / V_{max} (V_{max} = "
		<< normalization_string.str()
		<< " m/s)'\nset key horizontal outside bottom\n"
		   "set yrange [-1.05 : 1.05]\n"
		   "plot "
		     "'-' u 1:2 w l lw 2 t 'v_x', "
		     "'-' u 1:2 w l lw 2 t 'v_y', "
		     "'-' u 1:2 w l lw 2 t 'v_z'\n"
		;
	// vx
	for (size_t cell_i = 0; cell_i < std::tuple_size<Grid_T>::value; cell_i++) {
		const auto& cell_data = grid[cell_i];
		const auto& state = cell_data[MHD_T()];
		const auto& m = state[Momentum_Density_T()];
		const auto& rho = state[Mass_Density_T()];
		const auto v(m / rho);
		const double x = get_cell_center(cell_i) / (get_grid_end() - get_grid_start());
		gnuplot_file << x << " "
			<< v[0] / ranges[1] << "\n";
	}
	gnuplot_file << "end\n";
	// vy
	for (size_t cell_i = 0; cell_i < std::tuple_size<Grid_T>::value; cell_i++) {
		const auto& cell_data = grid[cell_i];
		const auto& state = cell_data[MHD_T()];
		const auto& m = state[Momentum_Density_T()];
		const auto& rho = state[Mass_Density_T()];
		const auto v(m / rho);
		const double x = get_cell_center(cell_i) / (get_grid_end() - get_grid_start());
		gnuplot_file << x << " "
			<< v[1] / ranges[1] << "\n";
	}
	gnuplot_file << "end\n";
	// vz
	for (size_t cell_i = 0; cell_i < std::tuple_size<Grid_T>::value; cell_i++) {
		const auto& cell_data = grid[cell_i];
		const auto& state = cell_data[MHD_T()];
		const auto& m = state[Momentum_Density_T()];
		const auto& rho = state[Mass_Density_T()];
		const auto v(m / rho);
		const double x = get_cell_center(cell_i) / (get_grid_end() - get_grid_start());
		gnuplot_file << x << " "
			<< v[2] / ranges[1] << "\n";
	}
	gnuplot_file << "end\n";

	normalization_string.str("");
	normalization_string << ranges[2];
	gnuplot_file
		<< "set term png enhanced\nset output '"
		<< B_plot_file_name
		<< "'\nset xlabel 'X / X_0 (X_0 = "
		<< boost::lexical_cast<std::string>(get_grid_end() - get_grid_start())
		<< " m)'\nset ylabel 'B / B_{max} (B_{max} = "
		<< normalization_string.str()
		<< " T)'\nset key horizontal outside bottom\n"
		   "set key horizontal outside bottom\n"
		   "set yrange [-1.05 : 1.05]\n"
		   "plot "
		     "'-' u 1:2 w l lw 2 t 'B_x', "
		     "'-' u 1:2 w l lw 2 t 'B_y', "
		     "'-' u 1:2 w l lw 2 t 'B_z'\n"
		     ;
	// Bx
	for (size_t cell_i = 0; cell_i < std::tuple_size<Grid_T>::value; cell_i++) {
		const double x
			= get_cell_center(cell_i)
			/ (get_grid_end() - get_grid_start());
		gnuplot_file << x << " "
			<< grid[cell_i][MHD_T()][Magnetic_Field_T()][0] / ranges[2] << "\n";
	}
	gnuplot_file << "end\n";
	// By
	for (size_t cell_i = 0; cell_i < std::tuple_size<Grid_T>::value; cell_i++) {
		const double x
			= get_cell_center(cell_i)
			/ (get_grid_end() - get_grid_start());
		gnuplot_file << x << " "
			<< grid[cell_i][MHD_T()][Magnetic_Field_T()][1] / ranges[2] << "\n";
	}
	gnuplot_file << "end\n";
	// Bz
	for (size_t cell_i = 0; cell_i < std::tuple_size<Grid_T>::value; cell_i++) {
		const double x
			= get_cell_center(cell_i)
			/ (get_grid_end() - get_grid_start());
		gnuplot_file << x << " "
			<< grid[cell_i][MHD_T()][Magnetic_Field_T()][2] / ranges[2] << "\n";
	}
	gnuplot_file << "end\n";

	return gnuplot_file_name;
}

/*!
Writes given advection simulation to an ascii file.

Returns the name of the file written which
is derived from given solver name.
*/
template <
	class MHD_T,
	class Mass_Density_T,
	class Momentum_Density_T,
	class Total_Energy_Density_T,
	class Magnetic_Field_T
> void save_mhd(
	const std::string& solver,
	const Grid_T& grid
) {
	const std::string file_name("mhd_" + solver + ".dat");

	std::ofstream outfile(file_name);
	outfile << std::setprecision(16) << std::scientific;

	for (const auto& cell_data: grid) {
		const auto& state = cell_data[MHD_T()];
		outfile
			<< state[Mass_Density_T()] << " "
			<< state[Momentum_Density_T()][0] << " "
			<< state[Momentum_Density_T()][1] << " "
			<< state[Momentum_Density_T()][2] << " "
			<< state[Total_Energy_Density_T()] << " "
			<< state[Magnetic_Field_T()][0] << " "
			<< state[Magnetic_Field_T()][1] << " "
			<< state[Magnetic_Field_T()][2] << "\n";
	}
	outfile << std::endl;
}


template <class T> T get_relative_error(const T a, const T b)
{
	if (a == T(0) && b == T(0)) {
		return {0};
	}

	return {std::fabs(a - b) / std::max(std::fabs(a), std::fabs(b))};
}


template <
	class MHD_T,
	class Mass_Density_T,
	class Momentum_Density_T,
	class Total_Energy_Density_T,
	class Magnetic_Field_T
> void verify_mhd(
	const Grid_T& grid,
	const std::string& file_name
) {
	std::ifstream infile(file_name);
	if (not infile.good()) {
		// try a subdirectory in case run from repository root
		infile.open("tests/mhd/" + file_name);

		if (not infile.good()) {
			std::cerr <<  __FILE__ << ":" << __LINE__
				<< " Couldn't open file " << file_name
				<< std::endl;
			abort();
		}
	}

	constexpr double max_error = 1e-4;

	for (size_t cell_i = 0; cell_i < std::tuple_size<Grid_T>::value; cell_i++) {
		typename Mass_Density_T::data_type ref_rho;
		typename Momentum_Density_T::data_type ref_mom;
		typename Total_Energy_Density_T::data_type ref_nrj;
		typename Magnetic_Field_T::data_type ref_mag;

		infile
			>> ref_rho
			>> ref_mom[0]
			>> ref_mom[1]
			>> ref_mom[2]
			>> ref_nrj
			>> ref_mag[0]
			>> ref_mag[1]
			>> ref_mag[2];

		const auto& state = grid[cell_i][MHD_T()];
		const auto rho = state[Mass_Density_T()];
		const auto mom = state[Momentum_Density_T()];
		const auto nrj = state[Total_Energy_Density_T()];
		const auto mag = state[Magnetic_Field_T()];

		if (get_relative_error(rho, ref_rho) > max_error) {
			std::cerr <<  __FILE__ << ":" << __LINE__
				<< " density " << rho << " != " << ref_rho
				<< " at cell " << cell_i
				<< std::endl;
			abort();
		}

		if (get_relative_error(mom[0], ref_mom[0]) > max_error) {
			std::cerr <<  __FILE__ << ":" << __LINE__
				<< " x momentum " << mom[0] << " != " << ref_mom[0]
				<< " at cell " << cell_i
				<< std::endl;
			abort();
		}
		if (get_relative_error(mom[1], ref_mom[1]) > max_error) {
			std::cerr <<  __FILE__ << ":" << __LINE__
				<< " y momentum " << mom[1] << " != " << ref_mom[1]
				<< " at cell " << cell_i
				<< std::endl;
			abort();
		}
		if (get_relative_error(mom[2], ref_mom[2]) > max_error) {
			std::cerr <<  __FILE__ << ":" << __LINE__
				<< " z momentum " << mom[2] << " != " << ref_mom[2]
				<< " at cell " << cell_i
				<< std::endl;
			abort();
		}

		if (get_relative_error(nrj, ref_nrj) > max_error) {
			std::cerr <<  __FILE__ << ":" << __LINE__
				<< " energy " << nrj << " != " << ref_nrj
				<< " at cell " << cell_i
				<< std::endl;
			abort();
		}

		if (get_relative_error(mag[0], ref_mag[0]) > max_error) {
			std::cerr <<  __FILE__ << ":" << __LINE__
				<< " x magnetic field " << mag[0] << " != " << ref_mag[0]
				<< " at cell " << cell_i
				<< std::endl;
			abort();
		}
		if (get_relative_error(mag[1], ref_mag[1]) > max_error) {
			std::cerr <<  __FILE__ << ":" << __LINE__
				<< " y magnetic field " << mag[1] << " != " << ref_mag[1]
				<< " at cell " << cell_i
				<< std::endl;
			abort();
		}
		if (get_relative_error(mag[2], ref_mag[2]) > max_error) {
			std::cerr <<  __FILE__ << ":" << __LINE__
				<< " z magnetic field " << mag[2] << " != " << ref_mag[2]
				<< " at cell " << cell_i
				<< std::endl;
			abort();
		}
	}
}


int main(int argc, char* argv[])
{
	bool save = false, plot = false, no_verify = false, verbose = false;
	std::string solver;
	boost::program_options::options_description options(
		"Usage: program_name [options], where options are:"
	);
	options.add_options()
		("help", "Print this help message")
		("solver",
			boost::program_options::value<std::string>(&solver)->default_value("hll_athena"),
			"Solver to use, available: hll_athena")
		("save", "Save end result to ascii file")
		("plot", "Plot results using gnuplot")
		("no-verify", "Do not verify against reference results")
		("verbose", "Print run time information");

	// read options from command line
	boost::program_options::variables_map option_variables;
	boost::program_options::store(boost::program_options::parse_command_line(argc, argv, options), option_variables);
	boost::program_options::notify(option_variables);

	if (option_variables.count("help") > 0) {
		cout << options << endl;
		return EXIT_SUCCESS;
	}

	if (option_variables.count("save") > 0) {
		save = true;
	}

	if (option_variables.count("plot") > 0) {
		plot = true;
	}

	if (option_variables.count("no-verify") > 0) {
		no_verify = true;
	}

	if (option_variables.count("verbose") > 0) {
		verbose = true;
	}


	Grid_T grid;

	if (verbose) {
		cout << "Initializing MHD" << endl;
	}
	initialize_mhd<
		pamhd::mhd::MHD_State_Conservative,
		pamhd::mhd::MHD_Flux_Conservative,
		pamhd::mhd::Mass_Density,
		pamhd::mhd::Momentum_Density,
		pamhd::mhd::Total_Energy_Density,
		pamhd::mhd::Magnetic_Field
	>(grid, adiabatic_index, vacuum_permeability);

	if (plot) {
		const std::string mhd_gnuplot_file_name
			= plot_mhd<
				pamhd::mhd::MHD_State_Conservative,
				pamhd::mhd::Mass_Density,
				pamhd::mhd::Momentum_Density,
				pamhd::mhd::Total_Energy_Density,
				pamhd::mhd::Magnetic_Field
			>(grid, 0, adiabatic_index, vacuum_permeability);
		system(("gnuplot " + mhd_gnuplot_file_name).c_str());
		if (verbose) {
			cout << "Initial state plotted from file " << mhd_gnuplot_file_name << endl;
		}
	}

	const double
		simulation_duration = 1,
		mhd_plot_interval = simulation_duration / 5;
	double
		max_dt = 0,
		simulation_time = 0,
		mhd_next_plot = mhd_plot_interval;
	while (simulation_time < simulation_duration) {

		const double
			CFL = 0.5,
			// don't step over the final simulation time
			until_end = simulation_duration - simulation_time,
			time_step = std::min(CFL * max_dt, until_end);

		max_dt = std::numeric_limits<double>::max();

		if (verbose) {
			cout << "Solving MHD at time " << simulation_time
				<< " s with time step " << time_step << " s" << endl;
		}
		max_dt = std::min(
			max_dt,
			solve_mhd<
				pamhd::mhd::MHD_State_Conservative,
				pamhd::mhd::MHD_Flux_Conservative,
				pamhd::mhd::Mass_Density,
				pamhd::mhd::Momentum_Density,
				pamhd::mhd::Total_Energy_Density,
				pamhd::mhd::Magnetic_Field
			>(
				"hll_athena",
				grid,
				time_step,
				adiabatic_index,
				vacuum_permeability
			)
		);

		simulation_time += time_step;

		if (plot && mhd_next_plot <= simulation_time) {
			mhd_next_plot += mhd_plot_interval;

			if (verbose) {
				cout << "Plotting simulation at time " << simulation_time << endl;
			}
			const std::string mhd_gnuplot_file_name
				= plot_mhd<
					pamhd::mhd::MHD_State_Conservative,
					pamhd::mhd::Mass_Density,
					pamhd::mhd::Momentum_Density,
					pamhd::mhd::Total_Energy_Density,
					pamhd::mhd::Magnetic_Field
				>(grid, simulation_time, adiabatic_index, vacuum_permeability);
			system(("gnuplot " + mhd_gnuplot_file_name).c_str());
		}
	}

	if (save) {
		if (verbose) {
			cout << "Saving MHD at time " << simulation_time << endl;
		}
		save_mhd<
			pamhd::mhd::MHD_State_Conservative,
			pamhd::mhd::Mass_Density,
			pamhd::mhd::Momentum_Density,
			pamhd::mhd::Total_Energy_Density,
			pamhd::mhd::Magnetic_Field
		>("hll_athena", grid);
	}

	if (not no_verify) {
		const std::string reference_name("mhd_" + solver + ".ref");

		if (verbose) {
			cout << "Verifying result against file " << reference_name << endl;
		}

		verify_mhd<
			pamhd::mhd::MHD_State_Conservative,
			pamhd::mhd::Mass_Density,
			pamhd::mhd::Momentum_Density,
			pamhd::mhd::Total_Energy_Density,
			pamhd::mhd::Magnetic_Field
		>(grid, reference_name);
	}

	if (verbose) {
		cout << "Simulation finished at time " << simulation_time << endl;
	}

	return EXIT_SUCCESS;
}
