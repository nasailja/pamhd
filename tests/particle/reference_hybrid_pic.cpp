/*
Reference test program for particle solver of PAMHD.

Copyright 2014, 2015 Ilja Honkonen
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
#include "chrono"
#include "cmath"
#include "cstdlib"
#include "fstream"
#include "iomanip"
#include "iostream"
#include "limits"
#include "random"
#include "sstream"
#include "string"
#include "type_traits"

#include "boost/lexical_cast.hpp"
#include "boost/program_options.hpp"
#include "Eigen/Geometry" 
#include "gensimcell.hpp"
#include "prettyprint.hpp"

#include "particle/accumulate.hpp"
#include "particle/common.hpp"
#include "particle/solve.hpp"
#include "particle/variables.hpp"
#include "poisson/solver.hpp"


using namespace std;


//! Initial positive charge required for quasi-neutrality
struct Ion_Charge
{
	using data_type = double;
};

//! Weighted average of electron charge density of nearby particles
struct Electron_Charge
{
	using data_type = double;
};

//! Total electric charge in a cell
struct Total_Charge
{
	using data_type = double;
};

//! Cell-centered electric field for propagating particles
struct Electric_Field
{
	using data_type = double;
};

//! Cell-centered electric potential obtained from Poisson's equation
struct Electric_Potential
{
	using data_type = double;
};

/*!
Serial program doesn't need target process variable in external particles
*/
struct Particles_Ext {
	using data_type = pamhd::particle::Particles_Storage;
};

using Cell = gensimcell::Cell<
	gensimcell::Never_Transfer,
	pamhd::particle::Particles_Int,
	Particles_Ext,
	Ion_Charge,
	Electron_Charge,
	Total_Charge,
	Electric_Field,
	Electric_Potential
>;

constexpr size_t grid_size = 30;
using Grid = std::array<Cell, grid_size>;


constexpr double
	//adiabatic_index = 5.0 / 3.0,    
	//vacuum_permeability = 4e-7 * M_PI,
	//proton_mass = 1.672621777e-27,
	//proton_charge_mass_ratio = elementary_charge / proton_mass,
	vacuum_permittivity = 1, //8.854187817e-12,
	electron_mass = 1, //9.10938291e-31,
	elementary_charge = 1, //1.602176565e-19,
	electron_charge_mass_ratio = -elementary_charge / electron_mass,
	particle_temp_nrj_ratio = 1; //1.3806488e-23;


/*!
Wraps given index + offset to index within grid_size [0, grid_size[.
*/
size_t wrap(const size_t index, const int offset, const size_t grid_size)
{
	if (offset < 0) {
		return (index + grid_size - (std::abs(offset) % grid_size)) % grid_size;
	} else {
		return (index + offset) % grid_size;
	}
}

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
	return 2 * M_PI;
}

/*!
Returns physical length of the grid.
*/
constexpr double get_grid_length()
{
	return get_grid_end() - get_grid_start();
}


/*!
Returns the size of a grid cell.
*/
constexpr double get_cell_size()
{
	static_assert(
		std::tuple_size<Grid>::value > 0,
		"Grid must have at least one cell"
	);
	return (get_grid_end() - get_grid_start()) / std::tuple_size<Grid>::value;
}


/*!
Returns minimum coordinate of a cell located at given index in given grid.

Index starts from 0.
Returns a quiet NaN in case of error.
*/
double get_cell_min(const size_t index)
{
	static_assert(
		std::tuple_size<Grid>::value > 0,
		"Grid must have at least one cell"
	);

	if (index >= std::tuple_size<Grid>::value) {
		return std::numeric_limits<double>::quiet_NaN();
	}

	return index * get_cell_size();
}


/*!
Returns maximum coordinate of a cell located at given index in given grid.

Index starts from 0.
Returns a quiet NaN in case of error.
*/
double get_cell_max(const size_t index)
{
	static_assert(
		std::tuple_size<Grid>::value > 0,
		"Grid must have at least one cell"
	);

	if (index >= std::tuple_size<Grid>::value) {
		return std::numeric_limits<double>::quiet_NaN();
	}

	return (index + 1) * get_cell_size();
}


/*!
Returns center of a cell located at given index in given grid.

Index starts from 0.
Returns a quiet NaN in case of error.
*/
double get_cell_center(const size_t index)
{
	static_assert(
		std::tuple_size<Grid>::value > 0,
		"Grid must have at least one cell"
	);

	if (index >= std::tuple_size<Grid>::value) {
		return std::numeric_limits<double>::quiet_NaN();
	}

	return (0.5 + index) * get_cell_size();
}


/*!
Returns center of a neighbor located at given offset from a cell.
*/
double get_neighbor_center(const size_t index, const int offset)
{
	return get_cell_center(index) + offset * get_cell_size();
}


/*!
Returns best fits of a_n * sin(x + b_n) for electric field.

First value is average whose offset is 0, second value is
best fit to E - avg of a0 * sin(x + b0), third value is
best fit to E - avg - a0*sin(x+b0) of a1*sin(2*x+b1), etc.
*/
template<
	class Electric_Field_T
> std::vector<
	std::pair<double, double>
> get_electric_field_modes(
	const size_t nr_of_modes,
	const Grid& grid
) {
	std::vector<std::pair<double, double>> ret_val;

	std::array<double, grid_size> electric_field;
	for (size_t i = 0; i < grid_size; i++) {
		electric_field[i] = grid[i][Electric_Field_T()];
	}

	double avg_E = 0, max_E = 0;
	for (const auto& ele: electric_field) {
		avg_E += ele;
		if (max_E < fabs(ele)) {
			max_E = fabs(ele);
		}
	}
	avg_E /= grid_size;

	ret_val.push_back(std::make_pair(avg_E, 0));

	for (auto& ele: electric_field) {
		ele -= avg_E;
	}

	std::mt19937_64 random_source;
	random_source.seed(
		std::chrono::duration_cast<
			std::chrono::nanoseconds
		>(
			std::chrono::high_resolution_clock::now()
			- std::chrono::time_point<std::chrono::high_resolution_clock>()
		).count()
	);
	std::uniform_real_distribution<> offset_generator(0, 2 * M_PI);
	std::uniform_real_distribution<> magnitude_generator(0, max_E);

	constexpr size_t iterations = 10000;

	for (size_t mode = 1; mode <= nr_of_modes; mode++) {
		double l2_norm = 0;

		const auto fit_func
			= [&mode](
				const double x,
				const double magnitude,
				const double offset
			) {
				return magnitude * sin(mode * 2 * M_PI * x / get_grid_length() + offset);
			};

		const auto get_l2_norm
			= [&fit_func, &electric_field](
				const double magnitude,
				const double offset
			) {
				double norm = 0;
				for (size_t i = 0; i < grid_size; i++) {
					const double x = get_cell_center(i);
					norm += pow(electric_field[i] - fit_func(x, magnitude, offset), 2);
				}
				return sqrt(norm);
			};

		double offset = 0;
		l2_norm = std::numeric_limits<double>::max();
		for (size_t iter = 0; iter < iterations; iter++) {
			const double
				new_offset = offset_generator(random_source),
				new_norm = get_l2_norm(1, new_offset);

			if (l2_norm > new_norm) {
				l2_norm = new_norm;
				offset = new_offset;
			}
		}

		double magnitude = 1;
		l2_norm = std::numeric_limits<double>::max();
		for (size_t iter = 0; iter < iterations; iter++) {
			const double
				new_magnitude = magnitude_generator(random_source),
				new_norm = get_l2_norm(new_magnitude, offset);

			if (l2_norm > new_norm) {
				l2_norm = new_norm;
				magnitude = new_magnitude;
			}
		}

		ret_val.push_back(std::make_pair(magnitude, offset));

		for (size_t i = 0; i < grid_size; i++) {
			electric_field[i] -= fit_func(get_cell_center(i), magnitude, offset);
		}
	}

	return ret_val;
}


/*!
Accumulates electron charge to grid cells and sets total charge.
*/
template<
	class Particles_Internal_T,
	class Particle_Mass_T,
	class Particle_Charge_Mass_Ratio_T,
	class Particle_Position_T,
	class Electron_Charge_T,
	class Ion_Charge_T,
	class Total_Charge_T
> void accumulate_charge(Grid& grid)
{
	const Particle_Mass_T Mas{};
	const Particle_Charge_Mass_Ratio_T C2M{};
	const Particle_Position_T Pos{};
	const Electron_Charge_T Ele{};
	const Ion_Charge_T Ion{};
	const Total_Charge_T Tot{};

	const auto
		cell_length = get_cell_size(),
		half_cell = cell_length / 2;

	for (auto& cell: grid) {
		cell[Ele] = 0;
	}

	for (size_t cell_i = 0; cell_i < grid_size; cell_i++) {

		const size_t
			neigh_neg_i = wrap(cell_i, -1, grid_size),
			neigh_pos_i = wrap(cell_i, +1, grid_size);

		const typename Particle_Position_T::data_type
			cell_min(get_cell_min(cell_i), 0, 0),
			cell_max(get_cell_max(cell_i), 1, 1),
			neigh_neg_min(cell_min[0] - get_cell_size(), 0, 0),
			neigh_neg_max(cell_min[0], 1, 1),
			neigh_pos_min(cell_max[0], 0, 0),
			neigh_pos_max(cell_max[0] + get_cell_size(), 1, 1);

		auto
			&cell_charge = grid[cell_i][Ele],
			&neigh_neg_charge = grid[neigh_neg_i][Ele],
			&neigh_pos_charge = grid[neigh_pos_i][Ele];

		const auto& particles = grid[cell_i][Particles_Internal_T()]();
		for (const auto& particle: particles) {

			const auto particle_charge = particle[C2M] * particle[Mas];
			const typename Particle_Position_T::data_type
				particle_min(particle[Pos][0] - half_cell, 0, 0),
				particle_max(particle[Pos][0] + half_cell, 1, 1);

			cell_charge += pamhd::particle::get_accumulated_value(
				particle_charge,
				particle_min,
				particle_max,
				cell_min,
				cell_max
			);

			neigh_neg_charge += pamhd::particle::get_accumulated_value(
				particle_charge,
				particle_min,
				particle_max,
				neigh_neg_min,
				neigh_neg_max
			);

			neigh_pos_charge += pamhd::particle::get_accumulated_value(
				particle_charge,
				particle_min,
				particle_max,
				neigh_pos_min,
				neigh_pos_max
			);
		}

		if (cell_i > 1) {
			auto& cell_data = grid[cell_i - 1];
			cell_data[Tot] = cell_data[Ion] + cell_data[Ele];
		}
	}

	grid[0][Tot] = grid[0][Ion] + grid[0][Ele];
	grid[grid_size - 1][Tot] = grid[grid_size - 1][Ion] + grid[grid_size - 1][Ele];
}


template<
	class Total_Charge_T,
	class Electric_Field_T,
	class Electric_Potential_T
> void set_electric_field(
	Grid& grid,
	const double vacuum_permittivity
) {
	const Electric_Field_T Ele{};
	const Electric_Potential_T Pot{};

	constexpr double cell_length = get_cell_size();

	// rhs is total charge density
	std::vector<double> rhs;
	for (const auto& cell: grid) {
		rhs.push_back(-cell[Total_Charge_T()] / vacuum_permittivity / cell_length);
	}

	const auto solution
		= pamhd::poisson::solve_failsafe(
			grid_size,
			1.0,
			1.0,
			cell_length,
			1.0,
			1.0,
			rhs
		);

	for (size_t cell_i = 0; cell_i < grid_size; cell_i++) {
		grid[cell_i][Pot] = solution[cell_i];

		const size_t
			neigh_neg_i = wrap(cell_i, -1, grid_size),
			neigh_pos_i = wrap(cell_i, +1, grid_size);

		grid[cell_i][Ele]
			= -(solution[neigh_pos_i] - solution[neigh_neg_i])
			/ (2 * cell_length);
	}
}


/*!
random_source is assumed to be std::mt19937 or similar.

total_mass is the total mass to assign
to created particles in each cell.

species_mass is the mass to use for
probability distribution of particle velocities.
*/
template <
	class Particles_T,
	class Particle_T,
	class Particle_Mass_T,
	class Particle_Charge_Mass_Ratio_T,
	class Particle_Position_T,
	class Particle_Velocity_T,
	class Electron_Charge_T,
	class Ion_Charge_T,
	class Total_Charge_T,
	class Random_Source
> void initialize_particles(
	const Eigen::Vector3d& bulk_velocity,
	const double temperature,
	const size_t particles_per_cell,
	const double charge_mass_ratio,
	const double total_mass_per_cell,
	const double species_mass,
	const double particle_temp_nrj_ratio,
	Grid& grid,
	Random_Source& random_source
) {
	for (size_t cell = 0; cell < grid_size; cell++) {
		const double
			cell_center = get_cell_center(cell),
			cell_min = get_cell_min(cell),
			cell_max = get_cell_max(cell),
			strength_of_nonlinearity = 0.01,
			final_mass
				= total_mass_per_cell
				+ total_mass_per_cell
					* strength_of_nonlinearity
					* cos(2 * M_PI * cell_center / get_grid_length());

		grid[cell][Particles_T()]()
			= pamhd::particle::create_particles<
				Particle_T,
				Particle_Mass_T,
				Particle_Charge_Mass_Ratio_T,
				Particle_Position_T,
				Particle_Velocity_T,
				Random_Source
			>(
				bulk_velocity,
				Eigen::Vector3d(cell_min, 0, 0),
				Eigen::Vector3d(cell_max, 1, 1),
				Eigen::Vector3d(temperature, temperature, temperature),
				particles_per_cell,
				charge_mass_ratio,
				final_mass,
				species_mass,
				particle_temp_nrj_ratio,
				random_source
			);

		// average total charge = 0
		grid[cell][Ion_Charge_T()]	= -total_mass_per_cell * charge_mass_ratio;
	}

	accumulate_charge<
		Particles_T,
		Particle_Mass_T,
		Particle_Charge_Mass_Ratio_T,
		Particle_Position_T,
		Electron_Charge_T,
		Ion_Charge_T,
		Total_Charge_T
	>(grid);
}


template<
	class Particles_Internal_T,
	class Particles_External_T,
	class Mass_T,
	class Charge_Mass_Ratio_T,
	class Position_T,
	class Velocity_T,
	class Electric_Field_T
> double solve_particles(Grid& grid, const double dt)
{
	static_assert(
		std::tuple_size<Grid>::value > 0,
		"Grid must have at least one cell"
	);

	using std::fabs;
	using std::max;
	using std::min;

	const Position_T Pos{};
	const Velocity_T Vel{};
	const Charge_Mass_Ratio_T C2M{};
	const Electric_Field_T Ele{};

	double max_dt = std::numeric_limits<double>::max();

	// propagate particles
	for (size_t cell_i = 0; cell_i < grid_size; cell_i++) {

		const double
			x_neg = get_neighbor_center(cell_i, -1),
			x_mid = get_cell_center(cell_i),
			x_pos = get_neighbor_center(cell_i, 1),
			cell_size = get_cell_size(),
			cell_start = x_mid - cell_size / 2,
			cell_end = x_mid + cell_size / 2;

		// neighborhood cell data
		auto
			&cell = grid[cell_i],
			&neigh_neg = grid[wrap(cell_i, -1, grid_size)],
			&neigh_pos = grid[wrap(cell_i, 1, grid_size)];

		// data to be interpolated to particle positions
		const Eigen::Array3d
			x_all{x_neg, x_mid, x_pos},
			E_all{
				neigh_neg[Ele],
				cell[Ele],
				neigh_pos[Ele]
			};

		auto& particles_int = cell[Particles_Internal_T()]();
		for (size_t particle_i = 0; particle_i < particles_int.size(); particle_i++) {
			auto& particle = particles_int[particle_i];

			/*
			Interpolate variables required for Lorentz force
			to particle position.

			Geometric coefficients for interpolating
			cell-centered values linearly

			               _0__
			             _/    \p_     particle between 2nd & 3rd cell
                       _/         0
			          0
			cells: |--X--|--X--|--X--|
			coeffs:   |   1st   | (max(0, < 0) == 0)
			                |2nd|
			                    | | 3rd
			*/
			const Eigen::Array3d geom_coeffs{
				max(0.0, 1 - (particle[Pos][0] - x_all[0]) / (x_all[1] - x_all[0])),
				1 - fabs(particle[Pos][0] - x_all[1]) / (x_all[1] - x_all[0]),
				max(0.0, 1 - (x_all[2] - particle[Pos][0]) / (x_all[2] - x_all[1]))
			};

			// variables at particle position
			const double E_interp = (E_all * geom_coeffs).sum();

			// propagate particle
			std::tie(
				particle[Pos],
				particle[Vel]
			) = pamhd::particle::propagate(
				particle[Pos],
				particle[Vel],
				{E_interp, 0, 0},
				{0, 0, 0},
				particle[C2M],
				dt
			);

			// return maximum allowed time step for next step
			const auto minmax_step
				= pamhd::particle::get_minmax_step(
					0,
					1.0 / 20.0,
					cell_size,
					particle[C2M],
					particle[Vel],
					{0, 0, 0}
				);
			//std::cout << minmax_step.second << std::endl;

			max_dt = min(max_dt, minmax_step.second);

			// assign particle to another cell
			if (particle[Pos][0] < cell_start or particle[Pos][0] > cell_end) {
				bool assign_to_neg_neigh = false;
				if (particle[Pos][0] < cell_start) {
					assign_to_neg_neigh = true;
				}

				// maybe wrap around particle position
				if (particle[Pos][0] < get_grid_start()) {
					particle[Pos][0] += get_grid_end() - get_grid_start();
				} else if (particle[Pos][0] > get_grid_end()) {
					particle[Pos][0] -= get_grid_end() - get_grid_start();
				}

				if (assign_to_neg_neigh) {
					neigh_neg[Particles_External_T()]().push_back(particle);
				} else {
					neigh_pos[Particles_External_T()]().push_back(particle);
				}

				// remove from original cell
				particles_int[particle_i] = particles_int[particles_int.size() - 1];
				particles_int.pop_back(); // particle no longer a valid reference
				particle_i--;
			}
		}
	}

	// assign particles to cell at their new position
	for (auto& cell: grid) {
		auto
			&particles_int = cell[Particles_Internal_T()](),
			&particle_ext = cell[Particles_External_T()]();

		particles_int.insert(
			particles_int.end(),
			particle_ext.cbegin(),
			particle_ext.cend()
		);
		particle_ext.clear();
	}

	return max_dt;
}


void plot_particles_bulk(
	const double simulation_time,
	const Grid& grid
) {
	std::ostringstream time_string;
	time_string
		<< std::setw(7)
		<< std::setfill('0')
		<< size_t(simulation_time * 1000);

	const std::string
		gnuplot_file_name(
			"particle_"
			+ time_string.str()
			+ "_ms.dat"
		),
		field_plot_file_name(
			"E_field_"
			+ time_string.str()
			+ "_ms.png"
		),
		electron_charge_plot_file_name(
			"electron_charge_"
			+ time_string.str()
			+ "_ms.png"
		),
		total_charge_plot_file_name(
			"E_charge_"
			+ time_string.str()
			+ "_ms.png"
		),
		potential_plot_file_name(
			"E_potential_"
			+ time_string.str()
			+ "_ms.png"
		);

	std::ofstream gnuplot_file(gnuplot_file_name);

	gnuplot_file
		<< std::setprecision(20)
		<<
			"set term png enhanced\n"
			"set xlabel 'X'\n"
			"set format x '%.1e'\n"
			"set format y '%.5e'\n"
			"set xtics "
		<< boost::lexical_cast<string>((get_grid_end() - get_grid_start()) / 5)
		<< "\n";

	// E
	gnuplot_file << "set output '" << field_plot_file_name <<
		"'\nset ylabel 'Electric field'\n"
		"plot '-' using 1:2 with line linewidth 2 title ''\n";
	for (size_t cell_i = 0; cell_i < grid_size; cell_i++) {
		gnuplot_file << get_cell_center(cell_i) << " "
			<< grid[cell_i][Electric_Field()] << "\n";
	}

	// electron charge
	gnuplot_file << "end\nset output '" << electron_charge_plot_file_name
		<< "'\nset ylabel 'Electron charge'\n"
		   "plot '-' u 1:2 w l lw 2 t ''\n";
	for (size_t cell_i = 0; cell_i < grid_size; cell_i++) {
		gnuplot_file << get_cell_center(cell_i) << " "
			<< grid[cell_i][Electron_Charge()] << "\n";
	}

	// total charge
	gnuplot_file << "end\nset output '" << total_charge_plot_file_name
		<< "'\nset ylabel 'Total electric charge'\n"
		   "plot '-' u 1:2 w l lw 2 t ''\n";
	for (size_t cell_i = 0; cell_i < grid_size; cell_i++) {
		gnuplot_file << get_cell_center(cell_i) << " "
			<< grid[cell_i][Total_Charge()] << "\n";
	}

	// potential
	gnuplot_file << "end\nset output '" << potential_plot_file_name
		<< "'\nset ylabel 'Electric potential'\n"
		   "plot '-' u 1:2 w l lw 2 t ''\n";
	for (size_t cell_i = 0; cell_i < grid_size; cell_i++) {
		gnuplot_file << get_cell_center(cell_i) << " "
			<< grid[cell_i][Electric_Potential()] << "\n";
	}
	gnuplot_file << "end\n";

	gnuplot_file.close();

	system(("gnuplot " + gnuplot_file_name).c_str());
}



/*!
Plots number of particles as a function of two of their variables.
*/
template<
	class Particles_T,
	class Horizontal_Variable,
	class Vertical_Variable,
	size_t Width,
	size_t Height
> void plot_pvec_pvec(
	const Grid& grid,
	const double simulation_time,
	size_t index_of_horizontal,
	size_t index_of_vertical
) {
	static_assert(
		std::tuple_size<Grid>::value > 0,
		"Grid must have at least one cell"
	);
	static_assert(
		Width > 0,
		"Horizontal direction must have at least one bin"
	);
	static_assert(
		Height > 0,
		"Vertical direction must have at least one bin"
	);

	std::array<std::array<size_t, Width>, Height> binned;
	for (auto& row: binned) {
		for (auto& number_of_particles: row) {
			number_of_particles = 0;
		}
	}

	double
		horiz_min = std::numeric_limits<double>::max(),
		horiz_max = std::numeric_limits<double>::lowest(),
		vert_min = std::numeric_limits<double>::max(),
		vert_max = std::numeric_limits<double>::lowest();

	for (const auto& cell: grid) {
		for (const auto& particle: cell[Particles_T()]()) {

			const auto& horizontal = particle[Horizontal_Variable()];
			if (horiz_min > horizontal[index_of_horizontal]) {
				horiz_min = horizontal[index_of_horizontal];
			}
			if (horiz_max < horizontal[index_of_horizontal]) {
				horiz_max = horizontal[index_of_horizontal];
			}

			const auto& vertical = particle[Vertical_Variable()];
			if (vert_min > vertical[index_of_vertical]) {
				vert_min = vertical[index_of_vertical];
			}
			if (vert_max < vertical[index_of_vertical]) {
				vert_max = vertical[index_of_vertical];
			}
		}
	}

	const double
		horiz_bin_size = (horiz_max - horiz_min) / Width,
		vert_bin_size = (vert_max - vert_min) / Height;

	for (const auto& cell: grid) {
		for (const auto& particle: cell[Particles_T()]()) {
			const auto& horizontal = particle[Horizontal_Variable()];
			const size_t horiz_bin
				= std::min(
					size_t(
						std::floor(
							(horizontal[index_of_horizontal] - horiz_min)
							/ horiz_bin_size
						)
					),
					Width - 1
				);

			const auto& vertical = particle[Vertical_Variable()];
			const size_t vert_bin
				= std::min(
					size_t(
						std::floor(
							(vertical[index_of_vertical] - vert_min)
							/ vert_bin_size
						)
					),
					Height - 1
				);

			binned[vert_bin][horiz_bin]++;
		}
	}

	std::ostringstream time_string;
	time_string
		<< std::setw(10)
		<< std::setfill('0')
		<< size_t(simulation_time * 1000000)
		<< "_us";

	const std::string
		separator("_"),
		file_name_prefix(
			"particle" + separator
			+ Horizontal_Variable::get_name() + "_"
			+ boost::lexical_cast<std::string>(index_of_horizontal) + "_"
			+ Vertical_Variable::get_name() + "_"
			+ boost::lexical_cast<std::string>(index_of_vertical) + "_"
			+ time_string.str()
			+ "."
		),
		gnuplot_file_name(file_name_prefix + "dat"),
		plot_file_name(file_name_prefix + "png");


	std::ofstream gnuplot_file(gnuplot_file_name);

	gnuplot_file
		<< "set term png enhanced\n"
		<< "set output '" << plot_file_name
		<< "'\nset format x '%.1e'\n"
		<< "set format y '%.1e'"
		<< "\nset xlabel '" << Horizontal_Variable::get_name() + " "
		<< boost::lexical_cast<std::string>(index_of_horizontal)
		<< "'\nset ylabel '" << Vertical_Variable::get_name() + " "
		<< boost::lexical_cast<std::string>(index_of_vertical)
		<< "'\nset cblabel 'Number of particles'\nx0 = "
		<< boost::lexical_cast<std::string>(horiz_min)
		<< "\nset xtics "
		<< boost::lexical_cast<string>((horiz_max - horiz_min) / 6)
		<< "\nset ytics "
		<< boost::lexical_cast<string>((vert_max - vert_min) / 6)
		<< "\nxincr = " << boost::lexical_cast<std::string>(horiz_bin_size)
		<< "\nyincr = " << boost::lexical_cast<std::string>(vert_bin_size)
		<< "\ny0 = " << boost::lexical_cast<std::string>(vert_min)
		<< "\nset xrange ["
		<< boost::lexical_cast<std::string>(horiz_min - horiz_bin_size / 2)
		<< " : " << boost::lexical_cast<std::string>(horiz_max - horiz_bin_size / 2)
		<< "]\nset yrange ["
		<< boost::lexical_cast<std::string>(vert_min - vert_bin_size / 2)
		<< " : " << boost::lexical_cast<std::string>(vert_max - vert_bin_size / 2)
		<< "]\nplot '-' u ($1*xincr + x0):($2*yincr + y0):3 matrix with image title ''\n";

	for (const auto& row: binned) {
		for (const auto& number_of_particles: row) {
			gnuplot_file << number_of_particles << " ";
		}
		gnuplot_file << "\n";
	}
	gnuplot_file << "end" << std::endl;

	system(("gnuplot " + gnuplot_file_name).c_str());
}


template <
	class Particles_T,
	class Particle_T,
	class Particle_Velocity_T,
	class Particle_Position_T
> std::string save_particles(
	const Grid& grid,
	const double simulation_time
) {
	std::ostringstream time_string;
	time_string
		<< std::setw(6)
		<< std::setfill('0')
		<< size_t(simulation_time * 1000);

	const std::string
		gnuplot_file_name(
			"particle_"
			+ time_string.str()
			+ "_ms.dat"
		),
		position_plot_file_name(
			"particle_pos_"
			+ time_string.str()
			+ "_ms.png"
		),
		velocity_plot_file_name(
			"particle_vel_"
			+ time_string.str()
			+ "_ms.png"
		),
		velocity_x_plot_file_name(
			"particle_vel_x_"
			+ time_string.str()
			+ "_ms.png"
		),
		velocity_y_plot_file_name(
			"particle_vel_y_"
			+ time_string.str()
			+ "_ms.png"
		),
		velocity_z_plot_file_name(
			"particle_vel_z_"
			+ time_string.str()
			+ "_ms.png"
		),
		vel_dist_plot_file_name(
			"particle_vdist_"
			+ time_string.str()
			+ "_ms.png"
		),
		vel_spectrum_plot_file_name(
			"particle_vspec_"
			+ time_string.str()
			+ "_ms.png"
		);

	std::ofstream gnuplot_file(gnuplot_file_name);

	gnuplot_file
		<< "set term png enhanced\n"
		   "set xlabel 'X'\n"
		   "set ylabel 'Cell index'\n";
	/*	   "set output '" << position_plot_file_name << "'\n"
		   "plot '-' using 1:2 with points title 'Particle position'\n";
	for (size_t cell_i = 0; cell_i < grid.size(); cell_i++) {
		const auto& cell = grid[cell_i];

		for (const auto& particle: cell[Particles_T()]()) {

			const auto& position = particle[Particle_Position_T()];
			gnuplot_file << position[0] << " " << cell_i << "\n";
		}
	}
	gnuplot_file << "end\n";

	gnuplot_file
		<< "set output '" << velocity_plot_file_name << "'\n"
		<< "set xlabel 'V_x'\n"
		   "set ylabel 'V_y'\n"
		   "set zlabel 'V_z'\n"
		   "splot '-' using 1:2:3 with points title 'Particle velocity'\n";
	for (const auto& cell: grid) {
		for (const auto& particle: cell[Particles_T()]()) {

			const auto& velocity = particle[Particle_Velocity_T()];

			// always plot 3 components for position
			for (size_t i = 0; i < 3; i++) {
				if (i < velocity.size()) {
					gnuplot_file << velocity[i] << " ";
				} else {
					gnuplot_file << "0 ";
				}
			}
			gnuplot_file << "\n";
		}
	}
	gnuplot_file << "end\n";

	gnuplot_file
		<< "set output '" << velocity_x_plot_file_name << "'\n"
		<< "set xlabel 'V_y'\n"
		   "set ylabel 'V_z'\n"
		   "plot '-' using 1:2 with points title 'Particle velocity'\n";
	for (const auto& cell: grid) {
		for (const auto& particle: cell[Particles_T()]()) {

			const auto& velocity = particle[Particle_Velocity_T()];
			for (size_t i = 1; i < 3; i++) {
				if (i < velocity.size()) {
					gnuplot_file << velocity[i] << " ";
				} else {
					gnuplot_file << "0 ";
				}
			}
			gnuplot_file << "\n";
		}
	}
	gnuplot_file << "end\n";

	gnuplot_file
		<< "set output '" << velocity_y_plot_file_name << "'\n"
		<< "set xlabel 'V_x'\n"
		   "set ylabel 'V_z'\n"
		   "plot '-' using 1:2 with points title 'Particle velocity'\n";
	for (const auto& cell: grid) {
		for (const auto& particle: cell[Particles_T()]()) {

			const auto& velocity = particle[Particle_Velocity_T()];

			// always plot 2 components for position
			for (size_t i = 0; i < 3; i++) {
				if (i == 1) {
					continue;
				}
				if (i < velocity.size()) {
					gnuplot_file << velocity[i] << " ";
				} else {
					gnuplot_file << "0 ";
				}
			}
			gnuplot_file << "\n";
		}
	}
	gnuplot_file << "end\n";

	gnuplot_file
		<< "set output '" << velocity_z_plot_file_name << "'\n"
		<< "set xlabel 'V_x'\n"
		   "set ylabel 'V_y'\n"
		   "plot '-' using 1:2 with points title 'Particle velocity'\n";
	for (const auto& cell: grid) {
		for (const auto& particle: cell[Particles_T()]()) {

			const auto& velocity = particle[Particle_Velocity_T()];

			// always plot 2 components for position
			for (size_t i = 0; i < 2; i++) {
				if (i < velocity.size()) {
					gnuplot_file << velocity[i] << " ";
				} else {
					gnuplot_file << "0 ";
				}
			}
			gnuplot_file << "\n";
		}
	}
	gnuplot_file << "end\n";

	gnuplot_file
		<< "set output '" << vel_dist_plot_file_name << "'\n"
		<< "binwidth = 5000\n"
		   "set boxwidth binwidth\n"
		   "bin(x,width)=width*floor(x/width) + binwidth/2.0\n"
		   "set xlabel 'Velocity (m/s)'\n"
		   "set ylabel 'Number of particles'\n"
		   "plot '-' using (bin($1,binwidth)):(1.0) smooth freq with boxes title ''\n";
	for (const auto& cell: grid) {
		for (const auto& particle: cell[Particles_T()]()) {

			const auto& velocity = particle[Particle_Velocity_T()];
			double v_magnitude = 0;
			for (size_t i = 0; i < 3; i++) {
				if (i < velocity.size()) {
					v_magnitude += std::pow(velocity[i], 2);
				}
			}
			v_magnitude = std::sqrt(v_magnitude);
			gnuplot_file << v_magnitude << "\n";
		}
	}
	gnuplot_file << "end" << std::endl;*/

	/*
	Prepare velocity spectrum data
	*/
	constexpr size_t velocity_bins = 50;

	// velocity spectrogram as a function of tube position
	std::array<
		std::array<
			// number of particles
			size_t,
			// cell index bins
			std::tuple_size<Grid>::value
		>,
		// particle velocity bins
		velocity_bins
	> vspec_data;
	for (auto& vbin: vspec_data) {
		for (auto& number_of_particles: vbin) {
			number_of_particles = 0;
		}
	}

	double max_v = -1;
	for (const auto& cell: grid) {
		for (const auto& particle: cell[Particles_T()]()) {
			const double velocity = particle[Particle_Velocity_T()].norm();
			max_v = std::max(max_v, velocity);
		}
	}

	const double vbin_size = max_v / velocity_bins;

	gnuplot_file
		<< "set output '" << vel_spectrum_plot_file_name
		<< "'\nset xlabel 'Position (X / X_0, X_0 = "
		<< boost::lexical_cast<std::string>(get_grid_end() - get_grid_start())
		<< ")\nset ylabel 'Velocity (m/s)'\n"
		   "set cblabel 'Number of particles'\nx0 = "
		<< boost::lexical_cast<std::string>(get_grid_start())
		<< "\nxincr = "
		<< boost::lexical_cast<std::string>(
			1.0 / grid_size
		)
		<< "\nyincr = " << boost::lexical_cast<std::string>(vbin_size)
		<< "\ny0 = 0\n"
		   "plot '-' u ($1*xincr - x0):($2*yincr - y0):3 matrix with image title ''\n";

	for (size_t cell_i = 0; cell_i < grid.size(); cell_i++) {
		const auto& cell = grid[cell_i];
		for (const auto& particle: cell[Particles_T()]()) {
			const double velocity = particle[Particle_Velocity_T()].norm();
			const size_t vbin
				= std::min(
					size_t(std::floor(velocity / vbin_size)),
					velocity_bins - 1
				);
			vspec_data[vbin][cell_i]++;
		}
	}

	for (const auto& vbin: vspec_data) {
		for (const auto& number_of_particles: vbin) {
			gnuplot_file << number_of_particles << " ";
		}
		gnuplot_file << "\n";
	}
	gnuplot_file << "end" << std::endl;

	return gnuplot_file_name;
}

template<
	class Horiz_Var,
	class Vert_Var,
	size_t horiz_cells,
	size_t vert_cells
> struct Plotter {
	const Grid& grid;

	Plotter(const Grid& given_grid) : grid(given_grid) {};

	void operator()(
		const double simulation_time,
		const size_t i1,
		const size_t i2
	) const {
		plot_pvec_pvec<
			pamhd::particle::Particles_Int,
			Horiz_Var,
			Vert_Var,
			horiz_cells,
			vert_cells
		>(this->grid, simulation_time, i1, i2);
	}
};


void plot_time_series(
	const std::vector<double>& time_stamps,
	const std::vector<double>& data,
	const std::string& plot_file_prefix,
	const std::string& ylabel
) {
	if (time_stamps.size() != data.size()) {
		std::cerr <<  __FILE__ << " (" << __LINE__ << "): "
			<< "Number of time stamps doens't match number of data points"
			<< std::endl;
		abort();
	}

	const std::string
		gnuplot_file_name(plot_file_prefix + ".dat"),
		plot_file_name(plot_file_prefix + ".png");

	std::ofstream gnuplot_file(gnuplot_file_name);

	gnuplot_file
		<< std::setprecision(20)
		<< "set term png enhanced\nset logscale y\nset xlabel 'Time'\nset ylabel '"
		<< ylabel
		<< "'\nset format x '%.1e'\nset format y '%.1e'\nset output '"
		<< plot_file_name
		<< "\nplot '-' using 1:2 with line linewidth 2 title ''\n";
	for (size_t i = 0; i < time_stamps.size(); i++) {
		gnuplot_file << time_stamps[i] << " " << data[i] << "\n";
	}
	gnuplot_file << "end\n";

	gnuplot_file.close();

	system(("gnuplot " + gnuplot_file_name).c_str());
}



int main(int argc, char* argv[])
{
	const size_t particles_per_cell = 2000;
	const double
		temperature = 0.4,
		electron_mass_density = M_PI,
		simulation_duration = 100;
	const Eigen::Vector3d
		bulk_velocity(0, 0, 0);


	std::mt19937 random_source;
	random_source.seed(1);

	bool save = false, plot = false, no_verify = false, verbose = false;
	const std::string solver_str("default");
	boost::program_options::options_description options(
		"Usage: program_name [options], where options are:"
	);
	options.add_options()
		("help", "Print this help message")
		("save", "Save end result to ascii file")
		("plot", "Plot results using gnuplot")
		("no-verify", "Do not verify against reference results")
		("verbose", "Print run time information");

	boost::program_options::variables_map option_variables;
	boost::program_options::store(
		boost::program_options::parse_command_line(argc, argv, options),
		option_variables
	);
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


	Grid grid;

	const Plotter<
		pamhd::particle::Position,
		pamhd::particle::Velocity,
		grid_size,
		40
	> r_v_plotter(grid);
	const Plotter<
		pamhd::particle::Velocity,
		pamhd::particle::Velocity,
		40,
		40
	> v_v_plotter(grid);

	std::vector<double> time_series_sim_time, time_series_E2;

	if (verbose) {
		cout << "Initializing particles" << endl;
	}

	initialize_particles<
		pamhd::particle::Particles_Int,
		pamhd::particle::Particle,
		pamhd::particle::Mass,
		pamhd::particle::Charge_Mass_Ratio,
		pamhd::particle::Position,
		pamhd::particle::Velocity,
		Electron_Charge,
		Ion_Charge,
		Total_Charge
	>(
		bulk_velocity,
		temperature,
		particles_per_cell,
		electron_charge_mass_ratio,
		electron_mass_density * get_cell_size(),
		electron_mass,
		particle_temp_nrj_ratio,
		grid,
		random_source
	);

	set_electric_field<
		Total_Charge,
		Electric_Field,
		Electric_Potential
	>(grid, vacuum_permittivity);

	std::ofstream electric_field_file;

	constexpr size_t nr_of_E_modes = 3;
	if (save) {
		if (verbose) {
			cout << "Saving things..." << endl;
		}

		electric_field_file.open("electric_field.dat");
		if (not electric_field_file.is_open()) {
			cerr << "Couldn't open electric field file." << endl;
			abort();
		}

		const auto E_modes
			= get_electric_field_modes<Electric_Field>(nr_of_E_modes, grid);

		electric_field_file <<
			"#set term png enhanced\n"
			"#set output 'E_field_normal_modes.png'\n"
			"#set xlabel 'X'\n"
			"#set format x '%.1e'\n"
			"#set format y '%.5e'\n0 ";
		for (size_t i = 0; i < E_modes.size(); i++) {
			electric_field_file << E_modes[i].first << " " << E_modes[i].second << " ";
		}
		electric_field_file << endl;
	}

	if (plot) {
		if (verbose) {
			cout << "Plotting initial state of particles" << endl;
		}

		r_v_plotter(0, 0, 0);
		r_v_plotter(0, 0, 1);
		r_v_plotter(0, 0, 2);
		v_v_plotter(0, 0, 1);
		v_v_plotter(0, 0, 2);
		v_v_plotter(0, 1, 2);

		plot_particles_bulk(0, grid);
	}


	const double
		particle_plot_interval = simulation_duration / 100,
		save_interval = simulation_duration / 10000;

	double
		max_dt = 0,
		simulation_time = 0,
		particle_next_plot = particle_plot_interval,
		next_save = save_interval;
	while (simulation_time < simulation_duration) {

		const double
			time_step_factor = 0.5, // applied to CFL(-like) limits
			// don't step over the final simulation time
			until_end = simulation_duration - simulation_time,
			time_step = std::min(time_step_factor * max_dt, until_end);

		max_dt = std::numeric_limits<double>::max();

		if (verbose) {
			cout << "Solving particles at time " << simulation_time
				<< " s with time step " << time_step << endl;
		}

		max_dt = std::min(
			max_dt,
			solve_particles<
				pamhd::particle::Particles_Int,
				Particles_Ext,
				pamhd::particle::Mass,
				pamhd::particle::Charge_Mass_Ratio,
				pamhd::particle::Position,
				pamhd::particle::Velocity,
				Electric_Field
			>(grid, time_step)
		);

		accumulate_charge<
			pamhd::particle::Particles_Int,
			pamhd::particle::Mass,
			pamhd::particle::Charge_Mass_Ratio,
			pamhd::particle::Position,
			Electron_Charge,
			Ion_Charge,
			Total_Charge
		>(grid);

		set_electric_field<
			Total_Charge,
			Electric_Field,
			Electric_Potential
		>(grid, vacuum_permittivity);

		simulation_time += time_step;

		time_series_sim_time.push_back(simulation_time);
		double new_E2 = 0;
		for (const auto& cell: grid) {
			new_E2 += vacuum_permittivity * pow(cell[Electric_Field()], 2);
		}
		time_series_E2.push_back(new_E2 / grid_size);

		if (save && next_save <= simulation_time) {
			next_save += save_interval;

			if (verbose) {
				cout << "Saving at time " << simulation_time << endl;
			}

			const auto E_modes
				= get_electric_field_modes<Electric_Field>(nr_of_E_modes, grid);

			electric_field_file << simulation_time << " ";
			for (size_t i = 0; i < E_modes.size(); i++) {
				electric_field_file
					<< E_modes[i].first << " "
					<< E_modes[i].second << " ";
			}
			electric_field_file << endl;
		}

		if (plot && particle_next_plot <= simulation_time) {
			particle_next_plot += particle_plot_interval;

			if (verbose) {
				cout << "Plotting particles at time " << simulation_time << endl;
			}

			r_v_plotter(simulation_time, 0, 0);
			r_v_plotter(simulation_time, 0, 1);
			r_v_plotter(simulation_time, 0, 2);
			v_v_plotter(simulation_time, 0, 1);
			v_v_plotter(simulation_time, 0, 2);
			v_v_plotter(simulation_time, 1, 2);

			plot_particles_bulk(simulation_time, grid);
		}
	}

	plot_time_series(
		time_series_sim_time,
		time_series_E2,
		"E2_time_series",
		"Total second power of electric field"
	);

	if (save) {
		electric_field_file.close();
	}

	if (verbose) {
		cout << "Simulation finished at time " << simulation_time << endl;
	}

	return EXIT_SUCCESS;
}
