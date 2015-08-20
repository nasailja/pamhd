/*
Program for propagating charged particles in 1+3d with analytic fields.

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

#include "cmath"
#include "cstdlib"
#include "fstream"
#include "iostream"
#include "random"

#include "boost/lexical_cast.hpp"
#include "boost/numeric/odeint.hpp"
#include "boost/program_options.hpp"
#include "Eigen/Core"
#include "Eigen/Geometry"
#include "mpi.h"

#include "particle/common.hpp"

constexpr double
	electron_mass = 9.10938291e-31,
	electron_charge = -1.602176565e-19,
	electron_c2m = electron_charge / electron_mass,
	vacuum_permittivity = 8.854187817e-12,
	light_speed = 299792458,
	light_speed2 = light_speed*light_speed;


using namespace std;
using namespace boost::numeric::odeint;
using namespace Eigen;


/*
Returns electric and magnetic fields at given time and position.

length == length of simulation box,
frequency == frequency of EM wave,
e1, b1 == amplitudes of EM wave perpendicular to sim box
b0 == component of magnetic field parallel to sim box
*/
struct Fields
{
	const double
		length,
		frequency,
		e1,
		b0,
		b1,
		inv_wavelength_mid,
		inv_wavelength_ends;

	Fields(
		const double given_length,
		const double given_frequency,
		const double given_e1,
		const double given_b0,
		const double given_b1,
		const double given_wavelength_mid,
		const double given_wavelength_ends
	) :
		length(given_length),
		frequency(given_frequency),
		e1(given_e1),
		b0(given_b0),
		b1(given_b1),
		inv_wavelength_mid(1.0 / given_wavelength_mid),
		inv_wavelength_ends(1.0 / given_wavelength_ends)
	{}

	std::string print() const
	{
		std::string ret_val("2*pi*x*");
		ret_val += boost::lexical_cast<std::string>(this->inv_wavelength_ends);
		ret_val += " - 2*pi*t*";
		ret_val += boost::lexical_cast<std::string>(frequency);
		ret_val += " + (1/";
		ret_val += boost::lexical_cast<std::string>(1.0/this->inv_wavelength_ends);
		ret_val += " - 1/";
		ret_val += boost::lexical_cast<std::string>(1.0/this->inv_wavelength_mid);
		ret_val += ")*2*";
		ret_val += boost::lexical_cast<std::string>(this->length);
		ret_val += "*cos(pi*x/";
		ret_val += boost::lexical_cast<std::string>(this->length);
		ret_val += ")";
		return ret_val;
	}

	std::array<Eigen::Vector3d, 2> operator()(
		const double time,
		const double position
	) const {
		const double
			arg
				= 2 * M_PI * position * this->inv_wavelength_ends
				- 2 * M_PI * frequency * time
				+ (this->inv_wavelength_ends - this->inv_wavelength_mid)
					* 2 * this->length
					* cos(M_PI * position / this->length);

		const Eigen::Vector3d
			E{
				0,
				this->e1 * sin(arg),
				this->e1 * sin(arg + M_PI_2)
			},
			B{
				this->b0,
				this->b1 * sin(arg - M_PI_2),
				this->b1 * sin(arg)
			};

		return {E, B};
	}
};


// [0] == position, [1] == velocity
using state_t = std::array<Vector3d, 2>;

struct Propagator
{
	const double c2m;
	const Fields fields;

	// c2m = charge to mass ratio of particle
	Propagator(
		const double given_c2m,
		const Fields given_fields
	) :
		c2m(given_c2m),
		fields(given_fields)
	{}

	void operator()(const state_t& state, state_t& change, const double time) const
	{
		const double relativity_factor = sqrt(1 - state[1].squaredNorm() / light_speed2);

		const auto EB = this->fields(time, state[0][0]);
		change[0] = state[1];
		change[1]
			= this->c2m
			* relativity_factor
			* (EB[0] + state[1].cross(EB[1]));
	}
};



int main(int argc, char* argv[])
{
	using std::pow;
	using std::sqrt;

	if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
		std::cerr << "Couldn't initialize MPI." << std::endl;
		abort();
	}

	MPI_Comm comm = MPI_COMM_WORLD;

	int rank_tmp = 0, comm_size_tmp = 0;
	if (MPI_Comm_rank(comm, &rank_tmp) != MPI_SUCCESS) {
		std::cerr << "Couldn't obtain MPI rank." << std::endl;
		abort();
	}
	if (MPI_Comm_size(comm, &comm_size_tmp) != MPI_SUCCESS) {
		std::cerr << "Couldn't obtain size of MPI communicator." << std::endl;
		abort();
	}

	if (rank_tmp < 0) {
		std::cerr << "Negative MPI rank not supported: " << rank_tmp << std::endl;
		abort();
	}
	const size_t rank = size_t(rank_tmp);
	if (comm_size_tmp < 0) {
		std::cerr << "Negative MPI communicator size not supported: " << comm_size_tmp << std::endl;
		abort();
	}
	const size_t comm_size = size_t(comm_size_tmp);

	bool verbose = true;
	std::string outfile_name("out");
	double
		length = 1e5,
		time_start = 0,
		time_length = 1e4,
		save_dt = 0,
		save_dx = 0,
		time_step_factor = 1,
		em_frequency = 1e4,
		density_mid = 1e7,
		density_ends = 1e8,
		b0 = 1e-6,
		b1 = 1e-7,
		e1 = 0.1,
		vx_min = 1e6,
		vx_max = 1e7,
		vy_min = -1e7,
		vy_max = 1e7,
		vz_min = -1e7,
		vz_max = 1e7;
	size_t particles = 1, seed = 1;

	boost::program_options::options_description
		options("Usage: program_name [options], where options are");

	options.add_options()
		("help", "Print this help message")
		("quiet", "Don't save additional information")
		("outfile",
			boost::program_options::value<std::string>(&outfile_name)
				->default_value(outfile_name),
			"Append results to file arg_N, where N is MPI rank)")
		("time-start",
			boost::program_options::value<double>(&time_start)
				->default_value(time_start),
			"Start time of simulation (s)")
		("time-length",
			boost::program_options::value<double>(&time_length)
				->default_value(time_length),
			"Length of simulation in units of initial electron gyro period")
		("length",
			boost::program_options::value<double>(&length)
				->default_value(length),
			"Length of simulation box in units of smallest EM wavelength")
		("frequency",
			boost::program_options::value<double>(&em_frequency)
				->default_value(em_frequency),
			"Frequency of EM waves")
		("number-density-mid",
			boost::program_options::value<double>(&density_mid)
				->default_value(density_mid),
			"Electron number density in middle of box")
		("number-density-ends",
			boost::program_options::value<double>(&density_ends)
				->default_value(density_ends),
			"Electron number density at ends of box")
		("B0",
			boost::program_options::value<double>(&b0)
				->default_value(b0),
			"Magnetic field parallel to simulation box")
		("B1",
			boost::program_options::value<double>(&b1)
				->default_value(b1),
			"Magnetic field magnitude perpendicular to simulation box")
		("E1",
			boost::program_options::value<double>(&e1)
				->default_value(e1),
			"Electric field magnitude perpendicular to simulation box")
		("particles",
			boost::program_options::value<size_t>(&particles)
				->default_value(particles),
			"Approximate total number of particles to propagate (minimum is 1 / MPI rank)")
		("vx-min",
			boost::program_options::value<double>(&vx_min)
				->default_value(vx_min),
			"Minimum bulk velocity x component of electrons")
		("vx-max",
			boost::program_options::value<double>(&vx_max)
				->default_value(vx_max),
			"Maximum bulk velocity x component of electrons")
		("vy-min",
			boost::program_options::value<double>(&vy_min)
				->default_value(vy_min),
			"Minimum bulk velocity y component of electrons")
		("vy-max",
			boost::program_options::value<double>(&vy_max)
				->default_value(vy_max),
			"Maximum bulk velocity y component of electrons")
		("vz-min",
			boost::program_options::value<double>(&vz_min)
				->default_value(vz_min),
			"Minimum bulk velocity z component of electrons")
		("vz-max",
			boost::program_options::value<double>(&vz_max)
				->default_value(vz_max),
			"Maximum bulk velocity z component of electrons")
		("save-dt",
			boost::program_options::value<double>(&save_dt)
				->default_value(save_dt),
			"Save results at least every arg seconds, 0 saves "
			"initial and final states, < 0 doesn't save")
		("save-dx",
			boost::program_options::value<double>(&save_dx)
				->default_value(save_dx),
			"Save results at least every arg lengths in units of minimum EM wavelength")
		("time-step-factor",
			boost::program_options::value<double>(&time_step_factor)
				->default_value(time_step_factor),
			"Multiply maximum allowed time step with factor arg")
		("seed",
			boost::program_options::value<size_t>(&seed)
				->default_value(seed),
			"Use arg (+ MPI rank) as seed for random number generator");

	boost::program_options::variables_map option_variables;
	try {
		boost::program_options::store(
			boost::program_options::command_line_parser(argc, argv)
				.options(options)
				.run(),
			option_variables
		);
	} catch (std::exception& e) {
		std::cerr <<  __FILE__ << "(" << __LINE__ << "): "
			<< "Couldn't parse command line options: " << e.what()
			<< std::endl;
		MPI_Finalize();
		return EXIT_FAILURE;
	}
	boost::program_options::notify(option_variables);

	if (option_variables.count("help") > 0) {
		if (rank == 0) {
			cout << options << endl;
		}
		MPI_Finalize();
		return EXIT_SUCCESS;
	}

	if (option_variables.count("quiet") > 0) {
		verbose = false;
	}

	seed += rank;
	if (particles > comm_size) {
		particles /= comm_size;
	} else {
		particles = 1;
	}
	outfile_name += "_" + boost::lexical_cast<std::string>(rank);

	std::ofstream outfile(outfile_name, std::ios::app | std::ios::ate);
	if (not outfile.good()) {
		cerr << "Couldn't open file " << outfile_name << " for writing." << endl;
		return EXIT_FAILURE;
	}

	// derived parameters
	const double
		e_gyro_freq
			= abs(electron_charge)
			* sqrt(b0*b0 + 2*b1*b1)
			/ (2 * M_PI * electron_mass),
		e_plasma_freq_mid
			= sqrt(
				density_mid
				* pow(electron_charge, 2)
				/ (4 * M_PI * M_PI * electron_mass * vacuum_permittivity)
			),
		e_plasma_freq_ends
			= sqrt(
				density_ends
				* pow(electron_charge, 2)
				/ (4 * M_PI * M_PI * electron_mass * vacuum_permittivity)
			),
		wavelength_mid
			= light_speed
			/ em_frequency
			/ sqrt(
				1
				+ pow(e_plasma_freq_mid, 2)
					/ em_frequency
					/ (e_gyro_freq - em_frequency)
			),
		wavelength_ends
			= light_speed
			/ em_frequency
			/ sqrt(
				1
				+ pow(e_plasma_freq_ends, 2)
					/ em_frequency
					/ (e_gyro_freq - em_frequency)
			);

	length *= min(wavelength_mid, wavelength_ends);
	save_dx *= min(wavelength_mid, wavelength_ends);
	time_length /= e_gyro_freq;

	const double min_wavelength = std::min(wavelength_mid, wavelength_ends);

	const Fields fields(
		length,
		em_frequency,
		e1,
		b0,
		b1,
		wavelength_mid,
		wavelength_ends
	);

	if (verbose) {
		outfile
			<< "# Parameters\n# Start time (s): " << time_start
			<< "\n# Simulation length: " << time_length << " s, "
			<< time_length * e_gyro_freq << " electron gyro periods (g_e)"
			<< "\n# Simulation box length (m): " << length
			<< "\n# Number of particles: " << particles
			<< "\n# Minimum electron bulk velocity (m/s): "
			<< vx_min << " " << vy_min << " " << vz_min
			<< "\n# Maximum electron bulk velocity (m/s): "
			<< vx_max << " " << vy_max << " " << vz_max
			<< "\n# EM wave frequency (f_w): " << em_frequency
			<< " Hz\n# Electron number density in middle (1/m^3): " << density_mid
			<< "\n# Electron number density at ends (1/m^3): " << density_ends
			<< "\n# B0 (T): " << b0
			<< "\n# B1 (T): " << b1
			<< "\n# E1 (V/m): " << e1
			<< "\n# Random seed: " << seed
			<< "\n# Derived parameters\n# Electron gyro frequency (f_g): "
			<< e_gyro_freq << " Hz, " << e_gyro_freq / em_frequency << " f_w"
			<< "\n# Electron plasma frequency in middle: " << e_plasma_freq_mid
			<< " Hz, " << e_plasma_freq_mid / e_gyro_freq << " f_g"
			<< "\n# Electron plasma frequency at ends: " << e_plasma_freq_ends
			<< " Hz, " << e_plasma_freq_ends / e_gyro_freq << " f_g"
			<< "\n# EM wavelength in middle: " << wavelength_mid << " m, "
			<< wavelength_mid / length << " simulation boxes"
			<< "\n# EM wavelength at ends: " << wavelength_ends << " m, "
			<< wavelength_ends / length << " simulation boxes"
			<< "\n# EM phase speed in middle: " << wavelength_mid * em_frequency
			<< " m/s, " << wavelength_mid * em_frequency / light_speed << " c"
			<< "\n# EM phase speed at ends: " << wavelength_ends * em_frequency
			<< " m/s, " << wavelength_ends * em_frequency / light_speed << " c"
			<< "\n# Argument to E (0, E1*sin(arg), E1*sin(arg+pi/2) "
				"and B (B0, B1*sin(arg-pi/2), B1*sin(arg): " << fields.print()
			<< "\n# Each particle on separate line, format: r_x, r_y, ..., v_y, v_z, r_x, ..."
			<< std::endl;
	}

	Propagator propagator(electron_c2m, fields);
	runge_kutta_fehlberg78<state_t> stepper;

	std::mt19937_64 random_source;
	random_source.seed(seed);
	std::uniform_real_distribution<>
		vx_gen(vx_min, vx_max),
		vy_gen(vy_min, vy_max),
		vz_gen(vz_min, vz_max);

	for (size_t particle = 0; particle < particles; particle++) {

		Eigen::Vector3d initial_vel{
			vx_gen(random_source),
			vy_gen(random_source),
			vz_gen(random_source)
		};

		while (initial_vel.squaredNorm() >= light_speed2) {
			initial_vel = {
				vx_gen(random_source),
				vy_gen(random_source),
				vz_gen(random_source)
			};
		}

		state_t state{
			Eigen::Vector3d{0, 0, 0},
			initial_vel
		};

		if (initial_vel[0] < 0) {
			state[0][0] = length;
		}

		std::vector<double> data_to_save;
		if (save_dt >= 0 or save_dx >= 0) {
			data_to_save.push_back(state[0][0]);
			data_to_save.push_back(state[0][1]);
			data_to_save.push_back(state[0][2]);
			data_to_save.push_back(state[1][0]);
			data_to_save.push_back(state[1][1]);
			data_to_save.push_back(state[1][2]);
		}

		double
			simulation_time = 0,
			dt = 0,
			next_save_t = time_start;

		if (save_dt < 0) {
			next_save_t += 2 * time_length;
		}

		if (save_dx < 0) {
			save_dx = 2 * length;
		} else if (save_dx == 0) {
			save_dx = length;
		}

		double previous_save_x = state[0][0];
		while (true) {

			if (
				simulation_time >= next_save_t
				or abs(previous_save_x - state[0][0]) >= save_dx
			) {
				data_to_save.push_back(state[0][0]);
				data_to_save.push_back(state[0][1]);
				data_to_save.push_back(state[0][2]);
				data_to_save.push_back(state[1][0]);
				data_to_save.push_back(state[1][1]);
				data_to_save.push_back(state[1][2]);

				if (save_dt == 0) {
					next_save_t += time_length;
				} else {
					next_save_t += save_dt;
				}
				previous_save_x = state[0][0];
			}

			double max_time_step = std::numeric_limits<double>::max();
			const auto EB = fields(simulation_time, state[0][0]);
			max_time_step = std::min(
				max_time_step,
				pamhd::particle::get_minmax_step(
					1.0,
					1.0 / 32.0,
					min_wavelength / 8,
					electron_c2m,
					state[1],
					EB[0],
					EB[1]
				).second
			);
			max_time_step = std::min(max_time_step, 1 / em_frequency / 8);
			dt = max_time_step * time_step_factor;

			stepper.do_step(propagator, state, simulation_time, dt);
			if (state[0][0] < -0.01 * length or state[0][0] > 1.01 * length) {
				break;
			}

			simulation_time += dt;

			if (simulation_time > time_start + time_length) {
				std::cerr
					<< "Simulation time exceeded for particle with initial velocity "
					<< initial_vel[0] << " "
					<< initial_vel[1] << " "
					<< initial_vel[2] << ", current velocity relative to initial "
					<< state[1].norm() / initial_vel.norm() << " and current position "
					<< state[0][0] << " "
					<< state[0][1] << " "
					<< state[0][2]
					<< std::endl;
				break;
			}
		}

		if (save_dt >= 0 or save_dx >= 0) {
			data_to_save.push_back(state[0][0]);
			data_to_save.push_back(state[0][1]);
			data_to_save.push_back(state[0][2]);
			data_to_save.push_back(state[1][0]);
			data_to_save.push_back(state[1][1]);
			data_to_save.push_back(state[1][2]);
		}

		for (const auto& item: data_to_save) {
			outfile << item << " ";
		}
		outfile << "\n";
	}
	if (not outfile.good()) {
		cerr << "Writing results to " << outfile_name << " probably failed." << endl;
	}
	outfile.close();

	MPI_Finalize();

	return EXIT_SUCCESS;
}
