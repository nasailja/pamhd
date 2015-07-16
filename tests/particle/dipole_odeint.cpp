/*
Test for odeint particle propagator in a dipole magnetic field.

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


Propagation test after http://dx.doi.org/10.1029/2005JA011382
*/

#include "cmath"
#include "cstdlib"
#include "iomanip"
#include "iostream"

#include "boost/numeric/odeint.hpp"
#include "Eigen/Core"
#include "Eigen/Geometry"

#include "particle/solve.hpp"


using namespace std;
using namespace boost::numeric::odeint;
using namespace Eigen;
using namespace pamhd::particle;


constexpr double
	proton_mass = 1.672621777e-27,
	elementary_charge = 1.602176565e-19,
	proton_charge_mass_ratio = elementary_charge / proton_mass,
	earth_radius = 6.371e6;

const Eigen::Vector3d dipole_moment{0, 0, -8e22};


// state[0] = particle position, state[1] = velocity
using state_t = std::array<Vector3d, 2>;


Eigen::Vector3d get_earth_dipole(const Eigen::Vector3d& position)
{
	const double distance = position.norm();
	const Eigen::Vector3d
		direction = position / distance,
		projected_dip_mom = dipole_moment.dot(direction) * direction;

	return 1e-7 * (3 * projected_dip_mom - dipole_moment) / std::pow(distance, 3);
}

struct Particle_Propagator
{
	const double c2m;

	Particle_Propagator(const double& given_c2m) : c2m(given_c2m) {}

	void operator()(const state_t& state, state_t& change, const double)
	{
		change[0] = state[1];
		change[1] = this->c2m * state[1].cross(get_earth_dipole(state[0]));
	};
};


/*!
Returns propagated state and max relative errors in
kinetic energy and magnetic moment of particle.
*/
template<class Stepper> std::tuple<state_t, double, double> trace_trajectory(
	state_t state,
	const double time_step
) {
	using std::fabs;
	using std::max;

	constexpr double propagation_time = 424.1; // seconds, 1/100 of doi above

	Particle_Propagator propagator(proton_charge_mass_ratio);

	const double
		initial_energy = 0.5 * proton_mass * state[1].squaredNorm(),
		initial_magnetic_moment
			= proton_mass
			* std::pow(state[1][1], 2)
			/ (2 * get_earth_dipole(state[0]).norm());

	Stepper stepper;
	double nrj_error = 0, mom_error = 0;
	for (double time = 0; time < propagation_time; time += time_step) {
		stepper.do_step(propagator, state, time, time_step);

		const Eigen::Vector3d
			v(state[1]),
			b(get_earth_dipole(state[0])),
			v_perp_b(v - (v.dot(b) / b.squaredNorm()) * b);

		const double
			current_energy
				= 0.5 * proton_mass * state[1].squaredNorm(),
			current_magnetic_moment
				= proton_mass * v_perp_b.squaredNorm() / (2 * b.norm());

		nrj_error = max(
			nrj_error,
			fabs(current_energy - initial_energy) / initial_energy
		);
		mom_error = max(
			mom_error,
			fabs(current_magnetic_moment - initial_magnetic_moment)
				/ initial_magnetic_moment
		);
	}

	return std::make_tuple(state, nrj_error, mom_error);
}


/*
The following return true on success and known failures, false otherwise.
*/

bool check_result(
	modified_midpoint<state_t> stepper,
	const double initial_energy,
	const double initial_angle,
	const size_t step_divisor,
	const double order_of_accuracy,
	const double nrj_error,
	const double mom_error,
	const Eigen::Vector3d& position
) {
	// particle escapes
	if (step_divisor <= 8) {
		return true;
	}

	if (initial_energy == 1e4) {

		if (initial_angle == 90) {

			if (order_of_accuracy < 4) {
				return false;
			}
			if (step_divisor >= 128 and nrj_error > 1e-3) {
				return false;
			}
			if (step_divisor >= 128 and mom_error > 2e-2) {
				return false;
			}

			return true;

		} else {

			if (step_divisor == 16) {
				return true;
			}
			if (order_of_accuracy < 4) {
				return false;
			}
			if (step_divisor >= 256 and nrj_error > 1e-3) {
				return false;
			}
			if (step_divisor >= 512 and mom_error > 4e-2) {
				return false;
			}

			return true;
		}

	} else {

		if (initial_angle == 90) {

			if (order_of_accuracy < 4) {
				return false;
			}
			if (step_divisor >= 64 and nrj_error > 1e-3) {
				return false;
			}
			if (step_divisor >= 64 and mom_error > 0.6) {
				return false;
			}

			return true;

		} else {

			if (step_divisor <= 32) {
				return true;
			}
			if (order_of_accuracy < 4) {
				return false;
			}
			if (step_divisor >= 128 and nrj_error > 1e-3) {
				return false;
			}
			if (step_divisor >= 128 and mom_error > 2.1) {
				return false;
			}

			return true;
		}
	}
}

bool check_result(
	runge_kutta4<state_t> stepper,
	const double initial_energy,
	const double initial_angle,
	const size_t step_divisor,
	const double order_of_accuracy,
	const double nrj_error,
	const double mom_error,
	const Eigen::Vector3d& position
) {
	// particle escapes
	if (step_divisor <= 1) {
		return true;
	}

	if (initial_energy == 1e4) {

		if (initial_angle == 90) {

			if (step_divisor >= 32 and order_of_accuracy < 4) {
				return false;
			}
			if (step_divisor >= 128 and nrj_error > 1e-3) {
				return false;
			}
			if (step_divisor >= 64 and mom_error > 2e-2) {
				return false;
			}

			return true;

		} else {

			if (step_divisor <= 32) {
				return true;
			}
			if (order_of_accuracy < 4) {
				return false;
			}
			if (step_divisor >= 256 and nrj_error > 1e-3) {
				return false;
			}
			if (step_divisor >= 128 and mom_error > 4e-2) {
				return false;
			}

			return true;
		}

	} else {

		if (initial_angle == 90) {

			if (step_divisor >= 32 and order_of_accuracy < 4) {
				return false;
			}
			if (step_divisor >= 64 and nrj_error > 1e-3) {
				return false;
			}
			if (step_divisor >= 32 and mom_error > 0.6) {
				return false;
			}

			return true;

		} else {

			if (step_divisor <= 16) {
				return true;
			}
			if (order_of_accuracy < 4.5) {
				return false;
			}
			if (step_divisor >= 128 and nrj_error > 1e-3) {
				return false;
			}
			if (step_divisor >= 32 and mom_error > 2.1) {
				return false;
			}

			return true;
		}

	}
}

bool check_result(
	runge_kutta_cash_karp54<state_t> stepper,
	const double initial_energy,
	const double initial_angle,
	const size_t step_divisor,
	const double order_of_accuracy,
	const double nrj_error,
	const double mom_error,
	const Eigen::Vector3d& position
) {
	// particle escapes
	if (step_divisor <= 4) {
		return true;
	}

	if (initial_energy == 1e4) {

		if (initial_angle == 90) {

			if (order_of_accuracy < 4.5) {
				return false;
			}
			if (step_divisor >= 32 and nrj_error > 1e-3) {
				return false;
			}
			if (step_divisor >= 32 and mom_error > 2e-2) {
				return false;
			}

			return true;

		} else {

			if (step_divisor <= 8) {
				return true;
			}
			if (order_of_accuracy < 4.5) {
				return false;
			}
			if (step_divisor >= 128 and nrj_error > 1e-3) {
				return false;
			}
			if (step_divisor >= 128 and mom_error > 5e-2) {
				return false;
			}

			return true;
		}

	} else {

		if (initial_angle == 90) {

			if (order_of_accuracy < 4.5) {
				return false;
			}
			if (step_divisor >= 32 and nrj_error > 1e-3) {
				return false;
			}
			if (step_divisor >= 16 and mom_error > 0.6) {
				return false;
			}

			return true;

		} else {

			if (step_divisor <= 8) {
				return true;
			}
			if (order_of_accuracy < 4.5) {
				return false;
			}
			if (step_divisor >= 64 and nrj_error > 1e-3) {
				return false;
			}
			if (step_divisor >= 64 and mom_error > 2.1) {
				return false;
			}

			return true;
		}
	}
}

bool check_result(
	runge_kutta_dopri5<state_t> stepper,
	const double initial_energy,
	const double initial_angle,
	const size_t step_divisor,
	const double order_of_accuracy,
	const double nrj_error,
	const double mom_error,
	const Eigen::Vector3d& position
) {
	// particle escapes
	if (step_divisor <= 1) {
		return true;
	}

	if (initial_energy == 1e4) {

		if (initial_angle == 90) {

			if (step_divisor >= 64 and order_of_accuracy < 4.5) {
				return false;
			}
			if (step_divisor >= 64 and nrj_error > 1e-3) {
				return false;
			}
			if (step_divisor >= 32 and mom_error > 2e-2) {
				return false;
			}

			return true;

		} else {

			if (step_divisor <= 8) {
				return true;
			}
			if (step_divisor >= 128 and order_of_accuracy < 4.5) {
				return false;
			}
			if (step_divisor >= 128 and nrj_error > 1e-3) {
				return false;
			}
			if (step_divisor >= 64 and mom_error > 4e-2) {
				return false;
			}

			return true;
		}

	} else {

		if (initial_angle == 90) {

			if (step_divisor >= 64 and order_of_accuracy < 4.5) {
				return false;
			}
			if (step_divisor >= 32 and nrj_error > 1e-3) {
				return false;
			}
			if (step_divisor >= 4 and mom_error > 0.6) {
				return false;
			}

			return true;

		} else {

			if (step_divisor <= 8) {
				return true;
			}
			if (step_divisor >= 64 and order_of_accuracy < 4.5) {
				return false;
			}
			if (step_divisor >= 128 and nrj_error > 1e-3) {
				return false;
			}
			if (step_divisor >= 16 and mom_error > 2.1) {
				return false;
			}

			return true;
		}
	}
}

bool check_result(
	runge_kutta_fehlberg78<state_t> stepper,
	const double initial_energy,
	const double initial_angle,
	const size_t step_divisor,
	const double order_of_accuracy,
	const double nrj_error,
	const double mom_error,
	const Eigen::Vector3d& position
) {
	// particle escapes
	if (step_divisor <= 2) {
		return true;
	}

	if (position.norm() > 7 * earth_radius) {
		if (initial_angle == 30) {
			if (initial_energy == 1e4) {
				if (step_divisor <= 8) {
					return true;
				}
			} else if (step_divisor <= 4) {
				return true;
			}
		}

		return false;
	}

	if (initial_energy == 1e4) {

		if (initial_angle == 90) {

			if (step_divisor < 256 and order_of_accuracy < 6) {
				return false;
			}
			if (step_divisor >= 8 and nrj_error > 1e-3) {
				return false;
			}
			if (step_divisor >= 8 and mom_error > 2e-2) {
				return false;
			}

			return true;

		} else {

			if (step_divisor <= 8) {
				return true;
			}
			if (step_divisor < 512 and order_of_accuracy < 7) {
				return false;
			}
			if (step_divisor >= 32 and nrj_error > 1e-3) {
				return false;
			}
			if (step_divisor >= 32 and mom_error > 4e-2) {
				return false;
			}

			return true;
		}

	} else {

		if (initial_angle == 90) {

			if (step_divisor < 128 and order_of_accuracy < 8) {
				return false;
			}
			if (step_divisor >= 8 and nrj_error > 1e-3) {
				return false;
			}
			if (step_divisor >= 4 and mom_error > 0.6) {
				return false;
			}

			return true;

		} else {

			if (step_divisor <= 4) {
				return true;
			}
			if (step_divisor < 256 and order_of_accuracy < 4.5) {
				return false;
			}
			if (step_divisor >= 16 and nrj_error > 1e-3) {
				return false;
			}
			if (step_divisor >= 8 and mom_error > 2.1) {
				return false;
			}

			return true;
		}

	}
}


template<class Stepper> bool test_stepper()
{
	const Eigen::Vector3d
		guiding_center_start{5 * earth_radius, 0, 0},
		field_at_start{get_earth_dipole(guiding_center_start)};

	const double gyroperiod
		= 2 * M_PI / fabs(proton_charge_mass_ratio) / field_at_start.norm();

	/*cout <<
		"Energy (eV), pitch angle (deg), steps/gyroperiod: "
		"order of convergence, max relative error in energy\n";*/
	for (auto kinetic_energy: {1e4, 1e7}) { // in eV
	for (auto pitch_angle: {90.0, 30.0}) { // V from B, degrees

		const double particle_speed
			= sqrt(2 * kinetic_energy * proton_charge_mass_ratio);

		const Eigen::Vector3d
			initial_velocity{
				0,
				sin(pitch_angle / 180 * M_PI) * particle_speed,
				cos(pitch_angle / 180 * M_PI) * particle_speed
			};

		const double gyroradius
			= initial_velocity[1]
			/ proton_charge_mass_ratio
			/ field_at_start.norm();

		const Eigen::Vector3d
			offset_from_gc{gyroradius, 0, 0},
			initial_position{guiding_center_start + offset_from_gc};

		double old_nrj_error = 0;
		size_t old_step_divisor = 0;

		for (size_t step_divisor = 1; step_divisor <= (1 << 10); step_divisor *= 2) {

			state_t state{{initial_position, initial_velocity}};
			const double time_step = gyroperiod / step_divisor;
			double
				max_nrj_error = -1,
				max_mom_error = -1,
				order_of_accuracy = 0;

			std::tie(
				state,
				max_nrj_error,
				max_mom_error
			) = trace_trajectory<Stepper>(state, time_step);

			/*cout
				<< kinetic_energy << " "
				<< pitch_angle << " "
				<< step_divisor << ": "
				<< max_nrj_error << " "
				<< max_mom_error
				<< endl;*/

			if (old_step_divisor > 0) {
				order_of_accuracy
						= -log(max_nrj_error / old_nrj_error)
						/ log(double(step_divisor) / old_step_divisor);
			}

			if (
				not check_result(
					Stepper(),
					kinetic_energy,
					pitch_angle,
					step_divisor,
					order_of_accuracy,
					max_nrj_error,
					max_mom_error,
					state[0]
				)
			) {
				return false;
			}

			old_nrj_error = max_nrj_error;
			old_step_divisor = step_divisor;
		}
	}}

	return true;
}


int main()
{
	//test_stepper<euler<state_t>>(); doesn't work almost at all

	if (not test_stepper<modified_midpoint<state_t>>()) {
		std::cerr <<  __FILE__ << " (" << __LINE__ << "): "
			<< "Unexpected failure in modified_midpoint stepper."
			<< std::endl;
		return EXIT_FAILURE;
	}

	//test_stepper<rosenbrock4<state_t>>(); doesn't compile

	if (not test_stepper<runge_kutta4<state_t>>()) {
		std::cerr <<  __FILE__ << " (" << __LINE__ << "): "
			<< "Unexpected failure in runge_kutta4 stepper."
			<< std::endl;
		return EXIT_FAILURE;
	}

	if (not test_stepper<runge_kutta_cash_karp54<state_t>>()) {
		std::cerr <<  __FILE__ << " (" << __LINE__ << "): "
			<< "Unexpected failure in runge_kutta_cash_karp54 stepper."
			<< std::endl;
		return EXIT_FAILURE;
	}

	if (not test_stepper<runge_kutta_dopri5<state_t>>()) {
		std::cerr <<  __FILE__ << " (" << __LINE__ << "): "
			<< "Unexpected failure in runge_kutta_dopri5 stepper."
			<< std::endl;
		return EXIT_FAILURE;
	}

	if (not test_stepper<runge_kutta_fehlberg78<state_t>>()) {
		std::cerr <<  __FILE__ << " (" << __LINE__ << "): "
			<< "Unexpected failure in runge_kutta_fehlberg78 stepper."
			<< std::endl;
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}
