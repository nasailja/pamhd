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
#include "Eigen/Core" // must be included before gensimcell.hpp
#include "Eigen/Geometry"

#include "particle/variables.hpp"
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


//! returns propagated state and max relative error in kinetic energy
template<class Stepper> std::pair<state_t, double> trace_trajectory(
	state_t state,
	const double time_step
) {
	using std::fabs;
	using std::max;

	constexpr double propagation_time = 424.1; // seconds, 1/100 of doi above

	Particle_Propagator propagator(proton_charge_mass_ratio);

	const double initial_energy
		= 0.5 * proton_mass * state[1].squaredNorm();

	Stepper stepper;
	double error = 0;
	for (double time = 0; time < propagation_time; time += time_step) {
		stepper.do_step(propagator, state, time, time_step);
		//cout << state[0][0] / earth_radius << " " << state[0][1] / earth_radius << " " << state[0][2] / earth_radius << endl;

		const double current_energy
			= 0.5 * proton_mass * state[1].squaredNorm();

		error = max(error, fabs(current_energy - initial_energy) / initial_energy);
	}

	return std::make_pair(state, error);
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
	const double error,
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
			if (step_divisor >= 128 and error > 1e-3) {
				return false;
			} else {
				return true;
			}

		} else {

			if (step_divisor == 16) {
				return true;
			}
			if (order_of_accuracy < 4) {
				return false;
			}
			if (step_divisor >= 256 and error > 1e-3) {
				return false;
			} else {
				return true;
			}
		}

	} else {

		if (initial_angle == 90) {

			if (order_of_accuracy < 4) {
				return false;
			}
			if (step_divisor >= 64 and error > 1e-3) {
				return false;
			} else {
				return true;
			}

		} else {

			if (step_divisor <= 32) {
				return true;
			}
			if (order_of_accuracy < 4) {
				return false;
			}
			if (step_divisor >= 128 and error > 1e-3) {
				return false;
			} else {
				return true;
			}
		}

	}
}

bool check_result(
	runge_kutta4<state_t> stepper,
	const double initial_energy,
	const double initial_angle,
	const size_t step_divisor,
	const double order_of_accuracy,
	const double error,
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
			if (step_divisor >= 128 and error > 1e-3) {
				return false;
			} else {
				return true;
			}

		} else {

			// particle escapes
			if (step_divisor <= 32) {
				return true;
			}
			if (order_of_accuracy < 4) {
				return false;
			}
			if (step_divisor >= 256 and error > 1e-3) {
				return false;
			} else {
				return true;
			}
		}

	} else {

		if (initial_angle == 90) {

			if (step_divisor >= 32 and order_of_accuracy < 4) {
				return false;
			}
			if (step_divisor >= 64 and error > 1e-3) {
				return false;
			} else {
				return true;
			}

		} else {

			if (step_divisor <= 16) {
				return true;
			}
			if (order_of_accuracy < 4.5) {
				return false;
			}
			if (step_divisor >= 128 and error > 1e-3) {
				return false;
			} else {
				return true;
			}
		}

	}
}

bool check_result(
	runge_kutta_cash_karp54<state_t> stepper,
	const double initial_energy,
	const double initial_angle,
	const size_t step_divisor,
	const double order_of_accuracy,
	const double error,
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
			if (step_divisor >= 32 and error > 1e-3) {
				return false;
			} else {
				return true;
			}

		} else {

			if (step_divisor <= 8) {
				return true;
			}
			if (order_of_accuracy < 4.5) {
				return false;
			}
			if (step_divisor >= 128 and error > 1e-3) {
				return false;
			} else {
				return true;
			}
		}

	} else {

		if (initial_angle == 90) {

			if (order_of_accuracy < 4.5) {
				return false;
			}
			if (step_divisor >= 32 and error > 1e-3) {
				return false;
			} else {
				return true;
			}

		} else {

			if (step_divisor <= 8) {
				return true;
			}
			if (order_of_accuracy < 4.5) {
				return false;
			}
			if (step_divisor >= 64 and error > 1e-3) {
				return false;
			} else {
				return true;
			}
		}

	}
}

bool check_result(
	runge_kutta_dopri5<state_t> stepper,
	const double initial_energy,
	const double initial_angle,
	const size_t step_divisor,
	const double order_of_accuracy,
	const double error,
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
			if (step_divisor >= 64 and error > 1e-3) {
				return false;
			} else {
				return true;
			}

		} else {

			if (step_divisor <= 8) {
				return true;
			}
			if (step_divisor >= 128 and order_of_accuracy < 4.5) {
				return false;
			}
			if (step_divisor >= 128 and error > 1e-3) {
				return false;
			} else {
				return true;
			}
		}

	} else {

		if (initial_angle == 90) {

			if (step_divisor >= 64 and order_of_accuracy < 4.5) {
				return false;
			}
			if (step_divisor >= 32 and error > 1e-3) {
				return false;
			} else {
				return true;
			}

		} else {

			if (step_divisor <= 8) {
				return true;
			}
			if (step_divisor >= 64 and order_of_accuracy < 4.5) {
				return false;
			}
			if (step_divisor >= 128 and error > 1e-3) {
				return false;
			} else {
				return true;
			}
		}

	}
}

bool check_result(
	runge_kutta_fehlberg78<state_t> stepper,
	const double initial_energy,
	const double initial_angle,
	const size_t step_divisor,
	const double order_of_accuracy,
	const double error,
	const Eigen::Vector3d& position
) {
	// particle escapes
	if (step_divisor <= 2) {
		return true;
	}

	if (initial_energy == 1e4) {

		if (initial_angle == 90) {

			if (step_divisor < 256 and order_of_accuracy < 6) {
				return false;
			}
			if (step_divisor >= 8 and error > 1e-3) {
				return false;
			} else {
				return true;
			}

		} else {

			if (step_divisor <= 8) {
				return true;
			}
			if (step_divisor < 512 and order_of_accuracy < 7) {
				return false;
			}
			if (step_divisor >= 32 and error > 1e-3) {
				return false;
			} else {
				return true;
			}
		}

	} else {

		if (initial_angle == 90) {

			if (step_divisor < 128 and order_of_accuracy < 8) {
				return false;
			}
			if (step_divisor >= 8 and error > 1e-3) {
				return false;
			} else {
				return true;
			}

		} else {

			if (step_divisor <= 4) {
				return true;
			}
			if (step_divisor < 256 and order_of_accuracy < 4.5) {
				return false;
			}
			if (step_divisor >= 16 and error > 1e-3) {
				return false;
			} else {
				return true;
			}
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

		const Eigen::Vector3d initial_velocity{
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

		double old_max_error = 0;
		size_t old_step_divisor = 0;

		/// smaller step divisor makes particle escape
		/// convergence of energy error stops at largest step divisor
		for (size_t step_divisor = 1; step_divisor <= (1 << 10); step_divisor *= 2) {

			//cout << kinetic_energy << " " << pitch_angle << " " << step_divisor << ": ";

			state_t state{{initial_position, initial_velocity}};
			const double time_step = gyroperiod / step_divisor;
			double max_relative_error = -1;

			std::tie(
				state,
				max_relative_error
			) = trace_trajectory<Stepper>(state, time_step);

			// check that position is sane
			if (state[0].norm() > 7 * earth_radius) {
				/*std::cerr <<  __FILE__ << " (" << __LINE__ << "): "
					<< "Particle escaped from earth with step divisor " << step_divisor
					<< std::endl;*/
			} else if (old_step_divisor > 0) {
				const double order_of_accuracy
						= -log(max_relative_error / old_max_error)
						/ log(double(step_divisor) / old_step_divisor);

				//cout << order_of_accuracy << " " << max_relative_error << endl;

				if (
					not check_result(
						Stepper(),
						kinetic_energy,
						pitch_angle,
						step_divisor,
						order_of_accuracy,
						max_relative_error,
						state[0]
					)
				) {
					return false;
				}
			}

			old_max_error = max_relative_error;
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
			<< "modified_midpoint stepper failed."
			<< std::endl;
		return EXIT_FAILURE;
	}

	//test_stepper<rosenbrock4<state_t>>(); doesn't compile

	if (not test_stepper<runge_kutta4<state_t>>()) {
		std::cerr <<  __FILE__ << " (" << __LINE__ << "): "
			<< "runge_kutta4 stepper failed."
			<< std::endl;
		return EXIT_FAILURE;
	}

	if (not test_stepper<runge_kutta_cash_karp54<state_t>>()) {
		std::cerr <<  __FILE__ << " (" << __LINE__ << "): "
			<< "runge_kutta_cash_karp54 stepper failed."
			<< std::endl;
		return EXIT_FAILURE;
	}

	if (not test_stepper<runge_kutta_dopri5<state_t>>()) {
		std::cerr <<  __FILE__ << " (" << __LINE__ << "): "
			<< "runge_kutta_dopri5 stepper failed."
			<< std::endl;
		return EXIT_FAILURE;
	}

	if (not test_stepper<runge_kutta_fehlberg78<state_t>>()) {
		std::cerr <<  __FILE__ << " (" << __LINE__ << "): "
			<< "runge_kutta_fehlberg78 stepper failed."
			<< std::endl;
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}
