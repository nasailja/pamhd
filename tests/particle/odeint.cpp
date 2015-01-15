/*
Tests odeint version of particle propagator.

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
#include "cmath"
#include "cstdlib"
#include "iostream"

#include "boost/numeric/odeint.hpp"
#include "Eigen/Core" // must be included before gensimcell.hpp
#include "Eigen/Geometry"
#include "gensimcell.hpp"

#include "particle/variables.hpp"

using namespace std;
using namespace boost::numeric::odeint;
using namespace Eigen;
using namespace pamhd::particle;

using state_t = std::array<Vector3d, 2>;

struct propagate_particle
{
	const Vector3d E, B;
	const double c2m;

	propagate_particle(
		const Vector3d& given_E,
		const Vector3d& given_B,
		const double& given_c2m
	) :
		E(given_E),
		B(given_B),
		c2m(given_c2m)
	{}

	void operator()(const state_t& state, state_t& change, const double)
	{
		change[0] = state[1];
		change[1] = this->c2m * (this->E + state[1].cross(this->B));
	};
};

int main()
{
	using std::fabs;
	using std::log;

	const Position Pos{};
	const Velocity Vel{};
	const Mass Mas{};
	const Charge_Mass_Ratio C2M{};

	double time, dt, gyro_radius;
	Particle particle;

	particle[C2M] = 1;
	const Vector3d
		electric_field = Vector3d(1, 0, 0),
		magnetic_field = Vector3d(0, 3, 0);
	const double
		vel_drift = electric_field[0] / magnetic_field[1],
		gyro_period = 2 * M_PI / particle[C2M] / magnetic_field[1];

	//euler<state_t> stepper;
	//modified_midpoint<state_t> stepper;
	// doesn't compile rosenbrock4<state_t> stepper;
	runge_kutta4<state_t> stepper;
	//runge_kutta_cash_karp54<state_t> stepper;
	//runge_kutta_dopri5<state_t> stepper;
	//runge_kutta_fehlberg78<state_t> stepper;

	// test convergence of E.cross(B) drift
	double
		old_distance_rx = std::numeric_limits<double>::max(),
		old_distance_rz = std::numeric_limits<double>::max(),
		old_distance_vx = std::numeric_limits<double>::max(),
		old_distance_vz = std::numeric_limits<double>::max();
	size_t old_nr_of_steps = 0;
	for (size_t steps = 2; steps <= 4096; steps *= 2) {

		particle[Pos] = Vector3d(0, 0, 0);
		particle[Vel] = Vector3d(0, 0, 1);
		dt = gyro_period / steps;

		propagate_particle propagator(electric_field, magnetic_field, particle[C2M]);

		state_t state{{particle[Pos], particle[Vel]}};
		integrate_const(stepper, propagator, state, 0.0, gyro_period, dt);
		particle[Pos] = state[0];
		particle[Vel] = state[1];

		const double
			distance_rx = fabs(particle[Pos][0]),
			distance_ry = fabs(particle[Pos][1]),
			distance_rz = fabs(particle[Pos][2] - vel_drift * gyro_period),
			distance_vx = fabs(particle[Vel][0]),
			distance_vy = fabs(particle[Vel][1]),
			distance_vz = fabs(particle[Vel][2] - 1);

		if (distance_ry > 0) {
			std::cerr <<  __FILE__ << ":" << __LINE__
				<< " Non-zero y coordinate for particle with " << steps
				<< " steps / gyro period: " << particle[Pos][1]
				<< std::endl;
			abort();
		}

		if (distance_vy > 0) {
			std::cerr <<  __FILE__ << ":" << __LINE__
				<< " Non-zero y velocity for particle with " << steps
				<< " steps / gyro period: " << particle[Vel][1]
				<< std::endl;
			abort();
		}

		if (old_nr_of_steps > 0) {
			const double
				order_of_accuracy_rx
					= -log(distance_rx / old_distance_rx)
					/ log(double(steps) / old_nr_of_steps),
				order_of_accuracy_rz
					= -log(distance_rz / old_distance_rz)
					/ log(double(steps) / old_nr_of_steps),
				order_of_accuracy_vx
					= -log(distance_vx / old_distance_vx)
					/ log(double(steps) / old_nr_of_steps),
				order_of_accuracy_vz
					= -log(distance_vz / old_distance_vz)
					/ log(double(steps) / old_nr_of_steps);

			if (steps < 4096 and order_of_accuracy_rx < 4) {
				std::cerr <<  __FILE__ << ":" << __LINE__
					<< " Bad convergence for x position of particle with " << steps
					<< " steps / gyro period: " << order_of_accuracy_rx
					<< std::endl;
				abort();
			}
			if (steps > 16 and order_of_accuracy_rz < 3.5) {
				std::cerr <<  __FILE__ << ":" << __LINE__
					<< " Bad convergence for z position of particle with " << steps
					<< " steps / gyro period: " << order_of_accuracy_rx
					<< std::endl;
				abort();
			}

			if (steps > 8 and order_of_accuracy_vx < 3.5) {
				std::cerr <<  __FILE__ << ":" << __LINE__
					<< " Bad convergence for x velocity of particle with " << steps
					<< " steps / gyro period: " << order_of_accuracy_vx
					<< std::endl;
				abort();
			}
			if (steps < 4096 and order_of_accuracy_vz < 4) {
				std::cerr <<  __FILE__ << ":" << __LINE__
					<< " Bad convergence for z velocity of particle with " << steps
					<< " steps / gyro period: " << order_of_accuracy_vz
					<< std::endl;
				abort();
			}

		}

		old_distance_rx = distance_rx;
		old_distance_rz = distance_rz;
		old_distance_vx = distance_vx;
		old_distance_vz = distance_vz;
		old_nr_of_steps = steps;
	}

	return EXIT_SUCCESS;
}
