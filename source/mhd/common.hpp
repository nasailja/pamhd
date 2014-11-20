/*
Common MHD functions of PAMHD.

Copyright 2014 Ilja Honkonen
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice, this
  list of conditions and the following disclaimer in the documentation and/or
  other materials provided with the distribution.

* Neither the name of copyright holders nor the names of their contributors
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

#ifndef PAMHD_MHD_COMMON_HPP
#define PAMHD_MHD_COMMON_HPP


#include "boost/lexical_cast.hpp"
#include "cmath"
#include "stdexcept"
#include "string"


namespace pamhd {
namespace mhd {


/*!
Returns pressure in given MHD state.

Throws std::domain_error if given a state with non-positive mass density.

Returns negative pressure if total energy in given state
is smaller than kinetic + magnetic energies.

Example:
struct Mass_Density {using data_type = double;};
...
using MHD = gensimcell::Cell<..., Mass_Density, ...>;
const MHD state;
double p = get_pressure<MHD, Mass_Density, ...>(state);
*/
template <
	class MHD_T,
	class Mass_Density_T,
	class Momentum_Density_T,
	class Total_Energy_Density_T,
	class Magnetic_Field_T
> double get_pressure(
	const MHD_T& state,
	const double adiabatic_index,
	const double vacuum_permeability
) throw(std::domain_error) {
	const auto rho = state[Mass_Density_T()];
	if (rho <= 0) {
		throw std::domain_error(
			std::string("Non-positive mass density given to ")
			+ __func__
			+ std::string(": ")
			+ boost::lexical_cast<std::string>(rho)
		);
	}

	const auto
		&mom = state[Momentum_Density_T()],
		&mag = state[Magnetic_Field_T()];
	const auto
		nrj = state[Total_Energy_Density_T()],
		kinetic_energy = 0.5 * mom.squaredNorm() / rho,
		magnetic_energy = 0.5 * mag.squaredNorm() / vacuum_permeability;

	return (nrj - kinetic_energy - magnetic_energy) * (adiabatic_index - 1);
}


/*!
Returns flux from given state into positive x direction.

Throws std::domain_error if given a state with non-positive mass density.
*/
template <
	class MHD_T,
	class Mass_Density_T,
	class Momentum_Density_T,
	class Total_Energy_Density_T,
	class Magnetic_Field_T
> MHD_T get_flux(
	const MHD_T& state,
	const double adiabatic_index,
	const double vacuum_permeability
) throw(std::domain_error) {
	// shorthand notation for variables
	const Mass_Density_T Mas{};
	const Momentum_Density_T Mom{};
	const Total_Energy_Density_T Nrj{};
	const Magnetic_Field_T Mag{};

	const auto rho = state[Mas];
	if (rho <= 0) {
		throw std::domain_error(
			std::string("Non-positive mass density given to ")
			+ __func__
			+ std::string(": ")
			+ boost::lexical_cast<std::string>(rho)
		);
	}

	const auto
		inv_permeability = 1.0 / vacuum_permeability,
		pressure_magnetic = inv_permeability * state[Mag].squaredNorm() / 2,
		pressure_thermal
			= get_pressure<
				MHD_T,
				Mass_Density_T,
				Momentum_Density_T,
				Total_Energy_Density_T,
				Magnetic_Field_T
			>(state, adiabatic_index, vacuum_permeability);

	const auto velocity = state[Mom] / rho;

	MHD_T flux;

	flux[Mas] = state[Mom][0];

	flux[Mom]
		= state[Mas] * velocity[0] * velocity
		- state[Mag][0] * state[Mag] * inv_permeability;
	flux[Mom][0] += pressure_thermal + pressure_magnetic;

	flux[Nrj]
		= velocity[0] * (state[Nrj] + pressure_thermal + pressure_magnetic)
		- state[Mag][0] * velocity.dot(state[Mag]) * inv_permeability;

	flux[Mag] = velocity[0] * state[Mag] - state[Mag][0] * velocity;
	flux[Mag][0] = 0;

	return flux;
}


/*!
Returns total energy density in given primitive MHD state.

Throws std::domain_error if given a state with non-positive mass density.

\see get_pressure()
*/
template <
	class MHD_Primitive_T,
	class Mass_Density_T,
	class Velocity_T,
	class Pressure_T,
	class Magnetic_Field_T
> double get_total_energy_density(
	const MHD_Primitive_T& state,
	const double adiabatic_index,
	const double vacuum_permeability
) throw(std::domain_error) {
	if (state[Mass_Density_T()] <= 0) {
		throw std::domain_error(
			std::string("Non-positive mass density given to ")
			+ __func__
			+ std::string(": ")
			+ boost::lexical_cast<std::string>(state[Mass_Density_T()])
		);
	}

	const auto
		kinetic_energy
			= 0.5 * state[Mass_Density_T()]
			* state[Velocity_T()].squaredNorm(),
		magnetic_energy
			= 0.5 * state[Magnetic_Field_T()].squaredNorm()
			/ vacuum_permeability;

	return
		state[Pressure_T()] / (adiabatic_index - 1)
		+ kinetic_energy
		+ magnetic_energy;
}


/*!
Returns speed of sound wave in given MHD state.

Throws std::domain_error if given a state with
non-positive mass density or pressure.

\see get_pressure() for help with the template arguments.
*/
template <
	class MHD_T,
	class Mass_Density_T,
	class Momentum_Density_T,
	class Total_Energy_Density_T,
	class Magnetic_Field_T
> double get_sound_speed(
	const MHD_T& state,
	const double adiabatic_index,
	const double vacuum_permeability
) throw(std::domain_error) {
	const auto rho = state[Mass_Density_T()];
	if (rho <= 0) {
		throw std::domain_error(
			std::string("Non-positive mass density given to ")
			+ __func__
			+ std::string(": ")
			+ boost::lexical_cast<std::string>(rho)
		);
	}

	const auto pressure
		= get_pressure<
			MHD_T,
			Mass_Density_T,
			Momentum_Density_T,
			Total_Energy_Density_T,
			Magnetic_Field_T
		>(state, adiabatic_index, vacuum_permeability);
	if (pressure <= 0) {
		throw std::domain_error(
			std::string("Non-positive pressure given to ")
			+ __func__
			+ std::string(": ")
			+ boost::lexical_cast<std::string>(pressure)
		);
	}

	return std::sqrt(adiabatic_index * pressure / rho);
}


/*!
Returns speed of Alfvén wave in given MHD state.

Throws std::domain_error if given a state with non-positive mass density.

\see get_pressure() for help with the template arguments.
*/
template <
	class MHD_T,
	class Mass_Density_T,
	class Magnetic_Field_T
> double get_alfven_speed(
	const MHD_T& state,
	const double vacuum_permeability
) throw(std::domain_error) {
	const auto rho = state[Mass_Density_T()];
	if (rho <= 0) {
		throw std::domain_error(
			std::string("Non-positive mass density given to ")
			+ __func__
			+ std::string(": ")
			+ boost::lexical_cast<std::string>(rho)
		);
	}

	return state[Magnetic_Field_T()].norm() / std::sqrt(vacuum_permeability * rho);
}


/*!
Returns speed of fast magnetosonic wave in x dimension of given MHD state.

Returns a non-negative value.

\see get_pressure() for help with the template arguments.
*/
template <
	class MHD_T,
	class Mass_Density_T,
	class Momentum_Density_T,
	class Total_Energy_Density_T,
	class Magnetic_Field_T
> double get_fast_magnetosonic_speed(
	const MHD_T& state,
	const double adiabatic_index,
	const double vacuum_permeability
) {
	const auto
		B_mag = state[Magnetic_Field_T()].norm(),
		sound = get_sound_speed<
			MHD_T,
			Mass_Density_T,
			Momentum_Density_T,
			Total_Energy_Density_T,
			Magnetic_Field_T
		>(state, adiabatic_index, vacuum_permeability);

	if (B_mag == 0) {
		return sound;
	}

	const auto
		B_in_dim = state[Magnetic_Field_T()][0],
		sound2 = sound * sound,
		alfven2
			= std::pow(
				get_alfven_speed<
					MHD_T,
					Mass_Density_T,
					Magnetic_Field_T
				>(state, vacuum_permeability),
				2
			),
		squared_sum = sound2 + alfven2;

	return std::sqrt(
		0.5 * (
			squared_sum
			+ std::sqrt(
				std::pow(squared_sum, 2)
				- 4 * sound2 * alfven2 * std::pow(B_in_dim / B_mag, 2)
			)
		)
	);
}


/*!
Throws std::domain_error if given a state with non-positive mass density.
*/
template <
	class MHD_Primitive_T,
	class Velocity_T,
	class Pressure_T,
	class MHD_Conservative_T,
	class Mass_Density_T,
	class Momentum_Density_T,
	class Total_Energy_Density_T,
	class Magnetic_Field_T
> MHD_Primitive_T get_primitive(
	const MHD_Conservative_T& state,
	const double adiabatic_index,
	const double vacuum_permeability
) throw(std::domain_error) {
	if (state[Mass_Density_T()] <= 0) {
		throw std::domain_error(
			std::string("Non-positive mass density given to ")
			+ __func__
			+ std::string(": ")
			+ boost::lexical_cast<std::string>(state[Mass_Density_T()])
		);
	}

	MHD_Primitive_T ret_val;

	ret_val[Mass_Density_T()] = state[Mass_Density_T()];
	ret_val[Velocity_T()] = state[Momentum_Density_T()] / state[Mass_Density_T()];
	ret_val[Pressure_T()]
		= get_pressure<
			MHD_Conservative_T,
			Mass_Density_T,
			Momentum_Density_T,
			Total_Energy_Density_T,
			Magnetic_Field_T
		>(state, adiabatic_index, vacuum_permeability);
	ret_val[Magnetic_Field_T()] = state[Magnetic_Field_T()];

	return ret_val;
}


template <
	class MHD_Conservative_T,
	class MHD_Primitive_T,
	class Momentum_Density_T,
	class Total_Energy_Density_T,
	class Mass_Density_T,
	class Velocity_T,
	class Pressure_T,
	class Magnetic_Field_T
> MHD_Conservative_T get_conservative(
	const MHD_Primitive_T& state,
	const double adiabatic_index,
	const double vacuum_permeability
) {
	MHD_Conservative_T ret_val;


	ret_val[Mass_Density_T()] = state[Mass_Density_T()];
	ret_val[Momentum_Density_T()] = state[Velocity_T()] * state[Mass_Density_T()];
	ret_val[Magnetic_Field_T()] = state[Magnetic_Field_T()];

	ret_val[Total_Energy_Density_T()] = get_total_energy_density<
		MHD_Primitive_T,
		Mass_Density_T,
		Velocity_T,
		Pressure_T,
		Magnetic_Field_T
	>(state, adiabatic_index, vacuum_permeability);

	return ret_val;
}


/*!
Rotates components of given vector for get_flux_* functions.

Positive direction rotates vector for get_flux_* (assume
two states have identical y and z coordinates), negative
rotates it back.
*/
template<class Vector_T> Vector_T get_rotated_vector(
	const Vector_T& v,
	const int direction
) {
	Vector_T ret_val;

	switch(direction) {
	case -3:
		ret_val[0] = v[1];
		ret_val[1] = v[2];
		ret_val[2] = v[0];
		return ret_val;
	case -2:
		ret_val[0] = v[2];
		ret_val[1] = v[0];
		ret_val[2] = v[1];
		return ret_val;
	case -1:
	case 1:
		return v;
	case 2:
		ret_val[0] = v[1];
		ret_val[1] = v[2];
		ret_val[2] = v[0];
		return ret_val;
	case 3:
		ret_val[0] = v[2];
		ret_val[1] = v[0];
		ret_val[2] = v[1];
		return ret_val;
	default:
		std::cerr <<  __FILE__ << "(" << __LINE__ << ") "
			<< "Invalid direction: " << direction
			<< std::endl;
		abort();
	}
}


/*!
Returns given MHD state with vector variables rotated.

Rotates vector variables so that given direction becomes
the first one.
Negative value reverses the operation, e.g.
state == get_rotated(get_rotated(state, 2), -2);

Direction starts from 1 and values of -1 and 1 do nothing.
*/
template <
	class MHD_T,
	class Mass_Density_T,
	class Momentum_Density_T,
	class Total_Energy_Density_T,
	class Magnetic_Field_T
> typename MHD_T::data_type get_rotated_state(
	const typename MHD_T::data_type& state,
	const int direction
) {
	typename MHD_T::data_type ret_val;

	ret_val[Mass_Density_T()] = state[Mass_Density_T()];
	ret_val[Total_Energy_Density_T()] = state[Total_Energy_Density_T()];
	ret_val[Momentum_Density_T()] = get_rotated_vector(state[Momentum_Density_T()], direction);
	ret_val[Magnetic_Field_T()] = get_rotated_vector(state[Magnetic_Field_T()], direction);

	return ret_val;
}


/*!
Positive flux adds to given cell multiplied with given factor.

Throws std::domain_error if new state has
non-positive mass density of pressure.
*/
template <
	class Cell_T,
	class MHD_T,
	class MHD_Flux_T,
	class Mass_Density_T,
	class Momentum_Density_T,
	class Total_Energy_Density_T,
	class Magnetic_Field_T
> void apply_fluxes(
	Cell_T& cell_data,
	const double factor,
	const double adiabatic_index,
	const double vacuum_permeability
) throw(std::domain_error) {
	auto& state = cell_data[MHD_T()];
	const auto& flux = cell_data[MHD_Flux_T()];

	state += flux * factor;

	if (state[Mass_Density_T()] <= 0) {
		throw std::domain_error(
			"New state has non-positive mass density: "
			+ boost::lexical_cast<std::string>(state[Mass_Density_T()])
		);
	}

	const auto pressure
		= get_pressure<
			typename MHD_T::data_type,
			Mass_Density_T,
			Momentum_Density_T,
			Total_Energy_Density_T,
			Magnetic_Field_T
		>(
			state,
			adiabatic_index,
			vacuum_permeability
		);
	if (pressure <= 0) {
		throw std::domain_error(
			"New state has non-positive pressure: "
			+ boost::lexical_cast<std::string>(pressure)
		);
	}
}


/*!
Zeroes fluxes of given variables.
*/
template <
	class Cell_T,
	class MHD_Flux_T,
	class Mass_Density_T,
	class Momentum_Density_T,
	class Total_Energy_Density_T,
	class Magnetic_Field_T
> void zero_fluxes(Cell_T& cell_data)
{
	auto& flux = cell_data[MHD_Flux_T()];
	flux[Mass_Density_T()]         =
	flux[Momentum_Density_T()][0]  =
	flux[Momentum_Density_T()][1]  =
	flux[Momentum_Density_T()][2]  =
	flux[Total_Energy_Density_T()] =
	flux[Magnetic_Field_T()][0]    =
	flux[Magnetic_Field_T()][1]    =
	flux[Magnetic_Field_T()][2]    = 0;
}


}} // namespaces

#endif // ifndef PAMHD_MHD_COMMON_HPP
