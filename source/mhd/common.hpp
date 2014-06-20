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

#include "cmath"

namespace pamhd {
namespace mhd {


/*!
Returns pressure in given MHD state.

Example:
struct Mass_Density {using data_type = double;}; ...
using MHD = gensimcell::Cell<Mass_Density, ...>;
const MHD mhd;
double p = get_pressure<MHD, Mass_Density, ...>(mhd);
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
) {
	const auto
		&mom = state[Momentum_Density_T()],
		&B = state[Magnetic_Field_T()];
	const auto
		rho = state[Mass_Density_T()],
		nrj = state[Total_Energy_Density_T()],
		kinetic_energy = 0.5 * mom.squaredNorm() / rho,
		magnetic_energy = 0.5 * B.squaredNorm() / vacuum_permeability;

	return (nrj - kinetic_energy - magnetic_energy) * (adiabatic_index - 1);
}


/*!
Returns total energy density in given primitive MHD state.

\see get_pressure() for help with the template arguments.
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
) {
	const auto
		&v = state[Velocity_T()],
		&B = state[Magnetic_Field_T()];
	const auto
		rho = state[Mass_Density_T()],
		p = state[Pressure_T()],
		kinetic_energy = 0.5 * rho * v.squaredNorm(),
		magnetic_energy = 0.5 * B.squaredNorm() / vacuum_permeability;

	return p / (adiabatic_index - 1) + kinetic_energy + magnetic_energy;
}


/*!
Returns speed of sound wave in given MHD state.

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
) {
	const auto
		pressure
			= get_pressure<
				MHD_T,
				Mass_Density_T,
				Momentum_Density_T,
				Total_Energy_Density_T,
				Magnetic_Field_T
			>(state, adiabatic_index, vacuum_permeability),
		rho = state[Mass_Density_T()];

	return std::sqrt(adiabatic_index * pressure / rho);
}


/*!
Returns speed of Alfvén wave in given MHD state.

\see get_pressure() for help with the template arguments.
*/
template <
	class MHD_T,
	class Mass_Density_T,
	class Magnetic_Field_T
> double get_alfven_speed(
	const MHD_T& state,
	const double vacuum_permeability
) {
	const auto
		rho = state[Mass_Density_T()],
		B_mag = state[Magnetic_Field_T()].norm();

	return B_mag / std::sqrt(vacuum_permeability * rho);
}


/*!
Returns speed of fast magnetosonic wave in x dimension of given MHD state.

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
		B_in_dim = state[Magnetic_Field_T()][0],
		B_mag = state[Magnetic_Field_T()].norm(),
		sound2
			= std::pow(
				get_sound_speed<
					MHD_T,
					Mass_Density_T,
					Momentum_Density_T,
					Total_Energy_Density_T,
					Magnetic_Field_T
				>(state, adiabatic_index, vacuum_permeability),
				2
			),
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
) {
	MHD_Primitive_T ret_val;

	ret_val[Mass_Density_T()] = state[Mass_Density_T()];
	ret_val[Magnetic_Field_T()] = state[Magnetic_Field_T()];

	ret_val[Pressure_T()] = get_pressure<
		MHD_Conservative_T,
		Mass_Density_T,
		Momentum_Density_T,
		Total_Energy_Density_T,
		Magnetic_Field_T
	>(state, adiabatic_index, vacuum_permeability);

	if (state[Mass_Density_T()] == 0) {
		ret_val[Velocity_T()] = 0;
	} else {
		ret_val[Velocity_T()] = state[Momentum_Density_T()] / state[Mass_Density_T()];
	}

	return ret_val;
}


template <
	class MHD_Conservative_T,
	class Momentum_Density_T,
	class Total_Energy_Density_T,
	class MHD_Primitive_T,
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
Positive flux adds to given cell multiplied with given factor.
Zeroes flux afterwards.
*/
template <
	class Cell_T,
	class MHD_T,
	class MHD_Flux_T,
	class Mass_Density_T,
	class Momentum_Density_T,
	class Total_Energy_Density_T,
	class Magnetic_Field_T
> void apply_flux(
	Cell_T& cell_data,
	const double factor
) {
	auto& state = cell_data[MHD_T()];
	auto& flux = cell_data[MHD_Flux_T()];

	state += flux * factor;

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
