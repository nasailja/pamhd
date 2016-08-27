/*
Particle-assisted version of hll_athena.hpp.
Copyright 2003 James M. Stone
Copyright 2003 Thomas A. Gardiner
Copyright 2003 Peter J. Teuben
Copyright 2003 John F. Hawley
Copyright 2014, 2015, 2016 Ilja Honkonen

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef PAMHD_MHD_N_HLL_ATHENA_HPP
#define PAMHD_MHD_N_HLL_ATHENA_HPP


#include "cmath"
#include "limits"
#include "tuple"

#include "boost/lexical_cast.hpp"
#include "gensimcell.hpp"

#include "mhd/common.hpp"


namespace pamhd {
namespace mhd {
namespace athena {


/*!
Multi-population version of get_flux_hll() in hll_athena.hpp.

Splits returned flux into contributions from state_neg and state_pos
to all populations based on their fraction of mass vs total mass.

Ignores background magnetic field.
*/
template <
	class MHD,
	class Vector,
	class Mass_Density_T,
	class Momentum_Density_T,
	class Total_Energy_Density_T,
	class Magnetic_Field_T
> std::tuple<MHD, MHD, double> get_flux_N_hll(
	const MHD& state_neg,
	const MHD& state_pos,
	const Vector& /*bg_face_magnetic_field*/,
	const double& area,
	const double& dt,
	const double& adiabatic_index,
	const double& vacuum_permeability
) {
	const Vector bg_face_magnetic_field{0, 0, 0};

	using std::isnormal;
	using std::isfinite;

	const Mass_Density_T Mas{};
	const Momentum_Density_T Mom{};
	const Total_Energy_Density_T Nrj{};
	const Magnetic_Field_T Mag{};

	const auto
		flow_v_neg(state_neg[Mom] / state_neg[Mas]),
		flow_v_pos(state_pos[Mom] / state_pos[Mas]);

	const auto
		pressure_thermal_neg
			= get_pressure(
				state_neg[Mas],
				state_neg[Mom],
				state_neg[Nrj],
				state_neg[Mag],
				adiabatic_index,
				vacuum_permeability
			),

		pressure_thermal_pos
			= get_pressure(
				state_pos[Mas],
				state_pos[Mom],
				state_pos[Nrj],
				state_pos[Mag],
				adiabatic_index,
				vacuum_permeability
			),

		pressure_magnetic_neg
			= state_neg[Mag].squaredNorm() / (2 * vacuum_permeability),

		pressure_magnetic_pos
			= state_pos[Mag].squaredNorm() / (2 * vacuum_permeability),

		fast_magnetosonic_neg
			= get_fast_magnetosonic_speed(
				state_neg[Mas],
				state_neg[Mom],
				state_neg[Nrj],
				state_neg[Mag],
				bg_face_magnetic_field,
				adiabatic_index,
				vacuum_permeability
			),

		fast_magnetosonic_pos
			= get_fast_magnetosonic_speed(
				state_pos[Mas],
				state_pos[Mom],
				state_pos[Nrj],
				state_pos[Mag],
				bg_face_magnetic_field,
				adiabatic_index,
				vacuum_permeability
			),

		max_signal = std::max(fast_magnetosonic_neg, fast_magnetosonic_pos),

		max_signal_neg
			= (flow_v_neg[0] <= flow_v_pos[0])
			? flow_v_neg[0] - max_signal
			: flow_v_pos[0] - max_signal,

		max_signal_pos
			= (flow_v_neg[0] <= flow_v_pos[0])
			? flow_v_pos[0] + max_signal
			: flow_v_neg[0] + max_signal,

		bm = std::min(max_signal_neg, 0.0),
		bp = std::max(max_signal_pos, 0.0);

	if (not isnormal(pressure_thermal_neg) or pressure_thermal_neg < 0) {
		throw std::domain_error(
			"Invalid thermal pressure in state_neg: "
			+ boost::lexical_cast<std::string>(pressure_thermal_neg)
		);
	}
	if (not isfinite(pressure_magnetic_neg) or pressure_magnetic_neg < 0) {
		throw std::domain_error(
			"Invalid magnetic pressure in state_neg: "
			+ boost::lexical_cast<std::string>(pressure_magnetic_neg)
		);
	}
	if (not isfinite(max_signal_neg)) {
		throw std::domain_error(
			"Invalid max signal speed in state_neg: "
			+ boost::lexical_cast<std::string>(max_signal_neg)
		);
	}

	if (not isnormal(pressure_thermal_pos) or pressure_thermal_pos < 0) {
		throw std::domain_error(
			"Invalid thermal pressure in state_pos: "
			+ boost::lexical_cast<std::string>(pressure_thermal_pos)
		);
	}
	if (not isfinite(pressure_magnetic_pos) or pressure_magnetic_pos < 0) {
		throw std::domain_error(
			"Invalid magnetic pressure in state_pos: "
			+ boost::lexical_cast<std::string>(pressure_magnetic_pos)
		);
	}
	if (not isfinite(max_signal_pos)) {
		throw std::domain_error(
			"Invalid max signal speed in state_pos: "
			+ boost::lexical_cast<std::string>(max_signal_pos)
		);
	}

	if (not isnormal(bp - bm) or bp - bm < 0) {
		MHD flux;

		flux[Mas]    =
		flux[Mom][0] =
		flux[Mom][1] =
		flux[Mom][2] =
		flux[Nrj]    =
		flux[Mag][0] =
		flux[Mag][1] =
		flux[Mag][2] = 0;

		return std::make_tuple(flux, flux, 0);
	}


	MHD flux_neg, flux_pos;

	// compute L/R fluxes along the lines bm/bp: F_{L}-S_{L}U_{L}; F_{R}-S_{R}U_{R}
	flux_neg[Mas]
		= state_neg[Mom][0]
		- bm * state_neg[Mas];

	flux_pos[Mas]
		= state_pos[Mom][0]
		- bp * state_pos[Mas];

	flux_neg[Mom]
		= state_neg[Mom]
		* (flow_v_neg[0] - bm);

	flux_pos[Mom]
		= state_pos[Mom]
		* (flow_v_pos[0] - bp);

	flux_neg[Mom][0] += pressure_thermal_neg;

	flux_pos[Mom][0] += pressure_thermal_pos;

	flux_neg[Nrj]
		= flow_v_neg[0] * (pressure_thermal_neg + state_neg[Nrj])
		- bm * state_neg[Nrj];

	flux_pos[Nrj]
		= flow_v_pos[0] * (pressure_thermal_pos + state_pos[Nrj])
		- bp * state_pos[Nrj];

	const auto
		&B_neg = state_neg[Mag],
		&B_pos = state_pos[Mag];

	flux_neg[Mom][0]
		-= 0.5
		* (B_neg[0]*B_neg[0] - B_neg[1]*B_neg[1] - B_neg[2]*B_neg[2])
		/ vacuum_permeability;
	flux_neg[Mom][1] -= B_neg[0] * B_neg[1] / vacuum_permeability;
	flux_neg[Mom][2] -= B_neg[0] * B_neg[2] / vacuum_permeability;

	flux_pos[Mom][0]
		-= 0.5
		* (B_pos[0]*B_pos[0] - B_pos[1]*B_pos[1] - B_pos[2]*B_pos[2])
		/ vacuum_permeability;
	flux_pos[Mom][1] -= B_pos[0] * B_pos[1] / vacuum_permeability;
	flux_pos[Mom][2] -= B_pos[0] * B_pos[2] / vacuum_permeability;

	flux_neg[Nrj]
		+= pressure_magnetic_neg * flow_v_neg[0]
		- B_neg[0] * flow_v_neg.dot(B_neg) / vacuum_permeability;

	flux_pos[Nrj]
		+= pressure_magnetic_pos * flow_v_pos[0]
		- B_pos[0] * flow_v_pos.dot(B_pos) / vacuum_permeability;

	flux_neg[Mag] = B_neg * (flow_v_neg[0] - bm) - B_neg[0] * flow_v_neg;

	flux_pos[Mag] = B_pos * (flow_v_pos[0] - bp) - B_pos[0] * flow_v_pos;

	flux_pos[Mag][0] =
	flux_neg[Mag][0] = 0;

	flux_neg *= bp / (bp - bm) * area * dt;
	flux_pos *= bm / (bm - bp) * area * dt;

	return std::make_tuple(flux_neg, flux_pos, std::max(std::fabs(bp), std::fabs(bm)));
}


}}} // namespaces

#endif // ifndef PAMHD_MHD_N_HLL_ATHENA_HPP
