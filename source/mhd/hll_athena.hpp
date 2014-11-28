/*
HLL MHD solver of PAMHD adapted from Athena (https://trac.princeton.edu/Athena/wiki)
Copyright 2003 James M. Stone
Copyright 2003 Thomas A. Gardiner
Copyright 2003 Peter J. Teuben
Copyright 2003 John F. Hawley
Copyright 2014 Ilja Honkonen

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

#ifndef PAMHD_MHD_HLL_ATHENA_HPP
#define PAMHD_MHD_HLL_ATHENA_HPP


#include "cmath"
#include "limits"

#include "boost/lexical_cast.hpp"
#include "gensimcell.hpp"

#include "mhd/common.hpp"


namespace pamhd {
namespace mhd {
namespace athena {


/*!
Returns the flux between states.

\param [Mass_Density_T] Used to access mass density in MHD states
\param [state_neg] MHD state in negative x direction from the face
\param [state_pos] MHD state in positive x direction from the face
\param [area] Area of the face shared by volumes of state_neg and state_pos
\param [dt] Length of time for which flux is calculated
*/
template <
	class MHD_T,
	class Mass_Density_T,
	class Momentum_Density_T,
	class Total_Energy_Density_T,
	class Magnetic_Field_T
> std::pair<MHD_T, double> get_flux_hll(
	MHD_T state_neg,
	MHD_T state_pos,
	const double area,
	const double dt,
	const double adiabatic_index,
	const double vacuum_permeability
) {
	using std::isnormal;
	using std::isfinite;

	// shorthand for referring to variables
	const Mass_Density_T Rho{};
	const Momentum_Density_T Mom{};
	const Total_Energy_Density_T Nrj{};
	const Magnetic_Field_T Mag{};

	if (not isnormal(state_neg[Rho]) or state_neg[Rho] < 0) {
		throw std::domain_error(
			"Invalid density in state_neg: "
			+ boost::lexical_cast<std::string>(state_neg[Rho])
		);
	}
	if (not isfinite(state_neg[Mom][0])) {
		throw std::domain_error(
			"Invalid x momentum in state_neg: "
			+ boost::lexical_cast<std::string>(state_neg[Mom][0])
		);
	}
	if (not isfinite(state_neg[Mom][1])) {
		throw std::domain_error(
			"Invalid y momentum in state_neg: "
			+ boost::lexical_cast<std::string>(state_neg[Mom][1])
		);
	}
	if (not isnormal(state_neg[Nrj]) or state_neg[Nrj] < 0) {
		throw std::domain_error(
			"Invalid total energy in state_neg: "
			+ boost::lexical_cast<std::string>(state_neg[Nrj])
		);
	}
	if (not isfinite(state_neg[Mom][2])) {
		throw std::domain_error(
			"Invalid z momentum in state_neg: "
			+ boost::lexical_cast<std::string>(state_neg[Mom][2])
		);
	}
	if (not isfinite(state_neg[Mag][0])) {
		throw std::domain_error(
			"Invalid x magnetic field in state_neg: "
			+ boost::lexical_cast<std::string>(state_neg[Mag][0])
		);
	}
	if (not isfinite(state_neg[Mag][1])) {
		throw std::domain_error(
			"Invalid y magnetic field in state_neg: "
			+ boost::lexical_cast<std::string>(state_neg[Mag][1])
		);
	}
	if (not isfinite(state_neg[Mag][2])) {
		throw std::domain_error(
			"Invalid z magnetic field in state_neg: "
			+ boost::lexical_cast<std::string>(state_neg[Mag][2])
		);
	}

	if (not isnormal(state_pos[Rho]) or state_pos[Rho] < 0) {
		throw std::domain_error(
			"Invalid density in state_pos: "
			+ boost::lexical_cast<std::string>(state_pos[Rho])
		);
	}
	if (not isfinite(state_pos[Mom][0])) {
		throw std::domain_error(
			"Invalid x momentum in state_pos: "
			+ boost::lexical_cast<std::string>(state_pos[Mom][0])
		);
	}
	if (not isfinite(state_pos[Mom][1])) {
		throw std::domain_error(
			"Invalid y momentum in state_pos: "
			+ boost::lexical_cast<std::string>(state_pos[Mom][1])
		);
	}
	if (not isnormal(state_pos[Nrj]) or state_pos[Nrj] < 0) {
		throw std::domain_error(
			"Invalid total energy in state_pos: "
			+ boost::lexical_cast<std::string>(state_pos[Nrj])
		);
	}
	if (not isfinite(state_pos[Mom][2])) {
		throw std::domain_error(
			"Invalid z momentum in state_pos: "
			+ boost::lexical_cast<std::string>(state_pos[Mom][2])
		);
	}
	if (not isfinite(state_pos[Mag][0])) {
		throw std::domain_error(
			"Invalid x magnetic field in state_pos: "
			+ boost::lexical_cast<std::string>(state_pos[Mag][0])
		);
	}
	if (not isfinite(state_pos[Mag][1])) {
		throw std::domain_error(
			"Invalid y magnetic field in state_pos: "
			+ boost::lexical_cast<std::string>(state_pos[Mag][1])
		);
	}
	if (not isfinite(state_pos[Mag][2])) {
		throw std::domain_error(
			"Invalid z magnetic field in state_pos: "
			+ boost::lexical_cast<std::string>(state_pos[Mag][2])
		);
	}

	const auto
		velocity_neg(state_neg[Mom] / state_neg[Rho]),
		velocity_pos(state_pos[Mom] / state_pos[Rho]);

	const auto
		pressure_thermal_neg
			= get_pressure<
				MHD_T,
				Mass_Density_T,
				Momentum_Density_T,
				Total_Energy_Density_T,
				Magnetic_Field_T
			>(state_neg, adiabatic_index, vacuum_permeability),

		pressure_thermal_pos
			= get_pressure<
				MHD_T,
				Mass_Density_T,
				Momentum_Density_T,
				Total_Energy_Density_T,
				Magnetic_Field_T
			>(state_pos, adiabatic_index, vacuum_permeability),

		pressure_magnetic_neg
			= state_neg[Mag].squaredNorm() / (2 * vacuum_permeability),

		pressure_magnetic_pos
			= state_pos[Mag].squaredNorm() / (2 * vacuum_permeability),

		max_signal_neg
			= get_fast_magnetosonic_speed<
				MHD_T,
				Mass_Density_T,
				Momentum_Density_T,
				Total_Energy_Density_T,
				Magnetic_Field_T
			>(state_neg, adiabatic_index, vacuum_permeability),

		max_signal_pos
			= get_fast_magnetosonic_speed<
				MHD_T,
				Mass_Density_T,
				Momentum_Density_T,
				Total_Energy_Density_T,
				Magnetic_Field_T
			>(state_pos, adiabatic_index, vacuum_permeability),

		bm = std::min(velocity_neg[0] - max_signal_neg, 0.0),
		bp = std::max(velocity_pos[0] + max_signal_pos, 0.0);

	if (not isnormal(pressure_thermal_neg) or pressure_thermal_neg < 0) {
		throw std::domain_error(
			"Invalid thermal pressure in state_neg: "
			+ boost::lexical_cast<std::string>(pressure_thermal_neg)
		);
	}
	if (not isnormal(pressure_magnetic_neg) or pressure_magnetic_neg < 0) {
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
	if (not isnormal(pressure_magnetic_pos) or pressure_magnetic_pos < 0) {
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

	MHD_T flux_neg, flux_pos;

	// compute L/R fluxes along the lines bm/bp: F_{L}-S_{L}U_{L}; F_{R}-S_{R}U_{R}
	flux_neg[Rho]
		= state_neg[Mom][0]
		- bm * state_neg[Rho];

	flux_pos[Rho]
		= state_pos[Mom][0]
		- bp * state_pos[Rho];

	flux_neg[Mom]
		= state_neg[Mom]
		* (velocity_neg[0] - bm);

	flux_pos[Mom]
		= state_pos[Mom]
		* (velocity_pos[0] - bp);

	flux_neg[Mom][0] += pressure_thermal_neg;

	flux_pos[Mom][0] += pressure_thermal_pos;

	flux_neg[Nrj]
		= velocity_neg[0] * (pressure_thermal_neg + state_neg[Nrj])
		- bm * state_neg[Nrj];

	flux_pos[Nrj]
		= velocity_pos[0] * (pressure_thermal_pos + state_pos[Nrj])
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
		+= pressure_magnetic_neg * velocity_neg[0]
		- B_neg[0] * velocity_neg.dot(B_neg) / vacuum_permeability;

	flux_pos[Nrj]
		+= pressure_magnetic_pos * velocity_pos[0]
		- B_pos[0] * velocity_pos.dot(B_pos) / vacuum_permeability;

	flux_neg[Mag] = B_neg * (velocity_neg[0] - bm) - B_neg[0] * velocity_neg;

	flux_pos[Mag] = B_pos * (velocity_pos[0] - bp) - B_pos[0] * velocity_pos;

	flux_pos[Mag][0] =
	flux_neg[Mag][0] = 0;


	MHD_T ret_val = (flux_neg + flux_pos) / 2;
	const auto divisor = 2 * (bp - bm);
	if (isnormal(divisor) and divisor > 0) {
		ret_val += (flux_neg - flux_pos) * (bp + bm) / divisor;
	}

	ret_val *= area * dt;

	return std::make_pair(ret_val, std::max(std::fabs(bp), std::fabs(bm)));
}


}}} // namespaces

#endif // ifndef PAMHD_MHD_HLL_ATHENA_HPP
