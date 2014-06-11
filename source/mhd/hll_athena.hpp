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

#include "gensimcell.hpp"
#include "common.hpp"

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
	const MHD_T& state_neg,
	const MHD_T& state_pos,
	const double area,
	const double dt,
	const double adiabatic_index,
	const double vacuum_permeability
) {
	// shorthand for referring to variables
	const Mass_Density_T Rho{};
	const Momentum_Density_T Mom {};
	const Total_Energy_Density_T Nrj{};
	const Magnetic_Field_T Mag{};

	const auto
		velocity_neg(state_neg[Mom] / state_neg[Rho]),
		velocity_pos(state_pos[Mom] / state_pos[Rho]);

	const auto
		pressure_neg
			= get_pressure<
				MHD_T,
				Mass_Density_T,
				Momentum_Density_T,
				Total_Energy_Density_T,
				Magnetic_Field_T
			>(state_neg, adiabatic_index, vacuum_permeability),

		pressure_pos
			= get_pressure<
				MHD_T,
				Mass_Density_T,
				Momentum_Density_T,
				Total_Energy_Density_T,
				Magnetic_Field_T
			>(state_pos, adiabatic_index, vacuum_permeability),

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

	flux_neg[Mom][0] += pressure_neg;

	flux_pos[Mom][0] += pressure_pos;

	flux_neg[Nrj]
		= velocity_neg[0] * (pressure_neg + state_neg[Nrj])
		- bm * state_neg[Nrj];

	flux_pos[Nrj]
		= velocity_pos[0] * (pressure_pos + state_pos[Nrj])
		- bp * state_pos[Nrj];

	const auto
		B_neg = state_neg[Mag],
		B_pos = state_pos[Mag];

	flux_neg[Mom][0]
		-= 0.5
		* (B_neg[0]*B_neg[0] - B_neg[1]*B_neg[1] - B_neg[2]*B_neg[2])
		/ vacuum_permeability;
	flux_neg[Mom][1] -= B_neg[0]*B_neg[1] / vacuum_permeability;
	flux_neg[Mom][2] -= B_neg[0]*B_neg[2] / vacuum_permeability;

	flux_pos[Mom][0]
		-= 0.5
		* (B_pos[0]*B_pos[0] - B_pos[1]*B_pos[1] - B_pos[2]*B_pos[2])
		/ vacuum_permeability;
	flux_pos[Mom][1] -= B_pos[0]*B_pos[1] / vacuum_permeability;
	flux_pos[Mom][2] -= B_pos[0]*B_pos[2] / vacuum_permeability;

	flux_neg[Nrj]
		+= pressure_neg * velocity_neg[0]
		- B_neg[0] * velocity_neg.dot(B_neg) / vacuum_permeability;

	flux_pos[Nrj]
		+= pressure_pos * velocity_pos[0]
		- B_pos[0] * velocity_pos.dot(B_pos) / vacuum_permeability;

	flux_neg[Mag] = B_neg * (velocity_neg[0] - bm) - B_neg[0] * velocity_neg;

	flux_pos[Mag] = B_pos * (velocity_pos[0] - bp) - B_pos[0] * velocity_pos;

	flux_pos[Mag][0] =
	flux_neg[Mag][0] = 0;


	MHD_T ret_val;

	const auto factor = 0.5 * (bp + bm) / (bp - bm);
	ret_val[Rho]
		= 0.5 * (flux_neg[Rho] + flux_pos[Rho])
		+ factor * (flux_neg[Rho] - flux_pos[Rho]);

	ret_val[Mom]
		= 0.5 * (flux_neg[Mom] + flux_pos[Mom])
		+ factor * (flux_neg[Mom] - flux_pos[Mom]);

	ret_val[Nrj]
		= 0.5 * (flux_neg[Nrj] + flux_pos[Nrj])
		+ factor * (flux_neg[Nrj] - flux_pos[Nrj]);

	ret_val[Mag]
		= 0.5 * (flux_neg[Mag] + flux_pos[Mag])
		+ factor * (flux_neg[Mag] - flux_pos[Mag]);

	ret_val[Rho] *= area * dt;
	ret_val[Mom] *= area * dt;
	ret_val[Nrj] *= area * dt;
	ret_val[Mag] *= area * dt;

	return std::make_pair(ret_val, std::max(std::fabs(bp), std::fabs(bm)));
}


}}} // namespaces

#endif // ifndef PAMHD_MHD_HLL_ATHENA_HPP
