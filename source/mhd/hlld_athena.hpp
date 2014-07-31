/*
HLLD MHD solver of PAMHD adapted from Athena (https://trac.princeton.edu/Athena/wiki)
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

#ifndef PAMHD_MHD_HLLD_ATHENA_HPP
#define PAMHD_MHD_HLLD_ATHENA_HPP


#include "cmath"
#include "limits"

#include "gensimcell.hpp"

#include "mhd/common.hpp"


namespace pamhd {
namespace mhd {
namespace athena {


/*!
Returns the flux between states and maximum signal speed from interface.

Returns negative signal speed in case of failure, for example,
if encounters negative pressure, etc.

Takahiro Miyoshi, Kanya Kusano (M&K):
A multi-state HLL approximate Riemann solver for ideal MHD,
Journal of Computational Physics, 208, 315-344, 2005.

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
> std::pair<MHD_T, double> get_flux_hlld(
	MHD_T state_neg,
	MHD_T state_pos,
	const double area,
	const double dt,
	const double adiabatic_index,
	const double vacuum_permeability
) {
	constexpr double epsilon = std::numeric_limits<double>::epsilon();

	// shorthand notation for simulation variables
	const Mass_Density_T Mas{};
	const Momentum_Density_T Mom{};
	const Total_Energy_Density_T Nrj{};
	const Magnetic_Field_T Mag{};

	// div B = 0 requires identical Bx between states
	state_neg[Mag][0] =
	state_pos[Mag][0] = (state_neg[Mag][0] + state_pos[Mag][0]) / 2;

	/*
	Figure 3 in M&K, Us are MHD states, Ss are signal speeds
	from interface between state_neg (Ul) and state_pos (Ur)

	time

	^      Sl*      Sm    Sr*       Sr
	|        \ Ul**  |Ur** / Ur*  __/
	|Sl_  Ul* \     /    _/    __/
	|   \__    \    |   /   __/
	| Ul   \__  \  /  _/ __/      Ur
	|         \__\ | /__/
	|            \\///
	-----------------------------------> space

	Variable names:
	Ul   -> state_neg
	Ur   -> state_pos
	Ul*  -> state_s_neg
	Ul** -> state_s2_neg
	Ur*  -> state_s_pos
	Ur** -> state_s2_pos
	Sl*  -> signal_s_neg
	Sr*  -> signal_s_pos
	Sm   -> signal_middle
	*/

	// compute left & right wave speeds (equation 67)
	const auto
		inv_permeability = 1.0 / vacuum_permeability,
		inv_permeability_sqrt = 1.0 / std::sqrt(vacuum_permeability),
		interface_Bx2 = std::pow(state_neg[Mag][0], 2),

		fast_magnetosonic_neg
			= get_fast_magnetosonic_speed<
				MHD_T,
				Mass_Density_T,
				Momentum_Density_T,
				Total_Energy_Density_T,
				Magnetic_Field_T
			>(state_neg, adiabatic_index, vacuum_permeability),

		fast_magnetosonic_pos
			= get_fast_magnetosonic_speed<
				MHD_T,
				Mass_Density_T,
				Momentum_Density_T,
				Total_Energy_Density_T,
				Magnetic_Field_T
			>(state_pos, adiabatic_index, vacuum_permeability),

		max_signal = std::max(fast_magnetosonic_neg, fast_magnetosonic_pos);

	const auto
		flow_v_neg
			= state_neg[Mom]
			/ state_neg[Mas],

		flow_v_pos
			= state_pos[Mom]
			/ state_pos[Mas];

	const auto
		max_signal_neg
			= (flow_v_neg[0] <= flow_v_pos[0])
			? flow_v_neg[0] - max_signal
			: flow_v_pos[0] - max_signal,

		max_signal_pos
			= (flow_v_neg[0] <= flow_v_pos[0])
			? flow_v_pos[0] + max_signal
			: flow_v_neg[0] + max_signal;

	const MHD_T
		flux_neg
			= get_flux<
				MHD_T,
				Mass_Density_T,
				Momentum_Density_T,
				Total_Energy_Density_T,
				Magnetic_Field_T
			>(state_neg, adiabatic_index, vacuum_permeability),

		flux_pos
			= get_flux<
				MHD_T,
				Mass_Density_T,
				Momentum_Density_T,
				Total_Energy_Density_T,
				Magnetic_Field_T
			>(state_pos, adiabatic_index, vacuum_permeability);

	// return upwind flux if flow is supermagnetosonic
	if (max_signal_neg >= 0.0) {
		return std::make_pair(flux_neg * area * dt, std::fabs(max_signal_pos));
	} else if (max_signal_pos <= 0.0) {
		return std::make_pair(flux_pos * area * dt, std::fabs(max_signal_neg));
	}

	const auto
		signal_flow_diff_neg = max_signal_neg - flow_v_neg[0],
		signal_flow_diff_pos = max_signal_pos - flow_v_pos[0],

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

		nrj_magnetic_neg = 0.5 * state_neg[Mag].squaredNorm() * inv_permeability,
		nrj_magnetic_pos = 0.5 * state_pos[Mag].squaredNorm() * inv_permeability,

		// eqn 38
		signal_middle
			= (
				signal_flow_diff_pos * state_pos[Mom][0]
				- signal_flow_diff_neg * state_neg[Mom][0]
				- (pressure_pos + nrj_magnetic_pos)
				+ (pressure_neg + nrj_magnetic_neg)
			) / (
				signal_flow_diff_pos * state_pos[Mas]
				- signal_flow_diff_neg * state_neg[Mas]
			),

		signal_max_middle_diff_neg = max_signal_neg - signal_middle,
		signal_max_middle_diff_pos = max_signal_pos - signal_middle,

		// eqn 43
		density_s_neg
			= state_neg[Mas]
			* signal_flow_diff_neg
			/ signal_max_middle_diff_neg,
		density_s_sqrt_neg = std::sqrt(density_s_neg),

		density_s_pos
			= state_pos[Mas]
			* signal_flow_diff_pos
			/ signal_max_middle_diff_pos,
		density_s_sqrt_pos = std::sqrt(density_s_pos),

		// eqn 51
		signal_s_neg
			= signal_middle
			- std::fabs(state_neg[Mag][0])
				/ std::sqrt(density_s_neg * vacuum_permeability),
		signal_s_pos
			= signal_middle
			+ std::fabs(state_pos[Mag][0])
				/ std::sqrt(density_s_pos * vacuum_permeability),


		/*
		Compute intermediate states
		*/

		ptst
			= pressure_neg
			+ nrj_magnetic_neg
			+ state_neg[Mas]
				* signal_flow_diff_neg
				* (signal_middle - flow_v_neg[0]);

	if (pressure_neg < 0 or pressure_pos < 0) {
		return std::make_pair(MHD_T(), -1);
	}

	MHD_T state_s_neg;
	const auto tmp_neg
		= state_neg[Mas] * signal_flow_diff_neg * signal_max_middle_diff_neg
		- interface_Bx2 * inv_permeability;

	if (std::fabs(tmp_neg) < epsilon * ptst) {
		// degenerate case
		state_s_neg[Mom] = density_s_neg * flow_v_neg;
		state_s_neg[Mag] = state_neg[Mag];
	} else {
		// eqns 44 and 46
		const auto tmp1
			= state_neg[Mag][0]
			* (signal_middle - flow_v_neg[0])
			* inv_permeability
			/ tmp_neg;
		state_s_neg[Mom] = density_s_neg * (flow_v_neg - state_neg[Mag] * tmp1);

		// eqns 45 and 47
		const auto tmp2
			= (
				state_neg[Mas] * std::pow(signal_flow_diff_neg, 2)
				- interface_Bx2 * inv_permeability
			) / tmp_neg;
		state_s_neg[Mag] = state_neg[Mag] * tmp2;
	}
	state_s_neg[Mas] = density_s_neg;
	state_s_neg[Mom][0] = density_s_neg * signal_middle;
	state_s_neg[Mag][0] = state_neg[Mag][0];

	const auto v_dot_b_s_neg = state_s_neg[Mom].dot(state_s_neg[Mag]) / density_s_neg;

	// eqn 48
	state_s_neg[Nrj]
		= (
			signal_flow_diff_neg * state_neg[Nrj]
			- (pressure_neg + nrj_magnetic_neg) * flow_v_neg[0]
			+ ptst * signal_middle
			+ state_neg[Mag][0] * inv_permeability
				* (flow_v_neg.dot(state_neg[Mag]) - v_dot_b_s_neg)
		) / signal_max_middle_diff_neg;

	MHD_T state_s_pos;
	const auto tmp_pos
		= state_pos[Mas] * signal_flow_diff_pos * signal_max_middle_diff_pos
		- interface_Bx2 * inv_permeability;

	if (std::fabs(tmp_pos) < epsilon * ptst) {
		state_s_pos[Mom] = density_s_pos * flow_v_pos;
		state_s_pos[Mag] = state_pos[Mag];
	} else {
		const auto tmp1
			= state_pos[Mag][0]
			* (signal_middle - flow_v_pos[0])
			* inv_permeability
			/ tmp_pos;
		state_s_pos[Mom] = density_s_pos * (flow_v_pos - state_pos[Mag] * tmp1);

		const auto tmp2
			= (
				state_pos[Mas] * std::pow(signal_flow_diff_pos, 2)
				- interface_Bx2 * inv_permeability
			) / tmp_pos;
		state_s_pos[Mag] = state_pos[Mag] * tmp2;
	}
	state_s_pos[Mas] = density_s_pos;
	state_s_pos[Mom][0] = density_s_pos * signal_middle;
	state_s_pos[Mag][0] = state_pos[Mag][0];

	const auto v_dot_b_s_pos = state_s_pos[Mom].dot(state_s_pos[Mag]) / density_s_pos;

	// eqn 48
	state_s_pos[Nrj]
		= (
			signal_flow_diff_pos * state_pos[Nrj]
			- (pressure_pos + nrj_magnetic_pos) * flow_v_pos[0]
			+ ptst * signal_middle
			+ state_pos[Mag][0] * inv_permeability
				* (flow_v_pos.dot(state_pos[Mag]) - v_dot_b_s_pos)
		) / signal_max_middle_diff_pos;

	const auto
		flow_v_s_neg = state_s_neg[Mom] / density_s_neg,
		flow_v_s_pos = state_s_pos[Mom] / density_s_pos;

	MHD_T state_s2_neg, state_s2_pos;
	if (interface_Bx2 / 2 < epsilon * ptst) {
		state_s2_neg = state_s_neg;
		state_s2_pos = state_s_pos;
	} else {
		const auto inv_sum_density_sqrt
			= 1.0 / (density_s_sqrt_neg + density_s_sqrt_pos);

		const auto Bx_sign = (state_neg[Mag][0] > 0) ? 1.0 : -1.0;

		state_s2_neg[Mas] = state_s_neg[Mas];
		state_s2_pos[Mas] = state_s_pos[Mas];

		const auto flow_v_s2 = [&](){
			typename Momentum_Density_T::data_type temp
				= inv_sum_density_sqrt
				* (
					density_s_sqrt_neg * flow_v_s_neg
					+ density_s_sqrt_pos * flow_v_s_pos
					+ Bx_sign * inv_permeability_sqrt
						* (state_s_pos[Mag] - state_s_neg[Mag])
				);

			temp[0] = signal_middle;

			return temp;
		}();

		// eqs 59 & 60
		state_s2_neg[Mom] = state_s_neg[Mas] * flow_v_s2;
		state_s2_neg[Mom][0] = state_s_neg[Mom][0];

		state_s2_pos[Mom] = state_s_pos[Mas] * flow_v_s2;
		state_s2_pos[Mom][0] = state_s_pos[Mom][0];

		// eqs 61 & 62
		state_s2_neg[Mag] =
		state_s2_pos[Mag]
			= inv_sum_density_sqrt
			* (
				density_s_sqrt_neg * state_s_pos[Mag]
				+ density_s_sqrt_pos * state_s_neg[Mag]
				+ Bx_sign * density_s_sqrt_neg * density_s_sqrt_pos
					* (flow_v_s_pos - flow_v_s_neg)
					* std::sqrt(vacuum_permeability)
			);
		state_s2_neg[Mag][0] =
		state_s2_pos[Mag][0] = state_neg[Mag][0];

		// eqn 63
		const auto v_dot_b_s2 = flow_v_s2.dot(state_s2_neg[Mag]);

		state_s2_neg[Nrj]
			= state_s_neg[Nrj]
			- density_s_sqrt_neg
				* Bx_sign
				* (v_dot_b_s_neg - v_dot_b_s2)
				* inv_permeability_sqrt;

		state_s2_pos[Nrj]
			= state_s_pos[Nrj]
			+ density_s_sqrt_pos
				* Bx_sign
				* (v_dot_b_s_pos - v_dot_b_s2)
				* inv_permeability_sqrt;
	}

	MHD_T flux;
	if (signal_s_neg >= 0) {

		flux = flux_neg + (state_s_neg - state_neg) * max_signal_neg;

	} else if (signal_middle >= 0) {

		flux
			= flux_neg
			+ state_s2_neg * signal_s_neg
			- state_s_neg * (signal_s_neg - max_signal_neg)
			- state_neg * max_signal_neg;

	} else if (signal_s_pos > 0) {

		flux
			= flux_pos
			+ state_s2_pos * signal_s_pos
			- state_s_pos * (signal_s_pos - max_signal_pos)
			- state_pos * max_signal_pos;

	} else {

		flux = flux_pos + (state_s_pos - state_pos) * max_signal_pos;

	}
	flux[Mag][0] = 0;

	flux *= area * dt;

	return std::make_pair(
		flux,
		std::max(
			std::fabs(max_signal_neg),
			std::fabs(max_signal_pos)
		)
	);
}


}}} // namespaces

#endif // ifndef PAMHD_MHD_HLLD_ATHENA_HPP
