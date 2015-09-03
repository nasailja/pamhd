/*
Solves the MHD part of PAMHD using an external flux function.

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

#ifndef PAMHD_MHD_SOLVE_HPP
#define PAMHD_MHD_SOLVE_HPP


#include "cmath"
#include "limits"
#include "tuple"
#include "vector"

#include "dccrg.hpp"

#include "mhd/N_solve.hpp"
#include "mhd/variables.hpp"


namespace pamhd {
namespace mhd {


/*!
Advances MHD solution for one time step of length dt with given solver.

Returns the maximum allowed length of time step for the next step on this process.
*/
template <
	class Solver,
	class Cell,
	class Geometry,
	class Mass_Density_Getter,
	class Momentum_Density_Getter,
	class Total_Energy_Density_Getter,
	class Magnetic_Field_Getter,
	class Mass_Density_Flux_Getter,
	class Momentum_Density_Flux_Getter,
	class Total_Energy_Density_Flux_Getter,
	class Magnetic_Field_Flux_Getter,
	class Cell_Type_Getter,
	class Cell_Type
> std::pair<double, size_t> solve(
	const Solver solver,
	const size_t solve_start_index,
	dccrg::Dccrg<Cell, Geometry>& grid,
	const double dt,
	const double adiabatic_index,
	const double vacuum_permeability,
	const Mass_Density_Getter Mas,
	const Momentum_Density_Getter Mom,
	const Total_Energy_Density_Getter Nrj,
	const Magnetic_Field_Getter Mag,
	const Mass_Density_Flux_Getter Mas_f,
	const Momentum_Density_Flux_Getter Mom_f,
	const Total_Energy_Density_Flux_Getter Nrj_f,
	const Magnetic_Field_Flux_Getter Mag_f,
	const Cell_Type_Getter Cell_t,
	const Cell_Type normal_cell,
	const Cell_Type dont_solve_cell
) {
	// Use N fluid solver with same total and fluid mass
	return
		N_solve(
			solver,
			solve_start_index,
			grid,
			dt,
			adiabatic_index,
			vacuum_permeability,
			Mas,
			Mas,
			Mom,
			Nrj,
			Mag,
			Mas_f,
			Mom_f,
			Nrj_f,
			Mag_f,
			Cell_t,
			normal_cell,
			dont_solve_cell
		);
}


/*!
Applies the MHD solution to given cells.

Returns 1 + last index where solution was applied.
*/
template <
	class Cell,
	class Geometry,
	class Mass_Density_Getter,
	class Momentum_Density_Getter,
	class Total_Energy_Density_Getter,
	class Magnetic_Field_Getter,
	class Mass_Density_Flux_Getter,
	class Momentum_Density_Flux_Getter,
	class Total_Energy_Density_Flux_Getter,
	class Magnetic_Field_Flux_Getter,
	class Cell_Type_Getter,
	class Cell_Type
> void apply_fluxes(
	dccrg::Dccrg<Cell, Geometry>& grid,
	const double adiabatic_index,
	const double vacuum_permeability,
	const Mass_Density_Getter Mas,
	const Momentum_Density_Getter Mom,
	const Total_Energy_Density_Getter Nrj,
	const Magnetic_Field_Getter Mag,
	const Mass_Density_Flux_Getter Mas_f,
	const Momentum_Density_Flux_Getter Mom_f,
	const Total_Energy_Density_Flux_Getter Nrj_f,
	const Magnetic_Field_Flux_Getter Mag_f,
	const Cell_Type_Getter Cell_t,
	const Cell_Type normal_cell,
	const bool check_new_state = true
) {
	using std::get;
	const auto& cell_data_pointers = grid.get_cell_data_pointers();

	for (size_t i = 0 ; i < cell_data_pointers.size(); i++) {
		const auto& cell_id = get<0>(cell_data_pointers[i]);

		// process only inner xor outer cells
		if (cell_id == dccrg::error_cell) {
			continue;
		}

		const auto& offset = get<2>(cell_data_pointers[i]);
		if (offset[0] != 0 or offset[1] != 0 or offset[2] != 0) {
			continue;
		}

		auto* const cell_data = get<1>(cell_data_pointers[i]);
		if (Cell_t(*cell_data) != normal_cell) {
			continue;
		}

		const auto length = grid.geometry.get_length(cell_id);
		const double inverse_volume = 1.0 / (length[0] * length[1] * length[2]);
		try {
			apply_fluxes(
				*cell_data,
				inverse_volume,
				adiabatic_index,
				vacuum_permeability,
				Mas, Mom, Nrj, Mag,
				Mas_f, Mom_f, Nrj_f, Mag_f,
				check_new_state
			);
		} catch (const std::domain_error& error) {
			std::cerr <<  __FILE__ << "(" << __LINE__
				<< ") New MHD state for cell " << cell_id
				<< " at " << grid.geometry.get_center(cell_id)
				<< " would be unphysical because (" << error.what()
				<< ") with flux (* " << inverse_volume << "):\n"
				<< Mas_f(*cell_data) * inverse_volume << ", "
				<< Mom_f(*cell_data) * inverse_volume << ", "
				<< Nrj_f(*cell_data) * inverse_volume << ", "
				<< Mag_f(*cell_data) * inverse_volume
				<< "\ngiving:\n"
				<< Mas(*cell_data) << ", "
				<< Mom(*cell_data) << ", "
				<< Nrj(*cell_data) << ", "
				<< Mag(*cell_data)
				<< "\nwith pressure: "
				<< get_pressure(
					Mas(*cell_data),
					Mom(*cell_data),
					Nrj(*cell_data),
					Mag(*cell_data),
					adiabatic_index,
					vacuum_permeability
				)
				<< std::endl;
			abort();
		}

		Mas_f(*cell_data)    =
		Mom_f(*cell_data)[0] =
		Mom_f(*cell_data)[1] =
		Mom_f(*cell_data)[2] =
		Nrj_f(*cell_data)    =
		Mag_f(*cell_data)[0] =
		Mag_f(*cell_data)[1] =
		Mag_f(*cell_data)[2] = 0;
	}
}


}} // namespaces


#endif // ifndef PAMHD_MHD_SOLVE_HPP
