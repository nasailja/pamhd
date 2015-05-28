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
	class Magnetic_Field_Flux_Getter
> double solve(
	const Solver solver,
	const std::vector<uint64_t>& cells,
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
	const Magnetic_Field_Flux_Getter Mag_f
) {
	// Use N fluid solver with same total and fluid mass
	return
		N_solve(
			solver,
			cells,
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
			Mag_f
		);
}


/*!
Applies the MHD solution to given cells.
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
	class Magnetic_Field_Flux_Getter
> void apply_fluxes(
	const std::vector<uint64_t>& cells,
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
	const bool check_new_state = true
) {
	for (const auto& cell_id: cells) {
		auto* const cell_data = grid[cell_id];
		if (cell_data == NULL) {
			std::cerr <<  __FILE__ << "(" << __LINE__ << "): "
				<< " No data for cell " << cell_id
				<< std::endl;
			abort();
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
