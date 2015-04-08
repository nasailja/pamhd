/*
Initializes the MHD solution of PAMHD.

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

#ifndef PAMHD_MHD_INITIALIZE_HPP
#define PAMHD_MHD_INITIALIZE_HPP


#include "cmath"
#include "iostream"
#include "limits"

#include "dccrg.hpp"

#include "mhd/common.hpp"
#include "mhd/variables.hpp"


namespace pamhd {
namespace mhd {


/*!
Sets the initial state of MHD simulation and zeroes fluxes.

\param [Grid_T] Assumed to provide a dccrg API
\param [MHD_T] Used to access the MHD state in grid cells
\param [MHD_Flux_T] Used to access the change in MHD state for next step
\param [Mass_Density_T] Used to access mass density in MHD_T and MHD_Flux_T
\param [Momentum_Density_T] Used to access momentum density
\param [Total_Energy_Density_T] Used to access total energy density
\param [Magnetic_Field_T] Used to access magnetic field
\param [grid] Grid storing the cells to initialize
\param [cells] List of cells to initialize
\param [adiabatic_index] https://en.wikipedia.org/wiki/Heat_capacity_ratio
\param [vacuum_permeability] https://en.wikipedia.org/wiki/Vacuum_permeability
*/
template <
	class MHD_T,
	class MHD_Flux_T,
	class Mass_Density_T,
	class Momentum_Density_T,
	class Total_Energy_Density_T,
	class Magnetic_Field_T,
	class Init_Cond,
	class Cell,
	class Geometry
> void initialize(
	Init_Cond& init_cond,
	dccrg::Dccrg<Cell, Geometry>& grid,
	const std::vector<uint64_t>& cells,
	const double time,
	const double adiabatic_index,
	const double vacuum_permeability,
	const double proton_mass,
	const bool verbose
) {
	static_assert(
		std::is_same<typename MHD_T::data_type, typename MHD_Flux_T::data_type>::value,
		"The data types of variables MHD_T and MHD_Flux_T must be equal"
	);

	const Mass_Density_T Rho{};
	const Momentum_Density_T M{};
	const Total_Energy_Density_T E{};
	const Magnetic_Field_T B{};
	const Velocity V{};
	const Pressure P{};
	const Number_Density N{};

	if (verbose and grid.get_rank() == 0) {
		std::cout << "Setting default MHD state... ";
		std::cout.flush();
	}
	for (const auto cell_id: cells) {
		const auto
			cell_start = grid.geometry.get_min(cell_id),
			cell_end = grid.geometry.get_max(cell_id),
			cell_center = grid.geometry.get_center(cell_id);

		init_cond.add_cell(
			cell_id,
			cell_start,
			cell_end
		);

		auto* cell_data = grid[cell_id];
		if (cell_data == NULL) {
			std::cerr <<  __FILE__ << "(" << __LINE__ << ") No data for cell: "
				<< cell_id
				<< std::endl;
			abort();
		}

		// zero fluxes
		auto& flux = (*cell_data)[MHD_Flux_T()];
		flux[Rho]  =
		flux[E]    =
		flux[M][0] =
		flux[M][1] =
		flux[M][2] =
		flux[B][0] =
		flux[B][1] =
		flux[B][2] = 0;

		// set default state
		MHD_Primitive temp;
		temp[Rho] = init_cond.default_data.get_data(N, cell_center, time) * proton_mass;
		temp[V] = init_cond.default_data.get_data(V, cell_center, time);
		temp[P] = init_cond.default_data.get_data(P, cell_center, time);
		temp[B] = init_cond.default_data.get_data(B, cell_center, time);

		auto& state = (*cell_data)[MHD_T()];
		state = get_conservative<
			typename MHD_T::data_type,
			MHD_Primitive,
			Momentum_Density_T,
			Total_Energy_Density_T,
			Mass_Density_T,
			Velocity,
			Pressure,
			Magnetic_Field_T
		>(temp, adiabatic_index, vacuum_permeability);
	}

	// set non-default initial conditions
	if (verbose and grid.get_rank() == 0) {
		std::cout << "done\nSetting non-default initial MHD state... ";
		std::cout.flush();
	}
	for (size_t bdy_i = 0; bdy_i < init_cond.get_number_of_boundaries(); bdy_i++) {
		const auto& boundary_cells = init_cond.get_cells(bdy_i);

		for (const auto& cell_id: boundary_cells) {
			const auto cell_center = grid.geometry.get_center(cell_id);

			auto* cell_data = grid[cell_id];
			if (cell_data == NULL) {
				std::cerr <<  __FILE__ << "(" << __LINE__ << ") No data for cell: "
					<< cell_id
					<< std::endl;
				abort();
			}

			pamhd::mhd::MHD_Primitive temp;
			temp[Rho] = init_cond.get_data(N, bdy_i, cell_center, time) * proton_mass;
			temp[V] = init_cond.get_data(V, bdy_i, cell_center, time);
			temp[P] = init_cond.get_data(P, bdy_i, cell_center, time);
			temp[B] = init_cond.get_data(B, bdy_i, cell_center, time);

			auto& state = (*cell_data)[MHD_T()];
			state = get_conservative<
				typename MHD_T::data_type,
				MHD_Primitive,
				Momentum_Density_T,
				Total_Energy_Density_T,
				Mass_Density_T,
				Velocity,
				Pressure,
				Magnetic_Field_T
			>(temp, adiabatic_index, vacuum_permeability);
		}
	}
	if (verbose and grid.get_rank() == 0) {
		std::cout << "done" << std::endl;
	}
}

}} // namespaces

#endif // ifndef PAMHD_MHD_INITIALIZE_HPP
