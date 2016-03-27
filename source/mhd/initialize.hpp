/*
Initializes the MHD solution of PAMHD.

Copyright 2014, 2015, 2016 Ilja Honkonen
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

Getters should return a reference to data of corresponding variable
when given a simulation cell's data.
\param [Init_Cond] Initial condition class defined in MHD test program
\param [grid] Grid storing the cells to initialize
\param [cells] List of cells to initialize
\param [adiabatic_index] https://en.wikipedia.org/wiki/Heat_capacity_ratio
\param [vacuum_permeability] https://en.wikipedia.org/wiki/Vacuum_permeability
*/
template <
	class Init_Cond,
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
> void initialize(
	Init_Cond& init_cond,
	dccrg::Dccrg<Cell, Geometry>& grid,
	const std::vector<uint64_t>& cells,
	const double time,
	const double adiabatic_index,
	const double vacuum_permeability,
	const double proton_mass,
	const bool verbose,
	const Mass_Density_Getter Mas,
	const Momentum_Density_Getter Mom,
	const Total_Energy_Density_Getter Nrj,
	const Magnetic_Field_Getter Mag,
	const Mass_Density_Flux_Getter Mas_f,
	const Momentum_Density_Flux_Getter Mom_f,
	const Total_Energy_Density_Flux_Getter Nrj_f,
	const Magnetic_Field_Flux_Getter Mag_f
) {
	if (verbose and grid.get_rank() == 0) {
		std::cout << "Setting default MHD state... ";
		std::cout.flush();
	}
	// set default state
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

		auto* const cell_data = grid[cell_id];
		if (cell_data == NULL) {
			std::cerr <<  __FILE__ << "(" << __LINE__ << ") No data for cell: "
				<< cell_id
				<< std::endl;
			abort();
		}

		// zero fluxes
		Mas_f(*cell_data)    =
		Nrj_f(*cell_data)    =
		Mom_f(*cell_data)[0] =
		Mom_f(*cell_data)[1] =
		Mom_f(*cell_data)[2] =
		Mag_f(*cell_data)[0] =
		Mag_f(*cell_data)[1] =
		Mag_f(*cell_data)[2] = 0;

		const auto mass_density
			= proton_mass
			* [&](){
				try {
					return init_cond.default_data.get_data(Number_Density(), cell_center, time);
				} catch (mup::ParserError& e) {
					std::cout << "Couldn't get number density for default initial condition." << std::endl;
					throw;
				}
			}();
		const auto velocity
			= [&](){
				try {
					return init_cond.default_data.get_data(Velocity(), cell_center, time);
				} catch (mup::ParserError& e) {
					std::cout << "Couldn't get velocity for default initial condition." << std::endl;
					throw;
				}
			}();
		const auto pressure
			= [&](){
				try {
					return init_cond.default_data.get_data(Pressure(), cell_center, time);
				} catch (mup::ParserError& e) {
					std::cout << "Couldn't get pressure for default initial condition." << std::endl;
					throw;
				}
			}();
		const auto magnetic_field
			= [&](){
				try {
					return init_cond.default_data.get_data(Magnetic_Field(), cell_center, time);
				} catch (mup::ParserError& e) {
					std::cout << "Couldn't get magnetic field for default initial condition." << std::endl;
				throw;
				}
			}();

		Mas(*cell_data) = mass_density;
		Mom(*cell_data) = mass_density * velocity;
		Mag(*cell_data) = magnetic_field;
		Nrj(*cell_data) = get_total_energy_density(
			mass_density,
			velocity,
			pressure,
			magnetic_field,
			adiabatic_index,
			vacuum_permeability
		);
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

			auto* const cell_data = grid[cell_id];
			if (cell_data == NULL) {
				std::cerr <<  __FILE__ << "(" << __LINE__ << ") No data for cell: "
					<< cell_id
					<< std::endl;
				abort();
			}

			const auto mass_density
				= proton_mass
				* [&](){
					try {
						return init_cond.get_data(Number_Density(), bdy_i, cell_center, time);
					} catch (mup::ParserError& e) {
						std::cout << "Couldn't get density for initial condition geometry " << bdy_i << std::endl;
						throw;
					}
				}();
			const auto velocity
				= [&](){
					try {
						return init_cond.get_data(Velocity(), bdy_i, cell_center, time);
					} catch (mup::ParserError& e) {
						std::cout << "Couldn't get velocity for initial condition geometry " << bdy_i << std::endl;
						throw;
					}
				}();
			const auto pressure
				= [&](){
					try {
						return init_cond.get_data(Pressure(), bdy_i, cell_center, time);
					} catch (mup::ParserError& e) {
						std::cout << "Couldn't get pressure for initial condition geometry " << bdy_i << std::endl;
						throw;
					}
				}();
			const auto magnetic_field
				= [&](){
					try {
						return init_cond.get_data(Magnetic_Field(), bdy_i, cell_center, time);
					} catch (mup::ParserError& e) {
						std::cout << "Couldn't get magnetic field for initial condition geometry " << bdy_i << std::endl;
					throw;
					}
				}();

			Mas(*cell_data) = mass_density;
			Mom(*cell_data) = mass_density * velocity;
			Mag(*cell_data) = magnetic_field;
			Nrj(*cell_data) = get_total_energy_density(
				mass_density,
				velocity,
				pressure,
				magnetic_field,
				adiabatic_index,
				vacuum_permeability
			);
		}
	}
	if (verbose and grid.get_rank() == 0) {
		std::cout << "done" << std::endl;
	}
}

}} // namespaces

#endif // ifndef PAMHD_MHD_INITIALIZE_HPP
