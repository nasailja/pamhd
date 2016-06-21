/*
Handles boundary cell classification logic of MHD part of PAMHD.

Copyright 2015, 2016 Ilja Honkonen
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

#ifndef PAMHD_MHD_BOUNDARIES_HPP
#define PAMHD_MHD_BOUNDARIES_HPP


#include "cmath"
#include "limits"
#include "map"
#include "utility"
#include "vector"

#include "dccrg.hpp"

#include "mhd/common.hpp"
#include "mhd/variables.hpp"


namespace pamhd {
namespace mhd {


/*!
Prepares boundary information needed for MHD solver about each simulation cell.
*/
template<
	class Solver_Info,
	class Cell_Data,
	class Geometry,
	class Boundaries,
	class Boundary_Geometries,
	class Solver_Info_Getter
> void set_solver_info(
	dccrg::Dccrg<Cell_Data, Geometry>& grid,
	const Boundaries& boundaries,
	const Boundary_Geometries& geometries,
	const Solver_Info_Getter& Sol_Info
) {
	for (const auto& cell: grid.cells) {
		Sol_Info(*cell.data) = 0;
	}

	// number density
	for (
		size_t i = 0;
		i < boundaries.get_number_of_value_boundaries(pamhd::mhd::Number_Density());
		i++
	) {
		const auto& value_bdy = boundaries.get_value_boundary(pamhd::mhd::Number_Density(), i);
		const auto& geometry_id = value_bdy.get_geometry_id();
		const auto& cells = geometries.get_cells(geometry_id);
		for (const auto& cell: cells) {
			auto* const cell_data = grid[cell];
			if (cell_data == NULL) {
				std::cerr <<  __FILE__ << ":" << __LINE__ << std::endl;
				abort();
			}
			Sol_Info(*cell_data) |= Solver_Info::mass_density_bdy;
		}
	}

	// velocity
	for (
		size_t i = 0;
		i < boundaries.get_number_of_value_boundaries(pamhd::mhd::Velocity());
		i++
	) {
		const auto& value_bdy = boundaries.get_value_boundary(pamhd::mhd::Velocity(), i);
		const auto& geometry_id = value_bdy.get_geometry_id();
		const auto& cells = geometries.get_cells(geometry_id);
		for (const auto& cell: cells) {
			auto* const cell_data = grid[cell];
			if (cell_data == NULL) {
				std::cerr <<  __FILE__ << ":" << __LINE__ << std::endl;
				abort();
			}
			Sol_Info(*cell_data) |= Solver_Info::velocity_bdy;
		}
	}

	// pressure
	for (
		size_t i = 0;
		i < boundaries.get_number_of_value_boundaries(pamhd::mhd::Pressure());
		i++
	) {
		const auto& value_bdy = boundaries.get_value_boundary(pamhd::mhd::Pressure(), i);
		const auto& geometry_id = value_bdy.get_geometry_id();
		const auto& cells = geometries.get_cells(geometry_id);
		for (const auto& cell: cells) {
			auto* const cell_data = grid[cell];
			if (cell_data == NULL) {
				std::cerr <<  __FILE__ << ":" << __LINE__ << std::endl;
				abort();
			}
			Sol_Info(*cell_data) |= Solver_Info::pressure_bdy;
		}
	}

	// magnetic field
	for (
		size_t i = 0;
		i < boundaries.get_number_of_value_boundaries(pamhd::mhd::Magnetic_Field());
		i++
	) {
		const auto& value_bdy = boundaries.get_value_boundary(pamhd::mhd::Magnetic_Field(), i);
		const auto& geometry_id = value_bdy.get_geometry_id();
		const auto& cells = geometries.get_cells(geometry_id);
		for (const auto& cell: cells) {
			auto* const cell_data = grid[cell];
			if (cell_data == NULL) {
				std::cerr <<  __FILE__ << ":" << __LINE__ << std::endl;
				abort();
			}
			Sol_Info(*cell_data) |= Solver_Info::magnetic_field_bdy;
		}
	}

	// dont solve cells in which no variable is solved
	std::map<uint64_t, unsigned int> dont_solves;
	for (const auto& cell: boundaries.get_dont_solve_cells(pamhd::mhd::Number_Density())) {
		if (dont_solves.count(cell) > 0) {
			dont_solves[cell]++;
		} else {
			dont_solves[cell] = 1;
		}
	}
	for (const auto& cell: boundaries.get_dont_solve_cells(pamhd::mhd::Velocity())) {
		if (dont_solves.count(cell) > 0) {
			dont_solves[cell]++;
		} else {
			dont_solves[cell] = 1;
		}
	}
	for (const auto& cell: boundaries.get_dont_solve_cells(pamhd::mhd::Pressure())) {
		if (dont_solves.count(cell) > 0) {
			dont_solves[cell]++;
		} else {
			dont_solves[cell] = 1;
		}
	}
	for (const auto& cell: boundaries.get_dont_solve_cells(pamhd::mhd::Magnetic_Field())) {
		if (dont_solves.count(cell) > 0) {
			dont_solves[cell]++;
		} else {
			dont_solves[cell] = 1;
		}
	}

	for (const auto& item: dont_solves) {
		if (item.second >= 4) {
			auto* const cell_data = grid[item.first];
			if (cell_data == NULL) {
				std::cerr <<  __FILE__ << ":" << __LINE__ << std::endl;
				abort();
			}
			Sol_Info(*cell_data) |= Solver_Info::dont_solve;
		}
	}
}


template<
	class Cell_Data,
	class Grid_Geometry,
	class Boundaries,
	class Boundary_Geometries,
	class Mass_Getter,
	class Momentum_Getter,
	class Energy_Getter,
	class Magnetic_Field_Getter
> void apply_boundaries(
	dccrg::Dccrg<Cell_Data, Grid_Geometry>& grid,
	Boundaries& boundaries,
	const Boundary_Geometries& bdy_geoms,
	const double simulation_time,
	const Mass_Getter& Mas,
	const Momentum_Getter& Mom,
	const Energy_Getter& Nrj,
	const Magnetic_Field_Getter& Mag,
	const double proton_mass,
	const double adiabatic_index,
	const double vacuum_permeability
) {
	constexpr pamhd::mhd::Number_Density N{};
	constexpr pamhd::mhd::Velocity V{};
	constexpr pamhd::mhd::Pressure P{};
	constexpr pamhd::mhd::Magnetic_Field B{};

	/*
	Value boundaries
	*/

	for (
		size_t i = 0;
		i < boundaries.get_number_of_value_boundaries(N);
		i++
	) {
		auto& value_bdy = boundaries.get_value_boundary(N, i);
		const auto& geometry_id = value_bdy.get_geometry_id();
		const auto& cells = bdy_geoms.get_cells(geometry_id);
		for (const auto& cell: cells) {
			const auto c = grid.geometry.get_center(cell);
			const auto r = sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2]);
			const auto
				lat = asin(c[2] / r),
				lon = atan2(c[1], c[0]);

			const auto mass_density
				= proton_mass
				* value_bdy.get_data(
					simulation_time,
					c[0], c[1], c[2],
					r, lat, lon
				);

			auto* const cell_data = grid[cell];
			if (cell_data == nullptr) {
				std::cerr <<  __FILE__ << "(" << __LINE__ << std::endl;
				abort();
			}

			Mas(*cell_data) = mass_density;
		}
	}

	for (
		size_t i = 0;
		i < boundaries.get_number_of_value_boundaries(V);
		i++
	) {
		auto& value_bdy = boundaries.get_value_boundary(V, i);
		const auto& geometry_id = value_bdy.get_geometry_id();
		const auto& cells = bdy_geoms.get_cells(geometry_id);
		for (const auto& cell: cells) {
			const auto c = grid.geometry.get_center(cell);
			const auto r = sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2]);
			const auto
				lat = asin(c[2] / r),
				lon = atan2(c[1], c[0]);

			const auto velocity = value_bdy.get_data(
				simulation_time,
				c[0], c[1], c[2],
				r, lat, lon
			);

			auto* const cell_data = grid[cell];
			if (cell_data == nullptr) {
				std::cerr <<  __FILE__ << "(" << __LINE__ << std::endl;
				abort();
			}

			Mom(*cell_data) = Mas(*cell_data) * velocity;
		}
	}

	for (
		size_t i = 0;
		i < boundaries.get_number_of_value_boundaries(B);
		i++
	) {
		auto& value_bdy = boundaries.get_value_boundary(B, i);
		const auto& geometry_id = value_bdy.get_geometry_id();
		const auto& cells = bdy_geoms.get_cells(geometry_id);
		for (const auto& cell: cells) {
			const auto c = grid.geometry.get_center(cell);
			const auto r = sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2]);
			const auto
				lat = asin(c[2] / r),
				lon = atan2(c[1], c[0]);

			const auto magnetic_field = value_bdy.get_data(
				simulation_time,
				c[0], c[1], c[2],
				r, lat, lon
			);

			auto* const cell_data = grid[cell];
			if (cell_data == nullptr) {
				std::cerr <<  __FILE__ << "(" << __LINE__ << std::endl;
				abort();
			}

			Mag(*cell_data) = magnetic_field;
		}
	}

	for (
		size_t i = 0;
		i < boundaries.get_number_of_value_boundaries(P);
		i++
	) {
		auto& value_bdy = boundaries.get_value_boundary(P, i);
		const auto& geometry_id = value_bdy.get_geometry_id();
		const auto& cells = bdy_geoms.get_cells(geometry_id);
		for (const auto& cell: cells) {
			const auto c = grid.geometry.get_center(cell);
			const auto r = sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2]);
			const auto
				lat = asin(c[2] / r),
				lon = atan2(c[1], c[0]);

			const auto pressure = value_bdy.get_data(
				simulation_time,
				c[0], c[1], c[2],
				r, lat, lon
			);

			auto* const cell_data = grid[cell];
			if (cell_data == nullptr) {
				std::cerr <<  __FILE__ << "(" << __LINE__ << std::endl;
				abort();
			}

			Nrj(*cell_data) = pamhd::mhd::get_total_energy_density(
				Mas(*cell_data),
				pamhd::mhd::get_velocity(
					Mom(*cell_data),
					Mas(*cell_data)
				),
				pressure,
				Mag(*cell_data),
				adiabatic_index,
				vacuum_permeability
			);
		}
	}

	// copy up-to-date data to copy boundaries
	Cell::set_transfer_all(true, pamhd::mhd::MHD_State_Conservative());
	grid.update_copies_of_remote_neighbors();
	Cell::set_transfer_all(false, pamhd::mhd::MHD_State_Conservative());

	/*
	Copy boundaries
	*/

	for (const auto& item: boundaries.get_copy_boundary_cells(N)) {
		const auto
			&target_id = item[0],
			&source_id = item[1];

		auto
			*target_data = grid[target_id],
			*source_data = grid[source_id];

		if (target_data == nullptr) {
			std::cerr <<  __FILE__ << ":" << __LINE__ << std::endl;
			abort();
		}
		if (source_data == nullptr) {
			std::cerr <<  __FILE__ << ":" << __LINE__ << std::endl;
			abort();
		}

		Mas(*target_data) = Mas(*source_data);
	}

	for (const auto& item: boundaries.get_copy_boundary_cells(V)) {
		const auto
			&target_id = item[0],
			&source_id = item[1];

		auto
			*const target_data = grid[target_id],
			*const source_data = grid[source_id];

		if (target_data == nullptr) {
			std::cerr <<  __FILE__ << ":" << __LINE__ << std::endl;
			abort();
		}
		if (source_data == nullptr) {
			std::cerr <<  __FILE__ << ":" << __LINE__ << std::endl;
			abort();
		}

		Mom(*target_data) = Mom(*source_data);
	}

	for (const auto& item: boundaries.get_copy_boundary_cells(P)) {
		const auto
			&target_id = item[0],
			&source_id = item[1];

		auto
			*target_data = grid[target_id],
			*source_data = grid[source_id];

		if (target_data == nullptr) {
			std::cerr <<  __FILE__ << ":" << __LINE__ << std::endl;
			abort();
		}
		if (source_data == nullptr) {
			std::cerr <<  __FILE__ << ":" << __LINE__ << std::endl;
			abort();
		}

		Nrj(*target_data) = Nrj(*source_data);
	}

	for (const auto& item: boundaries.get_copy_boundary_cells(B)) {
		const auto
			&target_id = item[0],
			&source_id = item[1];

		auto
			*target_data = grid[target_id],
			*source_data = grid[source_id];

		if (target_data == nullptr) {
			std::cerr <<  __FILE__ << ":" << __LINE__ << std::endl;
			abort();
		}
		if (source_data == nullptr) {
			std::cerr <<  __FILE__ << ":" << __LINE__ << std::endl;
			abort();
		}

		Mag(*target_data) = Mag(*source_data);
	}
}


}} // namespaces


#endif // ifndef PAMHD_MHD_BOUNDARIES_HPP
