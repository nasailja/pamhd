/*


Copyright 2015 Ilja Honkonen
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

#ifndef PAMHD_PARTICLE_SOLVE_FIELDS_HPP
#define PAMHD_PARTICLE_SOLVE_FIELDS_HPP


#include "divergence/remove.hpp"


namespace pamhd {
namespace particle {


/*!
Accumulated velocity is saved into reference returned by Bulk_Momentum_Getter.
*/
template<
	class Particle_List_Getter,
	class Particle_Position_Getter,
	class Particle_Mass_Getter,
	class Particle_Momentum_Getter,
	class Bulk_Mass_Getter,
	class Bulk_Momentum_Getter,
	class Bulk_Velocity_Getter,
	class Accu_List_Bulk_Mass_Getter,
	class Accu_List_Bulk_Momentum_Getter,
	class Accu_List_Target_Getter,
	class Accu_List_Length_Getter,
	class Accumulation_List_Getter,
	class Cell,
	class Geometry
> void accumulate_velocity(
	const std::vector<uint64_t>& cell_ids,
	dccrg::Dccrg<Cell, Geometry>& grid,
	const double vacuum_permeability,
	Particle_List_Getter Part_List,
	Particle_Position_Getter Part_Pos,
	Particle_Mass_Getter Part_Mass,
	Particle_Momentum_Getter Part_Mom,
	Bulk_Mass_Getter Bulk_Mass,
	Bulk_Momentum_Getter Bulk_Momentum,
	Bulk_Velocity_Getter Bulk_Velocity,
	Accu_List_Bulk_Mass_Getter List_Bulk_Mass,
	Accu_List_Bulk_Momentum_Getter List_Bulk_Momentum,
	Accu_List_Target_Getter List_Target,
	Accu_List_Length_Getter Accu_List_Len,
	Accumulation_List_Getter Accu_List
) {
	for (const auto& cell: cell_ids) {
		auto* const cell_data = grid[cell];
		if (cell_data == nullptr) {
			std::cerr <<  __FILE__ << "(" << __LINE__ << ")" << std::endl;
			abort();
		}
		Bulk_Mass(*cell_data) = 0;
		Bulk_Momentum(*cell_data) = {0, 0, 0};
	}
	pamhd::particle::accumulate(
		cell_ids,
		grid,
		Part_List,
		Part_Pos,
		Part_Mass,
		Bulk_Mass,
		List_Bulk_Mass,
		List_Target,
		Accu_List_Len,
		Accu_List
	);
	pamhd::particle::accumulate(
		cell_ids,
		grid,
		Part_List,
		Part_Pos,
		Part_Mom,
		Bulk_Momentum,
		List_Bulk_Momentum,
		List_Target,
		Accu_List_Len,
		Accu_List
	);
}


/*
...
Calculates electric field using generalized Ohm's law.

Included terms: E = JxB - VxB.

Getters should return a reference to data of corresponding
variable when given data of one cell.

Given variables should have correct values also in
neighbors of given cells
*/
template<
	class Electric_Field_Getter,
	class Magnetic_Field_Getter,
	class Current_Getter,
	class Particle_List_Getter,
	class Particle_Position_Getter,
	class Particle_Mass_Getter,
	class Particle_Momentum_Getter,
	class Bulk_Mass_Getter,
	class Bulk_Momentum_Getter,
	class Bulk_Velocity_Getter,
	class Accu_List_Bulk_Mass_Getter,
	class Accu_List_Bulk_Momentum_Getter,
	class Accu_List_Target_Getter,
	class Accu_List_Length_Getter,
	class Accumulation_List_Getter,
	class Cell,
	class Geometry
> void calculate_electric_field(
	const std::vector<uint64_t>& cell_ids,
	dccrg::Dccrg<Cell, Geometry>& grid,
	const double vacuum_permeability,
	Electric_Field_Getter Electric_Field,
	Magnetic_Field_Getter Magnetic_Field,
	Current_Getter Current,
	Particle_List_Getter Part_List,
	Particle_Position_Getter Part_Pos,
	Particle_Mass_Getter Part_Mass,
	Particle_Momentum_Getter Part_Mom,
	Bulk_Mass_Getter Bulk_Mass,
	Bulk_Momentum_Getter Bulk_Momentum,
	Bulk_Velocity_Getter Bulk_Velocity,
	Accu_List_Bulk_Mass_Getter List_Bulk_Mass,
	Accu_List_Bulk_Momentum_Getter List_Bulk_Momentum,
	Accu_List_Target_Getter List_Target,
	Accu_List_Length_Getter Accu_List_Len,
	Accumulation_List_Getter Accu_List
) {
	pamhd::divergence::get_curl(
		cell_ids,
		grid,
		Magnetic_Field,
		Current
	);
	for (const auto& cell: cell_ids) {
		auto* const cell_data = grid[cell];
		if (cell_data == nullptr) {
			std::cerr <<  __FILE__ << "(" << __LINE__ << ")" << std::endl;
			abort();
		}
		Current(*cell_data) /= vacuum_permeability;
	}

	accumulate_velocity(
		cell_ids,
		grid,
		vacuum_permeability,
		Part_List,
		Part_Pos,
		Part_Mass,
		Part_Mom,
		Bulk_Mass,
		Bulk_Momentum,
		Bulk_Velocity,
		List_Bulk_Mass,
		List_Bulk_Momentum,
		List_Target,
		Accu_List_Len,
		Accu_List
	);

	for (const auto& cell: cell_ids) {
		auto* const cell_data = grid[cell];
		if (cell_data == nullptr) {
			std::cerr <<  __FILE__ << "(" << __LINE__ << ")" << std::endl;
			abort();
		}
		Bulk_Velocity(*cell_data) = Bulk_Momentum(*cell_data) / Bulk_Mass(*cell_data);
	}

	for (const auto& cell: cell_ids) {
		auto* const cell_data = grid[cell];
		if (cell_data == nullptr) {
			std::cerr <<  __FILE__ << "(" << __LINE__ << ")" << std::endl;
			abort();
		}
		Electric_Field(*cell_data)
			= (Current(*cell_data) - Bulk_Velocity(*cell_data))
			.cross(Magnetic_Field(*cell_data));
	}
}


/*!
If cell doesn't have particles and no_particles_allowed == true
skips that cell.
*/
template<
	class Particle_Mass_T,
	class Particle_Velocity_T,
	class Particle_Species_Mass_T,
	class Particle_Bulk_Mass_Getter,
	class Particle_Bulk_Momentum_Getter,
	class Particle_List_Getter,
	class MHD_Mass_Getter,
	class MHD_Momentum_Getter,
	class MHD_Energy_Getter,
	class MHD_Magnetic_Field_Getter,
	class Cell,
	class Geometry
> void fill_mhd_fluid_values(
	const std::vector<uint64_t>& cell_ids,
	dccrg::Dccrg<Cell, Geometry>& grid,
	const double adiabatic_index,
	const double vacuum_permeability,
	const double particle_temp_nrj_ratio,
	Particle_Bulk_Mass_Getter Particle_Bulk_Mass,
	Particle_Bulk_Momentum_Getter Particle_Bulk_Momentum,
	Particle_List_Getter Particle_List,
	MHD_Mass_Getter MHD_Mass,
	MHD_Momentum_Getter MHD_Momentum,
	MHD_Energy_Getter MHD_Energy,
	MHD_Magnetic_Field_Getter MHD_Magnetic_Field,
	const bool no_particles_allowed = false
) {
	for (const auto& cell: cell_ids) {
		auto* const cell_data = grid[cell];
		if (cell_data == nullptr) {
			std::cerr <<  __FILE__ << "(" << __LINE__ << ")" << std::endl;
			abort();
		}

		if (Particle_List(*cell_data).size() == 0) {
			if (no_particles_allowed) {
				continue;
			}

			std::cerr <<  __FILE__ << "(" << __LINE__ << "): "
				<< "No particles in cell " << cell
				<< std::endl;
			abort();
		}

		const auto length = grid.geometry.get_length(cell);
		const auto volume = length[0] * length[1] * length[2];

		MHD_Mass(*cell_data) = Particle_Bulk_Mass(*cell_data) / volume;
		MHD_Momentum(*cell_data) = Particle_Bulk_Momentum(*cell_data) / volume;

		const auto pressure
			= pamhd::particle::get_pressure<
				Particle_Mass_T,
				Particle_Velocity_T,
				Particle_Species_Mass_T
			>(
				Particle_List(*cell_data),
				particle_temp_nrj_ratio,
				volume
			);

		MHD_Energy(*cell_data)
			= pressure / (adiabatic_index - 1)
			+ MHD_Momentum(*cell_data).squaredNorm() / MHD_Mass(*cell_data) / 2
			+ MHD_Magnetic_Field(*cell_data).squaredNorm() / vacuum_permeability / 2;
	}
}


}} // namespaces

#endif
