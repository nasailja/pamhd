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
				MHD_Mass(*cell_data) = 0;
				MHD_Momentum(*cell_data) = {0, 0, 0};
				MHD_Energy(*cell_data) = 0;
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
		if (MHD_Mass(*cell_data) <= 0) {
			MHD_Momentum(*cell_data) = {0, 0, 0};
			MHD_Energy(*cell_data) = 0;
			continue;
		}

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
