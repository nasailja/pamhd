/*
Saves particle solution of PAMHD.

Copyright 2015 Ilja Honkonen
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice, this
  list of conditions and the following disclaimer in the documentation and/or
  other materials provided with the distribution.

* Neither the name of copyright holders nor the names of their contributors
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

#ifndef PAMHD_PARTICLE_SAVE_HPP
#define PAMHD_PARTICLE_SAVE_HPP


#include "iomanip"

#include "dccrg.hpp"
#include "mpi.h"

#include "particle/variables.hpp"


namespace pamhd {
namespace particle {


/*!
Saves the particle solution into a file with name derived from simulation time.

file_name_prefix is added to the beginning of the file name.

The transfer of all first level variables must be switched
off before this function is called. After save returns the
transfer of all first level variables is switched off.

Cell must be compatible with gensimcell (github.com/nasailja/gensimcell)

Return true on success, false otherwise.
*/
template <
	class Electric_Field_T,
	class Magnetic_Field_T,
	class Nr_Particles_T,
	class Particles_T,
	class Cell,
	class Geometry
> bool save(
	const std::string& file_name_prefix,
	dccrg::Dccrg<Cell, Geometry>& grid,
	const double simulation_time,
	const double vacuum_permeability
) {
	std::tuple<void*, int, MPI_Datatype> header{
		(void*) &vacuum_permeability,
		1,
		MPI_DOUBLE
	};

	std::ostringstream time_string;
	time_string
		<< std::scientific
		<< std::setprecision(3)
		<< simulation_time;

	// update number of internal particles
	for (const auto& cell_id: grid.get_cells()) {
		auto* const cell_data = grid[cell_id];
		if (cell_data == nullptr) {
			std::cerr << __FILE__ << "(" << __LINE__ << ")" << std::endl;
			abort();
		}
		(*cell_data)[Nr_Particles_T()]
			= (*cell_data)[Particles_T()].size();
	}

	Cell::set_transfer_all(
		true,
		Electric_Field_T(),
		Magnetic_Field_T(),
		Nr_Particles_T(),
		Particles_T()
	);
	const bool ret_val = grid.save_grid_data(
		file_name_prefix + "particle_" + time_string.str() + "_s.dc",
		0,
		header
	);
	Cell::set_transfer_all(
		false,
		Electric_Field_T(),
		Magnetic_Field_T(),
		Nr_Particles_T(),
		Particles_T()
	);

	return ret_val;
}


}} // namespaces

#endif // ifndef PAMHD_PARTICLE_SAVE_HPP
