/*
Saves the MHD solution of PAMHD.

Copyright 2014 Ilja Honkonen
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

#ifndef PAMHD_MHD_SAVE_HPP
#define PAMHD_MHD_SAVE_HPP

#include "cmath"
#include "limits"

#include "gensimcell.hpp"

namespace pamhd {
namespace mhd {


/*!
Saves the MHD solution into a file with name derived from simulation time. 

MHD_T and subsequent arguments refer to the MHD solution to save.

The transfer of all first level variables must be switched
off before this function is called. After save returns the
transfer of all first level variables is switched off.
Transfer of variables in MHD_T must be switched on.

Grid_T is assumed to provide the dccrg API.
*/
template <
	class Grid_T,
	class Cell_T,
	class MHD_T
> void save(
	Grid_T& grid,
	const double simulation_time
) {
	Cell_T::set_transfer_all(true, MHD_T());

	// get the file name
	std::ostringstream time_string;
	time_string
		<< std::scientific
		<< std::setprecision(6)
		<< simulation_time;

	// use an empty header with a sane address
	char dummy;
	std::tuple<void*, int, MPI_Datatype> header{(void*) &dummy, 0, MPI_BYTE};

	grid.save_grid_data(
		"mhd_" + time_string.str() + "_s.dc",
		0,
		header
	);

	Cell_T::set_transfer_all(false, MHD_T());
}

}} // namespaces

#endif // ifndef PAMHD_MHD_SAVE_HPP
