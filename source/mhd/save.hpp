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

#include "iomanip"

#include "gensimcell.hpp"

namespace pamhd {
namespace mhd {


class Save
{
public:

	static std::string get_header_template()
	{
		return {"fluxes = "};
	}

	static size_t get_header_size()
	{
		return get_header_template().size() + 2;
	}


	/*!
	Saves the MHD solution into a file with name derived from simulation time.

	file_name_prefix is added to the beginning of the file name.

	MHD_T and subsequent arguments refer to the MHD solution to save.

	The transfer of all first level variables must be switched
	off before this function is called. After save returns the
	transfer of all first level variables is switched off.
	Transfer of variables in MHD_T must be switched on.

	Grid_T is assumed to provide the dccrg API.

	Return true on success, false otherwise.
	*/
	template <
		class Grid_T,
		class Cell_T,
		class MHD_T,
		class MHD_Flux_T
	> static bool save(
		const std::string& file_name_prefix,
		Grid_T& grid,
		const double simulation_time,
		const bool save_fluxes
	) {
		Cell_T::set_transfer_all(true, MHD_T());

		std::string header_data(get_header_template());
		if (save_fluxes) {
			header_data += "y\n";
			Cell_T::set_transfer_all(true, MHD_Flux_T());
		} else {
			header_data += "n\n";
		}
		if (header_data.size() != get_header_size()) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< ": Invalid size for header: " << header_data.size()
				<< ", should be " << get_header_size()
				<< std::endl;
			return false;
		}

		std::tuple<void*, int, MPI_Datatype> header{
			(void*) header_data.data(),
			get_header_size(),
			MPI_BYTE
		};

		std::ostringstream time_string;
		time_string
			<< std::scientific
			<< std::setprecision(3)
			<< simulation_time;

		const bool ret_val = grid.save_grid_data(
			file_name_prefix + "mhd_" + time_string.str() + "_s.dc",
			0,
			header
		);

		if (save_fluxes) {
			Cell_T::set_transfer_all(false, MHD_Flux_T());
		}
		Cell_T::set_transfer_all(false, MHD_T());

		return ret_val;
	}
};

}} // namespaces

#endif // ifndef PAMHD_MHD_SAVE_HPP
