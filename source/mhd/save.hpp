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

	static std::string get_header_string_template()
	{
		return {"fluxes = "};
	}

	static size_t get_header_string_size()
	{
		return get_header_string_template().size() + 2;
	}

	static constexpr size_t nr_header_doubles = 3;


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
		const double adiabatic_index,
		const double proton_mass,
		const double vacuum_permeability,
		const bool save_fluxes
	) {
		std::string header_string(get_header_string_template());
		if (save_fluxes) {
			header_string += "y\n";
			Cell_T::set_transfer_all(true, MHD_Flux_T());
		} else {
			header_string += "n\n";
		}
		if (header_string.size() != get_header_string_size()) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< ": Invalid size for header: " << header_string.size()
				<< ", should be " << get_header_string_size()
				<< std::endl;
			return false;
		}

		// write physical constants
		const std::array<double, nr_header_doubles> header_doubles{
			adiabatic_index,
			proton_mass,
			vacuum_permeability
		};


		std::array<int, 2> counts{
			int(header_string.size()),
			nr_header_doubles
		};
		std::array<MPI_Aint, 2> displacements{
			0,
			reinterpret_cast<char*>(const_cast<double*>(header_doubles.data()))
				- static_cast<const char*>(header_string.data())
		};
		std::array<MPI_Datatype, 2> datatypes{
			MPI_BYTE,
			MPI_DOUBLE
		};

		MPI_Datatype header_datatype;
		if (
			MPI_Type_create_struct(
				2,
				counts.data(),
				displacements.data(),
				datatypes.data(),
				&header_datatype
			) != MPI_SUCCESS
		) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Couldn't create header datatype"
				<< std::endl;
			abort();
		}

		std::tuple<void*, int, MPI_Datatype> header{
			(void*) header_string.data(),
			1,
			header_datatype
		};

		std::ostringstream time_string;
		time_string
			<< std::scientific
			<< std::setprecision(3)
			<< simulation_time;

		Cell_T::set_transfer_all(true, MHD_T());
		if (save_fluxes) {
			Cell_T::set_transfer_all(true, MHD_Flux_T());
		}
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
