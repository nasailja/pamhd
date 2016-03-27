/*
Class representing one boundary of PAMHD.

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


#ifndef PAMHD_BOUNDARIES_BOUNDARY_HPP
#define PAMHD_BOUNDARIES_BOUNDARY_HPP


#include "cstdlib"
#include "iostream"
#include "string"
#include "utility"

#include "boost/program_options.hpp"
#include "prettyprint.hpp"

#include "boundaries/box.hpp"
#include "boundaries/sphere.hpp"
#include "boundaries/variable_to_option.hpp"


namespace pamhd {
namespace boundaries {


template<
	class Geometry_T,
	class Cell_T,
	class... Variables
> class Boundary
{
public:

	// TODO: default constructor

	Geometry_T geometry;

	Variable_To_Option<Variables...> boundary_data;


	void add_options(
		const std::string& option_name_prefix,
		boost::program_options::options_description& options
	) {
		this->geometry.add_options(option_name_prefix, options);
		this->boundary_data.add_options(option_name_prefix, options);
	}


	/*!
	geometry_parameters is given to the overlaps function of Geometry_T.

	Returns true if given cell was added to cell list of this boundary,
	false otherwise.
	*/
	template<class... Geometry_Parameters> bool add_cell(
		const Cell_T& cell,
		Geometry_Parameters&&... geometry_parameters
	) {
		if (
			this->geometry.overlaps(
				std::forward<Geometry_Parameters>(geometry_parameters)...
			)
		) {
			this->cells.push_back(cell);
			return true;
		}

		return false;
	}

	void clear_cells()
	{
		this->cells.clear();
	}

	const std::vector<Cell_T>& get_cells() const
	{
		return this->cells;
	}


private:
	// cells overlapping Geometry_T
	std::vector<Cell_T> cells;
};

}} // namespaces

#endif // ifndef PAMHD_BOUNDARIES_BOUNDARY_HPP
