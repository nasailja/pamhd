/*
Class for handling grid options.

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

#ifndef PAMHD_GRID_OPTIONS_HPP
#define PAMHD_GRID_OPTIONS_HPP

#include "array"

#include "boundaries/variable_to_option.hpp"


namespace pamhd {
namespace grid {

struct Number_Of_Cells {
	using data_type = std::array<int, 3>;

	static const std::string get_name()
	{
		return {"number of cells"};
	}
	static const std::string get_option_name()
	{
		return {"nr-cells"};
	}
	static const std::string get_option_help()
	{
		return {"Number of grid cells of refinement level 0"};
	}
};

struct Volume {
	using data_type = std::array<double, 3>;

	static const std::string get_name()
	{
		return {"volume"};
	}
	static const std::string get_option_name()
	{
		return {"volume"};
	}
	static const std::string get_option_help()
	{
		return {"Volume of simulation grid"};
	}
};

struct Start {
	using data_type = std::array<double, 3>;

	static const std::string get_name()
	{
		return {"start"};
	}
	static const std::string get_option_name()
	{
		return {"start"};
	}
	static const std::string get_option_help()
	{
		return {"Starting coordinate of simulation grid"};
	}
};

struct Periodic {
	using data_type = std::array<bool, 3>;

	static const std::string get_name()
	{
		return {"periodic"};
	}
	static const std::string get_option_name()
	{
		return {"periodic"};
	}
	static const std::string get_option_help()
	{
		return {"Whether grid is periodic (true) or not (false)"};
	}
};

class Options
{
public:

	void add_options(
		const std::string& option_name_prefix,
		boost::program_options::options_description& options
	) {
		this->data.add_options(
			option_name_prefix,
			options
		);
	}

	boundaries::Variable_To_Option<
		Number_Of_Cells,
		Volume,
		Start,
		Periodic
	> data;
};

}} // namespaces

#endif // ifndef PAMHD_GRID_OPTIONS_HPP
