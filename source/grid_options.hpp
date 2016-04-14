/*
Class for handling grid options.

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

#ifndef PAMHD_GRID_OPTIONS_HPP
#define PAMHD_GRID_OPTIONS_HPP


#include "array"
#include "cstdint"
#include "stdexcept"

#include "rapidjson/document.h"

#include "boundaries/simulation_variable_expression.hpp"


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
private:

	boundaries::Simulation_Variable_Expression<
		Number_Of_Cells,
		Volume,
		Start,
		Periodic
	> data;


public:

	Options()
	{
		data.set_expression(Number_Of_Cells(), "{1, 1, 1}");
		data.set_expression(Volume(), "{1, 1, 1}");
		data.set_expression(Start(), "{0, 0, 0}");
		data.set_expression(Periodic(), "{false, false, false}");
	}


	template<class Variable_T> typename Variable_T::data_type get_data(const Variable_T& v) {
		return this->data.get_data(v);
	}


	// muparserx doesn't support uint64_t so convert from int
	std::array<uint64_t, 3> get_data(const Number_Of_Cells& v) {
		auto tmp = this->data.get_data(v);
		return {uint64_t(tmp[0]), uint64_t(tmp[1]), uint64_t(tmp[2])};
	}


	void set_expressions(const rapidjson::Document& document)
	{
		if (not document.IsObject()) {
			return;
		}

		const auto grid_iter = document.FindMember("grid");
		if (grid_iter == document.MemberEnd()) {
			return;
		}

		const auto& grid = document["grid"];

		const auto periodic = grid.FindMember(
			Periodic::get_option_name().c_str()
		);
		if (periodic != document.MemberEnd()) {
			if (periodic->value.IsString()) {
				data.set_expression(Periodic(), periodic->value.GetString());
			} else {
				throw std::invalid_argument(
					"Value of " + Periodic::get_option_name() + " variable is not a string."
				);
			}
		}

		const auto nr_cells = grid.FindMember(
			Number_Of_Cells::get_option_name().c_str()
		);
		if (nr_cells != document.MemberEnd()) {
			if (nr_cells->value.IsString()) {
				data.set_expression(Number_Of_Cells(), nr_cells->value.GetString());
			} else {
				throw std::invalid_argument(
					"Value of " + Number_Of_Cells::get_option_name() + " variable is not a string."
				);
			}
		}

		const auto volume = grid.FindMember(
			Volume::get_option_name().c_str()
		);
		if (volume != document.MemberEnd()) {
			if (volume->value.IsString()) {
				data.set_expression(
					Volume(),
					volume->value.GetString()
				);
			} else {
				throw std::invalid_argument(
					"Value of " + Volume::get_option_name() + " variable is not a string."
				);
			}
		}

		const auto start = grid.FindMember(
			Start::get_option_name().c_str()
		);
		if (start != document.MemberEnd()) {
			if (start->value.IsString()) {
				data.set_expression(
					Start(),
					start->value.GetString()
				);
			} else {
				throw std::invalid_argument(
					"Value of " + Start::get_option_name() + " variable is not a string."
				);
			}
		}
	}

};

}} // namespaces

#endif // ifndef PAMHD_GRID_OPTIONS_HPP
