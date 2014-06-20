/*
Class for setting the initial condition of a simulation.

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

#ifndef PAMHD_INITIAL_CONDITION_HPP
#define PAMHD_INITIAL_CONDITION_HPP


#include "boundaries/box.hpp"


namespace pamhd {
namespace boundaries {


template<
	class Cell_T,
	class Vector_T,
	class... Variables
> class Initial_Condition
{
public:

	/*!
	Invalidates previous boundary indices.
	*/
	template<class... Boundary_Data> void add_boundary_box(
		const Vector_T& geometry_start,
		const Vector_T& geometry_end,
		const Boundary_Data&... boundary_data
	) {
		this->boxes.add_boundary(geometry_start, geometry_end, boundary_data...);
	}


	void add_options(
		boost::program_options::options_description& options,
		const std::string& option_name_prefix
	) {
		this->default_.add_options(options, option_name_prefix + "default.");
		this->boxes.add_options(options, option_name_prefix + "box.");
	}


	size_t add_cell(
		const Cell_T& cell,
		const Vector_T& cell_start,
		const Vector_T& cell_end
	) {
		return this->boxes.add_cell(cell, cell_start, cell_end);
	}


	template<
		class Variable
	> const std::vector<
		typename Variable::data_type
	>& get_default_data(
		const Variable& variable
	) {
		return this->default_[variable];
	}


	//! Returns number of existing boundaries
	size_t get_number_of_boundaries() const
	{
		return this->boxes.get_number_of_boundaries();
	}


	//! Returns cells which belong to given boundary
	const std::vector<Cell_T>& get_boundary_cells(
		const size_t boundary_index
	) const {
		return this->boxes.get_boundary_cells(boundary_index);
	}


	//! Returns data of given variable in given boundary
	template<class Variable> const typename Variable::data_type& get_boundary_data(
		const Variable&,
		const size_t boundary_index
	) const {
		return this->boxes.get_boundary_data(Variable(), boundary_index);
	}



private:

	// default initial data of a cell unless it belongs to a boundary
	Variables_To_Options<Variables...> default_;

	Box<Cell_T, Vector_T, Variables...> boxes;
};

}} // namespaces

#endif // ifndef PAMHD_INITIAL_CONDITION_HPP
