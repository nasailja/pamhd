/*
Class for handling all copy boundaries of a simulation.

Copyright 2016 Ilja Honkonen
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

#ifndef PAMHD_BOUNDARIES_MULTIVARIABLE_COPY_BOUNDARIES_HPP
#define PAMHD_BOUNDARIES_MULTIVARIABLE_COPY_BOUNDARIES_HPP


#include "array"
#include "stdexcept"
#include "string"
#include "utility"
#include "vector"

#include "rapidjson/document.h"

#include "boundaries/copy_boundaries.hpp"


namespace pamhd {
namespace boundaries {


/*!
Collection of simulation copy boundaries for given variables creatable from json data.
*/
template<
	class Cell_Id,
	class Geometry_Id,
	class... Variables
> class Multivariable_Copy_Boundaries
{};

//! Stops recursion over simulation variables, does nothing.
template<
	class Cell_Id,
	class Geometry_Id
> class Multivariable_Copy_Boundaries<Cell_Id, Geometry_Id>
{
public:
	void set(const rapidjson::Value&) const {}
	void get_geometry_ids() const {}
	void get_copy_sources() const {}
	void push_back_source() const {}
};

template<
	class Cell_Id,
	class Geometry_Id,
	class Current_Variable,
	class... Rest
> class Multivariable_Copy_Boundaries<
	Cell_Id,
	Geometry_Id,
	Current_Variable,
	Rest...
> :
	Multivariable_Copy_Boundaries<Cell_Id, Geometry_Id, Rest...>
{
public:

	// expose inherited versions of functions
	using Multivariable_Copy_Boundaries<Cell_Id, Geometry_Id, Rest...>::set;
	using Multivariable_Copy_Boundaries<Cell_Id, Geometry_Id, Rest...>::get_geometry_ids;
	using Multivariable_Copy_Boundaries<Cell_Id, Geometry_Id, Rest...>::get_copy_sources;
	using Multivariable_Copy_Boundaries<Cell_Id, Geometry_Id, Rest...>::push_back_source;

	void set(const rapidjson::Value& object)
	{
		const std::string name = Current_Variable::get_option_name();
		if (not object.HasMember(name.c_str())) {
			throw std::invalid_argument(
				std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
				+ "Object doesn't have a " + name + " key."
			);
		}

		this->copy_boundaries.set(object[name.c_str()]);
		Multivariable_Copy_Boundaries<Cell_Id, Geometry_Id, Rest...>::set(object);
	}


	//! geometry id for each copy boundary of given simulation variable.
	const std::vector<Geometry_Id>& get_geometry_ids(
		const Current_Variable&
	) const {
		return this->copy_boundaries.geometry_ids;
	}

	/*!
	Returns data sources of copy boundary cells for given simulation variable.

	\see Copy_Boundaries::copy_sources
	*/
	const std::vector<std::array<Cell_Id, 2>>& get_copy_sources(
		const Current_Variable&
	) const {
		return this->copy_boundaries.copy_sources;
	}


	/*!
	Copy_Boundaries::push_back_source() of given simulation variable.
	*/
	void push_back_source(
		const Current_Variable&,
		const std::array<Cell_Id, 2>& source
	) {
		this->copy_boundaries.push_back_source(source);
	}



private:

	Copy_Boundaries<Cell_Id, Geometry_Id, Current_Variable> copy_boundaries;
};

}} // namespaces

#endif // ifndef PAMHD_BOUNDARIES_MULTIVARIABLE_COPY_BOUNDARIES_HPP
