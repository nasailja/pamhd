/*
Class for handling all value boundaries of a simulation.

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

#ifndef PAMHD_BOUNDARIES_MULTIVARIABLE_VALUE_BOUNDARIES_HPP
#define PAMHD_BOUNDARIES_MULTIVARIABLE_VALUE_BOUNDARIES_HPP


#include "stdexcept"
#include "string"
#include "type_traits"

#include "rapidjson/document.h"

#include "boundaries/value_boundaries.hpp"


namespace pamhd {
namespace boundaries {


/*!
Collection of value boundaries creatable json data.
*/
template<class Geometry_Id, class... Variables> class Multivariable_Value_Boundaries {};

template<class Geometry_Id> class Multivariable_Value_Boundaries<Geometry_Id>
{
public:
	void set(const rapidjson::Value&) const {}
	void get_number_of_boundaries() const {}
	void get_value_boundary() const {}
	void get_data() const {}
};

template<
	class Geometry_Id,
	class Current_Variable,
	class... Rest
> class Multivariable_Value_Boundaries<
	Geometry_Id,
	Current_Variable,
	Rest...
> :
	public Multivariable_Value_Boundaries<Geometry_Id, Rest...>
{
public:

	// expose inherited versions of functions
	using Multivariable_Value_Boundaries<Geometry_Id, Rest...>::set;
	using Multivariable_Value_Boundaries<Geometry_Id, Rest...>::get_number_of_boundaries;
	using Multivariable_Value_Boundaries<Geometry_Id, Rest...>::get_value_boundary;
	using Multivariable_Value_Boundaries<Geometry_Id, Rest...>::get_data;


	void set(const rapidjson::Value& object)
	{
		const std::string name = Current_Variable::get_option_name();
		if (not object.HasMember(name.c_str())) {
			throw std::invalid_argument(
				std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
				+ "Object doesn't have a " + name + " key."
			);
		}

		this->value_boundaries.set(object[name.c_str()]);
		Multivariable_Value_Boundaries<Geometry_Id, Rest...>::set(object);
	}


	size_t get_number_of_boundaries(const Current_Variable&) const
	{
		return this->value_boundaries.get_number_of_boundaries();
	}

	/*!
	Boundary id is in range 0..number of boundaries - 1
	*/
	const Value_Boundary<
		Geometry_Id,
		Current_Variable
	>& get_value_boundary(
		const Current_Variable&,
		const size_t& boundary_id) const
	{
		return this->value_boundaries.get_value_boundary(boundary_id);
	}


	typename Current_Variable::data_type get_data(
		const Current_Variable&,
		const size_t& boundary_index,
		const double& t,
		const double& x,
		const double& y,
		const double& z,
		const double& radius,
		const double& latitude,
		const double& longitude
	) {
		return this->value_boundaries.get_data(
			boundary_index, t, x, y, z, radius, latitude, longitude
		);
	}


private:

	Value_Boundaries<Geometry_Id, Current_Variable> value_boundaries;

};

}} // namespaces

#endif // ifndef PAMHD_BOUNDARIES_MULTIVARIABLE_VALUE_BOUNDARIES_HPP
