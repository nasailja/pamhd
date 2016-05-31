/*
Class for handling all initial conditions of a simulation.

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

#ifndef PAMHD_BOUNDARIES_MULTIVARIABLE_INITIAL_CONDITIONS_HPP
#define PAMHD_BOUNDARIES_MULTIVARIABLE_INITIAL_CONDITIONS_HPP


#include "algorithm"
#include "array"
#include "stdexcept"
#include "string"
#include "type_traits"
#include "utility"

#include "rapidjson/document.h"

#include "boundaries/initial_conditions.hpp"


namespace pamhd {
namespace boundaries {


template<class Geometry_Id, class... Variables> class Multivariable_Initial_Conditions;

template<class Geometry_Id> class Multivariable_Initial_Conditions<Geometry_Id>
{
public:
	void set(const rapidjson::Value&) {}
	void get_default_data() {}
	void get_number_of_regions() {}
	void get_initial_condition() {}
	void get_data() {}
};

template<
	class Geometry_Id,
	class Current_Variable,
	class... Rest
> class Multivariable_Initial_Conditions<
	Geometry_Id,
	Current_Variable,
	Rest...
> :
	public Multivariable_Initial_Conditions<Geometry_Id, Rest...>
{
public:

	// expose inherited versions of functions
	using Multivariable_Initial_Conditions<Geometry_Id, Rest...>::get_default_data;
	using Multivariable_Initial_Conditions<Geometry_Id, Rest...>::get_number_of_regions;
	using Multivariable_Initial_Conditions<Geometry_Id, Rest...>::get_initial_condition;
	using Multivariable_Initial_Conditions<Geometry_Id, Rest...>::get_data;


	void set(const rapidjson::Value& object)
	{
		const std::string name = Current_Variable::get_option_name();
		if (not object.HasMember(name.c_str())) {
			throw std::invalid_argument(
				std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
				+ "Object doesn't have a " + name + " key."
			);
		}

		this->initial_conditions.set(object[name.c_str()]);
		Multivariable_Initial_Conditions<Geometry_Id, Rest...>::set(object);
	}


	typename Current_Variable::data_type get_default_data(
		const Current_Variable&,
		const double& t,
		const double& x,
		const double& y,
		const double& z,
		const double& radius,
		const double& latitude,
		const double& longitude
	) {
		return this->initial_conditions.get_default_data(
			t, x, y, z, radius, latitude, longitude
		);
	}


	size_t get_number_of_regions(const Current_Variable&) const
	{
		return this->initial_conditions.get_number_of_regions();
	}

	/*!
	Region id is in range 0..number of regions - 1
	*/
	const Initial_Condition<
		Geometry_Id, Current_Variable
	>& get_initial_condition(
		const Current_Variable&,
		const size_t& region_id
	) const {
		return this->initial_conditions.get_initial_condition(region_id);
	}


	typename Current_Variable::data_type get_data(
		const Current_Variable&,
		const size_t& region_index,
		const double& t,
		const double& x,
		const double& y,
		const double& z,
		const double& radius,
		const double& latitude,
		const double& longitude
	) {
		return this->initial_conditions.get_data(
			region_index, t, x, y, z, radius, latitude, longitude
		);
	}


private:

	Initial_Conditions<Geometry_Id, Current_Variable> initial_conditions;

};

}} // namespaces

#endif // ifndef PAMHD_BOUNDARIES_MULTIVARIABLE_INITIAL_CONDITIONS_HPP
