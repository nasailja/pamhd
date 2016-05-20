/*
Class for setting the initial condition of a simulation.

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

#ifndef PAMHD_BOUNDARIES_INITIAL_CONDITION_HPP
#define PAMHD_BOUNDARIES_INITIAL_CONDITION_HPP


#include "algorithm"
#include "array"
#include "iterator"
#include "stdexcept"
#include "string"
#include "type_traits"
#include "utility"

#include "rapidjson/document.h"
#include "mpParser.h"


namespace pamhd {
namespace boundaries {


/*!
Collection of boundaries creatable from arguments given on command line.
*/
template<
	class Geometry_Id,
	class Variable
> class Initial_Condition
{
public:


	Initial_Condition() :
		t_var(&t_val),
		x_var(&x_val), y_var(&y_val), z_var(&z_val),
		radius_var(&radius_val), lat_var(&lat_val), lon_var(&lon_val)
	{
		this->parser.DefineVar("t", this->t_var);
		this->parser.DefineVar("x", this->x_var);
		this->parser.DefineVar("y", this->y_var);
		this->parser.DefineVar("z", this->z_var);
		this->parser.DefineVar("radius", this->radius_var);
		this->parser.DefineVar("lat", this->lat_var);
		this->parser.DefineVar("lon", this->lon_var);
	}

	Initial_Condition(const Initial_Condition& other) = delete;
	Initial_Condition(Initial_Condition&& other) = delete;


	/*!
	Prepares initial condition from given rapidjson object.

	Given object must have "geometry_id" key with integer
	value and "value"
	key with either integer, string or object as value. If
	object then it must have keys denoting coordinate planes
	where data points are located and a "data" key.
	Coordinates and data are arrays of numbers.
	Coordinate keys are either "x", "y" and "z" for
	cartesian geometry or "radius", "lat" and "lon" for
	spherical data. "data" is in order
	(x_1, y_1, z_1), (x_2, y_1, z_1), ...

	Geometry ids refer to geometries of
	source/boundaries/geometries.hpp

	Example jsons which could be given as object:
	\verbatim
	{"geometry_id": 0, "value": 123}
	\endverbatim
	\verbatim
	{"geometry_id": 1, "value": "sin(t)"}
	\endverbatim
	\verbatim
	{
		"geometry_id": 0,
		"value": {
			"x": [1, 2, 3, 4],
			"y": [-3, 2, 5, 8],
			"z": [0.1, 0.2, 40, 123],
			"data": [1, 2, 3, ..., 64]
		}
	}
	\endverbatim
	*/
	void set(const rapidjson::Value& object)
	{
		if (not object.HasMember("geometry_id")) {
			throw std::invalid_argument(__FILE__ ": object doesn't have a geometry_id key.");
		}
		const auto& geometry_id = object["geometry_id"];

		if (not geometry_id.IsUint()) {
			throw std::invalid_argument(__FILE__ ": geometry_id is not unsigned int.");
		}

		if (not object.HasMember("value")) {
			throw std::invalid_argument(__FILE__ ": object doesn't have a value key.");
		}
		const auto& value = object["value"];

		// one scalar or vector value
		if (value.IsNumber() or value.IsArray()) {

			this->number_value = this->get_json_value(value);
			this->init_cond_type = 1;

		// math expression
		} else if (value.IsString()) {

			this->parser.SetExpr(value.GetString());
			this->init_cond_type = 2;

		// 3d grid of scalar values
		} else if (value.IsObject()) {

			if (value.HasMember("x")) {
				this->coordinate_type = 1;
			} else if (value.HasMember("radius")) {
				this->coordinate_type = 2;
			} else {
				throw std::invalid_argument(__FILE__ ": value object doesn't have either x or radius key.");
			}

			std::array<const char*, 3> components{"x", "y", "z"};
			if (this->coordinate_type == 2) {
				components[0] = "radius";
				components[1] = "lat";
				components[2] = "lon";
			}

			// load coordinates of data points
			size_t component_i = 0;
			for (const auto& component: components) {
				if (not value.HasMember(component)) {
					throw std::invalid_argument(__FILE__ ": value object doesn't have required key.");
				}

				const auto& coord = value[component];
				if (not coord.IsArray()) {
					throw std::invalid_argument(__FILE__ ": coordinate isn't array.");
				}

				this->coordinates[component_i].clear();
				for (auto i = coord.Begin(); i != coord.End(); i++) {
					this->coordinates[component_i].push_back(i->GetDouble());
				}

				if (
					not std::is_sorted(
						this->coordinates[component_i].cbegin(), this->coordinates[component_i].cend()
					)
				) {
					throw std::invalid_argument(__FILE__ ": coordinates aren't sorted in non-descending order.");
				}

				component_i++;
			}

			// load data
			if (value.HasMember("data")) {
				const auto& data = value["data"];

				if (not data.IsArray()) {
					throw std::invalid_argument(__FILE__ ": data isn't array.");
				}

				this->data.clear();
				for (auto i = data.Begin(); i != data.End(); i++) {
					this->data.push_back(this->get_json_value(*i));
				}
			} else {
				throw std::invalid_argument(__FILE__ ": value object doesn't have a x or radius key.");
			}

			if (this->coordinates[0].size() * this->coordinates[1].size() * this->coordinates[2].size() != this->data.size()) {
				throw std::invalid_argument(__FILE__ ": value's data has invalid size.");
			}

			this->init_cond_type = 3;

		} else {
			throw std::invalid_argument(__FILE__ ": value's type isn't supported.");
		}
	}


	typename Variable::data_type get_data(
		const double& t,
		const double& x,
		const double& y,
		const double& z,
		const double& radius,
		const double& latitude,
		const double& longitude
	) {
		switch (this->init_cond_type) {
		case 1:
			return this->number_value;

		case 2:
			this->t_val = t;
			this->x_val = x;
			this->y_val = y;
			this->z_val = z;
			this->radius_val = radius;
			this->lat_val = latitude;
			this->lon_val = longitude;
			return this->get_parsed_value();

		// nearest neighbor interpolation in point cloud
		case 3: {
			std::array<double, 3> coord{x, y, z};

			switch (this->coordinate_type) {
			case 1:
				break;
			case 2:
				coord[0] = radius;
				coord[1] = latitude;
				coord[2] = longitude;
				break;
			default:
				throw std::out_of_range(__FILE__ ": Invalid coordinate type.");
			}

			std::array<size_t, 3> coordinate_index{0, 0, 0};
			for (size_t i = 0; i < coordinates.size(); i++) {
				// location of first value not smaller than data point coordinate
				const auto after = std::lower_bound(
					coordinates[i].cbegin(),
					coordinates[i].cend(),
					coord[i]
				);

				if (after == coordinates[i].cbegin()) {
					continue;
				} else if (after == coordinates[i].cend()) {
					coordinate_index[i] = coordinates[i].size() - 1;
				} else {
					const auto before = after - 1;
					if (std::abs(coord[i] - *before) <= std::abs(coord[i] - *after)) {
						coordinate_index[i] = std::distance(coordinates[i].cbegin(), before);
					} else {
						coordinate_index[i] = std::distance(coordinates[i].cbegin(), after);
					}
				}
			}

			const size_t data_i
				= coordinate_index[0]
				+ coordinate_index[1] * coordinates[0].size()
				+ coordinate_index[2] * coordinates[0].size() * coordinates[1].size();

			return this->data[data_i]; }

		default:
			throw std::out_of_range(__FILE__ ": Invalid initial condition type.");
		}
	}


	void set_geometry_id(const Geometry_Id& id)
	{
		this->geometry_id = id;
	}


private:

	// region where this initial condition is applied
	Geometry_Id geometry_id;

	// initial condition if a number was given in json file
	typename Variable::data_type number_value;

	// initial conditino if a string was given in json file
	mup::ParserX parser = mup::ParserX(mup::pckCOMMON | mup::pckNON_COMPLEX | mup::pckMATRIX | mup::pckUNIT);
	mup::Value t_val, x_val, y_val, z_val, radius_val, lat_val, lon_val;
	mup::Variable t_var, x_var, y_var, z_var, radius_var, lat_var, lon_var;

	// initial condition if an object was given in json file
	std::array<std::vector<double>, 3> coordinates; // e.g. x,y,z coordinates of points
	int coordinate_type = -1; // cartesian == 1, geographic == 2

	// data in order x[0],y[0],z[0]; x[1],y[0],z[0]; ...
	std::vector<typename Variable::data_type> data;

	int init_cond_type = -1; // number == 1, string == 2, object == 3


	/*!
	Returns value of type Variable::data_type when it's signed integral.

	Reads a signed integral number from given rapidjson object, throws
	if given object contains a value of different type.

	Object is interpreted as signed 64 bit integer, if it doesn't fit
	into Variable::data_type this function should be specialized for
	smaller integers (use https://github.com/nasailja/pamhd/pulls or
	submit https://github.com/nasailja/pamhd/issues).

	Example:
	\verbatim
	// Variable = struct {using data_type = int};
	const char* json = "[123]"
	rapidjson::Document doc; doc.Parse(json);
	typename Variable::data_type var = this->get_json_value(doc[0]);
	\endverbatim
	*/
	template<class V = Variable> typename std::enable_if<
		std::is_integral<typename V::data_type>::value
			&& std::is_signed<typename V::data_type>::value,
		typename V::data_type
	>::type get_json_value(const rapidjson::Value& object)
	{
		if (not object.IsInt64()) {
			throw std::invalid_argument(__FILE__ ": object isn't of signed integral type.");
		}
		return object.GetInt64();
	}

	//! returns value of type Variable::data_type when it's unsigned integral.
	template<class V = Variable> typename std::enable_if<
		std::is_integral<typename V::data_type>::value
			&& !std::is_signed<typename V::data_type>::value,
		typename V::data_type
	>::type get_json_value(const rapidjson::Value& object)
	{
		if (not object.IsUint64()) {
			throw std::invalid_argument(__FILE__ ": object isn't of unsigned integral type.");
		}
		return object.GetUint64();
	}

	//! returns value of type Variable::data_type when it's floating point.
	template<class V = Variable> typename std::enable_if<
		std::is_floating_point<typename V::data_type>::value,
		typename V::data_type
	>::type get_json_value(const rapidjson::Value& object)
	{
		if (not object.IsNumber()) {
			throw std::invalid_argument(__FILE__ ": object isn't of floating point type.");
		}
		return object.GetDouble();
	}

	//! returns value of type Variable::data_type when it's array of signed integral values.
	template<class V = Variable> typename std::enable_if<
		std::tuple_size<typename V::data_type>::value >= 0
			&& std::is_integral<typename V::data_type::value_type>::value
			&& std::is_signed<typename V::data_type::value_type>::value,
		typename V::data_type
	>::type get_json_value(const rapidjson::Value& object)
	{
		constexpr auto size = std::tuple_size<typename V::data_type>::value;

		if (not object.IsArray()) {
			throw std::invalid_argument(__FILE__ ": object isn't an array.");
		}
		if (object.Size() != size) {
			throw std::invalid_argument(__FILE__ ": array has wrong size.");
		}

		typename Variable::data_type ret_val;
		for (size_t i = 0; i < size; i++) {
			if (not object[i].IsInt64()) {
				throw std::invalid_argument(__FILE__ ": element of array isn't signed integral number.");
			}
			ret_val[i] = object[i].GetInt64();
		}
		return ret_val;
	}

	//! returns value of type Variable::data_type when it's array of unsigned integral values.
	template<class V = Variable> typename std::enable_if<
		std::tuple_size<typename V::data_type>::value >= 0
			&& std::is_integral<typename V::data_type::value_type>::value
			&& !std::is_signed<typename V::data_type::value_type>::value,
		typename V::data_type
	>::type get_json_value(const rapidjson::Value& object)
	{
		constexpr auto size = std::tuple_size<typename V::data_type>::value;

		if (not object.IsArray()) {
			throw std::invalid_argument(__FILE__ ": object isn't an array.");
		}
		if (object.Size() != size) {
			throw std::invalid_argument(__FILE__ ": array has wrong size.");
		}

		typename Variable::data_type ret_val;
		for (size_t i = 0; i < size; i++) {
			if (not object[i].IsUint64()) {
				throw std::invalid_argument(__FILE__ ": element of array isn't unsigned integral number.");
			}
			ret_val[i] = object[i].GetUint64();
		}
		return ret_val;
	}

	//! returns value of type Variable::data_type when it's array of floating point values.
	template<class V = Variable> typename std::enable_if<
		std::tuple_size<typename V::data_type>::value >= 0
			&& std::is_floating_point<typename V::data_type::value_type>::value,
		typename V::data_type
	>::type get_json_value(const rapidjson::Value& object)
	{
		constexpr auto size = std::tuple_size<typename V::data_type>::value;

		if (not object.IsArray()) {
			throw std::invalid_argument(__FILE__ ": object isn't an array.");
		}
		if (object.Size() != size) {
			throw std::invalid_argument(__FILE__ ": array has wrong size.");
		}

		typename Variable::data_type ret_val;
		for (size_t i = 0; i < size; i++) {
			if (not object[i].IsNumber()) {
				throw std::invalid_argument(__FILE__ ": element of array isn't a number.");
			}
			ret_val[i] = object[i].GetDouble();
		}
		return ret_val;
	}


	#ifdef EIGEN_WORLD_VERSION

	//! returns value of type Variable::data_type when it's Eigen::Vector of floating point values with compile time size.
	template<class V = Variable> typename std::enable_if<
		V::data_type::ColsAtCompileTime == 1
			&& V::data_type::RowsAtCompileTime >= 1
			&& std::is_floating_point<typename V::data_type::value_type>::value,
		typename V::data_type
	>::type get_json_value(const rapidjson::Value& object)
	{
		constexpr auto size = V::data_type::RowsAtCompileTime;

		if (not object.IsArray()) {
			throw std::invalid_argument(__FILE__ ": object isn't an array.");
		}
		if (object.Size() != size) {
			throw std::invalid_argument(__FILE__ ": array has wrong size.");
		}

		typename Variable::data_type ret_val;
		for (size_t i = 0; i < size; i++) {
			if (not object[i].IsNumber()) {
				throw std::invalid_argument(__FILE__ ": element of array isn't a number.");
			}
			ret_val(i) = object[i].GetDouble();
		}
		return ret_val;
	}

	#endif


	//! used to evaluate math expression if Variable::data_type is integral
	template<class V = Variable> typename std::enable_if<
		std::is_integral<typename V::data_type>::value,
		typename V::data_type
	>::type get_parsed_value()
	{
		try {
			const auto& evaluated = this->parser.Eval();
			if (evaluated.GetType() != 'i' and evaluated.GetType() != 'f') {
				throw std::invalid_argument(
					std::string("Expression \"")
					+ this->parser.GetExpr()
					+ std::string("\" is neither integer nor floating point.")
				);
			}

			return evaluated.GetInteger();

		} catch(mup::ParserError &error) {
			std::cerr <<  __FILE__ << "(" << __LINE__<< ") "
				<< "Could not evaluate expression \"" << this->parser.GetExpr()
				<< "\" for variable " << Variable::get_name()
				<< ": " << error.GetMsg()
				<< std::endl;
			throw;
		}
	}


	//! used to evaluate math expression if Variable::data_type is floating point
	template<class V = Variable> typename std::enable_if<
		std::is_floating_point<typename V::data_type>::value,
		typename V::data_type
	>::type get_parsed_value()
	{
		try {
			const auto& evaluated = this->parser.Eval();
			if (evaluated.GetType() != 'i' and evaluated.GetType() != 'f') {
				throw std::invalid_argument(
					std::string("Expression \"")
					+ this->parser.GetExpr()
					+ std::string("\" is neither integer nor floating point.")
				);
			}

			return evaluated.GetFloat();

		} catch(mup::ParserError &error) {
			std::cerr <<  __FILE__ << "(" << __LINE__<< ") "
				<< "Could not evaluate expression \"" << this->parser.GetExpr()
				<< "\" for variable " << Variable::get_name()
				<< ": " << error.GetMsg()
				<< std::endl;
			throw;
		}
	}


	//! used to evaluate math expression if Variable::data_type is array of floating point numbers
	template<class V = Variable> typename std::enable_if<
		std::tuple_size<typename V::data_type>::value >= 0
			&& std::is_floating_point<typename V::data_type::value_type>::value,
		typename V::data_type
	>::type get_parsed_value()
	{
		const auto& evaluated
			= [this](){
				try {
					const auto& temp = this->parser.Eval();
					if (temp.GetType() != 'm') {
						throw std::invalid_argument(
							std::string("Expression \"")
							+ this->parser.GetExpr()
							+ std::string("\" is not an array.")
						);
					}

					return temp.GetArray();

				} catch(mup::ParserError &error) {
					std::cerr <<  __FILE__ << "(" << __LINE__<< ") "
						<< "Could not evaluate expression \"" << this->parser.GetExpr()
						<< "\" for variable " << Variable::get_name()
						<< ": " << error.GetMsg()
						<< std::endl;
					throw;
				}
			}();

		// TODO: merge with identical code in int version
		if (evaluated.GetRows() != 1) {
			throw std::invalid_argument(
				std::string("Invalid number of rows in expression \"")
				+ parser.GetExpr()
				+ std::string("\" for variable ")
				+ Variable::get_name()
				+ std::string(", should be 1\n")
			);
		}

		constexpr size_t N = std::tuple_size<typename V::data_type>::value;

		if (evaluated.GetCols() != N) {
			throw std::invalid_argument(
				std::string("Invalid number of columns in expression \"")
				+ parser.GetExpr()
				+ std::string("\" for variable ")
				+ Variable::get_name()
				+ std::string("\n")
			);
		}

		typename Variable::data_type ret_val;
		for (size_t i = 0; i < N; i++) {
			ret_val[i] = evaluated.At(int(i)).GetFloat();
		}

		return ret_val;
	}


	//! used to evaluate math expression if Variable::data_type is array of integral numbers
	template<class V = Variable> typename std::enable_if<
		std::tuple_size<typename V::data_type>::value >= 0
			&& std::is_integral<typename V::data_type::value_type>::value,
		typename V::data_type
	>::type get_parsed_value()
	{
		const auto& evaluated
			= [this](){
				try {
					const auto& temp = this->parser.Eval();
					if (temp.GetType() != 'm') {
						throw std::invalid_argument(
							std::string("Expression \"")
							+ this->parser.GetExpr()
							+ std::string("\" is not an array.")
						);
					}

					return temp.GetArray();

				} catch(mup::ParserError &error) {
					std::cerr <<  __FILE__ << "(" << __LINE__<< ") "
						<< "Could not evaluate expression \""
						<< this->parser.GetExpr() << "\": " << error.GetMsg()
						<< " for variable " << Variable::get_name()
						<< std::endl;
					throw;
				}
			}();

		if (evaluated.GetRows() != 1) {
			throw std::invalid_argument(
				std::string("Invalid number of rows in expression \"")
				+ parser.GetExpr()
				+ std::string("\" for variable ")
				+ Variable::get_name()
				+ std::string(", should be 1\n")
			);
		}

		constexpr size_t N = std::tuple_size<typename Variable::data_type>::value;

		if (evaluated.GetCols() != N) {
			throw std::invalid_argument(
				std::string("Invalid number of columns in expression \"")
				+ parser.GetExpr()
				+ std::string("\" for variable ")
				+ Variable::get_name()
				+ std::string("\n")
			);
		}

		typename Variable::data_type ret_val;
		for (size_t i = 0; i < N; i++) {
			ret_val[i] = evaluated.At(int(i)).GetInteger();
		}

		return ret_val;
	}


	#ifdef EIGEN_WORLD_VERSION

	//! used to evaluate math expression if Variable::data_type is Eigen::Matrix<double or float, N, 1>
	template<class V = Variable> typename std::enable_if<
		V::data_type::ColsAtCompileTime == 1
			&& V::data_type::RowsAtCompileTime >= 1
			&& std::is_floating_point<typename V::data_type::value_type>::value,
		typename V::data_type
	>::type get_parsed_value()
	{
		const auto& evaluated
			= [this](){
				try {
					const auto& temp = this->parser.Eval();
					if (temp.GetType() != 'm') {
						throw std::invalid_argument(
							std::string("Expression ")
							+ this->parser.GetExpr()
							+ std::string(" is not an array.")
						);
					}

					return temp.GetArray();

				} catch(mup::ParserError &error) {
					std::cerr <<  __FILE__ << "(" << __LINE__<< ") "
						<< "Could not evaluate expression \""
						<< this->parser.GetExpr() << "\": " << error.GetMsg()
						<< " of variable " << Variable::get_name()
						<< std::endl;
					throw;
				}
			}();

		if (evaluated.GetRows() != 1) {
			throw std::invalid_argument(
				std::string("Invalid number of rows in expression ")
				+ parser.GetExpr()
				+ std::string(" for variable ")
				+ Variable::get_name()
				+ std::string(", should be 1\n")
			);
		}

		constexpr size_t N = V::data_type::RowsAtCompileTime;

		if (evaluated.GetCols() != N) {
			throw std::invalid_argument(
				std::string("Invalid number of columns in expression ")
				+ parser.GetExpr()
				+ std::string(" for variable ")
				+ Variable::get_name()
				+ std::string("\n")
			);
		}

		typename Variable::data_type ret_val;
		for (size_t i = 0; i < N; i++) {
			ret_val(i) = evaluated.At(int(i)).GetFloat();
		}

		return ret_val;
	}

	#endif // ifdef EIGEN_WORLD_VERSION
};

}} // namespaces

#endif // ifndef PAMHD_BOUNDARIES_INITIAL_CONDITION_HPP
