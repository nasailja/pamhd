/*
Class representing one time-dependent boundary of PAMHD.

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


#ifndef PAMHD_BOUNDARIES_BOUNDARY_TIME_DEPENDENT_HPP
#define PAMHD_BOUNDARIES_BOUNDARY_TIME_DEPENDENT_HPP


#include "cstdlib"
#include "iostream"
#include "stdexcept"
#include "string"
#include "utility"

#include "boost/program_options.hpp"
#include "prettyprint.hpp"

#include "boundaries/variable_to_option_vector.hpp"


namespace pamhd {
namespace boundaries {


template<
	class Geometry_T,
	class Cell_T,
	class Time_T,
	class... Variables
> class Boundary_Time_Dependent
{
public:

	static_assert(
		Geometry_T::scalar == false,
		"Invalid Geometry_T given to Boundary_Time_Dependent, use e.g. Boxes instead."
	);

	Geometry_T geometries;

	Variable_To_Option_Vector<Variables...> boundary_data;


	void add_options(
		const std::string& option_name_prefix,
		boost::program_options::options_description& options
	) {
		options.add_options()
			((option_name_prefix + "time").c_str(),
				boost::program_options::value<
					std::vector<Time_T>
				>(&this->time_stamps)
					->composing()
					->default_value(this->time_stamps),
				"Time stamp");

		this->geometries.add_options(option_name_prefix, options);
		this->boundary_data.add_options(option_name_prefix, options);
	}


	template<class Variable> typename Variable::data_type get_data(
		const Variable& variable,
		const std::array<double, 3>& given_position,
		const Time_T& given_time
	) {
		if (this->time_stamps.size() != boundary_data.get_number_of_expressions()) {
			std::cerr <<  __FILE__ << "(" << __LINE__<< ") "
				<< "Number of expressions for data and data time stamps differ: "
				<< this->time_stamps.size() << " and "
				<< boundary_data.get_number_of_expressions()
				<< std::endl;
			abort();
		}

		const size_t index
			= [&given_time, this](){
				try {
					return this->get_index(double(given_time));
				} catch (std::exception& e) {
					std::cerr << "Couldn't get data index for variable "
						<< Variable::get_name() << " at time " << given_time
						<< std::endl;
					throw;
				}
			}();

		return this->boundary_data.get_data(
			variable,
			index,
			given_position,
			double(given_time)
		);
	}


	/*!
	Classifies cell at first time stamp >= than given time.
	Doesn't remove previous classification (including other times),
	use clear_cells() for that.

	geometry_parameters is given to the overlaps function of Geometry_T.

	Returns true if given cell was added to cell list of this boundary,
	false otherwise.
	*/
	template<class... Geometry_Parameters> bool add_cell(
		const Time_T& time,
		const Cell_T& cell,
		Geometry_Parameters&&... geometry_parameters
	) {
		const size_t index
			= [&time, this](){
				try {
					return this->get_index(time);
				} catch (std::exception& e) {
					std::cerr << "Couldn't get data index for time " << time << std::endl;
					throw;
				}
			}();

		if (
			this->geometries.overlaps(
				index,
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


	/*!
	Translates given time to index in geometries and boundary_data.
	*/
	size_t get_index(const Time_T& given_time)
	{
		if (this->time_stamps.size() == 0) {
			throw std::out_of_range(
				"No time stamps set using set_time_stamp() "
				"in time dependent boundary."
			);
		}

		if (this->time_stamps.size() == 1) {
			return 0;
		}

		for (size_t i = 1; i < this->time_stamps.size(); i++) {
			if (this->time_stamps[i - 1] >= this->time_stamps[i]) {
				std::cerr <<  __FILE__ << "(" << __LINE__<< ") "
					<< "Time stamp at index " << i - 1
					<< " is >= than at index " << i
					<< std::endl;
				abort();
			}

			if (given_time < this->time_stamps[i]) {
				return i - 1;
			}
		}

		return this->time_stamps.size() - 1;
	}


	size_t get_number_of_instances() const
	{
		return this->time_stamps.size();
	}


	void set_number_of_instances(const size_t new_size)
	{
		this->time_stamps.resize(new_size);
		this->geometries.set_number_of_instances(new_size);
		this->boundary_data.set_number_of_expressions(new_size);
	}


	void set_time_stamp(
		const size_t index,
		const Time_T& given_time
	) {
		if (index >= this->time_stamps.size()) {
			std::cerr <<  __FILE__ << "(" << __LINE__<< ") "
				<< "Invalid index for time stamp: " << index
				<< ", should be < " << this->time_stamps.size()
				<< std::endl;
			abort();
		}

		this->time_stamps[index] = given_time;
	}



private:

	std::vector<Cell_T> cells;
	std::vector<Time_T> time_stamps;
};

}} // namespaces

#endif // ifndef PAMHD_BOUNDARIES_BOUNDARY_TIME_DEPENDENT_HPP
