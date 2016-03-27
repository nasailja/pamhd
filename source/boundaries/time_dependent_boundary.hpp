/*
Class for setting time-dependent boundary condition of a simulation.

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

#ifndef PAMHD_TIME_DEPENDENT_BOUNDARY_HPP
#define PAMHD_TIME_DEPENDENT_BOUNDARY_HPP


#include "boundaries/box.hpp"


namespace pamhd {
namespace boundaries {

namespace detail {

//! Variable for time stamping other variables
template <class Given_Time_T> struct Time {
	using data_type = Given_Time_T;
	static const std::string get_option_name()
	{
		return {"time"};
	}
	static const std::string get_option_help()
	{
		return {"Timestamp of bounday data"};
	}
};

} // namespace detail


/*!
\param Cell_T Type of cells to classify
\param Vector_T Type of coordinate to use for geometry
\param Time_T Time type to use for time stamping other data
\param Variables Time-dependent variables to store
*/
template<
	class Geometry_T,
	class Cell_T,
	class Vector_T,
	class Time_T,
	class... Variables
> class Time_Dependent_Boundary
{
public:


	template<class... Boundary_Data> void add_boundary_data(
		const Vector_T& geometry_start,
		const Vector_T& geometry_end,
		const Time_T& given_time,
		const Boundary_Data&... boundary_data
	) {
		this->boundaries.add_boundary(
			geometry_start,
			geometry_end,
			detail::Time<Time_T>(),
			given_time,
			boundary_data...
		);
		while (
			this->classified_at_index.size()
			< this->boundaries.get_number_of_boundaries()
		) {
			this->classified_at_index.push_back(false);
		}
	}


	/*
	Returns true if boundary has been set, false otherwise.
	*/
	bool exists() const
	{
		if (this->boundaries.get_number_of_boundaries() == 0) {
			return false;
		} else {
			return true;
		}
	}


	void add_options(
		boost::program_options::options_description& options,
		const std::string& option_name_prefix
	) {
		this->boundaries.add_options(options, option_name_prefix + "box.");
	}


	size_t add_cell(
		const Cell_T& cell,
		const Vector_T& cell_start,
		const Vector_T& cell_end,
		const Time_T& given_time
	) {
		while (
			this->classified_at_index.size()
			< this->boundaries.get_number_of_boundaries()
		) {
			this->classified_at_index.push_back(false);
		}

		const size_t bdy_i = this->get_boundary_index(given_time);

		// cells have now been classified so don't use previous cell lists
		this->classified_at_index[bdy_i] = true;

		return this->boundaries.add_cell(bdy_i, cell, cell_start, cell_end);
	}


	//! Returns cells which belong to given boundary
	const std::vector<Cell_T>& get_boundary_cells(
		const Time_T& given_time
	) const {
		size_t boundary_index = this->get_boundary_index(given_time);
		// return cells from last time they were classified
		while (
			not this->classified_at_index[boundary_index]
			and boundary_index > 0
		) {
			boundary_index--;
		}
		return this->boundaries.get_boundary_cells(boundary_index);
	}


	//! Returns data of given variable in given boundary
	template<class Variable> const typename Variable::data_type& get_boundary_data(
		const Variable&,
		const Time_T& given_time
	) const {
		return this->boundaries.get_boundary_data(
			Variable(),
			this->get_boundary_index(given_time)
		);
	}


	//! Removes all cells from this boundary.
	void clear_cells()
	{
		this->boundaries.clear_cells();
	}


private:
	Box<Cell_T, Vector_T, detail::Time<Time_T>, Variables...> boundaries;

	// whether any cells have been classified at given boundary index
	std::vector<bool> classified_at_index;

	/*!
	Returns boundary index corresponding to given time.
	*/
	size_t get_boundary_index(const Time_T& given_time) const
	{
		const auto& times
			= this->boundaries.get_all_boundary_data(detail::Time<Time_T>());
		if (times.size() == 0) {
			std::cerr <<  __FILE__ << "(" << __LINE__<< ") "
				<< "Should have at least one time stamped item of data."
				<< std::endl;
			abort();
		}

		if (times.size() == 1) {
			return 0;
		}

		// check that time values are ascending
		Time_T previous = times[0];
		for (size_t i = 1; i < times.size(); i++) {
			const Time_T& next = times[i];
			if (next < previous) {
				std::cerr <<  __FILE__ << "(" << __LINE__<< ") "
					<< "Time stamp at index " << i << "(" << next
					<< ") has lower value than at index " << i - 1
					<< " (" << previous << ")"
					<< std::endl;
				abort();
			}
			previous = next;
		}

		size_t bdy_i = 0;
		while (bdy_i < times.size() - 1) {
			if (given_time < times[bdy_i + 1]) {
				break;
			}
			bdy_i++;
		}

		return bdy_i;
	}
};

}} // namespaces

#endif // ifndef PAMHD_TIME_DEPENDENT_BOUNDARY_HPP
