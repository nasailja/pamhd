/*
Box boundary class of PAMHD.

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


#ifndef PAMHD_BOUNDARIES_BOX_HPP
#define PAMHD_BOUNDARIES_BOX_HPP


#include "cstdlib"
#include "iostream"
#include "string"

#include "boost/program_options.hpp"
#include "prettyprint.hpp"

#include "boundaries/variables_to_options.hpp"


namespace pamhd {
namespace boundaries {


namespace box_detail {

template <class Vector_T> struct Start {
	using data_type = Vector_T;
	static const std::string get_option_name()
	{
		return {"start"};
	}
	static const std::string get_option_help()
	{
		return {"Box start coordinate"};
	}
};

template <class Vector_T> struct End {
	using data_type = Vector_T;
	static const std::string get_option_name()
	{
		return {"end"};
	}
	static const std::string get_option_help()
	{
		return {"Box end coordinate"};
	}
};

} // namespace box_detail


/*!
A box boundary class of PAMHD.

Stores variable data and cells which overlap with one
or more boundary boxes. Boundaries can be added manually
or through boost::program_options.

\param Cell_T Type of cells to store (e.g. int)
\param Vector_T Type of coordinates to use (e.g. std::array)
\param Variables Variables to store in each boundary box
see pamhd::boundaries::Box::add_boundary() for the required API.

Add boundaries manually with
pamhd::boundaries::Box:: add_boundary() or from the command
line with pamhd::boundaries::Box::add_options().

Classify cells that overlap with one or more boundary with
pamhd::boundaries::Box::add_cell() and retrieve cells and
data of different boundaries with
pamhd::boundaries::Box::get_number_of_boundaries(),
pamhd::boundaries::Box::get_boundary_cells() and
pamhd::boundaries::Box::get_boundary_data().
*/
template<
	class Cell_T,
	class Vector_T,
	class... Variables
> class Box
{
private:

	/*!
	Adds data of given variables to the internal storage.
	*/
	template<
		class... Boundary_Data
	> void add_boundary_data(
		const Boundary_Data&... boundary_data
	) {}

	//! Start/continues recursion over variables and their data
	template<class Variable, class Data, class... Rest> void add_boundary_data(
		const Variable&,
		const Data& data,
		const Rest&... rest
	) {
		static_assert(
			std::is_same<typename Variable::data_type, Data>::value,
			"Variable and its data must be equal"
		);
		this->boundaries[Variable()].push_back(data);
		this->add_boundary_data(rest...);
	}

	//! Stops recursion over variables and their data
	template<class Last_Variable, class Last_Data> void add_boundary_data(
		const Last_Variable&,
		const Last_Data& data
	) {
		this->boundaries[Last_Variable()].push_back(data);
	}


	/*!
	Checks that correct number of parameters was given to variables.
	*/
	template<class... Variable_Data_T> bool check_variable_options(
		const size_t required_number
	) const {
		return true;
	}

	//! Starts/continues recursion over variables
	template<class First, class... Rest> bool check_variable_options(
		const size_t required_number
	) const {
		if (this->boundaries[First()].size() != required_number) {
			std::cout << "Wrong number of values for variable "
				<< First::get_name() << " of test boundary: "
				<< this->boundaries[First()].size()
				<< ", should be " << required_number
				<< std::endl;
			return false;
		}

		return this->check_variable_options<Rest...>(required_number);
	}

	//! Stop recursion over varibles
	template<class Last> bool check_variable_options(
		const size_t required_number
	) const {
		if (this->boundaries[Last()].size() != required_number) {
			std::cout << "Wrong number of values for variable "
				<< Last::get_name() << " of test boundary: "
				<< this->boundaries[Last()].size()
				<< ", should be " << required_number
				<< std::endl;
			return false;
		}

		return true;
	}


	/*!
	Checks that options were given correctly by the user.

	Returns true in case they were, false otherwise.
	*/
	bool check_options() const
	{
		if (
			this->boundaries[Start_T()].size()
			!= this->boundaries[End_T()].size()
		) {
			std::cout << "Number of parameters for test boundaries not equal."
				<< std::endl;
			return false;
		}

		for (size_t
			bdy_i = 0;
			bdy_i < this->boundaries[Start_T()].size();
			bdy_i++
		) {
			const auto
				&bdy_start = this->boundaries[Start_T()][bdy_i],
				&bdy_end = this->boundaries[End_T()][bdy_i];

			for (size_t coord_i = 0; coord_i < bdy_start.size(); coord_i++) {
				if (bdy_start[coord_i] >= bdy_end[coord_i]) {
					std::cerr <<  __FILE__ << "(" << __LINE__<< ") "
						<< "Starting coordinate of boundary at index " << coord_i
						<< " is not smaller than ending coordinate: "
						<< bdy_start[coord_i] << " >= " << bdy_end[coord_i]
						<< std::endl;
					return false;
				}
			}
		}

		return true;
	}



public:

	//! Type of the start coordinate of boundary geometry
	using Start_T = box_detail::Start<Vector_T>;
	//! Type of the end coordinate of boundary geometry
	using End_T = box_detail::End<Vector_T>;

	/*!
	Adds a box boundary with given data to the list of box boundaries.

	\param geometry_start Start coordinate of the added box
	\param geometry_end End coordinate of the added box
	\param boundary_data Variables and their data of the added box

	Example:
	\code
	struct Is_Alive {
		using data_type = int;
		static const std::string get_name() { return {"is alive"}; }
		static const std::string get_option_name() { return {"is-alive"}; }
		static const std::string get_option_help() { return {"TODO"}; }
	};
	Box_Boundary<
		int,
		std::array<double, 3>,
		Is_Alive
	> boundaries;
	boundaries.add_boundary(
		{{0, 0, 0}},
		{{1, 1, 1}},
		Is_Alive{},
		Is_Alive::data_type{1}
	);
	\endcode
	*/
	template<class... Boundary_Data> void add_boundary(
		const Vector_T& geometry_start,
		const Vector_T& geometry_end,
		const Boundary_Data&... boundary_data
	) {
		this->boundaries[Start_T()].push_back(geometry_start);
		this->boundaries[End_T()].push_back(geometry_end);
		this->add_boundary_data(boundary_data...);
	}


	//! Returns start coordinates of existing box boundaries
	const std::vector<Start_T>& get_geometry_start() const
	{
		return this->boundaries[Start_T()];
	}

	//! Returns end coordinates of existing box boundaries
	const std::vector<End_T>& get_geometry_end() const
	{
		return this->boundaries[End_T()];
	}


	/*!
	Adds options required add a box boundary to the list of box boundaries.

	\param options Destination where to add the necessary options
	\param option_name_prefix Prefix to add to each necessary option

	Example:
	\code
	boost::program_options::options_description options(
		"Usage: program_name [options], where options are:"
	);
	...
	Box_Boundary<
		int,
		std::array<double, 3>,
		Is_Alive
	> boundaries;
	boundaries.add_options(options, "box.");
	...
	boost::program_options::notify(option_variables);
	\endcode
	*/
	void add_options(
		boost::program_options::options_description& options,
		const std::string& option_name_prefix
	) {
		this->boundaries.add_options(options, option_name_prefix);
	}


	/*!
	Adds a cell to all existing boundaries with which the cell overlaps.

	\param cell Cell to add
	\param cell_start Start coordinate of cell
	\param cell_end End coordinate of cell

	Returns number of boundaries cell was added to.

	*/
	size_t add_cell(
		const Cell_T& cell,
		const Vector_T& cell_start,
		const Vector_T& cell_end
	) {
		for (size_t i = 0; i < cell_start.size(); i++) {
			if (cell_start[i] > cell_end[i]) {
				std::cerr <<  __FILE__ << "(" << __LINE__<< ") "
					<< "Starting coordinate at index " << i
					<< " is larger than ending coordinate: "
					<< cell_start[i] << " > " << cell_end[i]
					<< std::endl;
				abort();
			}
		}

		if (not this->check_options()) {
			std::cerr <<  __FILE__ << "(" << __LINE__<< ") "
				<< "Inconsistent boudary parameters."
				<< std::endl;
			abort();
		}

		this->boundary_cells.resize(this->boundaries[Start_T()].size());

		size_t ret_val = 0;

		for (size_t
			bdy_i = 0;
			bdy_i < this->boundaries[Start_T()].size();
			bdy_i++
		) {
			const auto
				&bdy_start = this->boundaries[Start_T()][bdy_i],
				&bdy_end = this->boundaries[End_T()][bdy_i];

			bool overlaps = true;
			for (size_t coord_i = 0; coord_i < cell_start.size(); coord_i++) {
				if (
					cell_start[coord_i] > bdy_end[coord_i]
					or cell_end[coord_i] < bdy_start[coord_i]
				) {
					overlaps = false;
					break;
				}
			}

			if (overlaps) {
				this->boundary_cells[bdy_i].push_back(cell);
				ret_val++;
			}
		}

		return ret_val;
	}


	//! Returns number of existing boundaries
	size_t get_number_of_boundaries() const
	{
		return this->boundaries[Start_T()].size();
	}

	//! Returns cells which belong to given boundary
	const std::vector<Cell_T>& get_boundary_cells(
		const size_t boundary_index
	) const {
		return this->boundary_cells[boundary_index];
	}

	//! Returns data of given variable in given boundary
	template<class Variable> const typename Variable::data_type& get_boundary_data(
		const Variable&,
		const size_t boundary_index
	) const {
		return this->boundaries[Variable()][boundary_index];
	}


private:
	Variables_To_Options<
		Start_T,
		End_T,
		Variables...
	> boundaries;

	/*
	boundary_cells[i] == cells (partially)
	overlapping boundaries[Start_T/End_T][i]
	*/
	std::vector<std::vector<Cell_T>> boundary_cells;
};

}} // namespaces

#endif // ifndef PAMHD_BOUNDARIES_BOX_HPP
