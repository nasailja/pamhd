/*
An experiment for a generic boundary class with several variables.

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

#include "array"
#include "cstdlib"
#include "iostream"
#include "string"
#include "type_traits"

#include "boost/algorithm/string/classification.hpp"
#include "boost/algorithm/string/constants.hpp"
#include "boost/algorithm/string/split.hpp"
#include "boost/lexical_cast.hpp"
#include "boost/program_options.hpp"
#include "prettyprint.hpp"

#include "boundaries/variables_to_options.hpp"

using namespace pamhd::boundaries;

struct Momentum_Density {
	using data_type = std::array<double, 3>;

	static const std::string get_name()
	{
		return {"momentum density"};
	}
	static const std::string get_option_name()
	{
		return {"momentum-density"};
	}
	static const std::string get_option_help()
	{
		return {"Momentum density (kg s / m^2)"};
	}
};

struct Magnetic_Field {
	using data_type = std::array<double, 3>;

	static const std::string get_name()
	{
		return {"magnetic field"};
	}
	static const std::string get_option_name()
	{
		return {"magnetic-field"};
	}
	static const std::string get_option_help()
	{
		return {"Magnetic field (T)"};
	}
};


// Variables used by geometry part of boundary class
template <class Vector_T> struct Start {
	using data_type = Vector_T;
	static const std::string get_option_name()
	{
		return {"start"};
	}
	static const std::string get_option_help()
	{
		return {"Starting coordinate of boundary"};
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
		return {"Ending coordinate of boundary"};
	}
};


template<
	class Cell_T,
	class Vector_T,
	class... Variables
> class Test_Boundary
{
private:

	template<
		class... Boundary_Data
	> void add_boundary_data (
		const Boundary_Data&... boundary_data
	) {}

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

	template<class Last_Variable, class Last_Data> void add_boundary_data(
		const Last_Variable&,
		const Last_Data& data
	) {
		this->boundaries[Last_Variable()].push_back(data);
	}


	template<class... Variable_Data_T> bool check_variable_options(
		const size_t required_number
	) const {
		return true;
	}

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



public:

	// type of the start coordinate of boundary geometry
	using Start_T = Start<Vector_T>;
	// type of the end coordinate of boundary geometry
	using End_T = End<Vector_T>;

	template<class... Boundary_Data> void add_boundary(
		const Vector_T& geometry_start,
		const Vector_T& geometry_end,
		const Boundary_Data&... boundary_data
	) {
		this->boundaries[Start_T()].push_back(geometry_start);
		this->boundaries[End_T()].push_back(geometry_end);
		this->add_boundary_data(boundary_data...);
	}


	const std::vector<Start_T>& get_geometry_start() const
	{
		return this->boundaries[Start_T()];
	}

	const std::vector<End_T>& get_geometry_end() const
	{
		return this->boundaries[End_T()];
	}


	void add_options(
		boost::program_options::options_description& options,
		const std::string& option_name_prefix
	) {
		this->boundaries.add_options(options, option_name_prefix);
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

	/*!
	Adds cell to all boundaries of this type.

	Returns number of boundaries cell was included in.
	*/
	size_t add(
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

	size_t get_number_of_boundaries() const
	{
		return this->boundaries[Start_T()].size();
	}

	const std::vector<Cell_T>& get_boundary_cells(
		const size_t boundary_index
	) const {
		return this->boundary_cells[boundary_index];
	}

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


int main(int argc, char* argv[])
{
	Test_Boundary<
		int,
		std::array<double, 3>,
		Momentum_Density,
		Magnetic_Field
	> boundaries;

	// set one default boundary
	boundaries.add_boundary(
		{{0, 0, 0}},
		{{1, 1, 1}},
		Momentum_Density{},
		Momentum_Density::data_type{1, 2, 3},
		Magnetic_Field{},
		Magnetic_Field::data_type{4, 5, 6}
	);

	boost::program_options::options_description options(
		"Usage: program_name [options], where options are:"
	);
	options.add_options()("help", "Print options and their descriptions");
	boundaries.add_options(options, "test.");

	boost::program_options::variables_map option_variables;
	boost::program_options::store(
		boost::program_options::parse_command_line(argc, argv, options),
		option_variables
	);
	boost::program_options::notify(option_variables);

	if (option_variables.count("help") > 0) {
		std::cout << options << std::endl;
		return EXIT_SUCCESS;
	}

	if (boundaries.get_number_of_boundaries() != 1) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
			<< "Invalid number of test boundaries."
			<< std::endl;
		return EXIT_FAILURE;
	}

	if (boundaries.add(-1, {{0.1, 0.2, 0.3}}, {{0.6, 0.5, 0.4}}) != 1) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
			<< "Cell added to invalid number of boundaries."
			<< std::endl;
		return EXIT_FAILURE;
	}

	if (boundaries.get_boundary_cells(0).size() != 1) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
			<< "Invalid number of cells in boundary 0."
			<< std::endl;
		return EXIT_FAILURE;
	}

	if (boundaries.get_boundary_cells(0)[0] != -1) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
			<< "Invalid first cell in boundary 0."
			<< std::endl;
		return EXIT_FAILURE;
	}

	if (boundaries.get_boundary_data(Momentum_Density(), 0)[0] != 1) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
			<< "Invalid boundary data in boundary 0."
			<< std::endl;
		return EXIT_FAILURE;
	}

	if (boundaries.get_boundary_data(Magnetic_Field(), 0)[2] != 6) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
			<< "Invalid boundary data in boundary 0."
			<< std::endl;
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}
