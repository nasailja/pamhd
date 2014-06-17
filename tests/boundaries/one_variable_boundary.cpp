/*
An experiment on how a generic boundary class could be implemented.

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

#include "boost/algorithm/string/classification.hpp"
#include "boost/algorithm/string/constants.hpp"
#include "boost/algorithm/string/split.hpp"
#include "boost/lexical_cast.hpp"
#include "boost/program_options.hpp"
#include "prettyprint.hpp"

#include "boundaries/variables_to_options.hpp"

using namespace pamhd::boundaries;

// Boundary variable
struct Is_Alive {
	/*
	Avoid returning temporary references for
	boundary data by not using bool.

	Boundary data is stored in vector, which
	has a special version for vector<bool>
	that cannot return real references to its
	data.
	*/
	using data_type = int;

	static const std::string get_name()
	{
		return {"is alive"};
	}
	static const std::string get_option_name()
	{
		return {"is-alive"};
	}
	static const std::string get_option_help()
	{
		return {"Whether cell is alive (> 0) or not (<= 0)"};
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


namespace std {
	template<class T, size_t N> void validate(
		boost::any& value,
		const std::vector<string>& all_parsed,
		array<T, N>*,
		long
	) {
		boost::program_options::validators::check_first_occurrence(value);

		const string& parsed
			= boost::program_options::validators::get_single_string(all_parsed);

		vector<string> components;
		components = boost::algorithm::split(
			components,
			parsed,
			boost::algorithm::is_any_of(" \t,;"),
			boost::algorithm::token_compress_on
		);

		if (components.size() != N) {
			std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
				<< "Option must have " << N << " values "
				<< "separated by e.g. a , or a space but given option ("
				<< parsed << ") has " << components.size()
				<< std::endl;
			abort();
		}

		array<T, N> final;
		for (size_t i = 0; i < N; i++) {
			final[i] = boost::lexical_cast<T>(components[i]);
		}

		value = final;
	}
}


template<
	class Cell_T,
	class Vector_T,
	class Variable
> class Test_Boundary
{
public:

	// type of the start coordinate of boundary geometry
	using Start_T = Start<Vector_T>;
	// type of the end coordinate of boundary geometry
	using End_T = End<Vector_T>;

	void add_boundary(
		const Vector_T& geometry_start,
		const Vector_T& geometry_end,
		const typename Variable::data_type& boundary_data
	) {
		this->boundaries[Start_T()].push_back(geometry_start);
		this->boundaries[End_T()].push_back(geometry_end);
		this->boundaries[Variable()].push_back(boundary_data);
	}


	const std::vector<Start_T>& get_geometry_start() const
	{
		return this->boundaries[Start_T()];
	}

	const std::vector<End_T>& get_geometry_end() const
	{
		return this->boundaries[End_T()];
	}

	const std::vector<Variable>& operator[](const Variable&) const
	{
		return this->boundaries[Variable()];
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
			this->boundaries[Start_T()].size() != this->boundaries[End_T()].size()
			or this->boundaries[End_T()].size() != this->boundaries[Variable()].size()
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

	const typename Variable::data_type& get_boundary_data(
		const Variable&,
		const size_t boundary_index
	) const {
		return this->boundaries[Variable()][boundary_index];
	}


private:
	Variables_To_Options<
		Start_T,
		End_T,
		Variable
	> boundaries;

	/*
	boundary_cells[i] == cells (partially)
	overlapping boundaries[Start_T/End_T][i]
	*/
	std::vector<std::vector<Cell_T>> boundary_cells;
};


int main(int argc, char* argv[])
{
	Test_Boundary<int, std::array<double, 3>, Is_Alive> boundaries;
	// set one default boundary
	boundaries.add_boundary({{0, 0, 0}}, {{1, 1, 1}}, 1);

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

	if (boundaries.get_boundary_data(Is_Alive(), 0) != 1) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
			<< "Invalid boundary data in boundary 0."
			<< std::endl;
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}
