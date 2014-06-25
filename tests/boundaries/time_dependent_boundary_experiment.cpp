/*
An experiment for a time-dependent generic boundary class.

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
#include "boost/optional.hpp"
#include "boost/program_options.hpp"
#include "prettyprint.hpp"
#include "string"

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

template <class Time_T> struct Time {
	using data_type = Time_T;
	static const std::string get_option_name()
	{
		return {"time"};
	}
	static const std::string get_option_help()
	{
		return {"Timestamp of data"};
	}
};


template<
	class Cell_T,
	class Vector_T,
	class Given_Time_T,
	class Variable
> class Time_Dependent_Boundary
{
public:

	// type of the start coordinate of boundary geometry
	using Start_T = Start<Vector_T>;
	// type of the end coordinate of boundary geometry
	using End_T = End<Vector_T>;
	// type of time stamps for data
	using Time_T = Time<Given_Time_T>;

	void add_boundary(
		const Vector_T& geometry_start,
		const Vector_T& geometry_end,
		const typename Time_T::data_type& time,
		const typename Variable::data_type& boundary_data
	) {
		this->boundaries[Start_T()].push_back(geometry_start);
		this->boundaries[End_T()].push_back(geometry_end);
		this->boundaries[Time_T()].push_back(time);
		this->boundaries[Variable()].push_back(boundary_data);

		this->boundary_cells.resize(this->boundaries[Start_T()].size());
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
	Classify cell at given time.

	Returns 0 if cell wasn't included in this boundary and 1 if it was.
	*/
	size_t add(
		const Cell_T& cell,
		const Vector_T& cell_start,
		const Vector_T& cell_end,
		const typename Time_T::data_type& given_time
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

		const size_t bdy_i = this->get_boundary_index(given_time);

		// cells have now been classified so don't use previous cell lists
		if (not this->boundary_cells[bdy_i]) {
			this->boundary_cells[bdy_i] = std::vector<Cell_T>();
		}

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
			this->boundary_cells[bdy_i]->push_back(cell);
			return 1;
		}

		return 0;
	}


	const std::vector<Cell_T>& get_boundary_cells(
		const typename Time_T::data_type& given_time
	) const {
		size_t boundary_index = this->get_boundary_index(given_time);
		// return cells from last time they were classified
		while (
			!this->boundary_cells[boundary_index]
			and boundary_index > 0
		) {
			boundary_index--;
		}
		return *this->boundary_cells[boundary_index];
	}


	const typename Variable::data_type& get_boundary_data(
		const Variable&,
		const typename Time_T::data_type& given_time
	) const {
		return this->boundaries[Variable()][
			this->get_boundary_index(given_time)
		];
	}


private:
	Variables_To_Options<
		Start_T,
		End_T,
		Time_T,
		Variable
	> boundaries;

	std::vector<boost::optional<std::vector<Cell_T>>> boundary_cells;

	/*!
	Returns boundary index corresponding to given time.
	*/
	size_t get_boundary_index(const typename Time_T::data_type& given_time) const
	{
		const auto& times = this->boundaries[Time_T()];
		if (times.size() == 0) {
			std::cerr <<  __FILE__ << "(" << __LINE__<< ") "
				<< "Should have at least one time stamped item of data."
				<< std::endl;
			abort();
		}

		if (times.size() == 1) {
			return 0;
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


int main(int argc, char* argv[])
{
	Time_Dependent_Boundary<int, std::array<double, 1>, double, Is_Alive> boundaries;
	// set default boundary
	boundaries.add_boundary({{0}}, {{1}}, 0, 1);
	boundaries.add_boundary({{0}}, {{1}}, 1, 0);
	boundaries.add_boundary({{1}}, {{2}}, 2, 1);

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

	if (boundaries.add(-1, {{0.1}}, {{0.6}}, -1) != 1) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
			<< "Cell wasn't added to boundary."
			<< std::endl;
		return EXIT_FAILURE;
	}
	if (boundaries.get_boundary_cells(-1).size() != 1) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
			<< "Invalid number of cells in boundary."
			<< std::endl;
		return EXIT_FAILURE;
	}
	if (boundaries.get_boundary_cells(0.9).size() != 1) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
			<< "Invalid number of cells in boundary."
			<< std::endl;
		return EXIT_FAILURE;
	}
	if (boundaries.get_boundary_cells(1.1).size() != 1) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
			<< "Invalid number of cells in boundary."
			<< std::endl;
		return EXIT_FAILURE;
	}
	if (boundaries.get_boundary_cells(2.1).size() != 1) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
			<< "Invalid number of cells in boundary."
			<< std::endl;
		return EXIT_FAILURE;
	}


	if (boundaries.add(0, {{0.1}}, {{0.6}}, 0.5) != 1) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
			<< "Cell wasn't added to boundary."
			<< std::endl;
		return EXIT_FAILURE;
	}
	if (boundaries.get_boundary_cells(-1).size() != 2) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
			<< "Invalid number of cells in boundary."
			<< std::endl;
		return EXIT_FAILURE;
	}
	if (boundaries.get_boundary_cells(0.9).size() != 2) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
			<< "Invalid number of cells in boundary."
			<< std::endl;
		return EXIT_FAILURE;
	}
	if (boundaries.get_boundary_cells(1.1).size() != 2) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
			<< "Invalid number of cells in boundary."
			<< std::endl;
		return EXIT_FAILURE;
	}
	if (boundaries.get_boundary_cells(2.1).size() != 2) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
			<< "Invalid number of cells in boundary."
			<< std::endl;
		return EXIT_FAILURE;
	}


	if (boundaries.add(1, {{0.1}}, {{0.6}}, 1.1) != 1) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
			<< "Cell wasn't added to boundary."
			<< std::endl;
		return EXIT_FAILURE;
	}
	if (boundaries.get_boundary_cells(-1).size() != 2) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
			<< "Invalid number of cells in boundary."
			<< std::endl;
		return EXIT_FAILURE;
	}
	if (boundaries.get_boundary_cells(0.9).size() != 2) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
			<< "Invalid number of cells in boundary."
			<< std::endl;
		return EXIT_FAILURE;
	}
	if (boundaries.get_boundary_cells(1.1).size() != 1) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
			<< "Invalid number of cells in boundary."
			<< std::endl;
		return EXIT_FAILURE;
	}
	if (boundaries.get_boundary_cells(2.1).size() != 1) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
			<< "Invalid number of cells in boundary."
			<< std::endl;
		return EXIT_FAILURE;
	}


	if (boundaries.add(2, {{0.1}}, {{0.6}}, 1.9) != 1) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
			<< "Cell wasn't added to boundary."
			<< std::endl;
		return EXIT_FAILURE;
	}
	if (boundaries.get_boundary_cells(-1).size() != 2) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
			<< "Invalid number of cells in boundary."
			<< std::endl;
		return EXIT_FAILURE;
	}
	if (boundaries.get_boundary_cells(0.9).size() != 2) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
			<< "Invalid number of cells in boundary."
			<< std::endl;
		return EXIT_FAILURE;
	}
	if (boundaries.get_boundary_cells(1.1).size() != 2) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
			<< "Invalid number of cells in boundary."
			<< std::endl;
		return EXIT_FAILURE;
	}
	if (boundaries.get_boundary_cells(2.1).size() != 2) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
			<< "Invalid number of cells in boundary."
			<< std::endl;
		return EXIT_FAILURE;
	}


	if (boundaries.add(3, {{0.1}}, {{0.6}}, 2.1) != 0) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
			<< "Cell was added to boundary."
			<< std::endl;
		return EXIT_FAILURE;
	}
	if (boundaries.get_boundary_cells(-1).size() != 2) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
			<< "Invalid number of cells in boundary."
			<< std::endl;
		return EXIT_FAILURE;
	}
	if (boundaries.get_boundary_cells(0.9).size() != 2) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
			<< "Invalid number of cells in boundary."
			<< std::endl;
		return EXIT_FAILURE;
	}
	if (boundaries.get_boundary_cells(1.1).size() != 2) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
			<< "Invalid number of cells in boundary."
			<< std::endl;
		return EXIT_FAILURE;
	}
	if (boundaries.get_boundary_cells(2.1).size() != 0) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
			<< "Invalid number of cells in boundary."
			<< std::endl;
		return EXIT_FAILURE;
	}


	if (boundaries.add(4, {{1.1}}, {{1.6}}, 2.1) != 1) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
			<< "Cell wasn't added to boundary."
			<< std::endl;
		return EXIT_FAILURE;
	}
	if (boundaries.get_boundary_cells(-1).size() != 2) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
			<< "Invalid number of cells in boundary."
			<< std::endl;
		return EXIT_FAILURE;
	}
	if (boundaries.get_boundary_cells(0.9).size() != 2) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
			<< "Invalid number of cells in boundary."
			<< std::endl;
		return EXIT_FAILURE;
	}
	if (boundaries.get_boundary_cells(1.1).size() != 2) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
			<< "Invalid number of cells in boundary."
			<< std::endl;
		return EXIT_FAILURE;
	}
	if (boundaries.get_boundary_cells(2.1).size() != 1) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
			<< "Invalid number of cells in boundary."
			<< std::endl;
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}
