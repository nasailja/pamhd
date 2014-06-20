/*
Tests initial condition class of PAMHD.

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
#include "set"
#include "string"
#include "vector"

#include "boost/program_options.hpp"

#include "boundaries/initial_condition.hpp"


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


int main(int argc, char* argv[])
{
	std::string initial_condition_file_name;

	boost::program_options::options_description options(
		"Usage: program_name [options], where options are:"
	);

	options.add_options()
		("help", "Print options and their descriptions")
		("initial-help", "Print help for options of initial condition")
		("initial-file",
			boost::program_options::value<std::string>(&initial_condition_file_name)
				->default_value(""),
			"Read initial condition from file arg (do not read if empty string)");

	Initial_Condition<
		int,
		std::array<double, 3>,
		Momentum_Density,
		Magnetic_Field
	> initial_condition;

	initial_condition.add_boundary_box(
		{{10, 10, 10}},
		{{20, 20, 20}},
		Momentum_Density{},
		Momentum_Density::data_type{1, 2, 3},
		Magnetic_Field{},
		Magnetic_Field::data_type{4, 5, 6}
	);
	initial_condition.add_boundary_box(
		{{-20, -20, -20}},
		{{-10, -10, -10}},
		Momentum_Density{},
		Momentum_Density::data_type{-1, -2, -3},
		Magnetic_Field{},
		Magnetic_Field::data_type{-4, -5, -6}
	);
	initial_condition.add_boundary_box(
		{{-18, -18, -18}},
		{{-11, -11, -11}},
		Momentum_Density{},
		Momentum_Density::data_type{3, 2, 1},
		Magnetic_Field{},
		Magnetic_Field::data_type{-6, -5, -4}
	);

	boost::program_options::options_description initial_condition_options(
		"Options for initial condition"
	);
	initial_condition.add_options(initial_condition_options, "initial.");

	boost::program_options::variables_map option_variables;
	boost::program_options::store(
		boost::program_options::parse_command_line(
			argc,
			argv,
			options
		),
		option_variables
	);
	boost::program_options::notify(option_variables);

	if (option_variables.count("help") > 0) {
		std::cout << options << std::endl;
		return EXIT_SUCCESS;
	}
	if (option_variables.count("initial-help") > 0) {
		std::cout << initial_condition_options << std::endl;
		return EXIT_SUCCESS;
	}

	if (initial_condition_file_name != "") {
		boost::program_options::store(
			boost::program_options::parse_config_file<char>(
				initial_condition_file_name.c_str(),
				initial_condition_options
			),
			option_variables
		);
		boost::program_options::notify(option_variables);
	}

	const std::vector<int> cells{
		0, 1, 2, 3, 4, 5,
		50, 51, 52, 53, 54, 55,
		100, 101, 102, 103, 104, 105
	};

	const std::vector<
		std::pair<
			std::array<double, 3>,
			std::array<double, 3>
		>
	> cell_coords{
		{{0.1, 0.1, 0.1}, {0.2, 0.2, 0.2}},
		{{1.1, 1.1, 1.1}, {1.2, 1.2, 1.2}},
		{{2.1, 2.1, 2.1}, {2.2, 2.2, 2.2}},
		{{3.1, 3.1, 3.1}, {3.2, 3.2, 3.2}},
		{{4.1, 4.1, 4.1}, {4.2, 4.2, 4.2}},
		{{5.1, 5.1, 5.1}, {5.2, 5.2, 5.2}},

		{{11.1, 11.1, 11.1}, {11.2, 11.2, 11.2}},
		{{12.1, 12.1, 12.1}, {12.2, 12.2, 12.2}},
		{{13.1, 13.1, 13.1}, {13.2, 13.2, 13.2}},
		{{14.1, 14.1, 14.1}, {14.2, 14.2, 14.2}},
		{{15.1, 15.1, 15.1}, {15.2, 15.2, 15.2}},
		{{16.1, 16.1, 16.1}, {16.2, 16.2, 16.2}},

		{{-11.3, -11.3, -11.3}, {-11.2, -11.2, -11.2}},
		{{-12.3, -12.3, -12.3}, {-12.2, -12.2, -12.2}},
		{{-13.3, -13.3, -13.3}, {-13.2, -13.2, -13.2}},
		{{-14.3, -14.3, -14.3}, {-14.2, -14.2, -14.2}},
		{{-15.3, -15.3, -15.3}, {-15.2, -15.2, -15.2}},
		{{-16.3, -16.3, -16.3}, {-16.2, -16.2, -16.2}},
	};


	for (size_t i = 0; i < cells.size(); i++) {
		const size_t ret_val = initial_condition.add_cell(
			cells[i],
			cell_coords[i].first,
			cell_coords[i].second
		);

		if (i <= 5 and ret_val > 0) {
			std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
				<< i << ", " << ret_val
				<< std::endl;
			abort();
		}

		if (i > 5 and i <= 11 and ret_val != 1) {
			std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
				<< i << ", " << ret_val
				<< std::endl;
			abort();
		}

		if (i > 11 and ret_val != 2) {
			std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
				<< i << ", " << ret_val
				<< std::endl;
			abort();
		}
	}

	std::set<int> non_boundary_cells(cells.cbegin(), cells.cend());

	for (size_t i = 0; i < initial_condition.get_number_of_boundaries(); i++) {
		const auto& boundary_cells = initial_condition.get_boundary_cells(i);
		for (const auto cell: boundary_cells) {
			non_boundary_cells.erase(cell);
		}
	}

	if (non_boundary_cells.size() != 6) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
			<< "Wrong number of non-boundary cells: " << non_boundary_cells.size()
			<< std::endl;
		abort();
	}


	std::multiset<int> boundary_cells_test;

	for (size_t i = 0; i < initial_condition.get_number_of_boundaries(); i++) {
		const auto& boundary_cells = initial_condition.get_boundary_cells(i);
		for (const auto cell: boundary_cells) {
			boundary_cells_test.insert(cell);
		}
	}

	for (size_t i = 0; i < cells.size(); i++) {
		if (i <= 5 and boundary_cells_test.count(cells[i]) != 0) {
			std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
				<< cells[i] << ", " << boundary_cells_test.count(cells[i])
				<< std::endl;
			abort();
		}

		if (i > 5 and i <= 11 and boundary_cells_test.count(cells[i]) != 1) {
			std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
				<< cells[i] << ", " << boundary_cells_test.count(cells[i])
				<< std::endl;
			abort();
		}

		if (i > 11 and boundary_cells_test.count(cells[i]) != 2) {
			std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
				<< cells[i] << ", " << boundary_cells_test.count(cells[i])
				<< std::endl;
			abort();
		}
	}


	return EXIT_SUCCESS;
}
