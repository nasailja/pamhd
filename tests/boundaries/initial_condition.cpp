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

#include "boost/optional.hpp"
#include "boost/program_options.hpp"
#ifdef HAVE_EIGEN
#include "Eigen/Core" // must be included before initial_condition.hpp
#endif

#include "boundaries/initial_condition.hpp"


using namespace pamhd::boundaries;


struct Mass_Density {
	using data_type = double;

	static const std::string get_name()
	{
		return {"mass density"};
	}
	static const std::string get_option_name()
	{
		return {"mass-density"};
	}
	static const std::string get_option_help()
	{
		return {"Mass density (kg / m^3)"};
	}
};

struct Momentum_Density {
	#ifdef HAVE_EIGEN
	using data_type = Eigen::Vector3d;
	#else
	using data_type = std::array<double, 3>;
	#endif

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


#ifdef HAVE_EIGEN
#define MAKE_VEC(x, y, z) Eigen::Vector3d(x, y, z)
#define MAKE_BOX(a, b, c, d, e, f) MAKE_VEC(a, b, c), MAKE_VEC(d, e, f)
#else
#define MAKE_VEC(x, y, z) std::array<double, 3>{{x, y, z}}
#define MAKE_BOX(a, b, c, d, e, f) MAKE_VEC(a, b, c), MAKE_VEC(d, e, f)
#endif


int main(int argc, char* argv[])
{
	Initial_Condition<
		int,
		double,
		#ifdef HAVE_EIGEN
		Eigen::Vector3d,
		#else
		std::array<double, 3>,
		#endif
		Mass_Density,
		Momentum_Density
	> initial_condition;

	//! Set default values
	initial_condition.default_data.set_expression(
		Mass_Density(),
		"r[0] + r[1] + r[2] + t"
	);
	initial_condition.default_data.set_expression(
		Momentum_Density(),
		"{r[0] + t, r[1] + t, r[2] + t}"
	);

	boost::optional<size_t> result;

	result = initial_condition.add_box(MAKE_BOX(1, 1, 1, 3, 3, 3));
	if (not result or *result != 0) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
			<< "Invalid boundary id for added box: " << *result
			<< std::endl;
		abort();
	}

	result = initial_condition.add_sphere(
		MAKE_VEC(2, 2, 2),
		0.5
	);
	if (not result or *result != 1) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
			<< "Invalid boundary id for added sphere: " << *result 
			<< std::endl;
		abort();
	}

	initial_condition.set_expression(
		Mass_Density(),
		0,
		"r[0] + r[1] + r[2] + t + 1"
	);
	initial_condition.set_expression(
		Momentum_Density(),
		0,
		"{r[0] + t + 1, r[1] + t + 1, r[2] + t + 1}"
	);
	initial_condition.set_expression(
		Mass_Density(),
		1,
		"r[0] + r[1] + r[2] + t + 2"
	);
	initial_condition.set_expression(
		Momentum_Density(),
		1,
		"{r[0] + t + 2, r[1] + t + 2, r[2] + t + 2}"
	);


	std::string initial_condition_file_name("");

	boost::program_options::options_description
		general_options(
			"Usage: program_name [options], where options are"
		),
		initial_condition_options(
			"Options related to initial condition"
		);

	general_options.add_options()
		("help", "Print options and their descriptions")
		("initial-help", "Print help for options of initial condition")
		("initial-file",
			boost::program_options::value<std::string>(&initial_condition_file_name)
				->default_value(initial_condition_file_name),
			"Read initial condition from file arg (do not read if empty string)");
	initial_condition.add_initialization_options("", general_options);

	boost::program_options::variables_map option_variables;
	boost::program_options::store(
		boost::program_options::parse_command_line(
			argc,
			argv,
			general_options
		),
		option_variables
	);
	boost::program_options::notify(option_variables);

	if (option_variables.count("help") > 0) {
		std::cout << general_options << std::endl;
		return EXIT_SUCCESS;
	}

	initial_condition.add_options("", initial_condition_options);

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

	const std::vector<int> cells{0, 1, 2, 3, 4};

	const std::vector<
		std::pair<
			#ifdef HAVE_EIGEN
			Eigen::Vector3d,
			Eigen::Vector3d
			#else
			std::array<double, 3>,
			std::array<double, 3>
			#endif
		>
	> cell_coords{
		{MAKE_BOX(0.1, 0.1, 0.1, 0.2, 0.2, 0.2)},
		{MAKE_BOX(1.1, 1.1, 1.1, 1.2, 1.2, 1.2)},
		{MAKE_BOX(2.1, 2.1, 2.1, 2.2, 2.2, 2.2)},
		{MAKE_BOX(3.1, 3.1, 3.1, 3.2, 3.2, 3.2)},
		{MAKE_BOX(4.1, 4.1, 4.1, 4.2, 4.2, 4.2)}
	};

	if (cells.size() != cell_coords.size()) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< ")" << std::endl;
		abort();
	}


	for (size_t i = 0; i < cells.size(); i++) {
		const size_t ret_val = initial_condition.add_cell(
			cells[i],
			cell_coords[i].first,
			cell_coords[i].second
		);

		if (i == 0 and ret_val != 0) {
			std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
				<< i << ", " << ret_val
				<< std::endl;
			abort();
		}

		if (i == 1 and ret_val != 1) {
			std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
				<< i << ", " << ret_val
				<< std::endl;
			abort();
		}

		if (i == 2 and ret_val != 2) {
			std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
				<< i << ", " << ret_val
				<< std::endl;
			abort();
		}

		if (i > 2 and ret_val != 0) {
			std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
				<< i << ", " << ret_val
				<< std::endl;
			abort();
		}
	}

	std::set<int> non_boundary_cells(cells.cbegin(), cells.cend());

	for (size_t i = 0; i < initial_condition.get_number_of_boundaries(); i++) {
		const auto& boundary_cells = initial_condition.get_cells(i);
		for (const auto cell: boundary_cells) {
			non_boundary_cells.erase(cell);
		}
	}

	if (non_boundary_cells.size() != 3) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
			<< "Wrong number of non-boundary cells: " << non_boundary_cells.size()
			<< std::endl;
		abort();
	}


	std::multiset<int> boundary_cells_test;

	for (size_t i = 0; i < initial_condition.get_number_of_boundaries(); i++) {
		const auto& boundary_cells = initial_condition.get_cells(i);
		for (const auto cell: boundary_cells) {
			boundary_cells_test.insert(cell);
		}
	}

	for (size_t i = 0; i < cells.size(); i++) {
		if (i == 0 and boundary_cells_test.count(cells[i]) != 0) {
			std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
				<< cells[i] << ", " << boundary_cells_test.count(cells[i])
				<< std::endl;
			abort();
		}

		if (i == 1 and boundary_cells_test.count(cells[i]) != 1) {
			std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
				<< cells[i] << ", " << boundary_cells_test.count(cells[i])
				<< std::endl;
			abort();
		}

		if (i == 2 and boundary_cells_test.count(cells[i]) != 2) {
			std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
				<< cells[i] << ", " << boundary_cells_test.count(cells[i])
				<< std::endl;
			abort();
		}

		if (i > 2 and boundary_cells_test.count(cells[i]) != 0) {
			std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
				<< cells[i] << ", " << boundary_cells_test.count(cells[i])
				<< std::endl;
			abort();
		}
	}

	return EXIT_SUCCESS;
}
