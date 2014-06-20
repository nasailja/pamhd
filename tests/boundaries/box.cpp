/*
Tests the box boundary of PAMHD.

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

#include "boost/program_options.hpp"

#include "boundaries/box.hpp"


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
	std::string boundary_file_name;

	boost::program_options::options_description options(
		"Usage: program_name [options], where options are:"
	);

	options.add_options()
		("help", "Print options and their descriptions")
		("boundary-help", "Print help for options of initial and boundary condition(s)")
		("boundary-file",
			boost::program_options::value<std::string>(&boundary_file_name)
				->default_value(""),
			"Read boundaries from file arg (do not read if empty string)");

	Box<
		int,
		std::array<double, 3>,
		Momentum_Density,
		Magnetic_Field
	> boundaries;

	boundaries.add_boundary(
		{{0, 0, 0}},
		{{1, 1, 1}},
		Momentum_Density{},
		Momentum_Density::data_type{1, 2, 3},
		Magnetic_Field{},
		Magnetic_Field::data_type{4, 5, 6}
	);
	boundaries.add_boundary(
		{{0.5, 0.5, 0.5}},
		{{2, 2, 2}},
		Momentum_Density{},
		Momentum_Density::data_type{-1, -2, -3},
		Magnetic_Field{},
		Magnetic_Field::data_type{-4, -5, -6}
	);

	boost::program_options::options_description boundary_options(
		"Options for initial and boundary conditions"
	);
	boundaries.add_options(boundary_options, "test.");

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
	if (option_variables.count("boundary-help") > 0) {
		std::cout << boundary_options << std::endl;
		return EXIT_SUCCESS;
	}

	if (boundary_file_name != "") {
		boost::program_options::store(
			boost::program_options::parse_config_file<char>(
				boundary_file_name.c_str(),
				boundary_options
			),
			option_variables
		);
		boost::program_options::notify(option_variables);
	}


	if (boundaries.get_number_of_boundaries() != 2) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
			<< "Invalid number of test boundaries."
			<< std::endl;
		return EXIT_FAILURE;
	}

	if (boundaries.add_cell(-2, {{-0.6, -0.5, -0.4}}, {{-0.1, -0.2, -0.3}}) != 0) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
			<< "Cell added to invalid number of boundaries."
			<< std::endl;
		return EXIT_FAILURE;
	}
	if (boundaries.add_cell(-1, {{0.1, 0.2, 0.3}}, {{0.6, 0.5, 0.4}}) != 1) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
			<< "Cell added to invalid number of boundaries."
			<< std::endl;
		return EXIT_FAILURE;
	}
	if (boundaries.add_cell(0, {{0.1, 0.2, 0.3}}, {{0.6, 0.7, 0.8}}) != 2) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
			<< "Cell added to invalid number of boundaries."
			<< std::endl;
		return EXIT_FAILURE;
	}
	if (boundaries.add_cell(1, {{1.1, 1.2, 1.3}}, {{1.6, 1.7, 1.8}}) != 1) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
			<< "Cell added to invalid number of boundaries."
			<< std::endl;
		return EXIT_FAILURE;
	}

	if (boundaries.get_boundary_cells(0).size() != 2) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
			<< "Invalid number of cells in boundary 0."
			<< std::endl;
		return EXIT_FAILURE;
	}

	if (boundaries.get_boundary_cells(1).size() != 2) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
			<< "Invalid number of cells in boundary 1."
			<< std::endl;
		return EXIT_FAILURE;
	}

	if (boundaries.get_boundary_cells(0)[0] != -1) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
			<< "Invalid first cell in boundary 0."
			<< std::endl;
		return EXIT_FAILURE;
	}
	if (boundaries.get_boundary_cells(0)[1] != 0) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
			<< "Invalid second cell in boundary 0."
			<< std::endl;
		return EXIT_FAILURE;
	}

	if (boundaries.get_boundary_cells(1)[0] != 0) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
			<< "Invalid first cell in boundary 1."
			<< std::endl;
		return EXIT_FAILURE;
	}
	if (boundaries.get_boundary_cells(1)[1] != 1) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
			<< "Invalid second cell in boundary 1."
			<< std::endl;
		return EXIT_FAILURE;
	}

	if (boundaries.get_boundary_data(Momentum_Density(), 0)[0] != 1) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
			<< "Invalid boundary data in boundary 0."
			<< std::endl;
		return EXIT_FAILURE;
	}
	if (boundaries.get_boundary_data(Magnetic_Field(), 1)[2] != -6) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
			<< "Invalid boundary data in boundary 1."
			<< std::endl;
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}
