/*
Tests copy boundary class of PAMHD.

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

#include "boost/optional.hpp"
#include "boost/program_options.hpp"

#include "boundaries/copy_boundary.hpp"


using namespace pamhd::boundaries;

#define MAKE_VEC(x, y, z) std::array<double, 3>{x, y, z}

int main(int argc, char* argv[])
{
	Copy_Boundary<
		int,
		double,
		std::array<double, 3>
	> boundary;

	//! Set default values
	boost::optional<size_t> result;

	result = boundary.add_box(MAKE_VEC(0, 0, 0), MAKE_VEC(1, 1, 1));
	if (not result or *result != 0) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
			<< "Invalid boundary id for added box: " << *result
			<< std::endl;
		abort();
	}

	result = boundary.add_sphere(MAKE_VEC(2, 2, 2), 2);
	if (not result or *result != 1) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
			<< "Invalid boundary id for added sphere: " << *result 
			<< std::endl;
		abort();
	}

	std::string copy_boundary_file_name("");

	boost::program_options::options_description
		general_options(
			"Usage: program_name [options], where options are"
		),
		copy_boundary_options(
			"Options related to copy boundary condition"
		);

	general_options.add_options()
		("help", "Print options and their descriptions")
		("copy-boundary-help", "Print help for options of copy boundary condition")
		("copy-boundary-file",
			boost::program_options::value<std::string>(&copy_boundary_file_name)
				->default_value(copy_boundary_file_name),
			"Read copy boundary condition from file arg (not read if empty string)");
	boundary.add_initialization_options("", general_options);

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

	boundary.add_options("", copy_boundary_options);

	if (option_variables.count("initial-help") > 0) {
		std::cout << copy_boundary_options << std::endl;
		return EXIT_SUCCESS;
	}

	if (copy_boundary_file_name != "") {
		boost::program_options::store(
			boost::program_options::parse_config_file<char>(
				copy_boundary_file_name.c_str(),
				copy_boundary_options
			),
			option_variables
		);
		boost::program_options::notify(option_variables);
	}

	if (boundary.add_cell(-1, MAKE_VEC(0.1, 0.1, 0.1), MAKE_VEC(0.9, 0.9, 0.9)) != 2) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< ")" << std::endl;
		abort();
	}
	if (boundary.add_cell(-2, MAKE_VEC(1.1, 1.1, 1.1), MAKE_VEC(1.9, 1.9, 1.9)) != 1) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< ")" << std::endl;
		abort();
	}
	if (boundary.add_cell(-3, MAKE_VEC(4, 4, 4), MAKE_VEC(5, 5, 5)) != 0) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< ")" << std::endl;
		abort();
	}
	if (boundary.add_cell(-4, MAKE_VEC(5, 5, 5), MAKE_VEC(6, 6, 6)) != 0) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< ")" << std::endl;
		abort();
	}

	boundary.add_as_other_boundary(-5);
	boundary.add_as_other_boundary(-6);

	if (boundary.add_cell(-7, MAKE_VEC(0, 0, 0), MAKE_VEC(1, 1, 1)) != 2) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< ")" << std::endl;
		abort();
	}
	if (boundary.set_neighbors_of(-1, {-2}) != 0) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< ")" << std::endl;
		abort();
	}
	if (boundary.set_neighbors_of(-2, {-1, -3}) != 1) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< ")" << std::endl;
		abort();
	}
	if (boundary.set_neighbors_of(-7, {-1, -2, -3, -4, -5, -6}) != 2) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< ")" << std::endl;
		abort();
	}

	return EXIT_SUCCESS;
}
