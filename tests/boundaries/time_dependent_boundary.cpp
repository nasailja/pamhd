/*
Test for time-dependent boundary condition class of PAMHD.

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

#include "boundaries/time_dependent_boundary.hpp"

using namespace pamhd::boundaries;

struct Is_Alive {
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


int main(int argc, char* argv[])
{
	Time_Dependent_Boundary<int, std::array<double, 1>, double, Is_Alive> boundaries;
	// set default boundary
	boundaries.add_boundary_data(
		{{0}},
		{{1}},
		0,
		Is_Alive{},
		Is_Alive::data_type{1}
	);
	boundaries.add_boundary_data(
		{{0}},
		{{1}},
		1,
		Is_Alive{},
		Is_Alive::data_type{0}
	);
	boundaries.add_boundary_data(
		{{1}},
		{{2}},
		2,
		Is_Alive{},
		Is_Alive::data_type{1}
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

	if (boundaries.add_cell(-1, {{0.1}}, {{0.6}}, -1) != 1) {
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


	if (boundaries.add_cell(0, {{0.1}}, {{0.6}}, 0.5) != 1) {
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


	if (boundaries.add_cell(1, {{0.1}}, {{0.6}}, 1.1) != 1) {
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


	if (boundaries.add_cell(2, {{0.1}}, {{0.6}}, 1.9) != 1) {
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


	if (boundaries.add_cell(3, {{0.1}}, {{0.6}}, 2.1) != 0) {
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


	if (boundaries.add_cell(4, {{1.1}}, {{1.6}}, 2.1) != 1) {
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
