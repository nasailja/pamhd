/*
Parse test for grid options of PAMHD.

Copyright 2016 Ilja Honkonen
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


#include "boundaries/simulation_variable_expression.hpp"
#include "grid_options.hpp"
#include "prettyprint.hpp"


int main()
{
	const char json[] = "{"
		"\"grid\": {"
			"\"periodic\": \"{false, true, false}\","
			"\"nr-cells\": \"{1+1, 2+2, 3+3}\","
			"\"volume\": \"{4+4, 5+5, 6+6}\","
			"\"start\": \"{7-7, 8-8, 9-9}\""
		"}"
	"}";
	rapidjson::Document document;
	document.Parse(json);
	if (document.HasParseError()) {
		return EXIT_FAILURE;
	}

	pamhd::grid::Options options;
	options.set_expressions(document);

	const auto periodic = options.get_data(pamhd::grid::Periodic());
	if (periodic[0] != false or periodic[1] != true or periodic[2] != false) {
		std::cerr << "Incorrect number of simulation cells: "
			<< periodic << ", should be [0, 1, 0]"
			<< std::endl;
		return EXIT_FAILURE;
	}

	const auto nr_cells = options.get_data(pamhd::grid::Number_Of_Cells());
	if (nr_cells[0] != 2 or nr_cells[1] != 4 or nr_cells[2] != 6) {
		std::cerr << "Incorrect number of simulation cells: "
			<< nr_cells << ", should be [2, 4, 6]"
			<< std::endl;
		return EXIT_FAILURE;
	}

	const auto volume = options.get_data(pamhd::grid::Volume());
	if (volume[0] != 8 or volume[1] != 10 or volume[2] != 12) {
		std::cerr << "Incorrect simulation volume: "
			<< nr_cells << ", should be [8, 10, 12]"
			<< std::endl;
		return EXIT_FAILURE;
	}

	const auto start = options.get_data(pamhd::grid::Start());
	if (start[0] != 0 or start[1] != 0 or start[2] != 0) {
		std::cerr << "Incorrect simulation volume: "
			<< nr_cells << ", should be [0, 0, 0]"
			<< std::endl;
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}
