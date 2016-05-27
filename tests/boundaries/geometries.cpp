/*
Boundary geometries tests of PAMHD.

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


#include "array"
#include "cstdlib"
#include "iostream"
#include "set"

#include "boundaries/geometries.hpp"


int main()
{
	const char json[] = "{"
		"\"geometries\": ["
			"{\"box\": {"
				"\"start\": [0, 0, 0],"
				"\"end\": [1, 1, 1]"
			"}},"
			"{\"sphere\": {"
				"\"center\": [-1, -2, -3],"
				"\"radius\": 0.01"
			"}},"
			"{\"box\": {"
				"\"start\": [-10, -20, -30],"
				"\"end\": [-3.5, -2.25, -1.125]"
			"}}"
		"]"
	"}";

	rapidjson::Document document;
	document.Parse(json);
	if (document.HasParseError()) {
		std::cerr << "Couldn't parse json data: " << json << std::endl;
		return EXIT_FAILURE;
	}


	pamhd::boundaries::Geometries<
		unsigned int,
		std::array<double, 3>,
		double,
		unsigned int
	> geometries;
	geometries.set(document);

	const auto ids = geometries.get_geometry_ids();
	if (ids.size() != 3) {
		std::cerr << "Incorrect number of geometry ids: " << ids.size()
			<< " but should be 3"
			<< std::endl;
		return EXIT_FAILURE;
	}
	const std::set<int> ids_set(ids.cbegin(), ids.cend());
	for (unsigned int id = 0; id < ids.size(); id++) {
		if (ids_set.count(id) == 0) {
			std::cerr << "geometry id " << id << " doesn't exist"
				<< std::endl;
			return EXIT_FAILURE;
		}
	}


	auto result = geometries.overlaps({0.1, 0.1, 0.1}, {0.2, 0.2, 0.2}, 1);
	if (not result.first) {
		std::cerr << "cell 1 didn't overlap with any geometry." << std::endl;
		return EXIT_FAILURE;
	}
	if (result.second != 0) {
		std::cerr << "cell 1 didn't overlap with geometry 0 but with: "
			<< result.second
			<< std::endl;
		return EXIT_FAILURE;
	}


	result = geometries.overlaps({-0.2, -0.2, -0.2}, {-0.1, -0.1, -0.1}, 2);
	if (result.first) {
		std::cerr << "cell 2 overlapped with a geometry." << std::endl;
		return EXIT_FAILURE;
	}


	result = geometries.overlaps({-1.001, -1.999, -3.005}, {-1, -1.9, -3}, 3);
	if (not result.first) {
		std::cerr << "cell 3 didn't overlap with any geometry." << std::endl;
		return EXIT_FAILURE;
	}
	if (result.second != 1) {
		std::cerr << "cell 3 didn't overlap with geometry 1 but with: "
			<< result.second
			<< std::endl;
		return EXIT_FAILURE;
	}

	geometries.overlaps({-1.001, -1.999, -3.005}, {-1, -1.9, -3}, 4);
	const auto& cells1 = geometries.get_cells(1);
	if (cells1.size() != 2) {
		std::cerr << "Wrong number of cells in geometry 1: "
			<< cells1.size() << ", should be 2."
			<< std::endl;
		return EXIT_FAILURE;
	}


	result = geometries.overlaps({-99, -99, -99}, {99, 99, 99}, 5);
	if (not result.first) {
		std::cerr << "cell 5 didn't overlap with any geometry." << std::endl;
		return EXIT_FAILURE;
	}
	if (result.second != 0) {
		std::cerr << "cell 5 didn't overlap with geometry 0 but with: "
			<< result.second
			<< std::endl;
		return EXIT_FAILURE;
	}
	const auto& cells0 = geometries.get_cells(0);
	if (cells0.size() != 2) {
		std::cerr << "Wrong number of cells in geometry 0: "
			<< cells0.size() << ", should be 2."
			<< std::endl;
		return EXIT_FAILURE;
	}


	return EXIT_SUCCESS;
}
