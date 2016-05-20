/*
Tests initial condition class of PAMHD with scalar simulation variable.

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
#include "string"
#include "vector"

#include "boundaries/initial_condition.hpp"


using namespace pamhd::boundaries;


struct Mass_Density {
	using data_type = double;
	static const std::string get_name(){ return {"mass density"}; }
};

int main()
{
	const char json[] = "["
		"{"
			"\"geometry_id\": 1,"
			"\"value\": 1"
		"},"
		"{"
			"\"geometry_id\": 3,"
			"\"value\": \"t\""
		"},"
		"{"
			"\"geometry_id\": 8,"
			"\"value\": {"
				"\"x\": [1],"
				"\"y\": [-2],"
				"\"z\": [0.1],"
				"\"data\": [-1]"
			"}"
		"},"
		"{"
			"\"geometry_id\": 8,"
			"\"value\": {"
				"\"x\": [1],"
				"\"y\": [-2, -1, 1],"
				"\"z\": [0.1],"
				"\"data\": [-1, 0, 1]"
			"}"
		"},"
		"{"
			"\"geometry_id\": 8,"
			"\"value\": {"
				"\"x\": [1, 2],"
				"\"y\": [-2, 3],"
				"\"z\": [0.1, 1.1],"
				"\"data\": [-1, 0, 1, 2, 3, 4, 5, 6]"
			"}"
		"}"
	"]";

	rapidjson::Document document;
	document.Parse(json);
	if (document.HasParseError()) {
		std::cerr << "Couldn't parse json data: " << json << std::endl;
		return EXIT_FAILURE;
	}

	typename Mass_Density::data_type mass;

	Initial_Condition<unsigned int, Mass_Density> init0, init1, init2, init3, init4;

	/*{
		"geometry_id": 1,
		"value": 1
	}*/
	init0.set(document[0]);
	mass = init0.get_data(0, 0, 0, 0, 0, 0, 0);
	if (mass != 1) {
		std::cerr << __FILE__ "(" << __LINE__ << "): Wrong value for variable: "
			<< mass << ", should be 1"
			<< std::endl;
		return EXIT_FAILURE;
	}


	/*{"
		"geometry_id": 3,
		"value": "t"
	}*/
	init1.set(document[1]);
	mass = init1.get_data(0.125, 0, 0, 0, 0, 0, 0);
	if (mass != 0.125) {
		std::cerr << __FILE__ "(" << __LINE__ << "): Wrong value for variable: "
			<< mass << ", should be 0.125"
			<< std::endl;
		return EXIT_FAILURE;
	}
	mass = init1.get_data(1, 2, 3, 4, 5, 6, 7);
	if (mass != 1) {
		std::cerr << __FILE__ "(" << __LINE__ << "): Wrong value for variable: "
			<< mass << ", should be 0.125"
			<< std::endl;
		return EXIT_FAILURE;
	}


	/*{
		"geometry_id": 8,
		"value": {
			"x": [1],
			"y": [-2],
			"z": [0.1],
			"data": [-1]
		}
	}*/
	init2.set(document[2]);
	mass = init2.get_data(0, 0, -3, 0, 0, 0, 0);
	if (mass != -1) {
		std::cerr << __FILE__ "(" << __LINE__ << "): Wrong value for variable: "
			<< mass << ", should be -1"
			<< std::endl;
		return EXIT_FAILURE;
	}
	mass = init2.get_data(0, 0, -1.51, 0, 0, 0, 0);
	if (mass != -1) {
		std::cerr << __FILE__ "(" << __LINE__ << "): Wrong value for variable: "
			<< mass << ", should be -1"
			<< std::endl;
		return EXIT_FAILURE;
	}


	init3.set(document[3]);
	mass = init3.get_data(0, 0, -1.51, 0, 0, 0, 0);
	if (mass != -1) {
		std::cerr << __FILE__ "(" << __LINE__ << "): Wrong value for variable: "
			<< mass << ", should be -1"
			<< std::endl;
		return EXIT_FAILURE;
	}
	mass = init3.get_data(0, 0, 4, 0, 0, 0, 0);
	if (mass != 1) {
		std::cerr << __FILE__ "(" << __LINE__ << "): Wrong value for variable: "
			<< mass << ", should be 1"
			<< std::endl;
		return EXIT_FAILURE;
	}
	mass = init3.get_data(0, 0, 0.01, 0, 0, 0, 0);
	if (mass != 1) {
		std::cerr << __FILE__ "(" << __LINE__ << "): Wrong value for variable: "
			<< mass << ", should be 1"
			<< std::endl;
		return EXIT_FAILURE;
	}
	mass = init3.get_data(0, 123, -0.1, 321, 0, 0, 0);
	if (mass != 0) {
		std::cerr << __FILE__ "(" << __LINE__ << "): Wrong value for variable: "
			<< mass << ", should be 0"
			<< std::endl;
		return EXIT_FAILURE;
	}


	/*{
		"geometry_id": 8,
		"value": {
			"x": [1, 2],
			"y": [-2, 3],
			"z": [0.1, 1.1],
			"data": [-1, 0, 1, 2, 3, 4, 5, 6]
		}
	}*/
	init4.set(document[4]);
	mass = init4.get_data(0, 0, 0, 0, 0, 0, 0);
	if (mass != -1) {
		std::cerr << __FILE__ "(" << __LINE__ << "): Wrong value for variable: "
			<< mass << ", should be -1"
			<< std::endl;
		return EXIT_FAILURE;
	}
	mass = init4.get_data(0, 5, 5, 5, 0, 0, 0);
	if (mass != 6) {
		std::cerr << __FILE__ "(" << __LINE__ << "): Wrong value for variable: "
			<< mass << ", should be 6"
			<< std::endl;
		return EXIT_FAILURE;
	}
	mass = init4.get_data(0, 1.6, 0, 0, 0, 0, 0);
	if (mass != 0) {
		std::cerr << __FILE__ "(" << __LINE__ << "): Wrong value for variable: "
			<< mass << ", should be 0"
			<< std::endl;
		return EXIT_FAILURE;
	}
	mass = init4.get_data(0, 0, 2, 0, 0, 0, 0);
	if (mass != 1) {
		std::cerr << __FILE__ "(" << __LINE__ << "): Wrong value for variable: "
			<< mass << ", should be 1"
			<< std::endl;
		return EXIT_FAILURE;
	}
	mass = init4.get_data(0, 0, 0, 1, 0, 0, 0);
	if (mass != 3) {
		std::cerr << __FILE__ "(" << __LINE__ << "): Wrong value for variable: "
			<< mass << ", should be 3"
			<< std::endl;
		return EXIT_FAILURE;
	}
	mass = init4.get_data(0, 2, 0, 1, 0, 0, 0);
	if (mass != 4) {
		std::cerr << __FILE__ "(" << __LINE__ << "): Wrong value for variable: "
			<< mass << ", should be 4"
			<< std::endl;
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}
