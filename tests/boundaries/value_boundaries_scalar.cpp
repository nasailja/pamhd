/*
Tests value boundaries class of PAMHD with scalar simulation variable.

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


#include "cstdlib"
#include "iostream"
#include "string"

#include "boundaries/value_boundaries.hpp"


using namespace pamhd::boundaries;


struct Mass_Density {
	using data_type = double;
	static const std::string get_name(){ return {"mass density"}; }
};

int main()
{
	const char json[] = "{\"value-boundaries\": ["
		"{"
			"\"geometry-id\": 1,"
			"\"time-stamps\": [-3, -2, -1],"
			"\"values\": [1, 2, 3]"
		"},"
		"{"
			"\"geometry-id\": 3,"
			"\"time-stamps\": [1, 4, 4.25],"
			"\"values\": [\"t\", \"t*t\", \"2*t\"]"
		"},"
		"{"
			"\"geometry-id\": 8,"
			"\"time-stamps\": [1, 3, 5],"
			"\"values\": {"
				"\"x\": [1],"
				"\"y\": [-2],"
				"\"z\": [0.1],"
				"\"data\": [-1, 8, -33]"
			"}"
		"},"
		"{"
			"\"geometry-id\": 2,"
			"\"time-stamps\": [-100, 100, 100100],"
			"\"values\": {"
				"\"x\": [1],"
				"\"y\": [-2, -1, 1],"
				"\"z\": [0.1],"
				"\"data\": [-1, 0, 1, 2, -2, -4, 10, 20, 30]"
			"}"
		"}"
	"]}";

	rapidjson::Document document;
	document.Parse(json);
	if (document.HasParseError()) {
		std::cerr << "Couldn't parse json data: " << json << std::endl;
		return EXIT_FAILURE;
	}

	typename Mass_Density::data_type mass;

	Value_Boundaries<unsigned int, Mass_Density> boundaries;
	boundaries.set(document);

	mass = boundaries.get_data(0, 0, 0, 0, 0, 0, 0, 0);
	if (mass != 3) {
		std::cerr << __FILE__ "(" << __LINE__ << "): Wrong value for variable: "
			<< mass << ", should be 3"
			<< std::endl;
		return EXIT_FAILURE;
	}
	mass = boundaries.get_data(0, -1.51, 0, 0, 0, 0, 0, 0);
	if (mass != 2) {
		std::cerr << __FILE__ "(" << __LINE__ << "): Wrong value for variable: "
			<< mass << ", should be 2"
			<< std::endl;
		return EXIT_FAILURE;
	}
	mass = boundaries.get_data(0, -123, 0, 0, 0, 0, 0, 0);
	if (mass != 1) {
		std::cerr << __FILE__ "(" << __LINE__ << "): Wrong value for variable: "
			<< mass << ", should be 1"
			<< std::endl;
		return EXIT_FAILURE;
	}


	mass = boundaries.get_data(1, 0.125, 0, 0, 0, 0, 0, 0);
	if (mass != 0.125) {
		std::cerr << __FILE__ "(" << __LINE__ << "): Wrong value for variable: "
			<< mass << ", should be 0.125"
			<< std::endl;
		return EXIT_FAILURE;
	}
	mass = boundaries.get_data(1, 2.75, 2, 3, 4, 5, 6, 7);
	if (mass != 7.5625) {
		std::cerr << __FILE__ "(" << __LINE__ << "): Wrong value for variable: "
			<< mass << ", should be 7.5625"
			<< std::endl;
		return EXIT_FAILURE;
	}
	mass = boundaries.get_data(1, 4.126, -2, -3, -4, -5, -6, -7);
	if (mass != 8.252) {
		std::cerr << __FILE__ "(" << __LINE__ << "): Wrong value for variable: "
			<< mass << ", should be 8.252"
			<< std::endl;
		return EXIT_FAILURE;
	}


	mass = boundaries.get_data(2, 1.99, 0, -3, 0, 0, 0, 0);
	if (mass != -1) {
		std::cerr << __FILE__ "(" << __LINE__ << "): Wrong value for variable: "
			<< mass << ", should be -1"
			<< std::endl;
		return EXIT_FAILURE;
	}
	mass = boundaries.get_data(2, 4.01, 0, -1.51, 0, 0, 0, 0);
	if (mass != -33) {
		std::cerr << __FILE__ "(" << __LINE__ << "): Wrong value for variable: "
			<< mass << ", should be -33"
			<< std::endl;
		return EXIT_FAILURE;
	}


	mass = boundaries.get_data(3, -1, 0, -1.51, 0, 0, 0, 0);
	if (mass != -1) {
		std::cerr << __FILE__ "(" << __LINE__ << "): Wrong value for variable: "
			<< mass << ", should be -1"
			<< std::endl;
		return EXIT_FAILURE;
	}
	mass = boundaries.get_data(3, 1, 0, -0.1, 0, 0, 0, 0);
	if (mass != -2) {
		std::cerr << __FILE__ "(" << __LINE__ << "): Wrong value for variable: "
			<< mass << ", should be -2"
			<< std::endl;
		return EXIT_FAILURE;
	}
	mass = boundaries.get_data(3, 100000, 0, -0.01, 0, 0, 0, 0);
	if (mass != 20) {
		std::cerr << __FILE__ "(" << __LINE__ << "): Wrong value for variable: "
			<< mass << ", should be 20"
			<< std::endl;
		return EXIT_FAILURE;
	}
	mass = boundaries.get_data(3, 100000, 0, 0.01, 0, 0, 0, 0);
	if (mass != 30) {
		std::cerr << __FILE__ "(" << __LINE__ << "): Wrong value for variable: "
			<< mass << ", should be 30"
			<< std::endl;
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}
