/*
Tests initial conditions class of PAMHD with vector simulation variable.

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
#include "string"
#include "vector"

#include "boundaries/initial_conditions.hpp"


using namespace pamhd::boundaries;


struct Momentum_Density {
	using data_type = std::array<double, 3>;
	static const std::string get_name(){ return {"momentum density"}; }
};

int main()
{
	typename Momentum_Density::data_type momentum;

	const char json1[] =
		"{\"default\": [-3, -4, -5],"
		"\"initial-conditions\": ["
			"{"
				"\"geometry-id\": 1,"
				"\"value\": [1, 3, 5]"
			"}"
		"]}";
	rapidjson::Document doc1;
	doc1.Parse(json1);
	if (doc1.HasParseError()) {
		std::cerr << "Couldn't parse json data: " << json1 << std::endl;
		return EXIT_FAILURE;
	}

	Initial_Conditions<unsigned int, Momentum_Density> init1;
	init1.set(doc1);
	momentum = init1.get_default_data(0, 0, 0, 0, 0, 0, 0);
	if (momentum[1] != -4) {
		std::cerr << __FILE__ "(" << __LINE__ << "): Wrong default value for variable: "
			<< momentum[1] << ", should be -3"
			<< std::endl;
		return EXIT_FAILURE;
	}
	if (init1.get_number_of_regions() != 1) {
		std::cerr << __FILE__ "(" << __LINE__
			<< "): Wrong number of non-default initial conditions: "
			<< init1.get_number_of_regions() << ", should be 1"
			<< std::endl;
		return EXIT_FAILURE;
	}
	momentum = init1.get_data(0, 0, 0, 0, 0, 0, 0, 0);
	if (momentum[1] != 3) {
		std::cerr << __FILE__ "(" << __LINE__
			<< "): Wrong value for variable in non-default initial condition: "
			<< momentum[1] << ", should be 3"
			<< std::endl;
		return EXIT_FAILURE;
	}


	const char json2[] =
		"{\"default\": \"{t+2, t+3, t+4}\","
		"\"initial-conditions\": ["
			"{"
				"\"geometry-id\": 1,"
				"\"value\": [1, 0, -1]"
			"}"
		"]}";
	rapidjson::Document doc2;
	doc2.Parse(json2);
	if (doc2.HasParseError()) {
		std::cerr << "Couldn't parse json data: " << json2 << std::endl;
		return EXIT_FAILURE;
	}

	Initial_Conditions<unsigned int, Momentum_Density> init2;
	init2.set(doc2);
	momentum = init2.get_default_data(3, 0, 0, 0, 0, 0, 0);
	if (momentum[1] != 6) {
		std::cerr << __FILE__ "(" << __LINE__ << "): Wrong default value for variable: "
			<< momentum[1] << ", should be 6"
			<< std::endl;
		return EXIT_FAILURE;
	}
	if (init2.get_number_of_regions() != 1) {
		std::cerr << __FILE__ "(" << __LINE__
			<< "): Wrong number of non-default initial conditions: "
			<< init2.get_number_of_regions() << ", should be 1"
			<< std::endl;
		return EXIT_FAILURE;
	}
	momentum = init2.get_data(0, 0, 0, 0, 0, 0, 0, 0);
	if (momentum[1] != 0) {
		std::cerr << __FILE__ "(" << __LINE__
			<< "): Wrong value for variable in non-default initial condition: "
			<< momentum[1] << ", should be 0"
			<< std::endl;
		return EXIT_FAILURE;
	}


	const char json3[] =
		"{\"default\": {"
			"\"x\": [1],"
			"\"y\": [-2],"
			"\"z\": [0.1],"
			"\"data\": [[-1, -2, -3]]"
		"},"
		"\"initial-conditions\": ["
			"{"
				"\"geometry-id\": 1,"
				"\"value\": [1, 2, 3]"
			"}"
		"]}";
	rapidjson::Document doc3;
	doc3.Parse(json3);
	if (doc3.HasParseError()) {
		std::cerr << "Couldn't parse json data: " << json3 << std::endl;
		return EXIT_FAILURE;
	}

	Initial_Conditions<unsigned int, Momentum_Density> init3;
	init3.set(doc3);
	momentum = init3.get_default_data(3, 0, 0, 0, 0, 0, 0);
	if (momentum[1] != -2) {
		std::cerr << __FILE__ "(" << __LINE__ << "): Wrong default value for variable: "
			<< momentum[1] << ", should be -2"
			<< std::endl;
		return EXIT_FAILURE;
	}
	if (init3.get_number_of_regions() != 1) {
		std::cerr << __FILE__ "(" << __LINE__
			<< "): Wrong number of non-default initial conditions: "
			<< init3.get_number_of_regions() << ", should be 1"
			<< std::endl;
		return EXIT_FAILURE;
	}
	momentum = init3.get_data(0, 0, 0, 0, 0, 0, 0, 0);
	if (momentum[1] != 2) {
		std::cerr << __FILE__ "(" << __LINE__
			<< "): Wrong value for variable in non-default initial condition: "
			<< momentum[1] << ", should be 2"
			<< std::endl;
		return EXIT_FAILURE;
	}


	return EXIT_SUCCESS;
}
