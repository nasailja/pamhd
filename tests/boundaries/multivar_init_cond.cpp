/*
Tests multivariable initial conditions class of PAMHD.

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

#include "rapidjson/document.h"

#include "boundaries/multivariable_initial_conditions.hpp"


using namespace pamhd::boundaries;


struct Mass_Density {
	using data_type = double;
	static const std::string get_name(){ return {"mass density"}; }
	static const std::string get_option_name(){ return {"mass_density"}; }
};
struct Momentum_Density {
	using data_type = std::array<int, 4>;
	static const std::string get_name(){ return {"momentum density"}; }
	static const std::string get_option_name(){ return {"momentum_density"}; }
};

int main()
{
	constexpr Mass_Density Mass{};
	constexpr Momentum_Density Momentum{};

	typename Mass_Density::data_type mass;
	typename Momentum_Density::data_type momentum;

	const char json1[] = "{"
		"\"mass_density\": {"
			"\"default\": -3,"
			"\"regions\": ["
				"{"
					"\"geometry_id\": 1,"
					"\"value\": 1"
				"}"
			"]"
		"},"
		"\"momentum_density\": {"
			"\"default\": [-3, -2, -1, 0],"
			"\"regions\": ["
				"{"
					"\"geometry_id\": 1,"
					"\"value\": [1, 2, 4, 6]"
				"}"
			"]"
		"}"
	"}";
	rapidjson::Document doc1;
	doc1.Parse(json1);
	if (doc1.HasParseError()) {
		std::cerr << "Couldn't parse json data: " << json1 << std::endl;
		return EXIT_FAILURE;
	}

	Multivariable_Initial_Conditions<
		unsigned int,
		Mass_Density,
		Momentum_Density
	> init1;

	init1.set(doc1);

	mass = init1.get_default_data(Mass, 0, 0, 0, 0, 0, 0, 0);
	if (mass != -3) {
		std::cerr << __FILE__ "(" << __LINE__ << "): Wrong default value for variable: "
			<< mass << ", should be -3"
			<< std::endl;
		return EXIT_FAILURE;
	}
	if (init1.get_number_of_regions(Mass) != 1) {
		std::cerr << __FILE__ "(" << __LINE__
			<< "): Wrong number of non-default initial conditions: "
			<< init1.get_number_of_regions(Mass) << ", should be 1"
			<< std::endl;
		return EXIT_FAILURE;
	}
	mass = init1.get_data(Mass, 0, 0, 0, 0, 0, 0, 0, 0);
	if (mass != 1) {
		std::cerr << __FILE__ "(" << __LINE__
			<< "): Wrong value for variable in non-default initial condition: "
			<< mass << ", should be 1"
			<< std::endl;
		return EXIT_FAILURE;
	}

	momentum = init1.get_default_data(Momentum, 0, 0, 0, 0, 0, 0, 0);
	if (momentum[2] != -1) {
		std::cerr << __FILE__ "(" << __LINE__ << "): Wrong default value for variable: "
			<< momentum[2] << ", should be -1"
			<< std::endl;
		return EXIT_FAILURE;
	}
	if (init1.get_number_of_regions(Momentum) != 1) {
		std::cerr << __FILE__ "(" << __LINE__
			<< "): Wrong number of non-default initial conditions: "
			<< init1.get_number_of_regions(Momentum) << ", should be 1"
			<< std::endl;
		return EXIT_FAILURE;
	}
	momentum = init1.get_data(Momentum, 0, 0, 0, 0, 0, 0, 0, 0);
	if (momentum[2] != 4) {
		std::cerr << __FILE__ "(" << __LINE__
			<< "): Wrong value for variable in non-default initial condition: "
			<< momentum[2] << ", should be 4"
			<< std::endl;
		return EXIT_FAILURE;
	}


	return EXIT_SUCCESS;
}
