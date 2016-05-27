/*
Tests value boundary class of PAMHD with vector simulation variable.

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

#include "boundaries/value_boundary.hpp"


using namespace pamhd::boundaries;


struct Momentum_Density {
	using data_type = std::array<double, 2>;
	static const std::string get_name(){ return {"momentum density"}; }
};

int main()
{
	const char json[] = "["
		"{"
			"\"geometry_id\": 1,"
			"\"time_stamps\": [-3, -2, -1],"
			"\"values\": [[1, -1], [2, -2], [3, -3]]"
		"},"
		"{"
			"\"geometry_id\": 3,"
			"\"time_stamps\": [1, 4, 4.25],"
			"\"values\": [\"{t, -t}\", \"{t*t, -t*t}\", \"{2*t, -2*t}\"]"
		"},"
		"{"
			"\"geometry_id\": 8,"
			"\"time_stamps\": [1, 3, 5],"
			"\"values\": {"
				"\"x\": [1],"
				"\"y\": [-2],"
				"\"z\": [0.1],"
				"\"data\": [[-1, 1], [8, -8], [-33, 33]]"
			"}"
		"},"
		"{"
			"\"geometry_id\": 8,"
			"\"time_stamps\": [-100, 100, 100100],"
			"\"values\": {"
				"\"x\": [1],"
				"\"y\": [-2, -1, 1],"
				"\"z\": [0.1],"
				"\"data\": [[-1, 1], [0, 0], [1, -1], [2, -2], "
					"[-2, 2], [-4, 4], [10, -10], [20, -20], [30, -30]]"
			"}"
		"}"
	"]";

	rapidjson::Document document;
	document.Parse(json);
	if (document.HasParseError()) {
		std::cerr << "Couldn't parse json data: " << json << std::endl;
		return EXIT_FAILURE;
	}

	typename Momentum_Density::data_type momentum;

	Value_Boundary<unsigned int, Momentum_Density> bdy0, bdy1, bdy2, bdy3;

	bdy0.set(document[0]);
	momentum = bdy0.get_data(0, 0, 0, 0, 0, 0, 0);
	if (momentum[1] != -3) {
		std::cerr << __FILE__ "(" << __LINE__ << "): Wrong value for variable: "
			<< momentum[1] << ", should be -3"
			<< std::endl;
		return EXIT_FAILURE;
	}
	momentum = bdy0.get_data(-1.51, 0, 0, 0, 0, 0, 0);
	if (momentum[1] != -2) {
		std::cerr << __FILE__ "(" << __LINE__ << "): Wrong value for variable: "
			<< momentum[1] << ", should be -2"
			<< std::endl;
		return EXIT_FAILURE;
	}
	momentum = bdy0.get_data(-123, 0, 0, 0, 0, 0, 0);
	if (momentum[1] != -1) {
		std::cerr << __FILE__ "(" << __LINE__ << "): Wrong value for variable: "
			<< momentum[1] << ", should be -1"
			<< std::endl;
		return EXIT_FAILURE;
	}


	bdy1.set(document[1]);
	momentum = bdy1.get_data(0.125, 0, 0, 0, 0, 0, 0);
	if (momentum[1] != -0.125) {
		std::cerr << __FILE__ "(" << __LINE__ << "): Wrong value for variable: "
			<< momentum[1] << ", should be -0.125"
			<< std::endl;
		return EXIT_FAILURE;
	}
	momentum = bdy1.get_data(2.75, 2, 3, 4, 5, 6, 7);
	if (momentum[1] != -7.5625) {
		std::cerr << __FILE__ "(" << __LINE__ << "): Wrong value for variable: "
			<< momentum[1] << ", should be -7.5625"
			<< std::endl;
		return EXIT_FAILURE;
	}
	momentum = bdy1.get_data(4.126, -2, -3, -4, -5, -6, -7);
	if (momentum[1] != -8.252) {
		std::cerr << __FILE__ "(" << __LINE__ << "): Wrong value for variable: "
			<< momentum[1] << ", should be -8.252"
			<< std::endl;
		return EXIT_FAILURE;
	}


	bdy2.set(document[2]);
	momentum = bdy2.get_data(1.99, 0, -3, 0, 0, 0, 0);
	if (momentum[1] != 1) {
		std::cerr << __FILE__ "(" << __LINE__ << "): Wrong value for variable: "
			<< momentum[1] << ", should be 1"
			<< std::endl;
		return EXIT_FAILURE;
	}
	momentum = bdy2.get_data(4.01, 0, -1.51, 0, 0, 0, 0);
	if (momentum[1] != 33) {
		std::cerr << __FILE__ "(" << __LINE__ << "): Wrong value for variable: "
			<< momentum[1] << ", should be 33"
			<< std::endl;
		return EXIT_FAILURE;
	}


	bdy3.set(document[3]);
	momentum = bdy3.get_data(-1, 0, -1.51, 0, 0, 0, 0);
	if (momentum[1] != 1) {
		std::cerr << __FILE__ "(" << __LINE__ << "): Wrong value for variable: "
			<< momentum[1] << ", should be 1"
			<< std::endl;
		return EXIT_FAILURE;
	}
	momentum = bdy3.get_data(1, 0, -0.1, 0, 0, 0, 0);
	if (momentum[1] != 2) {
		std::cerr << __FILE__ "(" << __LINE__ << "): Wrong value for variable: "
			<< momentum[1] << ", should be 2"
			<< std::endl;
		return EXIT_FAILURE;
	}
	momentum = bdy3.get_data(100000, 0, -0.01, 0, 0, 0, 0);
	if (momentum[1] != -20) {
		std::cerr << __FILE__ "(" << __LINE__ << "): Wrong value for variable: "
			<< momentum[1] << ", should be -20"
			<< std::endl;
		return EXIT_FAILURE;
	}
	momentum = bdy3.get_data(100000, 0, 0.01, 0, 0, 0, 0);
	if (momentum[1] != -30) {
		std::cerr << __FILE__ "(" << __LINE__ << "): Wrong value for variable: "
			<< momentum[1] << ", should be -30"
			<< std::endl;
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}
