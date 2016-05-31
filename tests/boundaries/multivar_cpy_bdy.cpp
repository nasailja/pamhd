/*
Tests multivariable copy boundaries class of PAMHD.

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

#include "boundaries/multivariable_copy_boundaries.hpp"


using namespace pamhd::boundaries;


struct Mass_Density {
	using data_type = double;
	static const std::string get_name(){ return {"mass density"}; }
	static const std::string get_option_name(){ return {"mass_density"}; }
};

struct Momentum_Density {
	using data_type = std::array<double, 4>;
	static const std::string get_name(){ return {"momentum density"}; }
	static const std::string get_option_name(){ return {"momentum_density"}; }
};

int main()
{
	const char json[] = "{"
		"\"mass_density\": {"
			"\"copy_boundaries\": ["
				"{\"geometry_id\": 0},"
				"{\"geometry_id\": 1}"
			"]"
		"},"
		"\"momentum_density\": {"
			"\"copy_boundaries\": ["
				"{\"geometry_id\": 2},"
				"{\"geometry_id\": 4},"
				"{\"geometry_id\": 6}"
			"]"
		"}"
	"}";

	rapidjson::Document document;
	document.Parse(json);
	if (document.HasParseError()) {
		std::cerr << "Couldn't parse json data: " << json << std::endl;
		return EXIT_FAILURE;
	}

	Multivariable_Copy_Boundaries<
		int,
		size_t,
		Mass_Density,
		Momentum_Density
	> boundary;
	boundary.set(document);

	constexpr Mass_Density Mass{};
	constexpr Momentum_Density Momentum{};

	if (boundary.get_geometry_ids(Mass).size() != 2) {
		std::cerr << __FILE__ "(" << __LINE__ << "): Wrong number of boundaries: "
			<< boundary.get_geometry_ids(Mass).size() << ", should be 2"
			<< std::endl;
		return EXIT_FAILURE;
	}

	if (boundary.get_geometry_ids(Momentum).size() != 3) {
		std::cerr << __FILE__ "(" << __LINE__ << "): Wrong number of boundaries: "
			<< boundary.get_geometry_ids(Momentum).size() << ", should be 3"
			<< std::endl;
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}
