/*
Tests creation of box volumes from json data.

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

#include "prettyprint.hpp"

#include "boundaries/boxes_from_json.hpp"


int main()
{
	const char json[] = "{"
		"\"boxes\": ["
			"{"
				"\"start\": [0, 0, 0],"
				"\"end\": [1, 1, 1]"
			"},"
			"{"
				"\"start\": [-1, -2, -3],"
				"\"end\": [3.5, 2.25, 1.125]"
			"}"
		"]"
	"}";

	rapidjson::Document document;
	document.Parse(json);
	if (document.HasParseError()) {
		std::cerr << "Couldn't parse json data: " << json << std::endl;
		return EXIT_FAILURE;
	}

	const auto boxes = pamhd::boundaries::boxes_from_json<std::array<double, 3>>(document);

	if (boxes.size() != 2) {
		std::cerr << "Incorrect number of boxes: " << boxes.size() << std::endl;
		return EXIT_FAILURE;
	}


	if (boxes[0].start[0] != 0 or boxes[0].start[1] != 0 or boxes[0].start[2] != 0) {
		std::cerr << "Incorrect start coordinate for 1st box: " << boxes[0].start
			<< ", should be at [0, 0, 0]"
			<< std::endl;
		return EXIT_FAILURE;
	}

	if (boxes[0].end[0] != 1 or boxes[0].end[1] != 1 or boxes[0].end[2] != 1) {
		std::cerr << "Incorrect start coordinate for 1st box: " << boxes[0].start
			<< ", should be at [0, 0, 0]"
			<< std::endl;
		return EXIT_FAILURE;
	}


	if (boxes[1].start[0] != -1 or boxes[1].start[1] != -2 or boxes[1].start[2] != -3) {
		std::cerr << "Incorrect start coordinate for 1st box: " << boxes[0].start
			<< ", should be at [-1, -2, -3]"
			<< std::endl;
		return EXIT_FAILURE;
	}

	if (boxes[1].end[0] != 3.5 or boxes[1].end[1] != 2.25 or boxes[1].end[2] != 1.125) {
		std::cerr << "Incorrect start coordinate for 1st box: " << boxes[0].start
			<< ", should be at [3.5, 2.25, 1.125]"
			<< std::endl;
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}
