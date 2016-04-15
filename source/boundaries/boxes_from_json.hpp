/*
Box geometry class of PAMHD.

Copyright 2014, 2015, 2016 Ilja Honkonen
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


#ifndef PAMHD_BOXES_FROM_JSON_HPP
#define PAMHD_BOXES_FROM_JSON_HPP


#include "cstdlib"
#include "iostream"
#include "string"

#include "rapidjson/document.h"

#include "boundaries/box.hpp"


namespace pamhd {
namespace boundaries {


template<
	class Vector_T
> std::vector<
	Box<Vector_T>
> boxes_from_json(
	const rapidjson::Document& document
) {
	std::vector<Box<Vector_T>> boxes;

	const auto boxes_it = document.FindMember("boxes");
	if (boxes_it != document.MemberEnd()) {
		if (not boxes_it->value.IsArray()) {
			throw std::invalid_argument("boxes item is not an array.");
		}

		for (rapidjson::SizeType i = 0; i < boxes_it->value.Size(); i++) {
			const auto start_it = boxes_it->value[i].FindMember("start");
			if (start_it == document.MemberEnd()) {
				continue;
			}

			if (not start_it->value.IsArray()) {
				throw std::invalid_argument("start of box item isn't an array.");
			}

			if (start_it->value.Size() != 3) {
				throw std::invalid_argument("start of box item doesn't have a length of 3.");
			}

			const Vector_T start{
				start_it->value[0].GetDouble(),
				start_it->value[1].GetDouble(),
				start_it->value[2].GetDouble()
			};

			const auto end_it = boxes_it->value[i].FindMember("end");
			if (end_it == document.MemberEnd()) {
				throw std::invalid_argument("box doesn't have an end coordinate.");
			}

			if (not end_it->value.IsArray()) {
				throw std::invalid_argument("end of box item isn't an array.");
			}

			if (end_it->value.Size() != 3) {
				throw std::invalid_argument("end of box item doesn't have a length of 3.");
			}

			const Vector_T end{
				end_it->value[0].GetDouble(),
				end_it->value[1].GetDouble(),
				end_it->value[2].GetDouble()
			};

			Box<Vector_T> new_box;
			if (not new_box.set_geometry(start, end)) {
				throw std::invalid_argument("Couldn't set box geometry.");
			}

			boxes.push_back(new_box);
		}
	}

	return boxes;
}

}} // namespaces

#endif // ifndef PAMHD_BOXES_FROM_JSON_HPP
