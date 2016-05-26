/*
Class for handling all copy boundaries of a simulation.

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

#ifndef PAMHD_BOUNDARIES_COPY_BOUNDARIES_HPP
#define PAMHD_BOUNDARIES_COPY_BOUNDARIES_HPP


#include "stdexcept"
#include "string"
#include "utility"
#include "vector"

#include "rapidjson/document.h"

#include "boundaries/copy_boundaries.hpp"


namespace pamhd {
namespace boundaries {


/*!
Collection of simulation copy boundaries creatable from json data.
*/
template<
	class Geometry_Id,
	class Variable
> class Copy_Boundaries
{
public:

	Copy_Boundaries<
		Geometry_Id,
		Variable
	>() : geometry_ids(geometry_ids_rw) {}

	Copy_Boundaries<
		Geometry_Id,
		Variable
	>(const Copy_Boundaries& other) :
		geometry_ids(geometry_ids_rw),
		geometry_ids_rw(other.geometry_ids_rw)
	{}

	Copy_Boundaries<
		Geometry_Id,
		Variable
	>(Copy_Boundaries&& other) :
		geometry_ids(geometry_ids_rw),
		geometry_ids_rw(std::move(other.geometry_ids_rw))
	{}


	/*!
	Prepares copy boundaries from given rapidjson object.

	Example json object:
	\verbatim
	{
		"copy_boundaries": [
			{"geometry_id": 0},
			{"geometry_id": 1}
		]
	}
	\endverbatim
	*/
	void set(const rapidjson::Value& object)
	{
		if (not object.HasMember("copy_boundaries")) {
			return;
		}
		const auto& json_bdys = object["copy_boundaries"];

		if (not json_bdys.IsArray()) {
			throw std::invalid_argument(__FILE__ ": copy_boundaries is not an array.");
		}

		this->geometry_ids_rw.clear();
		this->geometry_ids_rw.reserve(json_bdys.Size());

		for (size_t i = 0; i < json_bdys.Size(); i++) {
			const auto& json_bdy = json_bdys[i];
			if (not json_bdy.HasMember("geometry_id")) {
				throw std::invalid_argument(
					std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
					+ std::to_string(i + 1)
					+ "th object in copy_boundaries doesn't have a geometry_id item."
				);
			}

			this->geometry_ids_rw.push_back(json_bdy["geometry_id"].GetUint64());
		}
	}


	/*!
	Geometry id of each corresponding (by position in vector) copy boundary.
	*/
	const std::vector<Geometry_Id>& geometry_ids;

private:

	std::vector<Geometry_Id> geometry_ids_rw;

};

}} // namespaces

#endif // ifndef PAMHD_BOUNDARIES_COPY_BOUNDARIES_HPP
