/*
Boundary geometries class of PAMHD

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


#ifndef PAMHD_GEOMETRIES_HPP
#define PAMHD_GEOMETRIES_HPP


#include "cstdlib"
#include "exception"
#include "iostream"
#include "limits"
#include "map"
#include "string"
#include "utility"

#include "rapidjson/document.h"

#include "boundaries/box.hpp"
#include "boundaries/sphere.hpp"


namespace pamhd {
namespace boundaries {

/*!
Vector_T is assumed to be std::array or similar.
*/
template<
	class Vector,
	class Scalar,
	class Cell_Id
> class Geometries
{
public:

	Geometries<Vector, Scalar, Cell_Id>() = default;
	Geometries<Vector, Scalar, Cell_Id>(const Geometries<Vector, Scalar, Cell_Id>&) = default;
	Geometries<Vector, Scalar, Cell_Id>(Geometries<Vector, Scalar, Cell_Id>&&) = default;

	Geometries<Vector, Scalar, Cell_Id>(const rapidjson::Value& object)
	{
		this->set_geometries(object);
	}


	/*!
	Creates geometries from given rapidjson object.

	Given object must be an array of objects each of which
	contains a supported geometry object.

	Example json from which "geometries" would be given as object:
	\verbatim
	{
		...,
		"geometries": [
			{"box": {
				"start": [0, 0, 0],
				"end": [1, 1, 1]
			}},
			{"sphere": {
				"center": [-1, -2, -3],
				"radius": 0.01
			}}
		],
		...
	}
	\endverbatim
	*/
	void set_geometries(const rapidjson::Value& object)
	{
		if (not object.IsArray()) {
			throw std::invalid_argument("geometries item is not an array (\"geometries\": [...]).");
		}

		for (rapidjson::SizeType i = 0; i < object.Size(); i++) {
			this->boxes.erase(i);
			this->spheres.erase(i);

			if (object[i].HasMember("box")) {
				this->boxes[i] = Box<Vector, Cell_Id>(object[i]["box"]);
			} else if (object[i].HasMember("sphere")) {
				this->spheres[i] = Sphere<Vector, Scalar, Cell_Id>(object[i]["sphere"]);
			} else {
				throw std::invalid_argument(__FILE__ ": Item in geometries array is not of a supported type.");
			}

		}
	}


	/*!
	Calls overlaps function of each existing geometry.

	If cube spanning given volume overlaps one or more geometry
	returns true as well as id of that geometry, otherwise
	returns false.
	*/
	std::pair<bool, unsigned int> overlaps(
		const Vector& cell_start,
		const Vector& cell_end,
		const Cell_Id& cell_id
	) {
		for (auto& box: this->boxes) {
			if (box.second.overlaps(cell_start, cell_end, cell_id)) {
				return std::make_pair(true, box.first);
			}
		}

		for (auto& sphere: this->spheres) {
			if (sphere.second.overlaps(cell_start, cell_end, cell_id)) {
				return std::make_pair(true, sphere.first);
			}
		}

		return std::make_pair(false, std::numeric_limits<unsigned int>::max());
	}


	/*!
	Returns cells belonging to given boundary id.

	Throws an exception if given boundary id doesn't exist.
	*/
	const std::vector<Cell_Id>& get_cells(const unsigned int boundary_id) const
	{
		if (this->boxes.count(boundary_id) > 0) {
			return this->boxes.at(boundary_id).cells;
		} else if (this->spheres.count(boundary_id) > 0) {
			return this->spheres.at(boundary_id).cells;
		}

		throw std::invalid_argument(__FILE__ ": Given boundary id doesn't exist.");
	}


	void clear_cells(const unsigned int boundary_id) {
		if (this->boxes.count(boundary_id) > 0) {
			this->boxes.at(boundary_id).second.clear();
		} else if (this->spheres.count(boundary_id) > 0) {
			this->spheres.at(boundary_id).clear();
		}
	}


	std::vector<unsigned int> get_geometry_ids()
	{
		std::vector<unsigned int> ids;
		ids.reserve(this->boxes.size() + this->spheres.size());

		for (const auto& box: this->boxes) {
			ids.push_back(box.first);
		}

		for (const auto& sphere: this->spheres) {
			ids.push_back(sphere.first);
		}

		return ids;
	}


private:

	std::map<
		unsigned int, // geometry id shared between all geometries
		Box<Vector, Cell_Id>
	> boxes;

	std::map<
		unsigned int,
		Sphere<Vector, Scalar, Cell_Id>
	> spheres;

};

}} // namespaces

#endif // ifndef PAMHD_GEOMETRIES_HPP
