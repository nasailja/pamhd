/*
Sphere boundary geometry class of PAMHD.

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


#ifndef PAMHD_BOUNDARIES_SPHERE_HPP
#define PAMHD_BOUNDARIES_SPHERE_HPP


#include "cmath"
#include "vector"


namespace pamhd {
namespace boundaries {

/*!
Vector_T is assumed to be std::array or similar.
*/
template<class Vector, class Scalar, class Cell_Id> class Sphere
{
public:


	Sphere() :
		center(center_rw),
		radius(radius_rw),
		cells(cells_rw)
	{}

	Sphere(
		const Vector& sphere_center,
		const Scalar& sphere_radius
	) :
		center(center_rw),
		radius(radius_rw),
		cells(cells_rw)
	{
		this->set_geometry(sphere_center, sphere_radius);
	}

	Sphere(const rapidjson::Value& object) :
		center(center_rw),
		radius(radius_rw),
		cells(cells_rw)
	{
		this->set_geometry(object);
	}

	Sphere(const Sphere<Vector, Scalar, Cell_Id>& other) :
		center(center_rw),
		radius(radius_rw),
		cells(cells_rw),
		center_rw(other.center_rw),
		radius_rw(other.radius_rw),
		cells_rw(other.cells_rw)
	{}

	Sphere(Sphere<Vector, Scalar, Cell_Id>&& other) :
		center(center_rw),
		radius(radius_rw),
		cells(cells_rw),
		center_rw(std::move(other.center_rw)),
		radius_rw(std::move(other.radius_rw)),
		cells_rw(std::move(other.cells_rw))
	{}


	Sphere& operator =(const Sphere<Vector, Scalar, Cell_Id>& other)
	{
		if (this != &other) {
			this->center_rw = other.center_rw;
			this->radius_rw = other.radius_rw;
			this->cells_rw = other.cells_rw;
		}

		return *this;
	}


	Sphere& operator =(Sphere<Vector, Scalar, Cell_Id>&& other)
	{
		if (this != &other) {
			this->center_rw = std::move(other.center_rw);
			this->radius_rw = std::move(other.radius_rw);
			this->cells_rw = std::move(other.cells_rw);
		}

		return *this;
	}


	/*!
	Returns true if cell spanning given volume overlaps
	this sphere, false otherwise.
	*/
	bool overlaps(
		const Vector& cell_start,
		const Vector& cell_end,
		const Cell_Id& cell_id
	) {
		for (size_t i = 0; i < size_t(cell_start.size()); i++) {
			if (cell_start[i] > cell_end[i]) {
				throw std::invalid_argument(__FILE__ ": Given box ends before it starts.");
			}
		}

		// https://stackoverflow.com/questions/4578967
		Scalar distance2 = 0;
		for (size_t i = 0; i < size_t(cell_start.size()); i++) {
			if (this->center[i] < cell_start[i]) {
				distance2 += std::pow(this->center[i] - cell_start[i], 2);
			} else if (this->center[i] > cell_end[i]) {
				distance2 += std::pow(this->center[i] - cell_end[i], 2);
			}
		}

		const bool overlaps = distance2 < std::pow(this->radius, 2);

		if (overlaps) {
			this->cells_rw.push_back(cell_id);
		}

		return overlaps;
	}


	void set_geometry(
		const Vector& given_center,
		const Scalar& given_radius
	) {
		if (given_radius <= 0) {
			throw std::invalid_argument(__FILE__ ": Non-positive radius given.");
		}

		this->center_rw = given_center;
		this->radius_rw = given_radius;
	}


	/*!
	Sets sphere geometry from given rapidjson object.

	Given object must have "center" and "radius" items
	of array of 3 numbers and double respectively.
	*/
	void set_geometry(const rapidjson::Value& object)
	{
		const auto& obj_center_i = object.FindMember("center");
		if (obj_center_i == object.MemberEnd()) {
			throw std::invalid_argument(__FILE__ ": Given object doesn't have a center item.");
		}
		const auto& obj_center = obj_center_i->value;

		if (not obj_center.IsArray()) {
			throw std::invalid_argument(__FILE__ ": Center item isn't an array.");
		}

		if (obj_center.Size() != 3) {
			throw std::invalid_argument(__FILE__ ": Center item doesn't have a length of 3.");
		}

		this->center_rw[0] = obj_center[0].GetDouble();
		this->center_rw[1] = obj_center[1].GetDouble();
		this->center_rw[2] = obj_center[2].GetDouble();

		const auto& obj_radius_i = object.FindMember("radius");
		if (obj_radius_i == object.MemberEnd()) {
			throw std::invalid_argument(__FILE__ ": Given object doesn't have an radius item.");
		}
		const auto& obj_radius = obj_radius_i->value;

		if (not obj_radius.IsNumber()) {
			throw std::invalid_argument(__FILE__ ": Radius item isn't a number.");
		}

		this->radius_rw = obj_radius.GetDouble();

		if (this->radius_rw <= 0) {
			throw std::invalid_argument(__FILE__ ": Given object's radius isn't positive.");
		}
	}


	void clear_cells() {
		this->cells_rw.clear();
	}


	//! center coordinate of sphere
	const Vector &center;
	//! radius of sphere
	const Scalar& radius;

	//! simulation cells overlapping with this sphere
	const std::vector<Cell_Id>& cells;


private:

	Vector center_rw;
	Scalar radius_rw;

	std::vector<Cell_Id> cells_rw;
};

}} // namespaces

#endif // ifndef PAMHD_BOUNDARIES_SPHERE_HPP
