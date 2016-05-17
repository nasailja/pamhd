/*
Box boundary geometry class of PAMHD.

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


#ifndef PAMHD_BOUNDARIES_BOX_HPP
#define PAMHD_BOUNDARIES_BOX_HPP


#include "stdexcept"
#include "vector"

#include "rapidjson/document.h"


namespace pamhd {
namespace boundaries {


/*!
Vector is assumed to be std::array or similar.
*/
template<class Vector, class Cell_Id> class Box
{
public:

	Box() :
		start(start_rw),
		end(end_rw),
		cells(cells_rw)
	{}

	Box(
		const Vector& box_start,
		const Vector& box_end
	) :
		start(start_rw),
		end(end_rw),
		cells(cells_rw)
	{
		this->set_geometry(box_start, box_end);
	}

	Box(const rapidjson::Value& object) :
		start(start_rw),
		end(end_rw),
		cells(cells_rw)
	{
		this->set_geometry(object);
	}

	Box(const Box<Vector, Cell_Id>& other) :
		start(start_rw),
		end(end_rw),
		cells(cells_rw),
		start_rw(other.start_rw),
		end_rw(other.end_rw),
		cells_rw(other.cells_rw)
	{}

	Box(Box<Vector, Cell_Id>&& other) :
		start(start_rw),
		end(end_rw),
		cells(cells_rw),
		start_rw(std::move(other.start_rw)),
		end_rw(std::move(other.end_rw)),
		cells_rw(std::move(other.cells_rw))
	{}


	Box& operator =(const Box<Vector, Cell_Id>& other)
	{
		if (this != &other) {
			this->start_rw = other.start_rw;
			this->end_rw = other.end_rw;
			this->cells_rw = other.cells_rw;
		}

		return *this;
	}


	Box& operator =(Box<Vector, Cell_Id>&& other)
	{
		if (this != &other) {
			this->start_rw = std::move(other.start_rw);
			this->end_rw = std::move(other.end_rw);
			this->cells_rw = std::move(other.cells_rw);
		}

		return *this;
	}


	/*!
	Returns true if given simulation cell overlaps this box, false otherwise.

	If given cell overlaps it's id is added to this->cells.
	*/
	bool overlaps(
		const Vector& cell_start,
		const Vector& cell_end,
		const Cell_Id& cell_id
	) {
		for (size_t i = 0; i < size_t(cell_start.size()); i++) {
			if (cell_start[i] > cell_end[i]) {
				throw std::invalid_argument(__FILE__ ": Given cell ends before it starts.");
			}
		}

		bool overlaps = true;
		for (size_t i = 0; i < size_t(cell_start.size()); i++) {
			if (
				cell_start[i] >= this->end[i]
				or cell_end[i] <= this->start[i]
			) {
				overlaps = false;
				break;
			}
		}

		if (overlaps) {
			this->cells_rw.push_back(cell_id);
		}

		return overlaps;
	}


	void clear_cells() {
		this->cells_rw.clear();
	}


	void set_geometry(
		const Vector& given_start,
		const Vector& given_end
	) {
		for (size_t i = 0; i < size_t(given_start.size()); i++) {
			if (given_end[i] <= given_start[i]) {
				throw std::invalid_argument(__FILE__ ": Given box ends before it starts.");
			}
		}

		this->start_rw = given_start;
		this->end_rw = given_end;
	}


	/*!
	Sets box geometry from given rapidjson object.

	Given object must have "start" and "end" items
	and both must be arrays of 3 numbers.
	*/
	void set_geometry(const rapidjson::Value& object)
	{
		const auto& obj_start_i = object.FindMember("start");
		if (obj_start_i == object.MemberEnd()) {
			throw std::invalid_argument(__FILE__ ": Given object doesn't have a start item.");
		}
		const auto& obj_start = obj_start_i->value;

		if (not obj_start.IsArray()) {
			throw std::invalid_argument(__FILE__ ": Start item isn't an array.");
		}

		if (obj_start.Size() != 3) {
			throw std::invalid_argument(__FILE__ ": Start item doesn't have a length of 3.");
		}

		this->start_rw[0] = obj_start[0].GetDouble();
		this->start_rw[1] = obj_start[1].GetDouble();
		this->start_rw[2] = obj_start[2].GetDouble();

		const auto& obj_end_i = object.FindMember("end");
		if (obj_end_i == object.MemberEnd()) {
			throw std::invalid_argument(__FILE__ ": Given object doesn't have an end item.");
		}
		const auto& obj_end = obj_end_i->value;

		if (not obj_end.IsArray()) {
			throw std::invalid_argument(__FILE__ ": End item isn't an array.");
		}

		if (obj_end.Size() != 3) {
			throw std::invalid_argument(__FILE__ ": End item doesn't have a length of 3.");
		}

		this->end_rw[0] = obj_end[0].GetDouble();
		this->end_rw[1] = obj_end[1].GetDouble();
		this->end_rw[2] = obj_end[2].GetDouble();

		for (size_t i = 0; i < 3; i++) {
			if (this->end_rw[i] <= this->start_rw[i]) {
				throw std::invalid_argument(__FILE__ ": Given object ends before it starts.");
			}
		}
	}


	const Vector
		//! start coordinates of box
		&start,
		//! end coordinates of box
		&end;

	//! simulation cells overlapping with this box
	const std::vector<Cell_Id>& cells;


private:

	Vector start_rw, end_rw;

	std::vector<Cell_Id> cells_rw;
};

}} // namespaces

#endif // ifndef PAMHD_BOUNDARIES_BOX_HPP
