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


#ifndef PAMHD_BOUNDARIES_BOX_HPP
#define PAMHD_BOUNDARIES_BOX_HPP


#include "cstdlib"
#include "iostream"
#include "string"

#include "boost/program_options.hpp"
#include "prettyprint.hpp"

#include "program_options_validators.hpp"


namespace pamhd {
namespace boundaries {


/*!
Vector_T is assumed to be std::array or similar.
*/
template<class Vector_T> class Box
{
public:

	Box() :
		start(start_rw),
		end(end_rw)
	{}

	Box(const Box<Vector_T>& other) :
		start(start_rw),
		end(end_rw),
		start_rw(other.start),
		end_rw(other.end)
	{}

	Box(Box<Vector_T>&& other) :
		start(start_rw),
		end(end_rw),
		start_rw(std::move(other.start_rw)),
		end_rw(std::move(other.end_rw))
	{}


	/*!
	Returns true if cube spanning given volume overlaps
	this box, false otherwise.
	*/
	bool overlaps(
		const Vector_T& cell_start,
		const Vector_T& cell_end
	) {
		for (size_t i = 0; i < size_t(cell_start.size()); i++) {
			if (cell_start[i] > cell_end[i]) {
				std::cerr <<  __FILE__ << "(" << __LINE__<< ") "
					<< "Starting coordinate of box at index " << i
					<< " is larger than ending coordinate: "
					<< cell_start[i] << " > " << cell_end[i]
					<< std::endl;
				abort();
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

		return overlaps;
	}


	bool set_geometry(
		const Vector_T& given_start,
		const Vector_T& given_end
	) {
		for (size_t i = 0; i < size_t(given_start.size()); i++) {
			if (given_end[i] <= given_start[i]) {
				return false;
			}
		}

		this->start_rw = given_start;
		this->end_rw = given_end;
		return true;
	}


	const Vector_T
		//! start coordinates of box
		&start,
		//! end coordinates of box
		&end;


private:

	Vector_T start_rw, end_rw;
};

}} // namespaces

#endif // ifndef PAMHD_BOUNDARIES_BOX_HPP
