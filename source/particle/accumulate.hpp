/*
Functions for accumulating particle data to grid cells.

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

#ifndef PAMHD_PARTICLE_ACCUMULATE_HPP
#define PAMHD_PARTICLE_ACCUMULATE_HPP


#include "exception"


namespace pamhd {
namespace particle {


/*!
Returns portion of value spanning value_min/max that is within cell_min/max.

value is assumed to have a constant value within its volume.
*/
template<class Scalar, class Vector> Scalar get_accumulated_value(
	const Scalar& value,
	const Vector& value_min,
	const Vector& value_max,
	const Vector& cell_min,
	const Vector& cell_max
) {
	using std::max;
	using std::min;

	for (size_t dim = 0; dim < size_t(value_min.size()); dim++) {
		if (value_min[dim] > value_max[dim]) {
			throw std::out_of_range("Volume of given value ends before it starts");
		}
	}
	for (size_t dim = 0; dim < size_t(cell_min.size()); dim++) {
		if (cell_min[dim] > cell_max[dim]) {
			throw std::out_of_range("Volume of given cell ends before it starts");
		}
	}

	const Vector
		interval_min{
			max(value_min[0], cell_min[0]),
			max(value_min[1], cell_min[1]),
			max(value_min[2], cell_min[2])
		},
		interval_max{
			min(value_max[0], cell_max[0]),
			min(value_max[1], cell_max[1]),
			min(value_max[2], cell_max[2])
		},
		intersection{
			max(
				decltype(interval_max[0] - interval_min[0]){0},
				interval_max[0] - interval_min[0]
			),
			max(
				decltype(interval_max[1] - interval_min[1]){0},
				interval_max[1] - interval_min[1]
			),
			max(
				decltype(interval_max[2] - interval_min[2]){0},
				interval_max[2] - interval_min[2]
			)
		};

	const auto
		value_vol
			= (value_max[0] - value_min[0])
			* (value_max[1] - value_min[1])
			* (value_max[2] - value_min[2]),
		intersection_vol = intersection[0] * intersection[1] * intersection[2];

	return value * intersection_vol / value_vol;
}


}} // namespaces

#endif // ifndef PAMHD_PARTICLE_ACCUMULATE_HPP
