/*
Tests for volume_range.

Copyright 2014 Ilja Honkonen
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
#include "cmath"
#include "cstdlib"
#include "iostream"

#include "volume_range.hpp"

using namespace std;
using namespace pamhd;

template <
	class Value_T,
	size_t Dimensions
> void check_number_of_samples(const size_t samples_per_dim)
{
	std::array<Value_T, Dimensions> start, end;

	for (size_t d = 0; d < Dimensions; d++) {
		end[d] = Value_T(d) + Value_T(1);
		start[d] = -end[d];
	}

	size_t number_of_samples = 0;
	for (const auto& coord: volume_range(start, end, samples_per_dim)) {
		number_of_samples++;
		for (size_t d = 0; d < Dimensions; d++) {
			if (coord[d] < start[d] or coord[d] > end[d]) {
				std::cerr <<  __FILE__ << "(" << __LINE__ << "): "
					<< d << " " << coord[d] << " " << start[d] << " " << end[d]
					<< std::endl;
				abort();
			}
		}
	}

	if (number_of_samples > std::pow(samples_per_dim, Dimensions)) {
		std::cerr <<  __FILE__ << "(" << __LINE__ << "): "
			<< Dimensions << " " << number_of_samples << " " << samples_per_dim
			<< std::endl;
		abort();
	}
}

int main()
{
	check_number_of_samples<double, 1>(1);
	check_number_of_samples<double, 2>(1);
	check_number_of_samples<double, 3>(1);
	check_number_of_samples<double, 10>(1);

	check_number_of_samples<double, 1>(2);
	check_number_of_samples<double, 2>(2);
	check_number_of_samples<double, 3>(2);
	check_number_of_samples<double, 10>(2);

	check_number_of_samples<double, 1>(10);
	check_number_of_samples<double, 2>(10);
	check_number_of_samples<double, 3>(10);
	check_number_of_samples<double, 10>(3);

	return EXIT_SUCCESS;
}
