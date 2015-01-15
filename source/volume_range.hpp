/*
Functions for iterating over a volume.

Copyright 2014, 2015 Ilja Honkonen
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

#ifndef PAMHD_VOLUME_RANGE_HPP
#define PAMHD_VOLUME_RANGE_HPP


#include "boost/coroutine/all.hpp"
#include "functional"
#include "type_traits"


namespace pamhd {

namespace detail {

//! iterates over values in dimension Index
template <
	template <class Value_T, size_t Length> class Container_T,
	class Value_T,
	size_t Length,
	size_t Index
> class volume_range_impl
{
public:
	static void iterate(
		typename boost::coroutines::coroutine<
			Container_T<Value_T, Length>
		>::push_type& sink,
		Container_T<Value_T, Length>& current,
		const Container_T<Value_T, Length>& start,
		const Container_T<Value_T, Length>& end,
		const Container_T<Value_T, Length>& increment
	) {
		for(
			current[Index] = start[Index] + increment[Index] / Value_T(2);
			current[Index] < end[Index];
			current[Index] += increment[Index]
		) {
			volume_range_impl<
				Container_T,
				Value_T,
				Length,
				Index - 1
			>::iterate(
				sink,
				current,
				start,
				end,
				increment
			);
		}
	}
};

//! iterates over values of the last or only dimension
template <
	template <class Value_T, size_t Length> class Container_T,
	class Value_T,
	size_t Length
> class volume_range_impl<Container_T, Value_T, Length, 0>
{
public:
	static void iterate(
		typename boost::coroutines::coroutine<
			Container_T<Value_T, Length>
		>::push_type& sink,
		Container_T<Value_T, Length>& current,
		const Container_T<Value_T, Length>& start,
		const Container_T<Value_T, Length>& end,
		const Container_T<Value_T, Length>& increment
	) {
		for(
			current[0] = start[0] + increment[0] / Value_T(2);
			current[0] < end[0];
			current[0] += increment[0]
		) {
			sink(current);
		}
	}
};

} // namespace detail

/*!
Produces coordinates inside given volume.

Supports either std::array or boost::array.
Requires only O(D) memory where D in the number of
dimensions in start and end. Given number of samples
is produced in each dimension so the total number is
samples^D. For each dimension separately the samples
have equal amounts of space on both sides so the first
and last samples are at a distance of dx from start
and end respectively and the distance between each
sample is 2*dx.

Example usage:
\verbatim
for (const auto& coord:
	pamhd::volume_range<std::array, double, 2>(
		{-1.0, -2.0},
		{ 1.0,  2.0},
		2
	)
) {
	for (const auto value: coord) {
		std::cout << value << " ";
	}
	std::cout << std::endl;
}
\endverbatim
or
\verbatim
size_t samples = 2;
array<double, 2>
	start{-1.0, -2.0},
	end{1.0, 2.0};

for (const auto& coord: volume_range(start, end, samples)) {
	for (const auto value: coord) {
		cout << value << " ";
	}
	cout << endl;
}
\endverbatim
prints
\verbatim
-0.5 -1 
 0.5 -1 
-0.5  1 
 0.5  1 
\endverbatim

Linking requires boost_coroutine and boost_system
*/
template <
	template <
		class Value_T,
		size_t Length
	> class Container_T,
	class Value_T,
	size_t Length
> typename boost::coroutines::coroutine<
	Container_T<Value_T, Length>
>::pull_type volume_range(
	const Container_T<Value_T, Length>& start,
	const Container_T<Value_T, Length>& end,
	const size_t samples_per_dim
) {
	static_assert(
		std::is_floating_point<Value_T>::value,
		"Only floating point value types supported, "
		"but patches more than welcome"
	);

	Container_T<Value_T, Length> increment;

	// TODO: could loop at comple time
	for (size_t i = 0; i < Length; i++) {
		increment[i] = (end[i] - start[i]) / samples_per_dim;
	}

	return
		typename boost::coroutines::coroutine<
			Container_T<Value_T, Length>
		>::pull_type(
			std::bind(
				&detail::volume_range_impl<
					Container_T,
					Value_T,
					Length,
					Length - 1
				>::iterate,
				std::placeholders::_1,
				Container_T<Value_T, Length>(),
				start,
				end,
				increment
			)
		);
}

} // namespace pamhd

#endif
