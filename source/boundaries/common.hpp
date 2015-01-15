/*
Functionality used by different parts of PAMHD boundaries.

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

#ifndef PAMHD_BOUNDARIES_COMMON_HPP
#define PAMHD_BOUNDARIES_COMMON_HPP


#include "array"
#include "type_traits"

namespace pamhd {
namespace boundaries {
namespace detail {


template<class T> struct is_std_array_double : std::false_type {};
template<size_t N> struct is_std_array_double<
	std::array<double, N>
> : std::true_type {};

template<class T> struct is_std_array_int : std::false_type {};
template<size_t N> struct is_std_array_int<
	std::array<int, N>
> : std::true_type {};

template<class T> struct is_std_array_bool : std::false_type {};
template<size_t N> struct is_std_array_bool<
	std::array<bool, N>
> : std::true_type {};


/*
https://gcc.gnu.org/bugzilla/show_bug.cgi?id=61836
#ifdef EIGEN_WORLD_VERSION

template<class T> struct is_eigen_vector_double : std::false_type {};
template<size_t N> struct is_eigen_vector_double<
	Eigen::Matrix<double, N, 1>
> : std::true_type {};

#endif // ifdef EIGEN_WORLD_VERSION
*/


}}} // namespaces

#endif // ifndef PAMHD_BOUNDARIES_COMMON_HPP
