/*


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


#ifndef PAMHD_BOUNDARIES_VALIDATORS_HPP
#define PAMHD_BOUNDARIES_VALIDATORS_HPP


#include "array"
#include "iostream"

#include "boost/algorithm/string/classification.hpp"
#include "boost/algorithm/string/constants.hpp"
#include "boost/algorithm/string/split.hpp"
#include "boost/algorithm/string/trim.hpp"
#include "boost/lexical_cast.hpp"


namespace std {

template<class T, size_t N> void validate(
	boost::any& value,
	const std::vector<string>& all_parsed,
	array<T, N>*,
	long
) {
	boost::program_options::validators::check_first_occurrence(value);

	const string parsed
		= boost::algorithm::trim_copy_if(
			boost::algorithm::trim_copy(
				boost::program_options::validators::get_single_string(all_parsed)
			),
			boost::algorithm::is_any_of("[]{}")
		);

	vector<string> components;
	components = boost::algorithm::split(
		components,
		parsed,
		boost::algorithm::is_any_of(" \t,"),
		boost::algorithm::token_compress_on
	);

	if (components.size() != N) {
		std::cerr << __FILE__ << "(" << __LINE__<< "): "
			<< "Option must have " << N << " values "
			<< "separated by e.g. a , or a space but given option ("
			<< parsed << ") has " << components.size()
			<< std::endl;
		abort();
	}

	array<T, N> final;
	for (size_t i = 0; i < N; i++) {
		final[i] = boost::lexical_cast<T>(components[i]);
	}

	value = final;
}

} // namespace std


#ifdef EIGEN_WORLD_VERSION
namespace Eigen {

template <
	// TODO: variadic (https://gcc.gnu.org/bugzilla/show_bug.cgi?id=61539)
	class Scalar,
	int Rows,
	int Cols,
	int Options,
	int MaxRows,
	int MaxCols
> void validate(
	boost::any& value,
	const std::vector<std::string>& all_parsed,
	Matrix<
		Scalar,
		Rows,
		Cols,
		Options,
		MaxRows,
		MaxCols
	>*,
	long
) {
	boost::program_options::validators::check_first_occurrence(value);

	const std::string parsed
		= boost::algorithm::trim_copy_if(
			boost::algorithm::trim_copy(
				boost::program_options::validators::get_single_string(all_parsed)
			),
			boost::algorithm::is_any_of("[]{}")
		);

	std::vector<std::string> components;
	components = boost::algorithm::split(
		components,
		parsed,
		boost::algorithm::is_any_of(" \t,"),
		boost::algorithm::token_compress_on
	);

	using Data_T = Matrix<
		Scalar,
		Rows,
		Cols,
		Options,
		MaxRows,
		MaxCols
	>;

	Data_T final;

	size_t i = 0;
	for (typename Data_T::Index row = 0; row < Rows; row++) {
	for (typename Data_T::Index col = 0; col < Cols; col++) {
		final(row, col) = boost::lexical_cast<Scalar>(components[i]);
		i++;
	}}

	value = final;
}

} // namespace Eigen

#endif // ifdef EIGEN_WORLD_VERSION

#endif // ifndef PAMHD_BOUNDARIES_VALIDATORS_HPP
