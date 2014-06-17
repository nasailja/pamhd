/*
PAMHD test for parsing Eigen types using boost::program_options.

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


#include "boost/algorithm/string/classification.hpp"
#include "boost/algorithm/string/constants.hpp"
#include "boost/algorithm/string/split.hpp"
#include "boost/lexical_cast.hpp"
#include "boost/program_options.hpp"
#include "cmath"
#include "cstdlib"
#include "Eigen/Core"
#include "iostream"
#include "string"
#include "fstream"
#include "iomanip"
#include "iostream"
#include "sstream"
#include "string"

using namespace std;

namespace Eigen {
	template <
		class Scalar,
		int Rows,
		int Cols,
		int Options,
		int MaxRows,
		int MaxCols
	> void validate(
		boost::any& value,
		const std::vector<std::string>& all_parsed,
		Eigen::Matrix<
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

		const std::string& parsed
			= boost::program_options::validators::get_single_string(all_parsed);

		std::vector<std::string> components;
		components = boost::algorithm::split(
			components,
			parsed,
			boost::algorithm::is_any_of(" \t,;"),
			boost::algorithm::token_compress_on
		);

		if (components.size() != 3) {
			std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
				<< "Eigen::Vector3d option must have 3 values "
				<< "separated by e.g. a , or a space but given option ("
				<< parsed << ") has " << components.size()
				<< std::endl;
			abort();
		}

		const Eigen::Vector3d final(
			boost::lexical_cast<double>(components[0]),
			boost::lexical_cast<double>(components[1]),
			boost::lexical_cast<double>(components[2])
		);

		value = final;
	}
}

int main(int argc, char* argv[])
{
	Eigen::Vector3d test(-1, -2, -3);
	boost::program_options::options_description options(
		"Usage: program_name [options], where options are:"
	);
	options.add_options()(
		"test",
		boost::program_options::value<Eigen::Vector3d>(&test)
			->default_value(Eigen::Vector3d(1, 2, 3)),
		"test"
	);

	boost::program_options::variables_map option_variables;
	boost::program_options::store(
		boost::program_options::parse_command_line(argc, argv, options),
		option_variables
	);
	boost::program_options::notify(option_variables);

	if (test[0] == -1) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
			<< "First component of test vector was not initialized."
			<< std::endl;
		abort();
	}
	if (test[1] == -2) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
			<< "First component of test vector was not initialized."
			<< std::endl;
		abort();
	}
	if (test[2] == -3) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
			<< "First component of test vector was not initialized."
			<< std::endl;
		abort();
	}

	return EXIT_SUCCESS;
}
