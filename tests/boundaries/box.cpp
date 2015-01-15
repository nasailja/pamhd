/*
Box geometry test of PAMHD.

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

#include "array"
#include "cstdlib"
#include "iostream"
#include "string"

#include "boost/program_options.hpp"
#ifdef HAVE_EIGEN
#include "Eigen/Core" // must be included before box.hpp
#endif

#include "boundaries/box.hpp"


#ifdef HAVE_EIGEN
#define MAKE_BOX(a, b, c, d, e, f) {a, b, c}, {d, e, f}
#else
#define MAKE_BOX(a, b, c, d, e, f) {{a, b, c}}, {{d, e, f}}
#endif


using namespace pamhd::boundaries;

int main(int argc, char* argv[])
{
	#ifdef HAVE_EIGEN
	Box<Eigen::Vector3d> box;
	#else
	Box<std::array<double, 3>> box;
	#endif

	if (not box.set_geometry(MAKE_BOX(-1, -2, -3, 1, 2, 3))) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
			<< "Couldn't set geometry."
			<< std::endl;
		abort();
	}

	boost::program_options::options_description options(
		"Usage: program_name [options], where options are:"
	);
	options.add_options()
		("help", "Print options and their descriptions");
	box.add_options("", options);

	boost::program_options::variables_map option_variables;
	boost::program_options::store(
		boost::program_options::parse_command_line<char>(
			argc,
			argv,
			options
		),
		option_variables
	);
	boost::program_options::notify(option_variables);

	if (option_variables.count("help") > 0) {
		std::cout << options << std::endl;
		return EXIT_SUCCESS;
	}

	if (box.overlaps(MAKE_BOX(-2, -3, -4, -1.1, -2.1, -3.1))) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
			<< "Cell overlaps box."
			<< std::endl;
		abort();
	}

	if (box.overlaps(MAKE_BOX(1.1, 2.1, 3.1, 2, 3, 4))) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
			<< "Cell overlaps box."
			<< std::endl;
		abort();
	}

	if (not box.overlaps(MAKE_BOX(0, 1, 2, 5, 5, 5))) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
			<< "Cell does not overlap box."
			<< std::endl;
		abort();
	}


	return EXIT_SUCCESS;
}
