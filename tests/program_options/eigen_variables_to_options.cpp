/*
Tests simulation variables (Eigen Vectors) to boost::program_options translator.

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


#include "boost/program_options.hpp"
#include "cstdlib"
#include "Eigen/Core" // must be included before variables_to_options.hpp
#include "iostream"
#include "string"

#include "boundaries/variables_to_options.hpp"

using namespace pamhd::boundaries;

// MHD variables
struct Momentum_Density {
	using data_type = Eigen::Vector3d;
	static const std::string get_name()
	{
		return {"momentum density"};
	}
	static const std::string get_option_name()
	{
		return {"momentum-density"};
	}
	static const std::string get_option_help()
	{
		return {"Momentum density (kg / s / m^2)"};
	}
};

struct Magnetic_Field {
	using data_type = Eigen::Vector3d;
	static const std::string get_name()
	{
		return {"magnetic field"};
	}
	static const std::string get_option_name()
	{
		return {"magnetic-field"};
	}
	static const std::string get_option_help()
	{
		return {"Magnetic field (T)"};
	}
};


int main(int argc, char* argv[])
{
	Variables_To_Options<Momentum_Density, Magnetic_Field> mhd_var_opts;
	// set default values before adding options to boost
	mhd_var_opts[Momentum_Density()].push_back({3, 4, 5});
	mhd_var_opts[Momentum_Density()].push_back({6, 7, 8});
	mhd_var_opts[Magnetic_Field()].push_back({-3, -4, -5});
	mhd_var_opts[Magnetic_Field()].push_back({-6, -7, -8});

	boost::program_options::options_description options(
		"Usage: program_name [options], where options are:"
	);
	options.add_options()("help", "Print options and their descriptions");
	mhd_var_opts.add_options(options, "mhd.");

	boost::program_options::variables_map option_variables;
	boost::program_options::store(
		boost::program_options::parse_command_line(argc, argv, options),
		option_variables
	);
	boost::program_options::notify(option_variables);

	if (option_variables.count("help") > 0) {
		/* TODO: remove newlines from eigen types
		Eigen::IOFormat format(
			Eigen::StreamPrecision,
			0,
			", ",
			", ",
			"[",
			"]"
		);*/
		std::cout << options << std::endl;
		return EXIT_SUCCESS;
	}

	if (mhd_var_opts[Momentum_Density()].size() != 2) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< ")" << std::endl;
		abort();
	}
	if (mhd_var_opts[Momentum_Density()][0][0] != 3) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< ")" << std::endl;
		abort();
	}
	if (mhd_var_opts[Momentum_Density()][1][0] != 6) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< ")" << std::endl;
		abort();
	}

	if (mhd_var_opts[Magnetic_Field()].size() != 2) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< ")" << std::endl;
		abort();
	}
	if (mhd_var_opts[Magnetic_Field()][0][0] != -3) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< ")" << std::endl;
		abort();
	}
	if (mhd_var_opts[Magnetic_Field()][1][0] != -6) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< ")" << std::endl;
		abort();
	}

	return EXIT_SUCCESS;
}
