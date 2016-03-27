/*
Tests simulation variable to boost::program_options translator.

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


#include "cstdlib"
#include "iostream"
#include "string"

#include "boost/program_options.hpp"
#ifdef HAVE_EIGEN
#include "Eigen/Core" // must be included before variable_to_option.hpp
#endif

#include "boundaries/variable_to_option.hpp"

using namespace pamhd::boundaries;

// Game of life variable
struct Is_Alive {
	using data_type = bool;
	static const std::string get_name()
	{
		return {"is alive"};
	}
	static const std::string get_option_name()
	{
		return {"is-alive"};
	}
	static const std::string get_option_help()
	{
		return {"Whether cell is alive (1, true, etc.) or not (false)"};
	}
};

// MHD variables
struct Mass_Density {
	using data_type = double;
	static const std::string get_name()
	{
		return {"mass density"};
	}
	static const std::string get_option_name()
	{
		return {"mass-density"};
	}
	static const std::string get_option_help()
	{
		return {"Mass density (kg / m^3)"};
	}
};

struct Momentum_Density {
	#ifdef HAVE_EIGEN
	using data_type = Eigen::Vector3d;
	#else
	using data_type = std::array<double, 3>;
	#endif
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
		return {"Momentum density (J / m^3)"};
	}
};


int main(int argc, char* argv[])
{
	Variable_To_Option<Is_Alive> gol_var_opts;
	Variable_To_Option<Mass_Density, Momentum_Density> mhd_var_opts;

	// set default values before adding options to boost
	gol_var_opts.set_expression(Is_Alive(), "true");
	mhd_var_opts.set_expression(Mass_Density(), "r[0] + r[1] + r[2] + t");
	mhd_var_opts.set_expression(Momentum_Density(), "{r[0], r[1], r[2]}");

	boost::program_options::options_description options(
		"Usage: program_name [options], where options are:"
	);
	options.add_options()("help", "Print options and their descriptions");
	gol_var_opts.add_options("gol.", options);
	mhd_var_opts.add_options("mhd.", options);

	boost::program_options::variables_map option_variables;
	boost::program_options::store(
		boost::program_options::parse_command_line(argc, argv, options),
		option_variables
	);
	boost::program_options::notify(option_variables);

	if (option_variables.count("help") > 0) {
		std::cout << options << std::endl;
		return EXIT_SUCCESS;
	}

	if (gol_var_opts.get_data(Is_Alive(), {{1, 2, 3}}, 4) != true) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< ")" << std::endl;
		abort();
	}

	if (mhd_var_opts.get_data(Mass_Density(), {{1, 2, 3}}, 4) != 1+2+3+4) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< ")" << std::endl;
		abort();
	}

	if (mhd_var_opts.get_data(Momentum_Density(), {{1, 2, 3}}, 4)[0] != 1) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< ")" << std::endl;
		abort();
	}
	if (mhd_var_opts.get_data(Momentum_Density(), {{1, 2, 3}}, 4)[1] != 2) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< ")" << std::endl;
		abort();
	}
	if (mhd_var_opts.get_data(Momentum_Density(), {{1, 2, 3}}, 4)[2] != 3) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< ")" << std::endl;
		abort();
	}

	return EXIT_SUCCESS;
}
