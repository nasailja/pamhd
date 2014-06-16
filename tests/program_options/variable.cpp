/*
PAMHD test for parsing gensimcell variable using boost::program_options.

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
#include "iostream"
#include "string"

struct Mass_Density {
	using data_type = double;
	static const std::string get_name() { return {"mass density"}; }
	static const std::string get_option_name() { return {"mass-density"}; }
	static const std::string get_option_help() { return {"Mass density in kg / m^3"}; }
};


template<class Variable> class Prog_Opts_Test
{
public:
	Prog_Opts_Test(const typename Variable::data_type& given_value) :
		value(given_value)
	{}

	void add_options(boost::program_options::options_description& options)
	{
		options.add_options()
			(("test." + Variable::get_option_name()).c_str(), 
				boost::program_options::value<
					typename Variable::data_type
				>(&this->value),
				Variable::get_option_help().c_str());
	}

	typename Variable::data_type& operator[](const Variable&)
	{
		return this->value;
	}

	const typename Variable::data_type& operator[](const Variable&) const
	{
		return this->value;
	}


private:
	// stores (default) value for Variable
	typename Variable::data_type value{};
};


int main(int argc, char* argv[])
{
	Prog_Opts_Test<Mass_Density> prog_opts_test(3);

	boost::program_options::options_description options(
		"Usage: program_name [options], where options are:"
	);
	options.add_options()("help", "Print options and their descriptions");
	prog_opts_test.add_options(options);

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

	if (prog_opts_test[Mass_Density()] != 3) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
			<< "Incorrect value for mass density: "
			<< prog_opts_test[Mass_Density()]
			<< std::endl;
	}

	return EXIT_SUCCESS;
}
