/*
PAMHD experiment for parsing arbitrary number of variables with muparserx.

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


#include "array"
#include "cstdlib"
#include "iostream"
#include "string"

#include "boost/program_options.hpp"
#include "mpParser.h"


struct Mass_Density {
	using data_type = double;
	static const std::string get_name() { return {"mass density"}; }
	static const std::string get_option_name() { return {"mass-density"}; }
	static const std::string get_option_help() { return {"Expression for mass density in kg / m^3"}; }
	static const std::string get_expression_variable() { return {"rho"}; }
};

struct Momentum_Density {
	using data_type = std::array<double, 3>;
	static const std::string get_name() { return {"momentum density"}; }
	static const std::string get_option_name() { return {"momentum-density"}; }
	static const std::string get_option_help() { return {"Expression for momentum density in kg / m^2 / s"}; }
};


template<class T> struct is_std_array_double : std::false_type {};
template<size_t N> struct is_std_array_double<
	std::array<double, N>
> : std::true_type
{};



template<class... Variables> class Prog_Opts_Test {};

template<
	class Current_Variable,
	class... Rest_Of_Variables
> class Prog_Opts_Test<Current_Variable, Rest_Of_Variables...> :
	public Prog_Opts_Test<Rest_Of_Variables...>
{
public:
	using Prog_Opts_Test<Rest_Of_Variables...>::add_options;
	using Prog_Opts_Test<Rest_Of_Variables...>::set_expression;
	using Prog_Opts_Test<Rest_Of_Variables...>::get_data;


	Prog_Opts_Test<Current_Variable, Rest_Of_Variables...>() :
		Prog_Opts_Test<Rest_Of_Variables...>(),
		variable(&value)
	{}


	void add_options(boost::program_options::options_description& options)
	{
		options.add_options()
			(Current_Variable::get_option_name().c_str(), 
				boost::program_options::value<std::string>(&this->expression),
				Current_Variable::get_option_help().c_str());

		Prog_Opts_Test<Rest_Of_Variables...>::add_options(options);
	}

	void set_expression(
		const Current_Variable&,
		const std::string given_expression
	) {
		this->expression = given_expression;
	}


	typename Current_Variable::data_type get_data(
		const Current_Variable& variable,
		const std::array<double, 3>& position
	) {
		this->parser.SetExpr(this->expression);

		return this->get_parsed_value<
			typename Current_Variable::data_type
		>(variable, position);
	}


private:

	//! Used if Current_Variable::data_type is double
	template<
		class Out_T
	> typename std::enable_if<
		std::is_same<Out_T, double>::value,
		Out_T
	>::type get_parsed_value(
		const Current_Variable&,
		const std::array<double, 3>& position
	) {
		this->value = mup::Value(3, 0);
		for (size_t i = 0; i < 3; i++) {
			this->value.At(int(i)) = position[i];
		}
		try {
			return this->parser.Eval().GetFloat();
		} catch(mup::ParserError &error) {
			std::cerr <<  __FILE__ << "(" << __LINE__<< ") "
				<< "Could not evaluate expression "
				<< this->parser.GetExpr() << ": " << error.GetMsg()
				<< std::endl;
			abort();
		}
	}


	//! Used if Current_Variable::data_type is std::array<double, N>
	template<
		class Out_T
	> typename std::enable_if<
		is_std_array_double<Out_T>::value,
		Out_T
	>::type get_parsed_value(
		const Current_Variable&,
		const std::array<double, 3>& position
	) {
		constexpr size_t N = std::tuple_size<Out_T>::value;

		this->value = mup::Value(N, 0);
		for (size_t i = 0; i < N; i++) {
			this->value.At(int(i)) = position[i];
		}

		const auto& evaluated
			= [this](){
				try {
					return parser.Eval().GetArray();
				} catch(mup::ParserError &error) {
					std::cerr <<  __FILE__ << "(" << __LINE__<< ") "
						<< "Could not evaluate expression "
						<< this->parser.GetExpr() << ": " << error.GetMsg()
						<< " of variable " << Current_Variable::get_name()
						<< std::endl;
					abort();
				}
			}();

		if (evaluated.GetRows() != 1) {
			std::cerr <<  __FILE__ << "(" << __LINE__<< ") "
				<< "Invalid number of rows in expression " << parser.GetExpr()
				<< " for variable " << Current_Variable::get_name()
				<< ": " << evaluated.GetRows() << ", should be 1"
				<< std::endl;
			abort();
		}
		if (evaluated.GetCols() != N) {
			std::cerr <<  __FILE__ << "(" << __LINE__<< ") "
				<< "Invalid number of columns in expression " << parser.GetExpr()
				<< " for variable " << Current_Variable::get_name()
				<< ": " << evaluated.GetCols() << ", should be " << N
				<< std::endl;
			abort();
		}

		std::array<double, N> ret_val;
		for (size_t i = 0; i < N; i++) {
			ret_val[i] = evaluated.At(int(i)).GetFloat();
		}

		return ret_val;
	}


	std::string expression;
	mup::ParserX parser;
	mup::Value value;
	mup::Variable variable;
};


template<class Last_Variable> class Prog_Opts_Test<Last_Variable>
{
public:

	Prog_Opts_Test<Last_Variable>() : variable(&value) {}


	void set_expression(
		const Last_Variable&,
		const std::string given_expression
	) {
		this->expression = given_expression;
	}


	void add_options(boost::program_options::options_description& options)
	{
		options.add_options()
			(Last_Variable::get_option_name().c_str(), 
				boost::program_options::value<std::string>(&this->expression),
				Last_Variable::get_option_help().c_str());
	}


	typename Last_Variable::data_type get_data(
		const Last_Variable& variable,
		const std::array<double, 3>& position
	) {
		this->parser.SetExpr(this->expression);

		return this->get_parsed_value<
			typename Last_Variable::data_type
		>(variable, position);
	}



private:

	//! Used if Last_Variable::data_type is double
	template<
		class Out_T
	> typename std::enable_if<
		std::is_same<Out_T, double>::value,
		Out_T
	>::type get_parsed_value(
		const Last_Variable&,
		const std::array<double, 3>& position
	) {
		this->value = mup::Value(3, 0);
		for (size_t i = 0; i < 3; i++) {
			this->value.At(int(i)) = position[i];
		}
		try {
			return this->parser.Eval().GetFloat();
		} catch(mup::ParserError &error) {
			std::cerr <<  __FILE__ << "(" << __LINE__<< ") "
				<< "Could not evaluate expression "
				<< this->parser.GetExpr() << ": " << error.GetMsg()
				<< " for variable " << Last_Variable::get_name()
				<< std::endl;
			abort();
		}
	}

	//! Used if Last_Variable::data_type is std::array<double, N>
	template<
		class Out_T
	> typename std::enable_if<
		is_std_array_double<Out_T>::value,
		Out_T
	>::type get_parsed_value(
		const Last_Variable&,
		const std::array<double, 3>& position
	) {
		constexpr size_t N = std::tuple_size<Out_T>::value;

		this->value = mup::Value(N, 0);
		for (size_t i = 0; i < N; i++) {
			this->value.At(int(i)) = position[i];
		}

		const auto& evaluated
			= [this](){
				try {
					return parser.Eval().GetArray();
				} catch(mup::ParserError &error) {
					std::cerr <<  __FILE__ << "(" << __LINE__<< ") "
						<< "Could not evaluate expression "
						<< this->parser.GetExpr() << ": " << error.GetMsg()
						<< " for variable " << Last_Variable::get_name()
						<< std::endl;
					abort();
				}
			}();

		if (evaluated.GetRows() != 1) {
			std::cerr <<  __FILE__ << "(" << __LINE__<< ") "
				<< "Invalid number of rows in expression " << parser.GetExpr()
				<< " of variable " << Last_Variable::get_name()
				<< ": " << evaluated.GetRows() << ", should be 1"
				<< std::endl;
			abort();
		}
		if (evaluated.GetCols() != N) {
			std::cerr <<  __FILE__ << "(" << __LINE__<< ") "
				<< "Invalid number of columns in expression " << parser.GetExpr()
				<< " of variable " << Last_Variable::get_name()
				<< ": " << evaluated.GetCols() << ", should be " << N
				<< std::endl;
			abort();
		}

		std::array<double, N> ret_val;
		for (size_t i = 0; i < N; i++) {
			ret_val[i] = evaluated.At(int(i)).GetFloat();
		}

		return ret_val;
	}


	std::string expression;
	mup::ParserX parser;
	mup::Value value;
	mup::Variable variable;
};


int main(int argc, char* argv[])
{
	Prog_Opts_Test<Mass_Density, Momentum_Density> prog_opts_test;
	prog_opts_test.set_expression(Mass_Density(), "1+1+1");
	prog_opts_test.set_expression(Momentum_Density(), "{-1, -2, -3}");

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

	if (prog_opts_test.get_data(Mass_Density(), {{0, 0, 0}}) != 3) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
			<< "Incorrect value for mass density: "
			<< prog_opts_test.get_data(Mass_Density(), {{0, 0, 0}})
			<< ", should be 3."
			<< std::endl;
		abort();
	}

	if (prog_opts_test.get_data(Momentum_Density(), {{0, 0, 0}})[2] != -3) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
			<< "Incorrect value for momentum density: "
			<< prog_opts_test.get_data(Momentum_Density(), {{0, 0, 0}})[2]
			<< ", should be -3."
			<< std::endl;
		abort();
	}

	return EXIT_SUCCESS;
}
