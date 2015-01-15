/*
Class for evaluating an expression for a variable via boost::program_options.

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

#ifndef PAMHD_VARIABLE_TO_OPTION_HPP
#define PAMHD_VARIABLE_TO_OPTION_HPP


#include "array"
#include "cstddef"
#include "type_traits"

#include "boost/program_options.hpp"
#include "mpParser.h"

#include "boundaries/common.hpp"


namespace pamhd {
namespace boundaries {


template<class... Variables> class Variable_To_Option {};



template<
	class Current_Variable,
	class... Rest_Of_Variables
> class Variable_To_Option<Current_Variable, Rest_Of_Variables...> :
	public Variable_To_Option<Rest_Of_Variables...>
{
public:
	using Variable_To_Option<Rest_Of_Variables...>::set_expression;
	using Variable_To_Option<Rest_Of_Variables...>::get_data;


	Variable_To_Option<Current_Variable, Rest_Of_Variables...>() :
		Variable_To_Option<Rest_Of_Variables...>(),
		r_var(&r_val),
		t_var(&t_val)
	{
		this->parser.DefineVar("r", this->r_var);
		this->parser.DefineVar("t", this->t_var);
	}


	void add_options(
		const std::string& prefix,
		boost::program_options::options_description& options
	) {
		this->add_options_impl(prefix, options);
	}


	void set_expression(
		const Current_Variable&,
		const std::string given_expression
	) {
		this->expression = given_expression;
	}


	typename Current_Variable::data_type get_data(
		const Current_Variable&,
		const std::array<double, 3>& given_position,
		const double given_time
	) {
		this->parser.SetExpr(this->expression);

		return this->get_parsed_value<
			typename Current_Variable::data_type
		>(given_position, given_time);
	}

	typename Current_Variable::data_type get_data(
		const Current_Variable& variable
	) {
		return this->get_data(variable, {{0, 0, 0}}, 0);
	}



protected:

	using Variable_To_Option<Rest_Of_Variables...>::add_options_impl;

	void add_options_impl(
		const std::string& prefix,
		boost::program_options::options_description& options
	) {
		options.add_options()
			((prefix + Current_Variable::get_option_name()).c_str(),
				boost::program_options::value<std::string>(&this->expression)
					->default_value(this->expression),
				Current_Variable::get_option_help().c_str());

		Variable_To_Option<Rest_Of_Variables...>::add_options_impl(prefix, options);
	}



private:

	std::string expression;
	mup::ParserX parser;
	mup::Value r_val, t_val;
	mup::Variable r_var, t_var;


	void set_expression_variables(
		const std::array<double, 3>& given_position,
		const double given_time
	) {
		this->r_val = mup::Value(3, 0);
		for (size_t i = 0; i < 3; i++) {
			this->r_val.At(int(i)) = given_position[i];
		}

		this->t_val = given_time;
	}


	//! Used if Current_Variable::data_type is bool
	template<
		class Out_T
	> typename std::enable_if<
		std::is_same<Out_T, bool>::value,
		Out_T
	>::type get_parsed_value(
		const std::array<double, 3>& given_position,
		const double given_time
	) {
		this->set_expression_variables(given_position, given_time);

		try {

			const auto& evaluated = this->parser.Eval();
			if (
				evaluated.GetType() != 'b'
				and evaluated.GetType() != 'i'
				and evaluated.GetType() != 'f'
			) {
				std::cerr <<  __FILE__ << "(" << __LINE__<< ") "
					<< "Could not evaluate expression "
					<< this->parser.GetExpr() << " as a boolean."
					<< std::endl;
				abort();
			}

			return evaluated.GetBool();

		} catch(mup::ParserError &error) {
			std::cerr <<  __FILE__ << "(" << __LINE__<< ") "
				<< "Could not evaluate expression "
				<< this->parser.GetExpr() << " as a boolean."
				<< std::endl;
			abort();
		}
	}


	//! Used if Current_Variable::data_type is int
	template<
		class Out_T
	> typename std::enable_if<
		std::is_same<Out_T, int>::value,
		Out_T
	>::type get_parsed_value(
		const std::array<double, 3>& given_position,
		const double given_time
	) {
		this->set_expression_variables(given_position, given_time);

		try {

			const auto& evaluated = this->parser.Eval();
			if (
				evaluated.GetType() != 'b'
				and evaluated.GetType() != 'i'
				and evaluated.GetType() != 'f'
			) {
				std::cerr <<  __FILE__ << "(" << __LINE__<< ") "
					<< "Could not evaluate expression "
					<< this->parser.GetExpr() << " as an integer."
					<< std::endl;
				abort();
			}

			return evaluated.GetInteger();

		} catch(mup::ParserError &error) {
			std::cerr <<  __FILE__ << "(" << __LINE__<< ") "
				<< "Could not evaluate expression "
				<< this->parser.GetExpr() << " as an integer."
				<< std::endl;
			abort();
		}
	}


	//! Used if Current_Variable::data_type is double
	template<
		class Out_T
	> typename std::enable_if<
		std::is_same<Out_T, double>::value,
		Out_T
	>::type get_parsed_value(
		const std::array<double, 3>& given_position,
		const double given_time
	) {
		this->set_expression_variables(given_position, given_time);

		try {

			const auto& evaluated = this->parser.Eval();
			if (
				evaluated.GetType() != 'b'
				and evaluated.GetType() != 'i'
				and evaluated.GetType() != 'f'
			) {
				std::cerr <<  __FILE__ << "(" << __LINE__<< ") "
					<< "Could not evaluate expression "
					<< this->parser.GetExpr() << " as a double."
					<< std::endl;
				abort();
			}

			return evaluated.GetFloat();

		} catch(mup::ParserError &error) {
			std::cerr <<  __FILE__ << "(" << __LINE__<< ") "
				<< "Could not evaluate expression "
				<< this->parser.GetExpr() << " as a double."
				<< std::endl;
			abort();
		}
	}


	//! Used if Current_Variable::data_type is std::array<double, N>
	template<
		class Out_T
	> typename std::enable_if<
		detail::is_std_array_double<Out_T>::value,
		Out_T
	>::type get_parsed_value(
		const std::array<double, 3>& given_position,
		const double given_time
	) {
		this->set_expression_variables(given_position, given_time);

		const auto& evaluated
			= [this](){
				try {

					const auto& temp = this->parser.Eval();
					if (temp.GetType() != 'm') {
						std::cerr <<  __FILE__ << "(" << __LINE__<< ") "
							<< "Could not evaluate expression "
							<< this->parser.GetExpr() << " as an array."
							<< std::endl;
						abort();
					}

					return temp.GetArray();

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

		constexpr size_t N = std::tuple_size<Out_T>::value;

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


	//! Used if Current_Variable::data_type is std::array<int, N>
	template<
		class Out_T
	> typename std::enable_if<
		detail::is_std_array_int<Out_T>::value,
		Out_T
	>::type get_parsed_value(
		const std::array<double, 3>& given_position,
		const double given_time
	) {
		this->set_expression_variables(given_position, given_time);

		const auto& evaluated
			= [this](){
				try {

					const auto& temp = this->parser.Eval();
					if (temp.GetType() != 'm') {
						std::cerr <<  __FILE__ << "(" << __LINE__<< ") "
							<< "Could not evaluate expression "
							<< this->parser.GetExpr() << " as an array."
							<< std::endl;
						abort();
					}

					return temp.GetArray();

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

		constexpr size_t N = std::tuple_size<Out_T>::value;

		if (evaluated.GetCols() != N) {
			std::cerr <<  __FILE__ << "(" << __LINE__<< ") "
				<< "Invalid number of columns in expression " << parser.GetExpr()
				<< " for variable " << Current_Variable::get_name()
				<< ": " << evaluated.GetCols() << ", should be " << N
				<< std::endl;
			abort();
		}

		std::array<int, N> ret_val;
		for (size_t i = 0; i < N; i++) {
			ret_val[i] = evaluated.At(int(i)).GetInteger();
		}

		return ret_val;
	}


	//! Used if Current_Variable::data_type is std::array<bool, N>
	template<
		class Out_T
	> typename std::enable_if<
		detail::is_std_array_bool<Out_T>::value,
		Out_T
	>::type get_parsed_value(
		const std::array<double, 3>& given_position,
		const double given_time
	) {
		this->set_expression_variables(given_position, given_time);

		const auto& evaluated
			= [this](){
				try {

					const auto& temp = this->parser.Eval();
					if (temp.GetType() != 'm') {
						std::cerr <<  __FILE__ << "(" << __LINE__<< ") "
							<< "Could not evaluate expression "
							<< this->parser.GetExpr() << " as an array."
							<< std::endl;
						abort();
					}

					return temp.GetArray();

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

		constexpr size_t N = std::tuple_size<Out_T>::value;

		if (evaluated.GetCols() != N) {
			std::cerr <<  __FILE__ << "(" << __LINE__<< ") "
				<< "Invalid number of columns in expression " << parser.GetExpr()
				<< " for variable " << Current_Variable::get_name()
				<< ": " << evaluated.GetCols() << ", should be " << N
				<< std::endl;
			abort();
		}

		std::array<bool, N> ret_val;
		for (size_t i = 0; i < N; i++) {
			ret_val[i] = evaluated.At(int(i)).GetBool();
		}

		return ret_val;
	}


	#ifdef EIGEN_WORLD_VERSION

	//! Used if Current_Variable::data_type is Eigen::Matrix<double, N, 1>
	template<
		class Out_T
	> typename std::enable_if<
		//detail::is_eigen_vector_double<Out_T>::value,
		std::is_same<Out_T, Eigen::Matrix<double, 3, 1>>::value,
		Out_T
	>::type get_parsed_value(
		const std::array<double, 3>& given_position,
		const double given_time
	) {
		this->set_expression_variables(given_position, given_time);

		const auto& evaluated
			= [this](){
				try {

					const auto& temp = this->parser.Eval();
					if (temp.GetType() != 'm') {
						std::cerr <<  __FILE__ << "(" << __LINE__<< ") "
							<< "Could not evaluate expression "
							<< this->parser.GetExpr() << " as an array."
							<< std::endl;
						abort();
					}

					return temp.GetArray();

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

		constexpr size_t N = Out_T::RowsAtCompileTime;

		if (evaluated.GetCols() != N) {
			std::cerr <<  __FILE__ << "(" << __LINE__<< ") "
				<< "Invalid number of columns in expression " << parser.GetExpr()
				<< " for variable " << Current_Variable::get_name()
				<< ": " << evaluated.GetCols() << ", should be " << N
				<< std::endl;
			abort();
		}

		Eigen::Matrix<double, N, 1> ret_val;
		for (size_t i = 0; i < N; i++) {
			ret_val(i) = evaluated.At(int(i)).GetFloat();
		}

		return ret_val;
	}

	#endif // ifdef EIGEN_WORLD_VERSION
};



template<class Last_Variable> class Variable_To_Option<Last_Variable>
{
public:

	Variable_To_Option<Last_Variable>() :
		r_var(&r_val),
		t_var(&t_val)
	{
		this->parser.DefineVar("r", this->r_var);
		this->parser.DefineVar("t", this->t_var);
	}


	void set_expression(
		const Last_Variable&,
		const std::string given_expression
	) {
		this->expression = given_expression;
	}


	void add_options(
		const std::string& prefix,
		boost::program_options::options_description& options
	) {
		this->add_options_impl(prefix, options);
	}


	typename Last_Variable::data_type get_data(
		const Last_Variable&,
		const std::array<double, 3>& given_position,
		const double given_time
	) {
		this->parser.SetExpr(this->expression);

		return this->get_parsed_value<
			typename Last_Variable::data_type
		>(given_position, given_time);
	}

	typename Last_Variable::data_type get_data(
		const Last_Variable& variable
	) {
		return this->get_data(variable, {{0, 0, 0}}, 0);
	}



protected:

	void add_options_impl(
		const std::string& prefix,
		boost::program_options::options_description& options
	) {
		options.add_options()
			((prefix + Last_Variable::get_option_name()).c_str(),
				boost::program_options::value<std::string>(&this->expression)
					->default_value(this->expression),
				Last_Variable::get_option_help().c_str());
	}



private:

	std::string expression;
	mup::ParserX parser;
	mup::Value r_val, t_val;
	mup::Variable r_var, t_var;


	void set_expression_variables(
		const std::array<double, 3>& given_position,
		const double given_time
	) {
		this->r_val = mup::Value(3, 0);
		for (size_t i = 0; i < 3; i++) {
			this->r_val.At(int(i)) = given_position[i];
		}

		this->t_val = given_time;
	}


	//! Used if Current_Variable::data_type is bool
	template<
		class Out_T
	> typename std::enable_if<
		std::is_same<Out_T, bool>::value,
		Out_T
	>::type get_parsed_value(
		const std::array<double, 3>& given_position,
		const double given_time
	) {
		this->set_expression_variables(given_position, given_time);

		try {

			const auto& evaluated = this->parser.Eval();
			if (
				evaluated.GetType() != 'b'
				and evaluated.GetType() != 'i'
				and evaluated.GetType() != 'f'
			) {
				std::cerr <<  __FILE__ << "(" << __LINE__<< ") "
					<< "Could not evaluate expression "
					<< this->parser.GetExpr() << " as a boolean."
					<< std::endl;
				abort();
			}

			return evaluated.GetBool();

		} catch(mup::ParserError &error) {
			std::cerr <<  __FILE__ << "(" << __LINE__<< ") "
				<< "Could not evaluate expression "
				<< this->parser.GetExpr() << ": " << error.GetMsg()
				<< std::endl;
			abort();
		}
	}


	//! Used if Current_Variable::data_type is int
	template<
		class Out_T
	> typename std::enable_if<
		std::is_same<Out_T, int>::value,
		Out_T
	>::type get_parsed_value(
		const std::array<double, 3>& given_position,
		const double given_time
	) {
		this->set_expression_variables(given_position, given_time);

		try {

			const auto& evaluated = this->parser.Eval();
			if (
				evaluated.GetType() != 'b'
				and evaluated.GetType() != 'i'
				and evaluated.GetType() != 'f'
			) {
				std::cerr <<  __FILE__ << "(" << __LINE__<< ") "
					<< "Could not evaluate expression "
					<< this->parser.GetExpr() << " as an integer."
					<< std::endl;
				abort();
			}

			return evaluated.GetInteger();

		} catch(mup::ParserError &error) {
			std::cerr <<  __FILE__ << "(" << __LINE__<< ") "
				<< "Could not evaluate expression "
				<< this->parser.GetExpr() << ": " << error.GetMsg()
				<< std::endl;
			abort();
		}
	}


	//! Used if Last_Variable::data_type is double
	template<
		class Out_T
	> typename std::enable_if<
		std::is_same<Out_T, double>::value,
		Out_T
	>::type get_parsed_value(
		const std::array<double, 3>& given_position,
		const double given_time
	) {
		this->set_expression_variables(given_position, given_time);

		try {

			const auto& evaluated = this->parser.Eval();
			if (
				evaluated.GetType() != 'b'
				and evaluated.GetType() != 'i'
				and evaluated.GetType() != 'f'
			) {
				std::cerr <<  __FILE__ << "(" << __LINE__<< ") "
					<< "Could not evaluate expression "
					<< this->parser.GetExpr() << " as a floating point."
					<< std::endl;
				abort();
			}

			return evaluated.GetFloat();

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
		detail::is_std_array_double<Out_T>::value,
		Out_T
	>::type get_parsed_value(
		const std::array<double, 3>& given_position,
		const double given_time
	) {
		this->set_expression_variables(given_position, given_time);

		const auto& evaluated
			= [this](){
				try {

					const auto& temp = this->parser.Eval();
					if (temp.GetType() != 'm') {
						std::cerr <<  __FILE__ << "(" << __LINE__<< ") "
							<< "Could not evaluate expression "
							<< this->parser.GetExpr() << " as an array."
							<< std::endl;
						abort();
					}

					return temp.GetArray();

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

		constexpr size_t N = std::tuple_size<Out_T>::value;

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


	//! Used if Last_Variable::data_type is std::array<int, N>
	template<
		class Out_T
	> typename std::enable_if<
		detail::is_std_array_int<Out_T>::value,
		Out_T
	>::type get_parsed_value(
		const std::array<double, 3>& given_position,
		const double given_time
	) {
		this->set_expression_variables(given_position, given_time);

		const auto& evaluated
			= [this](){
				try {

					const auto& temp = this->parser.Eval();
					if (temp.GetType() != 'm') {
						std::cerr <<  __FILE__ << "(" << __LINE__<< ") "
							<< "Could not evaluate expression "
							<< this->parser.GetExpr() << " as an array."
							<< std::endl;
						abort();
					}

					return temp.GetArray();

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

		constexpr size_t N = std::tuple_size<Out_T>::value;

		if (evaluated.GetCols() != N) {
			std::cerr <<  __FILE__ << "(" << __LINE__<< ") "
				<< "Invalid number of columns in expression " << parser.GetExpr()
				<< " of variable " << Last_Variable::get_name()
				<< ": " << evaluated.GetCols() << ", should be " << N
				<< std::endl;
			abort();
		}

		std::array<int, N> ret_val;
		for (size_t i = 0; i < N; i++) {
			ret_val[i] = evaluated.At(int(i)).GetInteger();
		}

		return ret_val;
	}


	//! Used if Last_Variable::data_type is std::array<bool, N>
	template<
		class Out_T
	> typename std::enable_if<
		detail::is_std_array_bool<Out_T>::value,
		Out_T
	>::type get_parsed_value(
		const std::array<double, 3>& given_position,
		const double given_time
	) {
		this->set_expression_variables(given_position, given_time);

		const auto& evaluated
			= [this](){
				try {

					const auto& temp = this->parser.Eval();
					if (temp.GetType() != 'm') {
						std::cerr <<  __FILE__ << "(" << __LINE__<< ") "
							<< "Could not evaluate expression "
							<< this->parser.GetExpr() << " as an array."
							<< std::endl;
						abort();
					}

					return temp.GetArray();

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

		constexpr size_t N = std::tuple_size<Out_T>::value;

		if (evaluated.GetCols() != N) {
			std::cerr <<  __FILE__ << "(" << __LINE__<< ") "
				<< "Invalid number of columns in expression " << parser.GetExpr()
				<< " of variable " << Last_Variable::get_name()
				<< ": " << evaluated.GetCols() << ", should be " << N
				<< std::endl;
			abort();
		}

		std::array<bool, N> ret_val;
		for (size_t i = 0; i < N; i++) {
			ret_val[i] = evaluated.At(int(i)).GetBool();
		}

		return ret_val;
	}


	#ifdef EIGEN_WORLD_VERSION

	//! Used if Current_Variable::data_type is Eigen::Matrix<double, N, 1>
	template<
		class Out_T
	> typename std::enable_if<
		//detail::is_eigen_vector_double<Out_T>::value,
		std::is_same<Out_T, Eigen::Matrix<double, 3, 1>>::value,
		Out_T
	>::type get_parsed_value(
		const std::array<double, 3>& given_position,
		const double given_time
	) {
		this->set_expression_variables(given_position, given_time);

		const auto& evaluated
			= [this](){
				try {

					const auto& temp = this->parser.Eval();
					if (temp.GetType() != 'm') {
						std::cerr <<  __FILE__ << "(" << __LINE__<< ") "
							<< "Could not evaluate expression "
							<< this->parser.GetExpr() << " as an array."
							<< std::endl;
						abort();
					}

					return temp.GetArray();

				} catch(mup::ParserError &error) {
					std::cerr <<  __FILE__ << "(" << __LINE__<< ") "
						<< "Could not evaluate expression "
						<< this->parser.GetExpr() << ": " << error.GetMsg()
						<< " of variable " << Last_Variable::get_name()
						<< std::endl;
					abort();
				}
			}();

		if (evaluated.GetRows() != 1) {
			std::cerr <<  __FILE__ << "(" << __LINE__<< ") "
				<< "Invalid number of rows in expression " << parser.GetExpr()
				<< " for variable " << Last_Variable::get_name()
				<< ": " << evaluated.GetRows() << ", should be 1"
				<< std::endl;
			abort();
		}

		constexpr size_t N = Out_T::RowsAtCompileTime;

		if (evaluated.GetCols() != N) {
			std::cerr <<  __FILE__ << "(" << __LINE__<< ") "
				<< "Invalid number of columns in expression " << parser.GetExpr()
				<< " for variable " << Last_Variable::get_name()
				<< ": " << evaluated.GetCols() << ", should be " << N
				<< std::endl;
			abort();
		}

		Eigen::Matrix<double, N, 1> ret_val;
		for (size_t i = 0; i < N; i++) {
			ret_val(i) = evaluated.At(int(i)).GetFloat();
		}

		return ret_val;
	}

	#endif // ifdef EIGEN_WORLD_VERSION
};




}} // namespaces

#endif // ifndef PAMHD_VARIABLE_TO_OPTION_HPP
