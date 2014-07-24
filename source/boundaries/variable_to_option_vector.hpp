/*
Time-dependent version of variable_to_option.hpp.

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

#ifndef PAMHD_BOUNDARIES_VARIABLE_TO_OPTION_VECTOR_HPP
#define PAMHD_BOUNDARIES_VARIABLE_TO_OPTION_VECTOR_HPP


#include "array"
#include "cstddef"
#include "type_traits"
#include "utility"

#include "boost/program_options.hpp"
#include "mpParser.h"
#include "prettyprint.hpp"

#include "boundaries/common.hpp"


namespace pamhd {
namespace boundaries {


template<class... Variables> class Variable_To_Option_Vector {};



template<
	class Current_Variable,
	class... Rest_Of_Variables
> class Variable_To_Option_Vector<Current_Variable, Rest_Of_Variables...> :
	public Variable_To_Option_Vector<Rest_Of_Variables...>
{
public:
	using Variable_To_Option_Vector<Rest_Of_Variables...>::set_expression;
	using Variable_To_Option_Vector<Rest_Of_Variables...>::get_data;


	Variable_To_Option_Vector<Current_Variable, Rest_Of_Variables...>() :
		Variable_To_Option_Vector<Rest_Of_Variables...>(),
		r_var(&r_val),
		t_var(&t_val)
	{
		this->parser.DefineVar("r", this->r_var);
		this->parser.DefineVar("t", this->t_var);
	}

	Variable_To_Option_Vector<Current_Variable, Rest_Of_Variables...>(
		const Variable_To_Option_Vector<Current_Variable, Rest_Of_Variables...>& other
	) : Variable_To_Option_Vector<Current_Variable, Rest_Of_Variables...>()
	{
		this->expressions = other.expressions;
	}

	Variable_To_Option_Vector<Current_Variable, Rest_Of_Variables...>(
		Variable_To_Option_Vector<Current_Variable, Rest_Of_Variables...>&& other
	) : Variable_To_Option_Vector<Current_Variable, Rest_Of_Variables...>()
	{
		this->expressions = std::move(other.expressions);
	}

	Variable_To_Option_Vector<Current_Variable, Rest_Of_Variables...>& operator=(
		const Variable_To_Option_Vector<Current_Variable, Rest_Of_Variables...>& other
	) {
		if (this != &other) {
			this->expressions = other.expressions;
		}

		return *this;
	}


	void add_options(
		const std::string& prefix,
		boost::program_options::options_description& options
	) {
		this->add_options_impl(prefix, options);
	}


	size_t get_number_of_expressions() const
	{
		return this->get_number_of_expressions_impl();
	}


	void set_number_of_expressions(const size_t new_size)
	{
		this->set_number_of_expressions_impl(new_size);
	}


	void clear_expressions()
	{
		this->clear_expressions_impl();
	}


	void set_expression(
		const Current_Variable&,
		const size_t index,
		const std::string given_expression
	) {
		if (index >= this->expressions.size()) {
			std::cerr <<  __FILE__ << "(" << __LINE__<< ") "
				<< "Invalid index: " << index
				<< ", should be less than " << this->expressions.size()
				<< std::endl;
			abort();
		}

		this->expressions[index] = given_expression;
	}


	typename Current_Variable::data_type get_data(
		const Current_Variable&,
		const size_t index,
		const std::array<double, 3>& given_position,
		const double given_time
	) {
		if (index >= this->expressions.size()) {
			std::cerr <<  __FILE__ << "(" << __LINE__<< ") "
				<< "Invalid index: " << index
				<< ", should be less than " << this->expressions.size()
				<< std::endl;
			abort();
		}

		this->parser.SetExpr(this->expressions[index]);

		return this->get_parsed_value<
			typename Current_Variable::data_type
		>(given_position, given_time);
	}

	typename Current_Variable::data_type get_data(
		const Current_Variable& variable,
		const size_t index
	) {
		return this->get_data(variable, index, {0, 0, 0}, 0);
	}



protected:

	using Variable_To_Option_Vector<Rest_Of_Variables...>::add_options_impl;
	using Variable_To_Option_Vector<Rest_Of_Variables...>::get_number_of_expressions_impl;
	using Variable_To_Option_Vector<Rest_Of_Variables...>::set_number_of_expressions_impl;
	using Variable_To_Option_Vector<Rest_Of_Variables...>::clear_expressions_impl;

	void add_options_impl(
		const std::string& prefix,
		boost::program_options::options_description& options
	) {
		options.add_options()
			((prefix + Current_Variable::get_option_name()).c_str(),
				boost::program_options::value<
					std::vector<std::string>
				>(&this->expressions)
					->composing()
					->default_value(this->expressions),
				Current_Variable::get_option_help().c_str());

		Variable_To_Option_Vector<Rest_Of_Variables...>::add_options_impl(prefix, options);
	}


	size_t get_number_of_expressions_impl() const
	{
		const size_t
			current_nr = this->expressions.size(),
			next_nr
				= Variable_To_Option_Vector<
					Rest_Of_Variables...
				>::get_number_of_expressions_impl();

		if (current_nr != next_nr) {
			std::cerr <<  __FILE__ << "(" << __LINE__<< ") "
				<< "Number of expressions differ for variable "
				<< Current_Variable::get_name() << " and the next one: "
				<< current_nr << " != " << next_nr
				<< std::endl;
			abort();
		}

		return current_nr;
	}


	void set_number_of_expressions_impl(const size_t new_size)
	{
		this->expressions.resize(new_size);
		Variable_To_Option_Vector<
			Rest_Of_Variables...
		>::set_number_of_expressions_impl(new_size);
	}


	void clear_expressions_impl()
	{
		this->expressions.clear();
		Variable_To_Option_Vector<Rest_Of_Variables...>::clear_expressions_impl();
	}



private:

	std::vector<std::string> expressions;
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
					<< this->parser.GetExpr() << " as a boolean."
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
					<< this->parser.GetExpr() << " as a boolean."
					<< std::endl;
				abort();
			}

			return evaluated.GetFloat();

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



template<class Last_Variable> class Variable_To_Option_Vector<Last_Variable>
{
public:

	Variable_To_Option_Vector<Last_Variable>() :
		r_var(&r_val),
		t_var(&t_val)
	{
		this->parser.DefineVar("r", this->r_var);
		this->parser.DefineVar("t", this->t_var);
	}

	Variable_To_Option_Vector<Last_Variable>(
		const Variable_To_Option_Vector<Last_Variable>& other
	) : Variable_To_Option_Vector<Last_Variable>()
	{
		this->expressions = other.expressions;
	}

	Variable_To_Option_Vector<Last_Variable>(
		Variable_To_Option_Vector<Last_Variable>&& other
	) : Variable_To_Option_Vector<Last_Variable>()
	{
		this->expressions = std::move(other.expressions);
	}

	Variable_To_Option_Vector<Last_Variable>& operator=(
		const Variable_To_Option_Vector<Last_Variable>& other
	) {
		if (this != &other) {
			this->expressions = other.expressions;
		}

		return *this;
	}


	void set_expression(
		const Last_Variable&,
		const size_t index,
		const std::string given_expression
	) {
		if (index >= this->expressions.size()) {
			std::cerr <<  __FILE__ << "(" << __LINE__<< ") "
				<< "Invalid index: " << index
				<< ", should be less than " << this->expressions.size()
				<< std::endl;
			abort();
		}

		this->expressions[index] = given_expression;
	}


	void add_options(
		const std::string& prefix,
		boost::program_options::options_description& options
	) {
		this->add_options_impl(prefix, options);
	}


	typename Last_Variable::data_type get_data(
		const Last_Variable&,
		const size_t index,
		const std::array<double, 3>& given_position,
		const double given_time
	) {
		if (index >= this->expressions.size()) {
			std::cerr <<  __FILE__ << "(" << __LINE__<< ") "
				<< "Invalid index: " << index
				<< ", should be less than " << this->expressions.size()
				<< std::endl;
			abort();
		}

		this->parser.SetExpr(this->expressions[index]);

		return this->get_parsed_value<
			typename Last_Variable::data_type
		>(given_position, given_time);
	}

	typename Last_Variable::data_type get_data(
		const Last_Variable& variable,
		const size_t index
	) {
		return this->get_data(variable, index, {0, 0, 0}, 0);
	}


	size_t get_number_of_expressions() const
	{
		return this->get_number_of_expressions_impl();
	}


	void set_number_of_expressions(const size_t new_size)
	{
		this->set_number_of_expressions_impl(new_size);
	}


	void clear_expressions()
	{
		this->clear_expressions_impl();
	}



protected:

	void add_options_impl(
		const std::string& prefix,
		boost::program_options::options_description& options
	) {
		options.add_options()
			((prefix + Last_Variable::get_option_name()).c_str(),
				boost::program_options::value<
					std::vector<std::string>
				>(&this->expressions)
					->composing()
					->default_value(this->expressions),
				Last_Variable::get_option_help().c_str());
	}


	size_t get_number_of_expressions_impl() const
	{
		return this->expressions.size();
	}


	void set_number_of_expressions_impl(const size_t new_size)
	{
		this->expressions.resize(new_size);
	}


	void clear_expressions_impl()
	{
		this->expressions.clear();
	}



private:

	std::vector<std::string> expressions;
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

#endif // ifndef PAMHD_BOUNDARIES_VARIABLE_TO_OPTION_VECTOR_HPP
