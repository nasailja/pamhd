/*
Class that translates simulation variables to boost::program_options.

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

#ifndef PAMHD_VARIABLES_TO_OPTIONS_HPP
#define PAMHD_VARIABLES_TO_OPTIONS_HPP


namespace pamhd {
namespace boundaries {


template<class... Variables> class Variables_To_Options {};



template<
	class Current_Variable,
	class... Rest_Of_Variables
> class Variables_To_Options<Current_Variable, Rest_Of_Variables...> :
	public Variables_To_Options<Rest_Of_Variables...>
{
public:

	using Variables_To_Options<Rest_Of_Variables...>::add_options;
	using Variables_To_Options<Rest_Of_Variables...>::operator[];

	void add_options(
		boost::program_options::options_description& options,
		const std::string& option_name_prefix
	) {
		// make text representation of default value
		std::string values_str(" ");
		for (const auto& value: this->values) {
			values_str += boost::lexical_cast<std::string>(value) + " ";
		}

		options.add_options()
			((option_name_prefix + Current_Variable::get_option_name()).c_str(), 
				boost::program_options::value<
					std::vector<typename Current_Variable::data_type>
				>(&this->values)
					->multitoken()
					->default_value(
						this->values,
						values_str.c_str()
					),
				Current_Variable::get_option_help().c_str());

		Variables_To_Options<Rest_Of_Variables...>::add_options(
			options,
			option_name_prefix
		);
	}


	std::vector<
		typename Current_Variable::data_type
	>& operator[](
		const Current_Variable&
	) {
		return this->values;
	}

	const std::vector<
		typename Current_Variable::data_type
	>& operator[](
		const Current_Variable&
	) const {
		return this->values;
	}


protected:
	std::vector<typename Current_Variable::data_type> values;
};



template<class Last_Variable> class Variables_To_Options<Last_Variable>
{
public:

	void add_options(
		boost::program_options::options_description& options,
		const std::string& option_name_prefix
	) {
		// make text representation of default value
		std::string values_str(" ");
		for (const auto& value: this->values) {
			values_str += boost::lexical_cast<std::string>(value) + " ";
		}

		options.add_options()
			((option_name_prefix + Last_Variable::get_option_name()).c_str(), 
				boost::program_options::value<
					std::vector<typename Last_Variable::data_type>
				>(&this->values)
					->multitoken()
					->default_value(
						this->values,
						values_str.c_str()
					),
				Last_Variable::get_option_help().c_str());
	}


	std::vector<
		typename Last_Variable::data_type
	>& operator[](
		const Last_Variable&
	) {
		return this->values;
	}

	const std::vector<
		typename Last_Variable::data_type
	>& operator[](
		const Last_Variable&
	) const {
		return this->values;
	}


protected:
	std::vector<typename Last_Variable::data_type> values;
};




}} // namespaces

#endif // ifndef PAMHD_VARIABLES_TO_OPTIONS_HPP
