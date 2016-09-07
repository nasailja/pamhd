/*
Handles options of MHD part of PAMHD.

Copyright 2016 Ilja Honkonen
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

#ifndef PAMHD_MHD_OPTIONS_HPP
#define PAMHD_MHD_OPTIONS_HPP


#include "cmath"
#include "string"

#include "rapidjson/document.h"


namespace pamhd {
namespace mhd {


struct Options
{
	Options() = default;
	Options(const Options& other) = default;
	Options(Options&& other) = default;

	Options(const rapidjson::Value& object)
	{
		this->set(object);
	};


	std::string
		solver = "roe_athena",
		lb_name = "RCB",
		output_directory = "";

	size_t
		poisson_iterations_max = 1000,
		poisson_iterations_min = 0;

	double
		save_mhd_n = -1,
		time_start = 0,
		time_length = 1,
		time_step_factor = 0.5,
		remove_div_B_n = 0.1,
		poisson_norm_stop = 1e-10,
		poisson_norm_increase_max = 10,
		adiabatic_index = 5.0 / 3.0,    
		vacuum_permeability = 4e-7 * M_PI,
		proton_mass = 1.672621777e-27,
		min_pressure = 0;

	void set(const rapidjson::Value& object) {
		using std::isnormal;

		if (not object.HasMember("save-mhd-n")) {
			throw std::invalid_argument(
				std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
				+ "JSON data doesn't have a save-mhd-n key."
			);
		}
		save_mhd_n = object["save-mhd-n"].GetDouble();

		if (not object.HasMember("time-start")) {
			throw std::invalid_argument(
				std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
				+ "JSON data doesn't have a time-start key."
			);
		}
		time_start = object["time-start"].GetDouble();

		if (not object.HasMember("time-length")) {
			throw std::invalid_argument(
				std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
				+ "JSON data doesn't have a time-length key."
			);
		}
		time_length = object["time-length"].GetDouble();

		if (not object.HasMember("time-step-factor")) {
			throw std::invalid_argument(
				std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
				+ "JSON data doesn't have a time-step-factor key."
			);
		}
		time_step_factor = object["time-step-factor"].GetDouble();

		if (not object.HasMember("remove-div-B-n")) {
			throw std::invalid_argument(
				std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
				+ "JSON data doesn't have a remove-div-B-n key."
			);
		}
		remove_div_B_n = object["remove-div-B-n"].GetDouble();

		if (not object.HasMember("poisson-norm-stop")) {
			throw std::invalid_argument(
				std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
				+ "JSON data doesn't have a poisson-norm-stop key."
			);
		}
		poisson_norm_stop = object["poisson-norm-stop"].GetDouble();

		if (not object.HasMember("poisson-norm-increase-max")) {
			throw std::invalid_argument(
				std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
				+ "JSON data doesn't have a poisson-norm-increase-max key."
			);
		}
		poisson_norm_increase_max = object["poisson-norm-increase-max"].GetDouble();

		if (not object.HasMember("adiabatic-index")) {
			throw std::invalid_argument(
				std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
				+ "JSON data doesn't have a adiabatic-index key."
			);
		}
		adiabatic_index = object["adiabatic-index"].GetDouble();
		if (not isnormal(adiabatic_index) or adiabatic_index < 0) {
			throw std::invalid_argument(
				std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
				+ "Invalid adiabatic index: " + std::to_string(adiabatic_index)
				+ ", should be > 0"
			);
		}

		if (not object.HasMember("vacuum-permeability")) {
			throw std::invalid_argument(
				std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
				+ "JSON data doesn't have a vacuum-permeability key."
			);
		}
		vacuum_permeability = object["vacuum-permeability"].GetDouble();
		if (not isnormal(vacuum_permeability) or vacuum_permeability < 0) {
			throw std::invalid_argument(
				std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
				+ "Invalid vacuum permeability: " + std::to_string(vacuum_permeability)
				+ ", should be > 0"
			);
		}

		if (not object.HasMember("proton-mass")) {
			throw std::invalid_argument(
				std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
				+ "JSON data doesn't have a proton-mass key."
			);
		}
		proton_mass = object["proton-mass"].GetDouble();
		if (not isnormal(proton_mass) or proton_mass < 0) {
			throw std::invalid_argument(
				std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
				+ "Invalid proton_mass: " + std::to_string(proton_mass)
				+ ", should be > 0"
			);
		}

		if (not object.HasMember("solver-mhd")) {
			throw std::invalid_argument(
				std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
				+ "JSON data doesn't have a solver-mhd key."
			);
		}
		solver = object["solver-mhd"].GetString();
		if (
			solver != "hll-athena"
			and solver != "hlld-athena"
			and solver != "roe-athena"
		) {
			throw std::invalid_argument(
				std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
				+ "Invalid mhd solver: " + solver
				+ ", should be one of hll-athena, hlld-athena, roe-athena."
			);
		}

		if (object.HasMember("output-directory")) {
			output_directory = object["output-directory"].GetString();
		}

		if (not object.HasMember("load-balancer")) {
			throw std::invalid_argument(
				std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
				+ "JSON data doesn't have a load-balancer key."
			);
		}
		lb_name = object["load-balancer"].GetString();

		if (not object.HasMember("minimum-pressure")) {
			throw std::invalid_argument(
				std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
				+ "JSON data doesn't have a minimum-pressure key."
			);
		}
		const auto& min_pressure_json = object["minimum-pressure"];
		if (not min_pressure_json.IsNumber()) {
			throw std::invalid_argument(
				std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
				+ "JSON item minimum-pressure is not a number."
			);
		}
		min_pressure = min_pressure_json.GetDouble();
	}
};

}} // namespaces


#endif // ifndef PAMHD_MHD_OPTIONS_HPP
