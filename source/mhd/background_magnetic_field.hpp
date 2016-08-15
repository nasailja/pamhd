/*
Handles background magnetic field of MHD part of PAMHD.

Copyright 2014, 2016 Ilja Honkonen
All rights reserved.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef PAMHD_MHD_BACKGROUND_MAGNETIC_FIELD_HPP
#define PAMHD_MHD_BACKGROUND_MAGNETIC_FIELD_HPP


#include "cmath"
#include "string"

#include "rapidjson/document.h"


namespace pamhd {
namespace mhd {


//! Vector assumed to support Eigen API.
template<class Vector> class Background_Magnetic_Field
{
	static_assert(
		Vector::SizeAtCompileTime == 3,
		"Only 3 component Eigen vectors supported"
	);

private:
	Vector constant = {0, 0, 0};
	std::vector<std::pair<Vector, Vector>> dipole_moments_positions;
	double min_distance = 0;


public:
	Background_Magnetic_Field() = default;
	Background_Magnetic_Field(const Background_Magnetic_Field& other) = default;
	Background_Magnetic_Field(Background_Magnetic_Field&& other) = default;
	Background_Magnetic_Field(const rapidjson::Value& object)
	{
		this->set(object);
	};


	void set(const rapidjson::Value& object) {
		using std::isnormal;

		if (not object.HasMember("background-magnetic-field")) {
			return;
		}
		const auto& bg_B = object["background-magnetic-field"];

		if (bg_B.HasMember("value")) {
			const auto& value = bg_B["value"];
			if (not value.IsArray()) {
				throw std::invalid_argument(
					std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
					+ "Background magnetic field value is not an array."
				);
			}

			if (value.Size() != 3) {
				throw std::invalid_argument(
					std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
					+ "Background magnetic field value must have 3 components."
				);
			}

			this->constant[0] = value[0].GetDouble();
			this->constant[1] = value[1].GetDouble();
			this->constant[2] = value[2].GetDouble();
		}

		if (bg_B.HasMember("minimum-distance")) {
			const auto& min_distance_json = bg_B["minimum-distance"];
			if (not min_distance_json.IsDouble()) {
				throw std::invalid_argument(
					std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
					+ "Background magnetic field minimum distance is not a double."
				);
			}

			this->min_distance = min_distance_json.GetDouble();
		}

		if (bg_B.HasMember("dipoles")) {
			const auto& dipoles = bg_B["dipoles"];
			if (not dipoles.IsArray()) {
				throw std::invalid_argument(
					std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
					+ "Dipoles item must be an array."
				);
			}

			for (rapidjson::SizeType i = 0; i < dipoles.Size(); i++) {
				if (not dipoles[i].HasMember("moment")) {
					throw std::invalid_argument(
						std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
						+ std::to_string(i + 1) + "th dipole moment missing a moment member."
					);
				}
				const auto& moment_json = dipoles[i]["moment"];
				if (not moment_json.IsArray()) {
					throw std::invalid_argument(
						std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
						+ std::to_string(i + 1) + "th dipole moment must be an array."
					);
				}
				if (moment_json.Size() != 3) {
					throw std::invalid_argument(
						std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
						+ std::to_string(i + 1) + "th dipole moment array must have 3 components."
					);
				}

				const Vector moment{
					moment_json[0].GetDouble(),
					moment_json[1].GetDouble(),
					moment_json[2].GetDouble(),
				};

				if (not dipoles[i].HasMember("position")) {
					throw std::invalid_argument(
						std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
						+ std::to_string(i + 1) + "th dipole moment missing a position member."
					);
				}

				const auto& position_json = dipoles[i]["position"];
				if (not position_json.IsArray()) {
					throw std::invalid_argument(
						std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
						+ std::to_string(i + 1) + "th dipole position must be an array."
					);
				}
				if (position_json.Size() != 3) {
					throw std::invalid_argument(
						std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
						+ std::to_string(i + 1) + "th dipole position array must have 3 components."
					);
				}

				const Vector position{
					position_json[0].GetDouble(),
					position_json[1].GetDouble(),
					position_json[2].GetDouble(),
				};

				this->dipole_moments_positions.emplace_back(moment, position);
			}
		}
	}


	/*
	Modified from https://github.com/nasailja/background_B/blob/master/source/dipole.hpp

	Returns magnetic field from given sources at given position.
	If given position is closer to a dipole than this->min_distance,
	that dipole contributes 2/3*vacuum_permeability*moment instead.
	*/
	Vector get_background_field(
		const Vector& field_position,
		const double& vacuum_permeability
	) const {
		using std::pow;

		Vector ret_val{this->constant};

		for (const auto dip_mom_pos: this->dipole_moments_positions) {

			const Vector r = field_position - dip_mom_pos.second;

			if (r.squaredNorm() < pow(this->min_distance, 2)) {
				ret_val += 2.0 / 3.0 * vacuum_permeability * dip_mom_pos.first;
				continue;
			}

			const double r3 = r.norm() * r.squaredNorm();
			const Vector
				r_unit = r / r.norm(),
				projected_dip_mom = dip_mom_pos.first.dot(r_unit) * r_unit;

			ret_val
				+= vacuum_permeability / (4 * M_PI)
				* (3 * projected_dip_mom - dip_mom_pos.first) / r3;
		}

		return ret_val;
	}

};

}} // namespaces


#endif // ifndef PAMHD_MHD_BACKGROUND_MAGNETIC_FIELD_HPP
