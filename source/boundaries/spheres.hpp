/*
Sphere geometry class of PAMHD supporting multiple spheres.

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


#ifndef PAMHD_BOUNDARIES_SPHERES_HPP
#define PAMHD_BOUNDARIES_SPHERES_HPP


#include "cstdlib"
#include "iostream"
#include "string"

#include "boost/program_options.hpp"
#include "prettyprint.hpp"

#include "program_options_validators.hpp"


namespace pamhd {
namespace boundaries {

/*!
Vector_T is assumed to be std::array or similar.
*/
template<class Vector_T, class Scalar_T> class Spheres
{
public:

	//! Class can read more than one instance using boost::program_options
	static constexpr bool scalar = false;


	size_t get_number_of_instances() const
	{
		return this->centers.size();
	}


	void set_number_of_instances(const size_t new_size)
	{
		this->centers.resize(new_size);
		this->radii.resize(new_size);
	}


	/*!
	Returns true if cell spanning given volume overlaps
	this sphere, false otherwise.
	*/
	bool overlaps(
		const size_t index,
		const Vector_T& cell_start,
		const Vector_T& cell_end
	) {
		if (not this->check_index(index)) {
			abort();
		}

		for (size_t i = 0; i < size_t(cell_start.size()); i++) {
			if (cell_start[i] > cell_end[i]) {
				std::cerr <<  __FILE__ << "(" << __LINE__<< ") "
					<< "Starting coordinate of cell at index " << i
					<< " is larger than ending coordinate: "
					<< cell_start[i] << " > " << cell_end[i]
					<< std::endl;
				abort();
			}
		}

		// https://stackoverflow.com/questions/4578967
		Scalar_T distance2 = 0;
		for (size_t i = 0; i < size_t(cell_start.size()); i++) {
			if (this->centers[index][i] < cell_start[i]) {
				distance2 += std::pow(this->centers[index][i] - cell_start[i], 2);
			} else if (this->centers[index][i] > cell_end[i]) {
				distance2 += std::pow(this->centers[index][i] - cell_end[i], 2);
			}
		}

		return distance2 < std::pow(this->radii[index], 2);
	}


	void add_options(
		const std::string& option_name_prefix,
		boost::program_options::options_description& options
	) {
		options.add_options()
			((option_name_prefix + "center").c_str(),
				boost::program_options::value<
					std::vector<Vector_T>
				>(&this->centers)
					->composing()
					->default_value(this->centers),
				"Center of the sphere")
			((option_name_prefix + "radius").c_str(),
				boost::program_options::value<
					std::vector<Scalar_T>
				>(&this->radii)
					->composing()
					->default_value(this->radii),
				"Radius of the sphere");
	}


	bool set_geometry(
		const size_t index,
		const Vector_T& given_center,
		const Scalar_T& given_radius
	) {
		if (not this->check_index(index)) {
			return false;
		}

		if (given_radius <= 0) {
			return false;
		}

		this->centers[index] = given_center;
		this->radii[index] = given_radius;
		return true;
	}

	const Vector_T& get_center(const size_t index) const
	{
		if (not this->check_index(index)) {
			abort();
		}

		return this->centers[index];
	}

	const Scalar_T& get_radius(const size_t index) const
	{
		if (not this->check_index(index)) {
			abort();
		}

		return this->radii[index];
	}


private:

	std::vector<Vector_T> centers;
	std::vector<Scalar_T> radii;


	bool check_index(const size_t index)
	{
		if (this->centers.size() != this->radii.size()) {
			std::cerr <<  __FILE__ << "(" << __LINE__<< ") "
				<< "Number of centers and radii for spheres differ: "
				<< this->centers.size() << ", " << this->radii.size()
				<< std::endl;
			return false;
		}

		if (index >= this->centers.size()) {
			std::cerr <<  __FILE__ << "(" << __LINE__<< ") "
				<< "Invalid index: " << index
				<< ", should be less than " << this->centers.size()
				<< std::endl;
			return false;
		}

		return true;
	}
};

}} // namespaces

#endif // ifndef PAMHD_BOUNDARIES_SPHERES_HPP
