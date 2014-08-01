/*
Copy boundary class of PAMHD.

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


#ifndef PAMHD_BOUNDARIES_COPY_BOUNDARY_HPP
#define PAMHD_BOUNDARIES_COPY_BOUNDARY_HPP


#include "cstdlib"
#include "iostream"
#include "string"
#include "unordered_map"
#include "unordered_set"
#include "utility"
#include "vector"

#include "boost/optional.hpp"
#include "boost/program_options.hpp"
#include "prettyprint.hpp"

#include "boundaries/box.hpp"
#include "boundaries/sphere.hpp"


namespace pamhd {
namespace boundaries {


template<
	class Cell_T,
	class Scalar_T,
	class Vector_T
> class Copy_Boundary
{
public:

	using cell_type = Cell_T;
	using scalar_type = Scalar_T;
	using vector_type = Vector_T;


	/*!
	Adds options for setting the number of boundaries.

	Use to query the number of boundaries to create whose options
	can then be added with add_options().

	Must not be called after add_options().
	*/
	void add_initialization_options(
		const std::string& option_name_prefix,
		boost::program_options::options_description& options
	) {
		options.add_options()
			((option_name_prefix + "nr-boxes").c_str(),
				boost::program_options::value<size_t>(&this->number_of_boxes)
					->default_value(this->number_of_boxes),
				"Number of boxes in value boundary condition")
			((option_name_prefix + "nr-spheres").c_str(),
				boost::program_options::value<size_t>(&this->number_of_spheres)
					->default_value(this->number_of_spheres),
				"Number of spheres in value boundary condition");
	}


	/*!
	Add options for as many boundaries as given by the nr-*
	variables set by add_initialization_options.
	*/
	void add_options(
		const std::string& option_name_prefix,
		boost::program_options::options_description& options
	) {
		this->boxes.resize(this->number_of_boxes);
		this->spheres.resize(this->number_of_spheres);

		size_t counter;

		counter = 0;
		for (auto& box: this->boxes) {
			counter++;
			box.add_options(
				option_name_prefix
					+ "box"
					+ boost::lexical_cast<std::string>(counter)
					+ ".",
				options
			);
		}

		counter = 0;
		for (auto& sphere: this->spheres) {
			counter++;
			sphere.add_options(
				option_name_prefix
					+ "sphere"
					+ boost::lexical_cast<std::string>(counter)
					+ ".",
				options
			);
		}
	}


	/*!
	geometry_parameters is given to the overlaps function of each geometry
	object, see their documentation for the required parameters.

	Returns number of boundaries to which given cell was added.
	*/
	template<class... Geometry_Parameters> size_t add_cell(
		const Cell_T& cell,
		Geometry_Parameters&&... geometry_parameters
	) {
		size_t ret_val = 0;

		for (auto& box: this->boxes) {
			if (
				box.overlaps(
					std::forward<Geometry_Parameters>(geometry_parameters)...
				)
			) {
				ret_val++;
			}
		}

		for (auto& sphere: this->spheres) {
			if (
				sphere.overlaps(
					std::forward<Geometry_Parameters>(geometry_parameters)...
				)
			) {
				ret_val++;
			}
		}

		if (ret_val > 0) {
			this->sources[cell];
		}

		return ret_val;
	}


	/*!
	Excludes given cell as source for cell belonging to this boundary.
	*/
	void add_as_other_boundary(const Cell_T& cell)
	{
		this->other_boundary_cells.insert(cell);
	}


	/*!
	Sets source cells of given cell based on given neighbor list.

	Boundary cells are excluded from the given neighbor list.

	Must not be called before classifying all cells with add_cell()
	and excluding other boundary cells with add_as_other_boundary().

	Returns number of non-boundary source cells for given cell.

	Replaces previous sources.
	*/
	size_t set_neighbors_of(
		const Cell_T& cell,
		const std::vector<Cell_T>& neighbors_of
	) {
		if (
			this->sources.count(cell) == 0
			and this->other_boundary_cells.count(cell) == 0
		) {
			std::cerr <<  __FILE__ << "(" << __LINE__<< ") "
				<< "Given cell " << cell << " is not a boundary cell."
				<< std::endl;
			abort();
		}

		this->sources.at(cell).clear();

		for (const auto& neighbor_of: neighbors_of) {
			if (
				this->sources.count(neighbor_of) == 0
				and this->other_boundary_cells.count(neighbor_of) == 0
			) {
				this->sources.at(cell).push_back(neighbor_of);
			}
		}

		return this->sources.at(cell).size();
	}


	const std::vector<Cell_T>& get_source_cells(const Cell_T& cell)
	{
		return this->sources.at(cell);
	}


	//! Returns number of existing geometries
	size_t get_number_of_geometries() const
	{
		return this->boxes.size() + this->spheres.size();
	}


	/*!
	Returns cells which belong to this boundary.
	*/
	std::vector<Cell_T> get_cells() const
	{
		std::vector<Cell_T> ret_val;

		for (const auto& item: this->sources) {
			ret_val.push_back(item.first);
		};

		return ret_val;
	}


	void clear_temporary_data()
	{
		this->cells.clear();
		this->other_boundary_cells.clear();
	}


	void clear_cells()
	{
		this->clear_temporary_data();
		this->sources.clear();
	}


	void clear_boundaries()
	{
		this->boxes.clear();
		this->spheres.clear();
	}


	/*!
	Returns boundary index of added boundary if successfull.
	*/
	template<class... Geometry_Parameters> boost::optional<size_t> add_box(
		Geometry_Parameters&&... geometry_parameters
	) {
		Box<Vector_T> new_box;

		if (
			not new_box.set_geometry(
				std::forward<Geometry_Parameters>(geometry_parameters)...
			)
		) {
			return boost::optional<size_t>();
		}

		this->boxes.push_back(new_box);

		this->number_of_boxes = this->boxes.size();

		return boost::optional<size_t>(this->number_of_boxes - 1);
	}


	/*!
	Same as add_box() but adds a sphere.
	*/
	template<class... Geometry_Parameters> boost::optional<size_t> add_sphere(
		Geometry_Parameters&&... geometry_parameters
	) {
		Sphere<Vector_T, Scalar_T> new_sphere;

		if (
			not new_sphere.set_geometry(
				std::forward<Geometry_Parameters>(geometry_parameters)...
			)
		) {
			return boost::optional<size_t>();
		}

		this->spheres.push_back(new_sphere);

		this->number_of_spheres = this->spheres.size();

		return boost::optional<size_t>(
			this->number_of_boxes + this->number_of_spheres - 1
		);
	}



private:

	size_t
		number_of_boxes = 0,
		number_of_spheres = 0;

	std::vector<Box<Vector_T>> boxes;
	std::vector<Sphere<Vector_T, Scalar_T>> spheres;

	// source cells from which cells in this boundary should copy
	std::unordered_map<Cell_T, std::vector<Cell_T>> sources;

	std::unordered_set<Cell_T> other_boundary_cells;
};

}} // namespaces

#endif // ifndef PAMHD_BOUNDARIES_COPY_BOUNDARY_HPP
