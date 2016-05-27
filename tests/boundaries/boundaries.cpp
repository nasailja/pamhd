/*
Tests boundaries class for one simulation variable of PAMHD.

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


#include "cstdlib"
#include "iostream"
#include "set"
#include "utility"

#include "dccrg.hpp"
#include "mpi.h"

#include "boundaries/boundaries.hpp"
#include "boundaries/geometries.hpp"


using namespace pamhd::boundaries;


// storage for classifier algorithm
struct Cell_Type {
	using data_type = int;
	static const std::string get_name(){ return {"cell type"}; }
};

// simulation cell type used in this test
struct Cell_Data {
	typename Cell_Type::data_type type;

	std::tuple<void*, int, MPI_Datatype> get_mpi_datatype()
	{
		return std::make_tuple(static_cast<void*>(&this->type), 1, MPI_INT);
	}
};

// returns a reference to given cell's mass density data
const auto Type
	= [](Cell_Data& cell_data)
		->typename Cell_Type::data_type&
	{
		return cell_data.type;
	};


int main(int argc, char* argv[])
{
	const char json[] = "{"
		"\"geometries\": ["
			// 1 layer of cells on negative x side of simulation volume
			"{\"box\": {"
				"\"start\": [-9, -9, -9],"
				"\"end\": [0.8, 9, 9]"
			"}},"
			// close to center of sim volume
			"{\"sphere\": {"
				"\"center\": [2.5, 2.5, 2.5],"
				"\"radius\": 0.4"
			"}},"
			// positive z side of sim vol
			"{\"box\": {"
				"\"start\": [-9, -9, 3.2],"
				"\"end\": [9, 9, 9]"
			"}}"
		"],"
		"\"value_boundaries\": ["
			"{"
				"\"geometry_id\": 0,"
				"\"time_stamps\": [-3],"
				"\"values\": [2]"
			"},"
			"{"
				"\"geometry_id\": 1,"
				"\"time_stamps\": [1],"
				"\"values\": [\"t\"]"
			"}"
		"],"
		"\"copy_boundaries\": ["
			"{\"geometry_id\": 2}"
		"]"
	"}";

	rapidjson::Document document;
	document.Parse(json);
	if (document.HasParseError()) {
		std::cerr << "Couldn't parse json data: " << json << std::endl;
		return EXIT_FAILURE;
	}

	Geometries<int, std::array<double, 3>, double, uint64_t> geometries;
	geometries.set(document);

	Boundaries<uint64_t, int, Cell_Type> boundaries;
	boundaries.set(document);

	MPI_Init(&argc, &argv);

	dccrg::Dccrg<Cell_Data> grid;
	if (not grid.initialize(
		{4, 4, 4},
		MPI_COMM_WORLD,
		"RANDOM",
		0,
		0,
		false, false, false
	)) {
		std::cerr << __FILE__ << ":" << __LINE__
			<< ": Couldn't initialize grid."
			<< std::endl;
		return EXIT_FAILURE;
	}

	for (const auto& cell: grid.cells) {
		const auto
			start = grid.geometry.get_min(cell.id),
			end = grid.geometry.get_max(cell.id);
		geometries.overlaps(start, end, cell.id);
	}

	// check that cells have been assigned to correct geometries
	const std::set<uint64_t>
		ref_geom0_cells{1, 5, 9, 13, 17, 21, 25, 29, 33, 37, 41, 45},
		ref_geom1_cells{43},
		ref_geom2_cells{50, 51, 52, 54, 55, 56, 58, 59, 60, 62, 63, 64},
		ref_geom02_cells{49, 53, 57, 61}; // these might end up in either one

	const auto geom_ids = geometries.get_geometry_ids();
	if (geom_ids.size() != 3) {
		std::cerr << __FILE__ << ":" << __LINE__
			<< ": Wrong number of geometries: " << geom_ids.size()
			<< ", should be 3."
			<< std::endl;
		return EXIT_FAILURE;
	}


	const auto& geom0_cells = geometries.get_cells(0);
	unsigned int local_geom0_cells_size = geom0_cells.size(), geom0_cells_size = 0;
	MPI_Allreduce(
		&local_geom0_cells_size,
		&geom0_cells_size,
		1,
		MPI_UNSIGNED,
		MPI_SUM,
		MPI_COMM_WORLD
	);

	if (geom0_cells_size < ref_geom0_cells.size()) {
		std::cerr << __FILE__ << ":" << __LINE__
			<< ": Too few cells in geometry 0: " << geom0_cells_size
			<< ", should be at least " << ref_geom0_cells.size()
			<< std::endl;
		return EXIT_FAILURE;
	}

	for (const auto& cell: geom0_cells) {
		if (
			ref_geom0_cells.count(cell) == 0
			and ref_geom02_cells.count(cell) == 0
		) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< ": Cell " << cell << " not in geometries 0 or 2"
				<< std::endl;
			return EXIT_FAILURE;
		}
	}


	const auto& geom1_cells = geometries.get_cells(1);
	unsigned int local_geom1_cells_size = geom1_cells.size(), geom1_cells_size = 0;
	MPI_Allreduce(
		&local_geom1_cells_size,
		&geom1_cells_size,
		1,
		MPI_UNSIGNED,
		MPI_SUM,
		MPI_COMM_WORLD
	);

	if (geom1_cells_size != ref_geom1_cells.size()) {
		std::cerr << __FILE__ << ":" << __LINE__
			<< ": Wrong number of cells in geometry 1: " << geom1_cells_size
			<< ", should be " << ref_geom1_cells.size()
			<< std::endl;
		return EXIT_FAILURE;
	}

	for (const auto& cell: geom1_cells) {
		if (ref_geom1_cells.count(cell) == 0) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< ": Cell " << cell << " not in geometry 1"
				<< std::endl;
			return EXIT_FAILURE;
		}
	}


	const auto& geom2_cells = geometries.get_cells(2);
	unsigned int local_geom2_cells_size = geom2_cells.size(), geom2_cells_size = 0;
	MPI_Allreduce(
		&local_geom2_cells_size,
		&geom2_cells_size,
		1,
		MPI_UNSIGNED,
		MPI_SUM,
		MPI_COMM_WORLD
	);

	if (geom2_cells_size < ref_geom2_cells.size()) {
		std::cerr << __FILE__ << ":" << __LINE__
			<< ": Too few cells in geometry 2: " << geom2_cells_size
			<< ", should be at least " << ref_geom2_cells.size()
			<< std::endl;
		return EXIT_FAILURE;
	}

	for (const auto& cell: geom2_cells) {
		if (
			ref_geom2_cells.count(cell) == 0
			and ref_geom02_cells.count(cell) == 0
		) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< ": Cell " << cell << " not in geometries 0 or 2"
				<< std::endl;
			return EXIT_FAILURE;
		}
	}



	boundaries.classify(grid, Type, geometries);

	// check that cells have been classified correctly
	unsigned int
		local_normal_cells_size = boundaries.get_normal_cells().size(),
		normal_cells_size = 0;
	MPI_Allreduce(
		&local_normal_cells_size,
		&normal_cells_size,
		1,
		MPI_UNSIGNED,
		MPI_SUM,
		MPI_COMM_WORLD
	);

	if (normal_cells_size != 35) {
		std::cerr << __FILE__ << ":" << __LINE__
			<< ": Wrong number of normal cells: "
			<< normal_cells_size
			<< ", should be 35"
			<< std::endl;
		return EXIT_FAILURE;
	}

	for (const auto& cell: boundaries.get_normal_cells()) {
		if (
			ref_geom0_cells.count(cell) > 0
			or ref_geom1_cells.count(cell) > 0
			or ref_geom2_cells.count(cell) > 0
			or ref_geom02_cells.count(cell) > 0
		) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< ": Normal cell " << cell << " belongs to a boundary geometry."
				<< std::endl;
			return EXIT_FAILURE;
		}
	}


	unsigned int
		local_value_cells_size = boundaries.get_value_boundary_cells().size(),
		value_cells_size = 0;
	MPI_Allreduce(
		&local_value_cells_size,
		&value_cells_size,
		1,
		MPI_UNSIGNED,
		MPI_SUM,
		MPI_COMM_WORLD
	);

	if (value_cells_size < 13) {
		std::cerr << __FILE__ << ":" << __LINE__
			<< ": Wrong number of value boundary cells: "
			<< value_cells_size
			<< ", should be at least 13"
			<< std::endl;
		return EXIT_FAILURE;
	}

	for (const auto& cell: boundaries.get_value_boundary_cells()) {
		if (
			ref_geom0_cells.count(cell) == 0
			and ref_geom1_cells.count(cell) == 0
			and ref_geom02_cells.count(cell) == 0
		) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< ": Value boundary cell " << cell << " not in geometries 0 or 1."
				<< std::endl;
			return EXIT_FAILURE;
		}
	}


	unsigned int
		local_copy_cells_size = boundaries.get_copy_boundary_cells().size(),
		copy_cells_size = 0;
	MPI_Allreduce(
		&local_copy_cells_size,
		&copy_cells_size,
		1,
		MPI_UNSIGNED,
		MPI_SUM,
		MPI_COMM_WORLD
	);

	if (copy_cells_size < 11) {
		std::cerr << __FILE__ << ":" << __LINE__
			<< ": Wrong number of copy boundary cells: "
			<< copy_cells_size
			<< ", should be at least 11"
			<< std::endl;
		return EXIT_FAILURE;
	}

	const std::set<std::array<uint64_t, 2>> ref_copy_cells{
		{49, 33}, // target, source
		{50, 34}, {51, 35}, {52, 36}, {53, 37}, {54, 38},
		{55, 39}, {56, 40}, {57, 41}, {58, 42},
		//{59, 43}, 43 is value boundary so 59 should be dont_solve
		{60, 44}, {61, 45}, {62, 46}, {63, 47}, {64, 48}
	};

	for (const auto& item: boundaries.get_copy_boundary_cells()) {
		if (ref_copy_cells.count(item) == 0) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< ": Copy boundary item " << item[0] << "<-" << item[1]
				<< " not in copy boundary cells."
				<< std::endl;
			return EXIT_FAILURE;
		}
	}


	unsigned int
		local_dont_cells_size = boundaries.get_dont_solve_cells().size(),
		dont_cells_size = 0;
	MPI_Allreduce(
		&local_dont_cells_size,
		&dont_cells_size,
		1,
		MPI_UNSIGNED,
		MPI_SUM,
		MPI_COMM_WORLD
	);

	if (dont_cells_size != 5) {
		std::cerr << __FILE__ << ":" << __LINE__
			<< ": Wrong number of dont solve cells: "
			<< dont_cells_size
			<< ", should be 5"
			<< std::endl;
		return EXIT_FAILURE;
	}

	for (const auto& cell: boundaries.get_dont_solve_cells())
		if (ref_geom02_cells.count(cell) == 0 and cell != 59) {
		std::cerr << __FILE__ << ":" << __LINE__
			<< ": Dont solve cell not in geometry 02 and not 59"
			<< std::endl;
		return EXIT_FAILURE;
	}


	MPI_Finalize();
	return EXIT_SUCCESS;
}
