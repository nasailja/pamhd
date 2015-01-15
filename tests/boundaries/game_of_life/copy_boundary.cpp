/*
Game of life test for copy boundary class of PAMHD.

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

#include "array"
#include "cstdlib"
#include "iostream"
#include "set"
#include "string"
#include "vector"

#include "boost/program_options.hpp"

#include "boundaries/copy_boundary.hpp"


using namespace std;
using namespace pamhd::boundaries;


struct Cell {
	bool is_alive = false;
	int live_neighbors = 0;
};


constexpr size_t
	width = 12,
	height = 12;

using Grid = array<array<Cell, width>, height>; // Grid[row][column]
Grid grid;

// row == 0 is top edge, column == 0 left edge
array<double, 3> get_cell_start(
	const Grid& grid,
	const array<size_t, 2>& cell
) {
	if (
		cell[1] >= grid.size()
		or cell[0] >= grid[cell[1]].size()
	) {
		return {{
			std::numeric_limits<double>::quiet_NaN(),
			std::numeric_limits<double>::quiet_NaN(),
			std::numeric_limits<double>::quiet_NaN()
		}};
	}

	return {{double(cell[1]), double(cell[0]), 0}};
}


array<double, 3> get_cell_end(
	const Grid& grid,
	const array<size_t, 2>& cell
) {
	if (
		cell[1] >= grid.size()
		or cell[0] >= grid[cell[1]].size()
	) {
		return {{
			std::numeric_limits<double>::quiet_NaN(),
			std::numeric_limits<double>::quiet_NaN(),
			std::numeric_limits<double>::quiet_NaN()
		}};
	}

	return {{double(cell[1] + 1), double(cell[0] + 1), 1}};
}


/*
Prints given game of life.

Prints live cells as 0, dead as .,
if a cell is a copy boundary cell
prints O if alive and , if dead.
*/
template <class Copy_Boundary_T> void print_game(
	const Grid& grid,
	const Copy_Boundary_T& boundary
) {
	std::unordered_set<typename Copy_Boundary_T::cell_type> boundary_cells;
	for (const auto& cell: boundary.get_cells()) {
		boundary_cells.insert(cell);
	}

	for (size_t row_i = 0; row_i < grid.size(); row_i++) {
	for (size_t col_i = 0; col_i < grid.size(); col_i++) {
			const auto& cell = grid[row_i][col_i];

			if (boundary_cells.count({row_i, col_i}) > 0) {
				if (cell.is_alive > 0) {
					cout << "O";
				} else {
					cout << ",";
				}
			} else {
				if (cell.is_alive > 0) {
					cout << "0";
				} else {
					cout << ".";
				}
			}
		}
		cout << endl;
	}
	cout << endl;
}


size_t get_number_of_live_cells(const Grid& grid)
{
	size_t live_cells = 0;
	for (const auto& row: grid) {
		for (const auto& cell: row) {
			if (cell.is_alive == true) {
				live_cells++;
			}
		}
	}
	return live_cells;
}

#define MAKE_VEC(x, y, z) std::array<double, 3>{{x, y, z}}
#define MAKE_BOX(a, b, c, d, e, f) MAKE_VEC(a, b, c), MAKE_VEC(d, e, f)


/*
Make unordered_map<array<T, N>, ...> work:
https://stackoverflow.com/questions/8026890
*/
namespace std
{
	template<class T, size_t N> struct hash<array<T, N>>
	{
		using argument_type = array<T, N>;
		using result_type = size_t;

		result_type operator()(const argument_type& a) const
		{
			hash<T> hasher;
			result_type h = 0;
			for (const auto& item: a) {
				h = h * 31 + hasher(item);
			}
			return h;
		}
	};
}

/*
Parts of GoL patterns are copy boundaries.
*/
int main(int argc, char* argv[])
{
	std::array<size_t, 2> cell;
	std::vector<std::array<size_t, 2>> neighbors_of;

	/*
	Set initial states
	*/

	// upper half of block
	grid[1][1].is_alive =
	grid[1][2].is_alive = true;

	// middle and upper part of blinker
	grid[1][8].is_alive =
	grid[2][8].is_alive = true;

	// upper half of beacon
	grid[6][1].is_alive =
	grid[6][2].is_alive =
	grid[7][1].is_alive =
	grid[7][2].is_alive = true;


	Copy_Boundary<
		std::array<size_t, 2>,
		double,
		std::array<double, 3>
	> boundary;

	/*
	Set copy boundary geometries
	*/
	boost::optional<size_t> result;

	// lower half of block
	result = boundary.add_box(MAKE_BOX(1.1, 2.1, 0, 2.9, 2.9, 1));
	if (not result or *result != 0) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
			<< "Invalid boundary id for added box: " << *result
			<< std::endl;
		abort();
	}

	// lower and right part of blinker
	result = boundary.add_box(MAKE_BOX(8.1, 3.1, 0, 8.9, 3.9, 1));
	if (not result or *result != 1) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
			<< "Invalid boundary id for added box: " << *result
			<< std::endl;
		abort();
	}
	result = boundary.add_box(MAKE_BOX(9.1, 2.1, 0, 9.9, 2.9, 1));
	if (not result or *result != 2) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
			<< "Invalid boundary id for added box: " << *result
			<< std::endl;
		abort();
	}

	// lower half of beacon
	result = boundary.add_box(MAKE_BOX(3.1, 8.1, 0, 4.9, 9.9, 1));
	if (not result or *result != 3) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
			<< "Invalid boundary id for added box: " << *result
			<< std::endl;
		abort();
	}


	std::string copy_boundary_file_name("");

	boost::program_options::options_description
		general_options(
			"Usage: program_name [options], where options are"
		),
		copy_boundary_options(
			"Options related to copy boundary condition"
		);

	general_options.add_options()
		("help", "Print options and their descriptions")
		("copy-boundary-help", "Print help for options of copy boundary condition")
		("copy-boundary-file",
			boost::program_options::value<std::string>(&copy_boundary_file_name)
				->default_value(copy_boundary_file_name),
			"Read copy boundary condition from file arg (not read if empty string)");
	boundary.add_initialization_options("", general_options);

	boost::program_options::variables_map option_variables;
	boost::program_options::store(
		boost::program_options::parse_command_line(
			argc,
			argv,
			general_options
		),
		option_variables
	);
	boost::program_options::notify(option_variables);

	if (option_variables.count("help") > 0) {
		std::cout << general_options << std::endl;
		return EXIT_SUCCESS;
	}

	boundary.add_options("", copy_boundary_options);

	if (option_variables.count("initial-help") > 0) {
		std::cout << copy_boundary_options << std::endl;
		return EXIT_SUCCESS;
	}

	if (copy_boundary_file_name != "") {
		boost::program_options::store(
			boost::program_options::parse_config_file<char>(
				copy_boundary_file_name.c_str(),
				copy_boundary_options
			),
			option_variables
		);
		boost::program_options::notify(option_variables);
	}

	// cells in lower half of block
	cell = {{2, 1}};
	if (
		boundary.add_cell(
			cell,
			get_cell_start(grid, cell),
			get_cell_end(grid, cell)
		) != 1
	) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< ")" << std::endl;
		abort();
	}
	cell = {{2, 2}};
	if (
		boundary.add_cell(
			cell,
			get_cell_start(grid, cell),
			get_cell_end(grid, cell)
		) != 1
	) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< ")" << std::endl;
		abort();
	}

	// copy cells of blinker
	cell = {{3, 8}};
	if (
		boundary.add_cell(
			cell,
			get_cell_start(grid, cell),
			get_cell_end(grid, cell)
		) != 1
	) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< ")" << std::endl;
		abort();
	}
	cell = {{2, 9}};
	if (
		boundary.add_cell(
			cell,
			get_cell_start(grid, cell),
			get_cell_end(grid, cell)
		) != 1
	) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< ")" << std::endl;
		abort();
	}

	// beacon
	cell = {{8, 3}};
	if (
		boundary.add_cell(
			cell,
			get_cell_start(grid, cell),
			get_cell_end(grid, cell)
		) != 1
	) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< ")" << std::endl;
		abort();
	}
	cell = {{8, 4}};
	if (
		boundary.add_cell(
			cell,
			get_cell_start(grid, cell),
			get_cell_end(grid, cell)
		) != 1
	) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< ")" << std::endl;
		abort();
	}
	cell = {{9, 3}};
	if (
		boundary.add_cell(
			cell,
			get_cell_start(grid, cell),
			get_cell_end(grid, cell)
		) != 1
	) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< ")" << std::endl;
		abort();
	}
	cell = {{9, 4}};
	if (
		boundary.add_cell(
			cell,
			get_cell_start(grid, cell),
			get_cell_end(grid, cell)
		) != 1
	) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< ")" << std::endl;
		abort();
	}


	/*
	Set source cells for copy boundary cells
	*/

	// block
	cell = {{2, 1}};
	neighbors_of = {{{1, 1}}};
	if (boundary.set_neighbors_of(cell, neighbors_of) != 1) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< ")" << std::endl;
		abort();
	}
	cell = {{2, 2}};
	neighbors_of = {{{1, 2}}};
	if (boundary.set_neighbors_of(cell, neighbors_of) != 1) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< ")" << std::endl;
		abort();
	}

	// blinker
	cell = {{3, 8}};
	neighbors_of = {{{1, 8}}};
	if (boundary.set_neighbors_of(cell, neighbors_of) != 1) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< ")" << std::endl;
		abort();
	}
	cell = {{2, 9}};
	neighbors_of = {{{2, 7}}};
	if (boundary.set_neighbors_of(cell, neighbors_of) != 1) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< ")" << std::endl;
		abort();
	}

	// bacon
	cell = {{8, 3}};
	neighbors_of = {{{7, 2}}};
	if (boundary.set_neighbors_of(cell, neighbors_of) != 1) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< ")" << std::endl;
		abort();
	}
	cell = {{8, 4}};
	neighbors_of = {{{6, 2}}};
	if (boundary.set_neighbors_of(cell, neighbors_of) != 1) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< ")" << std::endl;
		abort();
	}
	cell = {{9, 3}};
	neighbors_of = {{{7, 1}}};
	if (boundary.set_neighbors_of(cell, neighbors_of) != 1) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< ")" << std::endl;
		abort();
	}
	cell = {{9, 4}};
	neighbors_of = {{{6, 1}}};
	if (boundary.set_neighbors_of(cell, neighbors_of) != 1) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< ")" << std::endl;
		abort();
	}


	// play
	constexpr size_t max_turns = 20;
	for (size_t turn = 0; turn < max_turns; turn++) {

		// copy is_alive into copy boundary cells from their source cell
		for (const auto& cell: boundary.get_cells()) {
			const auto& source_cells = boundary.get_source_cells(cell);
			if (source_cells.size() != 1) {
				std::cerr <<  __FILE__ << "(" << __LINE__<< ")" << std::endl;
				abort();
			}

			grid[cell[0]][cell[1]]
				= grid[source_cells[0][0]][source_cells[0][1]];
		}

		//print_game(grid, boundary);

		// collect live neighbor counts, use periodic boundaries
		for (size_t row_i = 0; row_i < grid.size(); row_i++)
		for (size_t cell_i = 0; cell_i < grid[row_i].size(); cell_i++) {

			auto& current_cell = grid[row_i][cell_i];

			for (auto row_offset: {size_t(1), size_t(0), width - 1}) 
			for (auto cell_offset: {size_t(1), size_t(0), height - 1}) {

				if (row_offset == 0 and cell_offset == 0) {
					continue;
				}

				const auto& neighbor
					= grid[
						(row_i + row_offset) % height
					][
						(cell_i + cell_offset) % width
					];

				if (neighbor.is_alive) {
					current_cell.live_neighbors++;
				}
			}
		}

		// set new state
		for (size_t row_i = 0; row_i < grid.size(); row_i++)
		for (size_t cell_i = 0; cell_i < grid[row_i].size(); cell_i++) {

			auto& cell = grid[row_i][cell_i];
			if (cell.live_neighbors == 3) {
				cell.is_alive = true;
			} else if (cell.live_neighbors != 2) {
				cell.is_alive = false;
			}
			cell.live_neighbors = 0;
		}

		if (get_number_of_live_cells(grid) != ((turn % 2 == 0) ? 13 : 15)) {
			std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
				<< "Wrong number of live cells after step: " << turn
				<< std::endl;
			abort();
		}
	}

	//print_game(grid, boundary);

	return EXIT_SUCCESS;
}
