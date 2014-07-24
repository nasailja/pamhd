/*
Game of life test for initial condition class of PAMHD.

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

#include "array"
#include "cstdlib"
#include "iostream"
#include "set"
#include "string"
#include "vector"

#include "boost/program_options.hpp"

#include "boundaries/initial_condition.hpp"


using namespace std;
using namespace pamhd::boundaries;


struct Is_Alive {
	using data_type = bool;

	static const std::string get_name()
	{
		return {"is alive"};
	}
	static const std::string get_option_name()
	{
		return {"is-alive"};
	}
	static const std::string get_option_help()
	{
		return {"Whether a cell is alive (true) or not (false)"};
	}
};

struct Live_Neighbors {
	using data_type = int;

	static const std::string get_name()
	{
		return {"live neighbors"};
	}
	static const std::string get_option_name()
	{
		return {"live-neighbors"};
	}
	static const std::string get_option_help()
	{
		return {"Number of cell's live neighbors"};
	}
};


struct Cell {
	bool is_alive = false;
	int live_neighbors = 0;
};


constexpr size_t
	width = 8,
	height = 8;

using Grid = array<array<Cell, width>, height>;
Grid grid;


array<double, 3> get_cell_start(
	const Grid& grid,
	const array<size_t, 2>& cell
) {
	if (
		cell[1] >= grid.size()
		or cell[0] >= grid[cell[1]].size()
	) {
		return {
			std::numeric_limits<double>::quiet_NaN(),
			std::numeric_limits<double>::quiet_NaN()
		};
	}

	return {
		-1.0 + cell[0] * 2.0 / grid[cell[1]].size(),
		1.0 - (1 + cell[1]) * 2.0 / grid.size(),
		-1
	};
}


array<double, 3> get_cell_end(
	const Grid& grid,
	const array<size_t, 2>& cell
) {
	if (
		cell[1] >= grid.size()
		or cell[0] >= grid[cell[1]].size()
	) {
		return {
			std::numeric_limits<double>::quiet_NaN(),
			std::numeric_limits<double>::quiet_NaN()
		};
	}

	return {
		-1.0 + (1 + cell[0]) * 2.0 / grid[cell[1]].size(),
		1.0 - cell[1] * 2.0 / grid.size(),
		1
	};
}


/*
Prints given game of life to standard output,
using 0 for live cells and . for dead cells.
*/
void print_game(const Grid& grid)
{
	for (const auto& row: grid) {
		for (const auto& cell: row) {
			if (cell.is_alive > 0) {
				cout << "0";
			} else {
				cout << ".";
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
			if (cell.is_alive > 0) {
				live_cells++;
			}
		}
	}
	return live_cells;
}

#define MAKE_VEC(x, y, z) std::array<double, 3>{x, y, z}

int main(int argc, char* argv[])
{
	Initial_Condition<
		std::array<size_t, 2>,
		double,
		std::array<double, 3>,
		Is_Alive
	> initial_condition;

	// set defaults
	initial_condition.default_data.set_expression(
		Is_Alive(),
		"false"
	);

	boost::optional<size_t> result;

	result = initial_condition.add_box(
		MAKE_VEC(-0.5, 0.5, -1),
		MAKE_VEC(-0.25, 0.75, 1)
	);
	if (not result) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
			<< "Couldn't add boundary box."
			<< std::endl;
		abort();
	}

	result = initial_condition.add_box(
		MAKE_VEC(-0.25, 0.25, -1),
		MAKE_VEC(0, 0.5, 1)
	);
	if (not result) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
			<< "Couldn't add boundary box."
			<< std::endl;
		abort();
	}

	result = initial_condition.add_box(
		MAKE_VEC(-0.25, 0, -1),
		MAKE_VEC(0, 0.25, 1)
	);
	if (not result) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
			<< "Couldn't add boundary box."
			<< std::endl;
		abort();
	}

	result = initial_condition.add_box(
		MAKE_VEC(-0.5, 0, -1),
		MAKE_VEC(-0.25, 0.25, 1)
	);
	if (not result) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
			<< "Couldn't add boundary box."
			<< std::endl;
		abort();
	}

	result = initial_condition.add_box(
		MAKE_VEC(-0.75, 0, -1),
		MAKE_VEC(-0.5, 0.25, 1)
	);
	if (not result) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
			<< "Couldn't add boundary box."
			<< std::endl;
		abort();
	}

	for (size_t i = 0; i < initial_condition.get_number_of_boundaries(); i++) {
		initial_condition.set_expression(Is_Alive(), i, "true");
	}


	std::string initial_condition_file_name("");

	boost::program_options::options_description
		options(
			"Usage: program_name [options], where options are"
		),
		initial_condition_options(
			"Options for initial condition (must be read from file)"
		);

	options.add_options()
		("help", "Print options and their descriptions")
		("initial-help", "Print help for options of initial condition")
		("initial-file",
			boost::program_options::value<std::string>(&initial_condition_file_name)
				->default_value(initial_condition_file_name),
			"Read initial condition from file arg (do not read if empty string)");
	initial_condition.add_initialization_options("", options);

	boost::program_options::variables_map option_variables;
	boost::program_options::store(
		boost::program_options::parse_command_line(
			argc,
			argv,
			options
		),
		option_variables
	);
	boost::program_options::notify(option_variables);

	if (option_variables.count("help") > 0) {
		std::cout << options << std::endl;
		return EXIT_SUCCESS;
	}

	initial_condition.add_options("", initial_condition_options);

	if (option_variables.count("initial-help") > 0) {
		std::cout << initial_condition_options << std::endl;
		return EXIT_SUCCESS;
	}

	if (initial_condition_file_name != "") {
		boost::program_options::store(
			boost::program_options::parse_config_file<char>(
				initial_condition_file_name.c_str(),
				initial_condition_options
			),
			option_variables
		);
		boost::program_options::notify(option_variables);
	}

	// set default initial state
	for (size_t row_i = 0; row_i < grid.size(); row_i++)
	for (size_t cell_i = 0; cell_i < grid[row_i].size(); cell_i++) {
		grid[row_i][cell_i].is_alive
			= initial_condition.default_data.get_data(Is_Alive(), {0, 0, 0}, 0);
	}

	// classify cells into regions given as initial condition
	for (size_t row_i = 0; row_i < grid.size(); row_i++)
	for (size_t cell_i = 0; cell_i < grid[row_i].size(); cell_i++) {

		const std::array<size_t, 2> cell{row_i, cell_i};

		const size_t result = initial_condition.add_cell(
			cell,
			get_cell_start(grid, cell),
			get_cell_end(grid, cell)
		);

		if (result > 1) {
			std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
				<< "Cell " << row_i << ", " << cell_i
				<< " assigned to more than one region."
				<< std::endl;
			abort();
		}
	}

	// go through all boundaries and set their cells' data to given value
	for (size_t
		bdy_i = 0;
		bdy_i < initial_condition.get_number_of_boundaries();
		bdy_i++
	) {
		const auto& boundary_cells = initial_condition.get_cells(bdy_i);

		for (const auto& cell: boundary_cells) {
			const size_t
				row_i = cell[1],
				cell_i = cell[0];

			grid[row_i][cell_i].is_alive = initial_condition.get_data(
				Is_Alive(),
				bdy_i,
				{0, 0, 0},
				0
			);
		}
	}

	const size_t live_cells = get_number_of_live_cells(grid);

	// play
	constexpr size_t max_turns = 32;
	for (size_t turn = 0; turn < max_turns; turn++) {

		//print_game(grid);

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

		const size_t new_live_cells = get_number_of_live_cells(grid);
		if (new_live_cells != live_cells) {
			std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
				<< "Wrong number of live cells after step " << turn
				<< ": " << new_live_cells
				<< ", should be " << live_cells
				<< std::endl;
			abort();
		}
	}

	//print_game(grid);

	return EXIT_SUCCESS;
}
