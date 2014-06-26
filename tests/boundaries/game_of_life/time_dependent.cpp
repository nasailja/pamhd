/*
Game of life test for time-dependent boundary class of PAMHD.

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

#include "boost/lexical_cast.hpp"
#include "boost/program_options.hpp"

#include "boundaries/time_dependent_boundary.hpp"


using namespace std;
using namespace pamhd::boundaries;


struct Is_Alive {
	using data_type = int;

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
		return {"Whether a cell is alive (> 0) or not (<= 0)"};
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
	int is_alive = 0, live_neighbors = 0;
};


constexpr size_t
	width = 8,
	height = 8;

using Grid = array<array<Cell, width>, height>;
Grid grid;


array<double, 2> get_cell_start(
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
		1.0 - (1 + cell[1]) * 2.0 / grid.size()
	};
}


array<double, 2> get_cell_end(
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
		1.0 - cell[1] * 2.0 / grid.size()
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


int main(int argc, char* argv[])
{
	std::string boundary_file_name;

	boost::program_options::options_description
		options(
			"Usage: program_name [options], where options are:"
		),
		boundary_options(
			"Options for time-dependent boundary condition"
		);

	options.add_options()
		("help", "Print options and their descriptions")
		("boundary-help", "Print help for options of boundary condition")
		("boundary-file",
			boost::program_options::value<std::string>(&boundary_file_name)
				->default_value("tests/boundaries/game_of_life/gol_boundary.cfg"),
			"Read time-dependent \"boundary\" conditions from file arg");


	// create gliders from 3 time-dependent boundaries
	Time_Dependent_Boundary<
		std::array<size_t, 2>,
		std::array<double, 2>,
		size_t,
		Is_Alive
	> bdy1, bdy2, bdy3;
	bdy1.add_options(boundary_options, "bdy1.");
	bdy2.add_options(boundary_options, "bdy2.");
	bdy3.add_options(boundary_options, "bdy3.");


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
	if (option_variables.count("boundary-help") > 0) {
		std::cout << boundary_options << std::endl;
		return EXIT_SUCCESS;
	}

	boost::program_options::store(
		boost::program_options::parse_config_file<char>(
			boundary_file_name.c_str(),
			boundary_options
		),
		option_variables
	);
	boost::program_options::notify(option_variables);

	// set default initial state
	for (size_t row_i = 0; row_i < grid.size(); row_i++) {
	for (size_t cell_i = 0; cell_i < grid[row_i].size(); cell_i++) {
		grid[row_i][cell_i].is_alive = false;
	}}

	// play
	constexpr size_t max_turns = 27;
	for (size_t turn = 0; turn < max_turns; turn++) {

		// classify and set time-dependent boundary values twice during the game
		if (turn == 0 or turn == 1 or turn == 16) {
			for (size_t row_i = 0; row_i < grid.size(); row_i++) {
			for (size_t cell_i = 0; cell_i < grid[row_i].size(); cell_i++) {

				const std::array<size_t, 2> cell{row_i, cell_i};

				bdy1.add_cell(
					cell,
					get_cell_start(grid, cell),
					get_cell_end(grid, cell),
					turn
				);
				bdy2.add_cell(
					cell,
					get_cell_start(grid, cell),
					get_cell_end(grid, cell),
					turn
				);
				bdy3.add_cell(
					cell,
					get_cell_start(grid, cell),
					get_cell_end(grid, cell),
					turn
				);
			}}

			for (const auto cell: bdy1.get_boundary_cells(turn)) {
				const size_t row_i = cell[1], cell_i = cell[0];
				grid[row_i][cell_i].is_alive
					= bdy1.get_boundary_data(Is_Alive(), turn);
			}

			for (const auto cell: bdy2.get_boundary_cells(turn)) {
				const size_t row_i = cell[1], cell_i = cell[0];
				grid[row_i][cell_i].is_alive
					= bdy2.get_boundary_data(Is_Alive(), turn);
			}

			for (const auto cell: bdy3.get_boundary_cells(turn)) {
				const size_t row_i = cell[1], cell_i = cell[0];
				grid[row_i][cell_i].is_alive
					= bdy3.get_boundary_data(Is_Alive(), turn);
			}
		}

		const size_t live_cells = get_number_of_live_cells(grid);
		if (turn == 0 and live_cells != 0) {
			std::cerr <<  __FILE__ << "(" << __LINE__<< ") "
				<< "Wrong number of live cells at turn " << turn << ": " << live_cells
				<< std::endl;
			abort();
		}
		if (turn > 0 and turn < 16 and live_cells != 5) {
			std::cerr <<  __FILE__ << "(" << __LINE__<< ") "
				<< "Wrong number of live cells at turn " << turn << ": " << live_cells
				<< std::endl;
			abort();
		}
		if (turn == 21 and live_cells != 9) {
			std::cerr <<  __FILE__ << "(" << __LINE__<< ") "
				<< "Wrong number of live cells at turn " << turn << ": " << live_cells
				<< std::endl;
			abort();
		}
		if (turn == 25 and live_cells != 3) {
			std::cerr <<  __FILE__ << "(" << __LINE__<< ") "
				<< "Wrong number of live cells at turn " << turn << ": " << live_cells
				<< std::endl;
			abort();
		}
		if (turn == 26 and live_cells != 0) {
			std::cerr <<  __FILE__ << "(" << __LINE__<< ") "
				<< "Wrong number of live cells at turn " << turn << ": " << live_cells
				<< std::endl;
			abort();
		}

		//print_game(grid);

		// collect live neighbor counts, use periodic boundaries
		for (size_t row_i = 0; row_i < grid.size(); row_i++)
		for (size_t cell_i = 0; cell_i < grid[row_i].size(); cell_i++) {

			auto& current_cell = grid[row_i][cell_i];

			for (int row_offset: {1, 0, -1}) 
			for (int cell_offset: {1, 0, -1}) {

				if (row_offset == 0 and cell_offset == 0) {
					continue;
				}
				if (row_i == 0 and row_offset < 0) {
					continue;
				}
				if (row_i == grid.size() - 1 and row_offset > 0) {
					continue;
				}
				if (cell_i == 0 and cell_offset < 0) {
					continue;
				}
				if (cell_i == grid[row_i].size() - 1 and cell_offset > 0) {
					continue;
				}

				const auto& neighbor = grid[row_i + row_offset][cell_i + cell_offset];

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
	}

	//print_game(grid);

	return EXIT_SUCCESS;
}
