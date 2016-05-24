/*
Game of life test for initial condition class of PAMHD, expression version.

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


Uses one instance of initial condition class which fills entire game board
with data parsed from a mathematical expression.
*/

#include "array"
#include "cstdlib"
#include "iostream"
#include "set"
#include "string"
#include "vector"

#include "boundaries/initial_condition.hpp"
#include "rapidjson/document.h"


using namespace std;
using namespace pamhd::boundaries;


struct Is_Alive {
	using data_type = bool;
	static const std::string get_name(){ return {"is alive"}; }
};

struct Cell {
	typename Is_Alive::data_type is_alive = 0;
	int live_neighbors = 0;
};


constexpr size_t
	width = 8,
	height = 8;

/*!
Grid[y position][x position]

On screen x and y position and index increases
from left to right and top to bottom respectively
*/
using Grid = array<array<Cell, width>, height>;
Grid grid;


/*!
cell[0] == x , cell[1] == y
*/
array<double, 3> get_cell_center(
	const Grid& grid,
	const array<size_t, 2>& cell
) {
	if (
		cell[1] >= grid.size()
		or cell[0] >= grid[cell[1]].size()
	) {
		return {
			std::numeric_limits<double>::quiet_NaN(),
			std::numeric_limits<double>::quiet_NaN(),
			std::numeric_limits<double>::quiet_NaN()
		};
	}

	return {
		-1.0 + (0.5 + cell[0]) * 2.0 / grid[cell[1]].size(),
		-1.0 + (0.5 + cell[1]) * 2.0 / grid.size(),
		0
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


int main()
{
	const char json[] = "{"
		"\"geometry_id\": 2147483647,"
		"\"value\": \"abs(y - 0.1) < 0.1 ? 1 : 0\""
	"}";

	rapidjson::Document document;
	document.Parse(json);
	if (document.HasParseError()) {
		std::cerr << "Couldn't parse json data: " << json << std::endl;
		return EXIT_FAILURE;
	}


	Initial_Condition<int, Is_Alive> initial_condition;
	initial_condition.set(document.GetObject());
	for (size_t row_i = 0; row_i < grid.size(); row_i++) {
	for (size_t cell_i = 0; cell_i < grid[row_i].size(); cell_i++) {
		const auto center = get_cell_center(grid, {cell_i, row_i});

		grid[row_i][cell_i].is_alive
			= initial_condition.get_data(
				0,
				center[0], center[1], center[2],
				0, 0, 0
			);
	}}

	// play
	constexpr size_t max_turns = 5;
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

		const size_t live_cells = get_number_of_live_cells(grid);
		switch (turn) {
		case 0:
			if (live_cells != 24) {
				std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
					<< "Wrong number of live cells after step " << turn
					<< ": " << live_cells << ", should be 24"
					<< std::endl;
				abort();
			}
			break;
		case 1:
			if (live_cells != 16) {
				std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
					<< "Wrong number of live cells after step " << turn
					<< ": " << live_cells << ", should be 16"
					<< std::endl;
				abort();
			}
			break;
		case 2:
			if (live_cells != 48) {
				std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
					<< "Wrong number of live cells after step " << turn
					<< ": " << live_cells << ", should be 48"
					<< std::endl;
				abort();
			}
			break;
		default:
			if (live_cells != 0) {
				std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
					<< "Wrong number of live cells after step " << turn
					<< ": " << live_cells << ", should be 0"
					<< std::endl;
				abort();
			}
			break;
		}
	}

	//print_game(grid);

	return EXIT_SUCCESS;
}
