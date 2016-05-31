/*
Game of life test for copy boundaries class of PAMHD.

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


With uncommented print_games prints:
............
.00.....0...
.OO.....0,..
........O...
............
............
.00.........
.00.........
...OO.......
...OO.......
............
............

............
.00.........
.OO....00O..
........,...
............
............
.00.........
.0..........
...,O.......
...OO.......
............
............

...

where 0 and . are normal cells, and O and , are copy boundary cells.
*/

#include "array"
#include "cstdlib"
#include "iostream"
#include "limits"
#include "string"
#include "unordered_set"
#include "vector"

#include "rapidjson/document.h"

#include "boundaries/copy_boundaries.hpp"


using namespace std;
using namespace pamhd::boundaries;


struct Is_Alive {
	using data_type = bool;
};

struct Cell {
	typename Is_Alive::data_type is_alive = false;
	int live_neighbors = 0;
};

constexpr size_t
	width = 12,
	height = 12;

using Grid = array<array<Cell, width>, height>; // Grid[row][column]
Grid grid;

using Cell_Id = array<size_t, 2>;


/*
Make unordered_set<array<T, N>> work:
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
Prints given game of life.

Prints non-boundary live cells as 0 and dead as .
and copy boundary live cells as O and dead as ,
*/
template <class Boundary> void print_game(
	const Grid& grid,
	const Boundary& boundary
) {
	std::unordered_set<Cell_Id> boundary_cells;
	for (const auto& item: boundary.copy_sources) {
		boundary_cells.insert(item[0]);
	}

	for (size_t row_i = 0; row_i < grid.size(); row_i++) {
		for (size_t col_i = 0; col_i < grid.size(); col_i++) {
			const Cell_Id cell_id{row_i, col_i};
			const auto& cell_data = grid[row_i][col_i];

			if (boundary_cells.count(cell_id) > 0) {
				if (cell_data.is_alive > 0) {
					cout << "O";
				} else {
					cout << ",";
				}
			} else {
				if (cell_data.is_alive > 0) {
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


int main()
{
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


	Copy_Boundaries<Cell_Id, int, Is_Alive> boundaries;

	/*
	Set source cells for copy boundary cells
	*/

	// block
	boundaries.push_back_source({Cell_Id{2, 1}, Cell_Id{1, 1}});
	boundaries.push_back_source({Cell_Id{2, 2}, Cell_Id{1, 2}});

	// blinker
	boundaries.push_back_source({Cell_Id{3, 8}, Cell_Id{1, 8}});
	boundaries.push_back_source({Cell_Id{2, 9}, Cell_Id{2, 7}});

	// bacon
	boundaries.push_back_source({Cell_Id{8, 3}, Cell_Id{7, 2}});
	boundaries.push_back_source({Cell_Id{8, 4}, Cell_Id{6, 2}});
	boundaries.push_back_source({Cell_Id{9, 3}, Cell_Id{7, 1}});
	boundaries.push_back_source({Cell_Id{9, 4}, Cell_Id{6, 1}});


	// play
	constexpr size_t max_turns = 20;
	for (size_t turn = 0; turn < max_turns; turn++) {

		// copy is_alive into copy boundary cells from their source cell
		for (const auto& item: boundaries.copy_sources) {
			const auto
				&target = item[0],
				&source = item[1];
			grid[target[0]][target[1]] = grid[source[0]][source[1]];
		}

		//print_game(grid, boundaries);

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

	//print_game(grid, boundaries);

	return EXIT_SUCCESS;
}
