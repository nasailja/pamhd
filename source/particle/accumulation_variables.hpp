/*
Variables used in particle accumulation part of PAMHD.

Copyright 2015 Ilja Honkonen
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

#ifndef PAMHD_PARTICLE_ACCUMULATION_VARIABLES_HPP
#define PAMHD_PARTICLE_ACCUMULATION_VARIABLES_HPP


#include "vector"

#include "mpi.h" // must be included before gensimcell.hpp
#include "gensimcell.hpp"


namespace pamhd {
namespace particle {


//! average number of particles in a cell
struct Count {
	using data_type = double;
};

//! cell id when accumulating to cells of other processes
struct Target {
	using data_type = unsigned long long int;
};

//! particle data accumulated to one cell of another process
using Accumulated_To_Cell = gensimcell::Cell<
	gensimcell::Always_Transfer,
	Target,
	Count
>;

//! particle data accumulated to remote neighbors of local cell
struct Accumulated_To_Cells {
	using data_type = std::vector<Accumulated_To_Cell>;
};

//! number of cell into which particle data has been accumulated
struct Nr_Accumulated_To_Cells {
	using data_type = unsigned long long int;
};

//! Example cell type for the accumulation model without particles
using Accumulation_Cell = gensimcell::Cell<
	gensimcell::Optional_Transfer,
	// values accumulated to cell itself
	Count,
	// values accumulated to cell's neighbors
	Nr_Accumulated_To_Cells,
	Accumulated_To_Cells
>;


}} // namespaces

#endif // ifndef PAMHD_PARTICLE_ACCUMULATION_VARIABLES_HPP
