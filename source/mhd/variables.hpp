/*
MHD variables and cell class of PAMHD.

Copyright 2014 Ilja Honkonen
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice, this
  list of conditions and the following disclaimer in the documentation and/or
  other materials provided with the distribution.

* Neither the name of copyright holders nor the names of their contributors
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

#ifndef PAMHD_MHD_VARIABLES_HPP
#define PAMHD_MHD_VARIABLES_HPP

#include "Eigen/Core"

#include "gensimcell.hpp"

namespace pamhd {
namespace mhd {


/*
Variables used by magnetohydrodynamic (MHD) solver
*/

struct Mass_Density {
	using data_type = double;
	static constexpr char name[] = "mass";
};

struct Momentum_Density {
	using data_type = Eigen::Vector3d;
	static constexpr char name[] = "momentum";
};

struct Total_Energy_Density {
	using data_type = double;
	static constexpr char name[] = "energy";
};

struct Magnetic_Field {
	using data_type = Eigen::Vector3d;
	static constexpr char name[] = "magnetic field";
};

struct Velocity {
	using data_type = Eigen::Vector3d;
	static constexpr char name[] = "velocity";
};

struct Pressure {
	using data_type = double;
	static constexpr char name[] = "pressure";
};


//! Set of conservative MHD variables
using MHD_Conservative = gensimcell::Cell<
	Mass_Density,
	Momentum_Density,
	Total_Energy_Density,
	Magnetic_Field
>;

//! Set of primitive MHD variables
using MHD_Primitive = gensimcell::Cell<
	Mass_Density,
	Velocity,
	Pressure,
	Magnetic_Field
>;

//! Represents current state of MHD in simulation cell
struct MHD_State_Conservative {
	using data_type = MHD_Conservative;
};

//! Represents change of MHD in simulation cell for next step
struct MHD_Flux_Conservative {
	using data_type = MHD_Conservative;
};

//! Cell type used by MHD simulations of pamhd
using Cell_T = gensimcell::Cell<
	MHD_State_Conservative,
	MHD_Flux_Conservative
>;


}} // namespaces

#endif // ifndef PAMHD_MHD_VARIABLES_HPP
