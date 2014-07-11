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
	static const std::string get_name() { return {"mass density"}; }
	static const std::string get_option_name() { return {"mass-density"}; }
	static const std::string get_option_help() { return {"Plasma mass density (kg / m^3)"}; }
};

struct Momentum_Density {
	using data_type = Eigen::Vector3d;
	static const std::string get_name() { return {"momentum density"}; }
	static const std::string get_option_name() { return {"momentum-density"}; }
	static const std::string get_option_help() { return {"Plasma momentum density (kg / m^3 * m / s)"}; }
};

struct Total_Energy_Density {
	using data_type = double;
	static const std::string get_name() { return {"total energy density"}; }
	static const std::string get_option_name() { return {"total-energy-density"}; }
	static const std::string get_option_help() { return {"Plasma total energy density (J / m^3)"}; }
};

struct Magnetic_Field {
	using data_type = Eigen::Vector3d;
	static const std::string get_name() { return {"magnetic field"}; }
	static const std::string get_option_name() { return {"magnetic-field"}; }
	static const std::string get_option_help() { return {"Plasma magnetic field (T)"}; }
};

struct Velocity {
	using data_type = Eigen::Vector3d;
	static const std::string get_name() { return {"velocity"}; }
	static const std::string get_option_name() { return {"velocity"}; }
	static const std::string get_option_help() { return {"Plasma velocity (m / s)"}; }
};

struct Pressure {
	using data_type = double;
	static const std::string get_name() { return {"pressure"}; }
	static const std::string get_option_name() { return {"pressure"}; }
	static const std::string get_option_help() { return {"Plasma thermal pressure (Pa)"}; }
};

struct Number_Density {
	using data_type = double;
	static const std::string get_name() { return {"number density"}; }
	static const std::string get_option_name() { return {"number-density"}; }
	static const std::string get_option_help() { return {"Plasma number density (protons / m^3)"}; }
};

struct MPI_Rank {
	using data_type = int;
	static const std::string get_name() { return {"MPI rank"}; }
	static const std::string get_option_name() { return {"mpi-rank"}; }
	static const std::string get_option_help() { return {"Owner (MPI process) of cell"}; }
};

struct Cell_Type {
	using data_type = int;
	static const std::string get_name() { return {"cell type"}; }
	static const std::string get_option_name() { return {"cell-type"}; }
	static const std::string get_option_help() { return {"Owner (MPI process) of cell"}; }
};

//! Set of conservative MHD variables
using MHD_Conservative = gensimcell::Cell<
	gensimcell::Always_Transfer,
	Mass_Density,
	Momentum_Density,
	Total_Energy_Density,
	Magnetic_Field
>;

//! Set of primitive MHD variables
using MHD_Primitive = gensimcell::Cell<
	gensimcell::Never_Transfer,
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

using Cell = gensimcell::Cell<
	gensimcell::Optional_Transfer,
	MHD_State_Conservative,
	Cell_Type,
	MPI_Rank,
	MHD_Flux_Conservative
>;

}} // namespaces

#endif // ifndef PAMHD_MHD_VARIABLES_HPP
