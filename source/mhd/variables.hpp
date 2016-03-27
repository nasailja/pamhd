/*
MHD variables and cell class of PAMHD.

Copyright 2014, 2015 Ilja Honkonen
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

#include "Eigen/Core" // must be included before gensimcell.hpp
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

struct Magnetic_Field_Flux {
	using data_type = Eigen::Vector3d;
	static const std::string get_name() { return {"magnetic field flux"}; }
	static const std::string get_option_name() { return {"magnetic-field-flux"}; }
	static const std::string get_option_help() { return {"Flux of magnetic field (T)"}; }
};

//! stores B before divergence removal so B can be restored after failed removal
struct Magnetic_Field_Temp {
	using data_type = Eigen::Vector3d;
	static const std::string get_name() { return {"temporary magnetic field"}; }
	static const std::string get_option_name() { return {"temporary-magnetic-field"}; }
	static const std::string get_option_help() { return {"Temporary value of magnetic field in plasma (T)"}; }
};

//! stores change in B due to resistivity
struct Magnetic_Field_Resistive {
	using data_type = Eigen::Vector3d;
	static const std::string get_name() { return {"resistive magnetic field"}; }
	static const std::string get_option_name() { return {"resistive-magnetic-field"}; }
	static const std::string get_option_help() { return {"Change in magnetic field due to resistivity"}; }
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

// TODO: separate into own file
struct Cell_Type {
	using data_type = int;
	static const std::string get_name() { return {"cell type"}; }
	static const std::string get_option_name() { return {"cell-type"}; }
	static const std::string get_option_help() {
		return {"Cell type (normal, do-not-solve, value / copy boundary)"};
	}
};

struct Value_Boundary_Id {
	using data_type = int;
	static const std::string get_name() { return {"value boundary id"}; }
	static const std::string get_option_name() { return {"value-boundary-id"}; }
	static const std::string get_option_help() {
		return {"Value boundary cell's boundary id"};
	}
};

struct Copy_Source {
	using data_type = unsigned long long int;
	static const std::string get_name() { return {"source of copy cell"}; }
	static const std::string get_option_name() { return {"source-of-copy-cell"}; }
	static const std::string get_option_help() {
		return {"Id of cell that is source of data for copy boundary cell."};
	}
};

struct Magnetic_Field_Divergence {
	using data_type = double;
	static const std::string get_name() { return {"magnetic field divergence"}; }
	static const std::string get_option_name() { return {"magnetic-field-divergence"}; }
	static const std::string get_option_help() { return {"Divergence of plasma magnetic field (T/m)"}; }
};

struct Scalar_Potential_Gradient {
	using data_type = Eigen::Vector3d;
	static const std::string get_name() { return {"scalar potential gradient"}; }
	static const std::string get_option_name() { return {"scalar-potential-gradient"}; }
	static const std::string get_option_help() { return {"Gradient of scalar potential from Poisson's equation"}; }
};

struct Electric_Current_Density {
	using data_type = Eigen::Vector3d;
	static const std::string get_name() { return {"current density"}; }
	static const std::string get_option_name() { return {"current-density"}; }
	static const std::string get_option_help() { return {"Density of electric current"}; }
};

//! Electrical resistivity of plasma (n in E = -VxB + nJ + ...)
struct Resistivity {
	using data_type = double;
	static const std::string get_name() { return {"electrical resistivity"}; }
	static const std::string get_option_name() { return {"electrical-resistivity"}; }
	static const std::string get_option_help() { return {"Electrical resistivity"}; }
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
	static const std::string get_name() { return {"conservative MHD variables"}; }
	static const std::string get_option_name() { return {"MHD-conservative"}; }
	static const std::string get_option_help() {
		return {"Conservative MHD variables"};
	}
};

//! Represents change of MHD in simulation cell for next step
struct MHD_Flux_Conservative {
	using data_type = MHD_Conservative;
	static const std::string get_name() { return {"conservative MHD flux variables"}; }
	static const std::string get_option_name() { return {"MHD-flux-conservative"}; }
	static const std::string get_option_help() {
		return {"Conservative MHD flux variables"};
	}
};

// cell type for MHD test program
using Cell = gensimcell::Cell<
	gensimcell::Optional_Transfer,
	MHD_State_Conservative,
	Electric_Current_Density,
	Cell_Type,
	Copy_Source,
	Value_Boundary_Id,
	MPI_Rank,
	Resistivity,
	Magnetic_Field_Resistive,
	Magnetic_Field_Temp,
	Magnetic_Field_Divergence,
	Scalar_Potential_Gradient,
	MHD_Flux_Conservative
>;


// hydrodynamic type and variables for multifluid test program
using HD_Conservative = gensimcell::Cell<
	gensimcell::Always_Transfer,
	pamhd::mhd::Mass_Density,
	pamhd::mhd::Momentum_Density,
	pamhd::mhd::Total_Energy_Density
>;
struct HD1_State { using data_type = HD_Conservative; };
struct HD2_State { using data_type = HD_Conservative; };
struct HD1_Flux { using data_type = HD_Conservative; };
struct HD2_Flux { using data_type = HD_Conservative; };

// separate variables for fluids and fields in boundary code
struct Cell_Type_Fluid2 { using data_type = int; };
struct Cell_Type_Field { using data_type = int; };
struct Copy_Source_Fluid2 { using data_type = unsigned long long int; };
struct Copy_Source_Field { using data_type = unsigned long long int; };
struct Value_Boundary_Id_Fluid2 { using data_type = int; };
struct Value_Boundary_Id_Field { using data_type = int; };

// cell type for multipopulation MHD test program
using Cell2 = gensimcell::Cell<
	gensimcell::Optional_Transfer,
	pamhd::mhd::HD1_State,
	pamhd::mhd::HD2_State,
	pamhd::mhd::Magnetic_Field,
	pamhd::mhd::Electric_Current_Density,
	pamhd::mhd::Cell_Type,
	pamhd::mhd::Copy_Source,
	pamhd::mhd::Value_Boundary_Id,
	pamhd::mhd::Cell_Type_Fluid2,
	pamhd::mhd::Copy_Source_Fluid2,
	pamhd::mhd::Value_Boundary_Id_Fluid2,
	pamhd::mhd::Cell_Type_Field,
	pamhd::mhd::Value_Boundary_Id_Field,
	pamhd::mhd::Copy_Source_Field,
	pamhd::mhd::MPI_Rank,
	pamhd::mhd::Resistivity,
	pamhd::mhd::Magnetic_Field_Resistive,
	pamhd::mhd::Magnetic_Field_Temp,
	pamhd::mhd::Magnetic_Field_Divergence,
	pamhd::mhd::Scalar_Potential_Gradient,
	pamhd::mhd::HD1_Flux,
	pamhd::mhd::HD2_Flux,
	pamhd::mhd::Magnetic_Field_Flux
>;


}} // namespaces

#endif // ifndef PAMHD_MHD_VARIABLES_HPP
