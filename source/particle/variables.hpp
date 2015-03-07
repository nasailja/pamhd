/*
Particle variables and cell class of PAMHD.

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

#ifndef PAMHD_PARTICLE_VARIABLES_HPP
#define PAMHD_PARTICLE_VARIABLES_HPP


#include "array"
#include "vector"

#include "mpi.h" // must be included before gensimcell.hpp
#include "Eigen/Core" // must be included before gensimcell.hpp
#include "gensimcell.hpp"


namespace pamhd {
namespace particle {


/*
Variables used by particle solver
*/

struct Position {
	using data_type = Eigen::Vector3d;
	static const std::string get_name() { return {"position"}; }
	static const std::string get_option_name() { return {"position"}; }
	static const std::string get_option_help() { return {"Particle position"}; }
};

struct Velocity {
	using data_type = Eigen::Vector3d;
	static const std::string get_name() { return {"velocity"}; }
	static const std::string get_option_name() { return {"velocity"}; }
	static const std::string get_option_help() { return {"Particle velocity"}; }
};

//! Test particles have 0 mass.
struct Mass {
	using data_type = double;
	static const std::string get_name() { return {"mass"}; }
	static const std::string get_option_name() { return {"mass"}; }
	static const std::string get_option_help() { return {"Particle mass)"}; }
};

//! Represents the ratio of particle mass and particle charge
struct Charge_Mass_Ratio {
	using data_type = double;
	static const std::string get_name() { return {"charge mass ratio"}; }
	static const std::string get_option_name() { return {"charge-mass-ratio"}; }
	static const std::string get_option_help() { return {"Particle charge to mass ratio"}; }
};


/*!
Template for a particle type that stores basic particle data
and given optional variables.
*/
template<
	class... Extra_Variables
> using Particle_T = gensimcell::Cell<
	gensimcell::Always_Transfer,
	Position,
	Velocity,
	Mass,
	Charge_Mass_Ratio,
	Extra_Variables...
>;


//! Represents particles not moving between cells
struct Particles_Internal {
	using data_type = std::vector<Particle_T<>>;
};

//! Represents number of particles not moving between cells.
struct Nr_Particles_Internal {
	using data_type = unsigned long long int;
};


//! Represents destination cell of particles moving between processes
struct Destination_Cell {
	using data_type = unsigned long long int;
};

//! Represents particles moving between cells
struct Particles_External {
	using data_type = std::vector<Particle_T<Destination_Cell>>;
};

//! Represents number of particles moving between cells.
struct Nr_Particles_External {
	using data_type = unsigned long long int;
};


struct Electric_Field {
	using data_type = Eigen::Vector3d;
	static const std::string get_name() { return {"electric field"}; }
	static const std::string get_option_name() { return {"electric-field"}; }
	static const std::string get_option_help() { return {"Electric field used for propagating particles"}; }
};

struct Magnetic_Field {
	using data_type = Eigen::Vector3d;
	static const std::string get_name() { return {"magnetic field"}; }
	static const std::string get_option_name() { return {"magnetic-field"}; }
	static const std::string get_option_help() { return {"Magnetic field used for propagating particles"}; }
};


/*!
Cell type for particle model of PAMHD.

Fixed size variables are first so that they
are saved at a known position in the file by dccrg.
gensimcell puts the variables in an MPI datatype
in the same order as they are given here, which
dccrg uses to save the file.
*/
using Cell = gensimcell::Cell<
	gensimcell::Optional_Transfer,
	Electric_Field,
	Magnetic_Field,
	Nr_Particles_Internal,
	Nr_Particles_External,
	Particles_Internal,
	Particles_External
>;


}} // namespaces

#endif // ifndef PAMHD_PARTICLE_VARIABLES_HPP
