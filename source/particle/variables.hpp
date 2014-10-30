/*
Particle variables and cell class of PAMHD.

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

#ifndef PAMHD_PARTICLE_VARIABLES_HPP
#define PAMHD_PARTICLE_VARIABLES_HPP


#include "array"
#include "stdexcept"
#include "tuple"
#include "vector"

//#include "mpi.h" // must be included before gensimcell.hpp
#include "Eigen/Core" // must be included before gensimcell.hpp
#include "gensimcell.hpp" // mpi.h (if used) must be included before this


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
Template for a particle type that stores data
relevant to physics and given optional variables.
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


/*!
Basic particle type that only stores data relevant to physics.
*/
using Particle = Particle_T<>;


/*!
Stores a collection of particles that stay on the current process.
*/
struct Particles_Storage
{
	std::vector<Particle> particles;


	//! Returns a const reference to particle data
	const decltype(particles)& operator()() const
	{
		return this->particles;
	}

	//! Returns a reference to particle data
	decltype(particles)& operator()()
	{
		return const_cast<decltype(particles)&>(
			static_cast<const Particles_Storage&>(*this)()
		);
	}


	#if defined(MPI_VERSION) && (MPI_VERSION >= 2)
	//! Returns 0 bytes.
	std::tuple<void*, int, MPI_Datatype> get_mpi_datatype() const
	{
		return std::make_tuple((void*) this->particles.data(), 0, MPI_BYTE);
	}
	#endif
};

//! Represents particles not moving between processes
struct Particles_Int {
	using data_type = Particles_Storage;
};

//! Represents number of local particles, used for file I/O.
struct Nr_Particles_Int {
	using data_type = unsigned long long int;
	static const std::string get_name() { return {"number of particles"}; }
	static const std::string get_option_name() { return {"number-of-particles"}; }
	static const std::string get_option_help() { return {"Number of particles "}; }
};

//! Represents destination process of particles moving between processes
struct Destination_Process {
	using data_type = int;
	static const std::string get_name() { return {"destination process"}; }
	static const std::string get_option_name() { return {"destination-process"}; }
	static const std::string get_option_help() { return {"Process to whom particle should be transferred"}; }
};


/*!
Particle type for transferring particles between processes.

Particles that move from a cell to another cell owned by a
different process are changed to this type before updating
particle data between processes.
*/
using Particle_Transfer = Particle_T<Destination_Process>;


/*!
Stores a collection of particles that move between cells/processes.

Provides their data via get_mpi_datatype().
*/
class Particles_Transfer_Storage
{
public:
	std::vector<Particle_Transfer> particles;


	//! Returns a const reference to particle data
	const decltype(particles)& operator()() const
	{
		return this->particles;
	}

	//! Returns a reference to particle data
	decltype(particles)& operator()()
	{
		return const_cast<decltype(particles)&>(
			static_cast<const Particles_Transfer_Storage&>(*this)()
		);
	}


	#if defined(MPI_VERSION) && (MPI_VERSION >= 2)
	/*!
	Returns MPI Datatype corresponding to data of stored particles.

	Takes data of each particle from its get_mpi_datatype() method.
	*/
	std::tuple<
		void*,
		int,
		MPI_Datatype
	> get_mpi_datatype() throw(std::length_error)
	{
		std::vector<
			std::tuple<
				void*,
				int,
				MPI_Datatype
			>
		> transfer_info;
		transfer_info.reserve(this->particles.size());

		for (auto& particle: this->particles) {
			transfer_info.push_back(particle.get_mpi_datatype());
		}

		const auto nr_of_particles = transfer_info.size();
		if (nr_of_particles > std::numeric_limits<int>::max()) {
			throw std::domain_error(
				__func__ + std::string(": Number of particles exceeds int")
			);
		}

		std::vector<int> counts(nr_of_particles, 0);
		for (size_t i = 0; i < nr_of_particles; i++) {
			counts[i] = std::get<1>(transfer_info[i]);
		}

		std::vector<MPI_Aint> displacements(nr_of_particles, 0);
		for (size_t i = 0; i < nr_of_particles; i++) {
			displacements[i]
				= static_cast<char*>(std::get<0>(transfer_info[i]))
				- static_cast<char*>(std::get<0>(transfer_info[0]));
		}

		std::vector<MPI_Datatype> datatypes(nr_of_particles, MPI_DATATYPE_NULL);
		for (size_t i = 0; i < nr_of_particles; i++) {
			datatypes[i] = std::get<2>(transfer_info[i]);
		}

		MPI_Datatype final_datatype = MPI_DATATYPE_NULL;
		if (
			MPI_Type_create_struct(
				int(nr_of_particles),
				counts.data(),
				displacements.data(),
				datatypes.data(),
				&final_datatype
			) != MPI_SUCCESS
		) {
			return std::make_tuple(nullptr, -1, MPI_DATATYPE_NULL);
		}

		// free particle datatypes
		for (auto& datatype: datatypes) {
			if (datatype == MPI_DATATYPE_NULL) {
				continue;
			}
			int combiner = -1, tmp1 = -1, tmp2 = -1, tmp3 = -1;
			MPI_Type_get_envelope(datatype, &tmp1, &tmp2, &tmp3, &combiner);
			if (combiner != MPI_COMBINER_NAMED) {
				MPI_Type_free(&datatype);
			}
		}

		return std::make_tuple(
			std::get<0>(transfer_info[0]),
			1,
			final_datatype
		);
	}
	#endif
};


/*!
Records the number of particles to be transferred between processes.
*/
struct Nr_Particles_Ext {
	using data_type = unsigned long long int;
	static const std::string get_name() { return {"number of particles"}; }
	static const std::string get_option_name() { return {"number-of-particles"}; }
	static const std::string get_option_help() { return {"Number of particles "}; }
};

//! Variable used to refer to particles moving between cells/processes
struct Particles_Ext {
	using data_type = Particles_Transfer_Storage;
};


/*!


Fixed size variables are first so that they
are saved at a known position in the file by dccrg.
gensimcell puts the variables in an MPI datatype
in the same order as they are given here, which
dccrg uses to save the file.
*/
using Cell = gensimcell::Cell<
	gensimcell::Optional_Transfer,
	Nr_Particles_Int,
	Nr_Particles_Ext,
	Particles_Int,
	Particles_Ext
>;


}} // namespaces

#endif // ifndef PAMHD_PARTICLE_VARIABLES_HPP
