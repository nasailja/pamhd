/*
Common MHD functions of PAMHD.

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

#ifndef PAMHD_MHD_COMMON_HPP
#define PAMHD_MHD_COMMON_HPP


#include "cmath"
#include "stdexcept"
#include "string"

#include "boost/lexical_cast.hpp"


namespace pamhd {
namespace mhd {


/*!
Returns pressure.

Throws std::domain_error if given a state with non-positive mass density.

Returns negative pressure if total energy in given state
is smaller than kinetic + magnetic energies.

Mom and mag must be compatible with std::array<double, 3>.
*/
template <
	class Momentum_Density,
	class Magnetic_Field
> double get_pressure(
	const double& mass_density,
	const Momentum_Density& mom,
	const double& total_energy_density,
	const Magnetic_Field& mag,
	const double adiabatic_index,
	const double vacuum_permeability
) {
	if (mass_density <= 0) {
		throw std::domain_error(
			std::string("Non-positive mass density given to ")
			+ __func__
			+ std::string(": ")
			+ boost::lexical_cast<std::string>(mass_density)
		);
	}

	const auto
		kinetic_energy
			= 0.5 / mass_density
			* (
				mom[0] * mom[0]
				+ mom[1] * mom[1]
				+ mom[2] * mom[2]
			),
		magnetic_energy
			= 0.5 / vacuum_permeability
			* (
				mag[0] * mag[0]
				+ mag[1] * mag[1]
				+ mag[2] * mag[2]
			);

	return
		(total_energy_density - kinetic_energy - magnetic_energy)
		* (adiabatic_index - 1);
}


template<
	class Momentum,
	class Mass
> Momentum get_velocity(
	Momentum& mom,
	Mass& mas
) {
	if (mas > 0) {
		return mom / mas;
	} else {
		return {0, 0, 0};
	}
}


/*!
Returns flux from given state into positive x direction.

Throws std::domain_error if given a state with non-positive mass density.
*/
template <
	class Container,
	class Mass_Density_Getter,
	class Momentum_Density_Getter,
	class Total_Energy_Density_Getter,
	class Magnetic_Field_Getter
> Container get_flux(
	Container& data,
	const double adiabatic_index,
	const double vacuum_permeability,
	const Mass_Density_Getter Mas,
	const Momentum_Density_Getter Mom,
	const Total_Energy_Density_Getter Nrj,
	const Magnetic_Field_Getter Mag
) {
	if (Mas(data) <= 0) {
		throw std::domain_error(
			std::string("Non-positive mass density given to ")
			+ __func__
			+ std::string(": ")
			+ boost::lexical_cast<std::string>(Mas(data))
		);
	}

	const auto
		inv_permeability = 1.0 / vacuum_permeability,
		pressure_magnetic
			= 0.5 * inv_permeability
			* (
				Mag(data)[0] * Mag(data)[0]
				+ Mag(data)[1] * Mag(data)[1]
				+ Mag(data)[2] * Mag(data)[2]
			),
		pressure_thermal
			= get_pressure(
				Mas(data),
				Mom(data),
				Nrj(data),
				Mag(data),
				adiabatic_index,
				vacuum_permeability
			);

	const auto velocity = get_velocity(Mom(data), Mas(data));

	Container flux;

	Mas(flux) = Mom(data)[0];

	Mom(flux)
		= Mom(data) * velocity[0]
		- Mag(data)[0] * Mag(data) * inv_permeability;
	Mom(flux)[0] += pressure_thermal + pressure_magnetic;

	Nrj(flux)
		= velocity[0] * (Nrj(data) + pressure_thermal + pressure_magnetic)
		- Mag(data)[0] * velocity.dot(Mag(data)) * inv_permeability;

	Mag(flux) = velocity[0] * Mag(data) - Mag(data)[0] * velocity;
	Mag(flux)[0] = 0;

	return flux;
}


/*!
Returns total energy density.

Throws std::domain_error if given a state with non-positive mass density.

Vel and mag must be compatible with std::array<double, 3>.
*/
template <
	class Velocity,
	class Magnetic_Field
> double get_total_energy_density(
	const double& mass_density,
	const Velocity& vel,
	const double& pressure,
	const Magnetic_Field& mag,
	const double& adiabatic_index,
	const double& vacuum_permeability
) {
	if (mass_density <= 0) {
		throw std::domain_error(
			std::string("Non-positive mass density given to ")
			+ __func__
			+ std::string(": ")
			+ boost::lexical_cast<std::string>(mass_density)
		);
	}

	if (adiabatic_index <= 1) {
		throw std::domain_error(
			__func__
			+ std::string("Adiabatic index must be > 1: ")
			+ boost::lexical_cast<std::string>(adiabatic_index)
		);
	}

	if (vacuum_permeability <= 0) {
		throw std::domain_error(
			__func__
			+ std::string("Vacuum permeability must be > 0: ")
			+ boost::lexical_cast<std::string>(vacuum_permeability)
		);
	}

	const auto
		kinetic_energy
			= 0.5 * mass_density
			* (
				vel[0] * vel[0]
				+ vel[1] * vel[1]
				+ vel[2] * vel[2]
			),
		magnetic_energy
			= 0.5 / vacuum_permeability
			* (
				mag[0] * mag[0]
				+ mag[1] * mag[1]
				+ mag[2] * mag[2]
			);

	return
		pressure / (adiabatic_index - 1)
		+ kinetic_energy
		+ magnetic_energy;
}


/*!
Returns speed of sound wave.

Throws std::domain_error if given a state with
non-positive mass density or pressure.

Mom and mag must be compatible with std::array<double, 3>.
*/
template <
	class Momentum_Density,
	class Magnetic_Field
> double get_sound_speed(
	const double& mass_density,
	const Momentum_Density& mom,
	const double& total_energy_density,
	const Magnetic_Field& mag,
	const double adiabatic_index,
	const double vacuum_permeability
) {
	if (mass_density <= 0) {
		throw std::domain_error(
			std::string("Non-positive mass density given to ")
			+ __func__
			+ std::string(": ")
			+ boost::lexical_cast<std::string>(mass_density)
		);
	}

	const auto pressure
		= get_pressure(
			mass_density,
			mom,
			total_energy_density,
			mag,
			adiabatic_index,
			vacuum_permeability
		);
	if (pressure <= 0) {
		throw std::domain_error(
			std::string("Non-positive pressure given to ")
			+ __func__
			+ std::string(": ")
			+ boost::lexical_cast<std::string>(pressure)
		);
	}

	return std::sqrt(adiabatic_index * pressure / mass_density);
}


/*!
Returns speed of Alfvén wave.

Throws std::domain_error if given a state with non-positive mass density.

Mag must be compatible with std::array<double, 3>.
*/
template <
	class Magnetic_Field
> double get_alfven_speed(
	const double& mass_density,
	const Magnetic_Field& mag,
	const double vacuum_permeability
) {
	using std::sqrt;

	if (mass_density <= 0) {
		throw std::domain_error(
			std::string("Non-positive mass density given to ")
			+ __func__
			+ std::string(": ")
			+ boost::lexical_cast<std::string>(mass_density)
		);
	}

	const auto mag_mag = sqrt(mag[0]*mag[0] + mag[1]*mag[1] + mag[2]*mag[2]);
	return mag_mag / sqrt(vacuum_permeability * mass_density);
}


/*!
Returns speed of fast magnetosonic wave in first dimension.

Returns a non-negative value.

Mom and mag must be compatible with std::array<double, 3>.
*/
template <
	class Momentum_Density,
	class Magnetic_Field
> double get_fast_magnetosonic_speed(
	const double& mass_density,
	const Momentum_Density& mom,
	const double& total_energy_density,
	const Magnetic_Field& mag,
	const double adiabatic_index,
	const double vacuum_permeability
) {
	using std::pow;
	using std::sqrt;

	const auto
		mag_mag = sqrt(mag[0]*mag[0] + mag[1]*mag[1] + mag[2]*mag[2]),
		sound = get_sound_speed(
			mass_density,
			mom,
			total_energy_density,
			mag,
			adiabatic_index,
			vacuum_permeability
		);

	if (mag_mag == 0) {
		return sound;
	}

	const auto
		sound2 = sound * sound,
		alfven2
			= pow(
				get_alfven_speed(mass_density, mag, vacuum_permeability),
				2
			),
		squared_sum = sound2 + alfven2;

	return
		sqrt(
			0.5 * (
				squared_sum
				+ sqrt(
					pow(squared_sum, 2)
					- 4 * sound2 * alfven2 * pow(mag[0] / mag_mag, 2)
				)
			)
		);
}


/*!
Throws std::domain_error if given a state with non-positive mass density.
*/
template <
	class Primitive,
	class P_Mass_Density_Getter,
	class Velocity_Getter,
	class Pressure_Getter,
	class P_Magnetic_Field_Getter,
	class Conservative,
	class C_Mass_Density_Getter,
	class Momentum_Density_Getter,
	class Total_Energy_Density_Getter,
	class C_Magnetic_Field_Getter
> Primitive get_primitive(
	Conservative data,
	const double adiabatic_index,
	const double vacuum_permeability,
	const P_Mass_Density_Getter Mas_p,
	const Velocity_Getter Vel,
	const Pressure_Getter Pre,
	const P_Magnetic_Field_Getter Mag_p,
	const C_Mass_Density_Getter Mas_c,
	const Momentum_Density_Getter Mom,
	const Total_Energy_Density_Getter Nrj,
	const C_Magnetic_Field_Getter Mag_c
) {
	if (Mas_c(data) <= 0) {
		throw std::domain_error(
			std::string("Non-positive mass density given to ")
			+ __func__
			+ std::string(": ")
			+ boost::lexical_cast<std::string>(Mas_c(data))
		);
	}

	Primitive ret_val;

	Mas_p(ret_val) = Mas_c(data);
	Vel(ret_val) = get_velocity(Mom(data), Mas_c(data));
	Pre(ret_val)
		= get_pressure(
			Mas_c(data),
			Mom(data),
			Nrj(data),
			Mag_c(data),
			adiabatic_index,
			vacuum_permeability
		);
	Mag_p(ret_val) = Mag_c(data);

	return ret_val;
}


/*!
Rotates components of given vector for get_flux_* functions.

Positive direction rotates vector for get_flux_* (assume
two states have identical y and z coordinates), negative
rotates it back.
*/
template<class Vector> Vector get_rotated_vector(
	const Vector& v,
	const int direction
) {
	Vector ret_val;

	switch(direction) {
	case -3:
		ret_val[0] = v[1];
		ret_val[1] = v[2];
		ret_val[2] = v[0];
		return ret_val;
	case -2:
		ret_val[0] = v[2];
		ret_val[1] = v[0];
		ret_val[2] = v[1];
		return ret_val;
	case -1:
	case 1:
		return v;
	case 2:
		ret_val[0] = v[1];
		ret_val[1] = v[2];
		ret_val[2] = v[0];
		return ret_val;
	case 3:
		ret_val[0] = v[2];
		ret_val[1] = v[0];
		ret_val[2] = v[1];
		return ret_val;
	default:
		std::cerr <<  __FILE__ << "(" << __LINE__ << ") "
			<< "Invalid direction: " << direction
			<< std::endl;
		abort();
	}
}


/*!
Function call signature of all MHD flux functions.

Input:
    - Conservative MHD variables in two cells that share
      a face and are neighbors in the x dimension
    - Area shared between given cells
    - Length of time step for which to calculate flux

state_neg represents the MHD variables in the cell in the
negative x direction from the shared face, state_pos in
the cell in positive x direction from the face.

Output:
    - Flux of conservative MHD variables over time dt
      through area shared_area
    - Absolute value of maximum signal speed from shared face

See for example hll_athena.hpp for an implementation.
*/
template <
	class MHD
> using solver_t = std::function<
	std::tuple<
		MHD,
		double
	>(
		MHD& /* state_neg */,
		MHD& /* state_pos */,
		const double& /* shared_area */,
		const double& /* dt */,
		const double& /* adiabatic_index */,
		const double& /* vacuum_permeability */
	)
>;


/*!
Positive flux adds to given cell multiplied with given factor.

Throws std::domain_error if new state has non-positive
mass density or pressure check_new_state == true.
*/
template <
	class Container,
	class Mass_Density_Getter,
	class Momentum_Density_Getter,
	class Total_Energy_Density_Getter,
	class Magnetic_Field_Getter,
	class Mass_Density_Flux_Getter,
	class Momentum_Density_Flux_Getter,
	class Total_Energy_Density_Flux_Getter,
	class Magnetic_Field_Flux_Getter
> void apply_fluxes(
	Container& data,
	const double factor,
	const double adiabatic_index,
	const double vacuum_permeability,
	const Mass_Density_Getter Mas,
	const Momentum_Density_Getter Mom,
	const Total_Energy_Density_Getter Nrj,
	const Magnetic_Field_Getter Mag,
	const Mass_Density_Flux_Getter Mas_f,
	const Momentum_Density_Flux_Getter Mom_f,
	const Total_Energy_Density_Flux_Getter Nrj_f,
	const Magnetic_Field_Flux_Getter Mag_f
) {
	Mas(data) += Mas_f(data) * factor;
	Mom(data) += Mom_f(data) * factor;
	Nrj(data) += Nrj_f(data) * factor;
	Mag(data) += Mag_f(data) * factor;

	if (Mas(data) <= 0) {
		throw std::domain_error(
			"New state has negative mass density: "
			+ boost::lexical_cast<std::string>(Mas(data))
		);
	}

	const auto pressure
		= get_pressure(
			Mas(data),
			Mom(data),
			Nrj(data),
			Mag(data),
			adiabatic_index,
			vacuum_permeability
		);
	if (pressure <= 0) {
		throw std::domain_error(
			"New state has negative pressure: "
			+ boost::lexical_cast<std::string>(pressure)
		);
	}
}

template <
	class Container,
	class Mass_Density_Getters,
	class Momentum_Density_Getters,
	class Total_Energy_Density_Getters,
	class Magnetic_Field_Getter,
	class Mass_Density_Flux_Getters,
	class Momentum_Density_Flux_Getters,
	class Total_Energy_Density_Flux_Getters,
	class Magnetic_Field_Flux_Getter
> void apply_fluxes_N(
	Container& data,
	const double factor,
	const Mass_Density_Getters Mas,
	const Momentum_Density_Getters Mom,
	const Total_Energy_Density_Getters Nrj,
	const Magnetic_Field_Getter Mag,
	const Mass_Density_Flux_Getters Mas_f,
	const Momentum_Density_Flux_Getters Mom_f,
	const Total_Energy_Density_Flux_Getters Nrj_f,
	const Magnetic_Field_Flux_Getter Mag_f
) {
	Mag(data) += Mag_f(data) * factor;

	Mas.first(data) += Mas_f.first(data) * factor;
	Mom.first(data) += Mom_f.first(data) * factor;
	Nrj.first(data) += Nrj_f.first(data) * factor;

	Mas.second(data) += Mas_f.second(data) * factor;
	Mom.second(data) += Mom_f.second(data) * factor;
	Nrj.second(data) += Nrj_f.second(data) * factor;
}


}} // namespaces

#endif // ifndef PAMHD_MHD_COMMON_HPP
