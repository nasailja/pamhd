/*
Solves the MHD part of PAMHD using an external flux function.

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

#ifndef PAMHD_MHD_SOLVE_HPP
#define PAMHD_MHD_SOLVE_HPP


#include "cmath"
#include "functional"
#include "limits"
#include "vector"

#include "dccrg.hpp"

#include "mhd/variables.hpp"


namespace pamhd {
namespace mhd {


/*!
Advances MHD solution for one time step of length dt with given solver.

Returns the maximum allowed length of time step for the next step on this process.
*/
template <
	class Solver,
	class Cell,
	class Geometry,
	class Mass_Density_Getter,
	class Momentum_Density_Getter,
	class Total_Energy_Density_Getter,
	class Magnetic_Field_Getter,
	class Mass_Density_Flux_Getter,
	class Momentum_Density_Flux_Getter,
	class Total_Energy_Density_Flux_Getter,
	class Magnetic_Field_Flux_Getter
> double solve(
	const Solver solver,
	const std::vector<uint64_t>& cells,
	dccrg::Dccrg<Cell, Geometry>& grid,
	const double dt,
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
	if (not std::isfinite(dt) or dt < 0) {
		throw std::domain_error(
			"Invalid time step: "
			+ boost::lexical_cast<std::string>(dt)
		);
	}

	// shorthand for referring to variables of internal MHD data type
	const Mass_Density mas_int{};
	const Momentum_Density mom_int{};
	const Total_Energy_Density nrj_int{};
	const Magnetic_Field mag_int{};

	// maximum allowed next time step for cells of this process
	double max_dt = std::numeric_limits<double>::max();

	for (const auto& cell_id: cells) {

		auto* const cell_data = grid[cell_id];
		if (cell_data == NULL) {
			std::cerr <<  __FILE__ << "(" << __LINE__ << "): "
				<< " No data for cell " << cell_id
				<< " on process " << grid.get_rank()
				<< std::endl;
			abort();
		}

		const std::array<double, 3>
			cell_length = grid.geometry.get_length(cell_id),
			// area of cell perpendicular to each dimension
			cell_area{{
				cell_length[1] * cell_length[2],
				cell_length[0] * cell_length[2],
				cell_length[0] * cell_length[1]
			}};

		const auto face_neighbors = grid.get_face_neighbors_of(cell_id);
		for (const auto& item: face_neighbors) {
			const uint64_t neighbor_id = item.first;

			auto* const neighbor_data = grid[neighbor_id];
			if (neighbor_data == NULL) {
				std::cerr <<  __FILE__ << "(" << __LINE__ << "): "
					<< " No data for neighbor " << neighbor_id
					<< " owned by " << grid.get_process(neighbor_id)
					<< " on process " << grid.get_rank()
					<< std::endl;
				abort();
			}

			// direction of neighbor from cell
			const int neighbor_dir = item.second;
			if (
				neighbor_dir == 0
				or neighbor_dir > 3
				or neighbor_dir < -3
			) {
				std::cerr <<  __FILE__ << "(" << __LINE__ << "): "
					<< " Invalid neighbor direction " << neighbor_dir
					<< std::endl;
				abort();
			}

			if (grid.is_local(neighbor_id) and neighbor_dir < 0) {
				/*
				This case is handled when neighbor is the current cell
				and current cell is neighbor in positive direction
				*/
				continue;
			}

			const size_t neighbor_dim = size_t(abs(neighbor_dir) - 1);

			const std::array<double, 3>
				neighbor_length = grid.geometry.get_length(neighbor_id),
				neighbor_area{{
					neighbor_length[1] * neighbor_length[2],
					neighbor_length[0] * neighbor_length[2],
					neighbor_length[0] * neighbor_length[1]
				}};

			const double shared_area
				= std::min(cell_area[neighbor_dim], neighbor_area[neighbor_dim]);

			if (not std::isnormal(shared_area) or shared_area < 0) {
				throw std::domain_error(
					"Invalid area between cells "
					+ boost::lexical_cast<std::string>(cell_id)
					+ " and "
					+ boost::lexical_cast<std::string>(neighbor_id)
					+ ": "
					+ boost::lexical_cast<std::string>(shared_area)
				);
			}

			// take into account direction of neighbor from cell
			MHD_Conservative state_neg, state_pos;
			if (neighbor_dir > 0) {
				state_neg[mas_int] = Mas(*cell_data);
				state_neg[mom_int] = get_rotated_vector(Mom(*cell_data), abs(neighbor_dir));
				state_neg[Total_Energy_Density()] = Nrj(*cell_data);
				state_neg[mag_int] = get_rotated_vector(Mag(*cell_data), abs(neighbor_dir));

				state_pos[mas_int] = Mas(*neighbor_data);
				state_pos[mom_int] = get_rotated_vector(Mom(*neighbor_data), abs(neighbor_dir));
				state_pos[Total_Energy_Density()] = Nrj(*neighbor_data);
				state_pos[mag_int] = get_rotated_vector(Mag(*neighbor_data), abs(neighbor_dir));
			} else {
				state_pos[mas_int] = Mas(*cell_data);
				state_pos[mom_int] = get_rotated_vector(Mom(*cell_data), abs(neighbor_dir));
				state_pos[Total_Energy_Density()] = Nrj(*cell_data);
				state_pos[mag_int] = get_rotated_vector(Mag(*cell_data), abs(neighbor_dir));

				state_neg[mas_int] = Mas(*neighbor_data);
				state_neg[mom_int] = get_rotated_vector(Mom(*neighbor_data), abs(neighbor_dir));
				state_neg[Total_Energy_Density()] = Nrj(*neighbor_data);
				state_neg[mag_int] = get_rotated_vector(Mag(*neighbor_data), abs(neighbor_dir));
			}

			MHD_Conservative flux_neg, flux_pos;
			double max_vel;
			try {
				std::tie(
					flux_neg,
					flux_pos,
					max_vel
				) = solver(
					state_neg,
					state_pos,
					shared_area,
					dt,
					adiabatic_index,
					vacuum_permeability
				);
			} catch (const std::domain_error& error) {
				std::cerr <<  __FILE__ << "(" << __LINE__ << ") "
					<< "Solution failed between cells " << cell_id
					<< " and " << neighbor_id
					<< " at " << grid.geometry.get_center(cell_id)
					<< " and " << grid.geometry.get_center(neighbor_id)
					<< " in direction " << neighbor_dir
					<< " with states (mass, momentum, total energy, magnetic field): "
					<< Mas(*cell_data) << ", "
					<< Mom(*cell_data) << ", "
					<< Nrj(*cell_data) << ", "
					<< Mag(*cell_data) << " and "
					<< Mas(*neighbor_data) << ", "
					<< Mom(*neighbor_data) << ", "
					<< Nrj(*neighbor_data) << ", "
					<< Mag(*neighbor_data)
					<< " because: " << error.what()
					<< std::endl;
				abort();
			}

			max_dt = std::min(max_dt, cell_length[neighbor_dim] / max_vel);

			// rotate flux back
			flux_neg[mom_int] = get_rotated_vector(flux_neg[mom_int], -abs(neighbor_dir));
			flux_neg[mag_int] = get_rotated_vector(flux_neg[mag_int], -abs(neighbor_dir));
			flux_pos[mom_int] = get_rotated_vector(flux_pos[mom_int], -abs(neighbor_dir));
			flux_pos[mag_int] = get_rotated_vector(flux_pos[mag_int], -abs(neighbor_dir));

			if (neighbor_dir > 0) {
				Mas_f(*cell_data) -= flux_neg[mas_int] + flux_pos[mas_int];
				Mom_f(*cell_data) -= flux_neg[mom_int] + flux_pos[mom_int];
				Nrj_f(*cell_data) -= flux_neg[nrj_int] + flux_pos[nrj_int];
				Mag_f(*cell_data) -= flux_neg[mag_int] + flux_pos[mag_int];

				if (grid.is_local(neighbor_id)) {
					Mas_f(*neighbor_data) += flux_neg[mas_int] + flux_pos[mas_int];
					Mom_f(*neighbor_data) += flux_neg[mom_int] + flux_pos[mom_int];
					Nrj_f(*neighbor_data) += flux_neg[nrj_int] + flux_pos[nrj_int];
					Mag_f(*neighbor_data) += flux_neg[mag_int] + flux_pos[mag_int];
				}
			} else {
				Mas_f(*cell_data) += flux_neg[mas_int] + flux_pos[mas_int];
				Mom_f(*cell_data) += flux_neg[mom_int] + flux_pos[mom_int];
				Nrj_f(*cell_data) += flux_neg[nrj_int] + flux_pos[nrj_int];
				Mag_f(*cell_data) += flux_neg[mag_int] + flux_pos[mag_int];

				if (grid.is_local(neighbor_id)) {
					std::cerr <<  __FILE__ << "(" << __LINE__ << ") "
						"Invalid direction for adding flux to local neighbor."
						<< std::endl;
					abort();
				}
			}
		}
	}

	return max_dt;
}


/*!
Applies the MHD solution to given cells.
*/
template <
	class Cell,
	class Geometry,
	class Mass_Density_Getter,
	class Momentum_Density_Getter,
	class Total_Energy_Density_Getter,
	class Magnetic_Field_Getter,
	class Mass_Density_Flux_Getter,
	class Momentum_Density_Flux_Getter,
	class Total_Energy_Density_Flux_Getter,
	class Magnetic_Field_Flux_Getter
> void apply_fluxes(
	const std::vector<uint64_t>& cells,
	dccrg::Dccrg<Cell, Geometry>& grid,
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
	for (const auto& cell_id: cells) {
		auto* const cell_data = grid[cell_id];
		if (cell_data == NULL) {
			std::cerr <<  __FILE__ << "(" << __LINE__ << "): "
				<< " No data for cell " << cell_id
				<< std::endl;
			abort();
		}

		const auto length = grid.geometry.get_length(cell_id);
		const double inverse_volume = 1.0 / (length[0] * length[1] * length[2]);

		try {
			apply_fluxes(
				*cell_data,
				inverse_volume,
				adiabatic_index,
				vacuum_permeability,
				Mas, Mom, Nrj, Mag,
				Mas_f, Mom_f, Nrj_f, Mag_f
			);
		} catch (const std::domain_error& error) {
			std::cerr <<  __FILE__ << "(" << __LINE__
				<< ") New MHD state for cell " << cell_id
				<< " at " << grid.geometry.get_center(cell_id)
				<< " would be unphysical because (" << error.what()
				<< ") with flux (* " << inverse_volume << "):\n"
				<< Mas_f(*cell_data) * inverse_volume << ", "
				<< Mom_f(*cell_data) * inverse_volume << ", "
				<< Nrj_f(*cell_data) * inverse_volume << ", "
				<< Mag_f(*cell_data) * inverse_volume
				<< "\ngiving:\n"
				<< Mas(*cell_data) << ", "
				<< Mom(*cell_data) << ", "
				<< Nrj(*cell_data) << ", "
				<< Mag(*cell_data)
				<< "\nwith pressure: "
				<< get_pressure(
					Mas(*cell_data),
					Mom(*cell_data),
					Nrj(*cell_data),
					Mag(*cell_data),
					adiabatic_index,
					vacuum_permeability
				)
				<< std::endl;
			abort();
		}
	}
}

/*!
Zeros fluxes in given cells.
*/
template <
	class Cell,
	class Geometry,
	class Mass_Density_Getter,
	class Momentum_Density_Getter,
	class Total_Energy_Density_Getter,
	class Magnetic_Field_Getter
> void zero_fluxes(
	const std::vector<uint64_t>& cells,
	dccrg::Dccrg<Cell, Geometry>& grid,
	Mass_Density_Getter Mas,
	Momentum_Density_Getter Mom,
	Total_Energy_Density_Getter Nrj,
	Magnetic_Field_Getter Mag
) {
	for (const auto& cell_id: cells) {
		auto* const cell_data = grid[cell_id];
		if (cell_data == NULL) {
			std::cerr <<  __FILE__ << "(" << __LINE__ << "): "
				<< " No data for cell " << cell_id
				<< std::endl;
			abort();
		}

		Mas(*cell_data)    =
		Mom(*cell_data)[0] =
		Mom(*cell_data)[1] =
		Mom(*cell_data)[2] =
		Nrj(*cell_data)    =
		Mag(*cell_data)[0] =
		Mag(*cell_data)[1] =
		Mag(*cell_data)[2] = 0;
	}
}


}} // namespaces


#endif // ifndef PAMHD_MHD_SOLVE_HPP
