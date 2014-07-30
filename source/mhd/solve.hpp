/*
Solves the MHD part of PAMHD.

Copyright 2014 Ilja Honkonen
All rights reserved.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef PAMHD_MHD_SOLVE_HPP
#define PAMHD_MHD_SOLVE_HPP

#include "cmath"
#include "limits"
#include "vector"

#include "gensimcell.hpp"

#include "mhd/hll_athena.hpp"
#include "mhd/hlld_athena.hpp"

namespace pamhd {
namespace mhd {


/*!
Advances the MHD solution for one time step of length dt with given solver.

Returns the maximum allowed length of time step for the next step on this
process or a negative value in case of error.
*/
template <
	class Grid_T,
	class MHD_T,
	class MHD_Flux_T,
	class Mass_Density_T,
	class Momentum_Density_T,
	class Total_Energy_Density_T,
	class Magnetic_Field_T
> double solve(
	const std::string& solver,
	Grid_T& grid,
	const std::vector<uint64_t>& cells,
	const double dt,
	const double adiabatic_index,
	const double vacuum_permeability
) {
	static_assert(
		std::is_same<typename MHD_T::data_type, typename MHD_Flux_T::data_type>::value,
		"The data types of variables MHD_T and MHD_Flux_T must be equal"
	);

	const Mass_Density_T Rho{};
	const Momentum_Density_T M{};
	const Total_Energy_Density_T E{};
	const Magnetic_Field_T B{};

	// maximum allowed next time step for cells of this process
	double max_dt = std::numeric_limits<double>::max();

	for (const auto& cell_id: cells) {

		auto* cell_data = grid[cell_id];
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
			cell_area = {
				cell_length[1] * cell_length[2],
				cell_length[0] * cell_length[2],
				cell_length[0] * cell_length[1]
			};

		const auto face_neighbors = grid.get_face_neighbors_of(cell_id);
		for (const auto& item: face_neighbors) {
			const uint64_t neighbor_id = item.first;

			auto* neighbor_data = grid[neighbor_id];
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
				neighbor_area = {
					neighbor_length[1] * neighbor_length[2],
					neighbor_length[0] * neighbor_length[2],
					neighbor_length[0] * neighbor_length[1]
				};

			const double shared_area
				= std::min(cell_area[neighbor_dim], neighbor_area[neighbor_dim]);

			// take into account direction of neighbor from cell
			const typename MHD_T::data_type
				state_neg
					= [&](){
						if (neighbor_dir > 0) {
							return get_rotated_state<
								MHD_T,
								Mass_Density_T,
								Momentum_Density_T,
								Total_Energy_Density_T,
								Magnetic_Field_T
							>(
								(*cell_data)[MHD_T()],
								abs(neighbor_dir)
							);
						} else {
							return get_rotated_state<
								MHD_T,
								Mass_Density_T,
								Momentum_Density_T,
								Total_Energy_Density_T,
								Magnetic_Field_T
							>(
								(*neighbor_data)[MHD_T()],
								abs(neighbor_dir)
							);
						}
					}(),

				state_pos
					= [&](){
						if (neighbor_dir > 0) {
							return get_rotated_state<
								MHD_T,
								Mass_Density_T,
								Momentum_Density_T,
								Total_Energy_Density_T,
								Magnetic_Field_T
							>(
								(*neighbor_data)[MHD_T()],
								abs(neighbor_dir)
							);
						} else {
							return get_rotated_state<
								MHD_T,
								Mass_Density_T,
								Momentum_Density_T,
								Total_Energy_Density_T,
								Magnetic_Field_T
							>(
								(*cell_data)[MHD_T()],
								abs(neighbor_dir)
							);
						}
					}();


			typename MHD_Flux_T::data_type flux;
			double max_vel;
			if (solver == "hll_athena") {
				std::tie(
					flux,
					max_vel
				) = pamhd::mhd::athena::get_flux_hll<
					typename MHD_T::data_type,
					Mass_Density_T,
					Momentum_Density_T,
					Total_Energy_Density_T,
					Magnetic_Field_T
				>(
					state_neg,
					state_pos,
					shared_area,
					dt,
					adiabatic_index,
					vacuum_permeability
				);
			} else if (solver == "hlld_athena") {
				std::tie(
					flux,
					max_vel
				) = pamhd::mhd::athena::get_flux_hlld<
					typename MHD_T::data_type,
					Mass_Density_T,
					Momentum_Density_T,
					Total_Energy_Density_T,
					Magnetic_Field_T
				>(
					state_neg,
					state_pos,
					shared_area,
					dt,
					adiabatic_index,
					vacuum_permeability
				);
			} else {
				std::cerr <<  __FILE__ << "(" << __LINE__ << ") Invalid solver: "
					<< solver << ", use --help to list available solvers"
					<< std::endl;
				abort();
			}
			if (max_vel < 0) {
				std::cerr <<  __FILE__ << "(" << __LINE__ << ") "
					<< "Solution failed between cells " << cell_id
					<< " and " << neighbor_id
					<< " at " << grid.geometry.get_center(cell_id)
					<< " and " << grid.geometry.get_center(neighbor_id)
					<< " with states (mass, momentum, total energy, magnetic field): "
					<< (*cell_data)[MHD_T()][Mass_Density_T()] << ", "
					<< (*cell_data)[MHD_T()][Momentum_Density_T()] << ", "
					<< (*cell_data)[MHD_T()][Total_Energy_Density_T()] << ", "
					<< (*cell_data)[MHD_T()][Magnetic_Field_T()] << " and "
					<< (*neighbor_data)[MHD_T()][Mass_Density_T()] << ", "
					<< (*neighbor_data)[MHD_T()][Momentum_Density_T()] << ", "
					<< (*neighbor_data)[MHD_T()][Total_Energy_Density_T()] << ", "
					<< (*neighbor_data)[MHD_T()][Magnetic_Field_T()]
					<< std::endl;
				abort();
			}

			max_dt = std::min(max_dt, cell_length[neighbor_dim] / max_vel);

			// rotate flux back
			flux = get_rotated_state<
				MHD_T,
				Mass_Density_T,
				Momentum_Density_T,
				Total_Energy_Density_T,
				Magnetic_Field_T
			>(flux, -abs(neighbor_dir));


			const MHD_Flux_T Flux{};
			if (neighbor_dir > 0) {
				(*cell_data)[Flux][Rho] -= flux[Rho];
				(*cell_data)[Flux][M] -= flux[M];
				(*cell_data)[Flux][E] -= flux[E];
				(*cell_data)[Flux][B] -= flux[B];

				if (grid.is_local(neighbor_id)) {
					(*neighbor_data)[Flux][Rho] += flux[Rho];
					(*neighbor_data)[Flux][M] += flux[M];
					(*neighbor_data)[Flux][E] += flux[E];
					(*neighbor_data)[Flux][B] += flux[B];
				}
			} else {
				(*cell_data)[Flux][Rho] += flux[Rho];
				(*cell_data)[Flux][M] += flux[M];
				(*cell_data)[Flux][E] += flux[E];
				(*cell_data)[Flux][B] += flux[B];

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
	class Grid_T,
	class Cell_T,
	class MHD_T,
	class MHD_Flux_T
> void apply_fluxes(
	Grid_T& grid,
	const std::vector<uint64_t>& cells
) {
	static_assert(
		std::is_same<typename MHD_T::data_type, typename MHD_Flux_T::data_type>::value,
		"The data types of variables MHD_T and MHD_Flux_T must be equal"
	);

	for (const auto& cell_id: cells) {
		auto* cell_data = grid[cell_id];
		if (cell_data == NULL) {
			std::cerr <<  __FILE__ << "(" << __LINE__ << "): "
				<< " No data for cell " << cell_id
				<< std::endl;
			abort();
		}

		const auto length = grid.geometry.get_length(cell_id);
		const double inverse_volume = 1.0 / (length[0] * length[1] * length[2]);

		apply_fluxes<Cell_T, MHD_T, MHD_Flux_T>((*cell_data), inverse_volume);
	}
}

/*!
Applies the MHD solution to given cells.

Zeros fluxes.
*/
template <
	class Grid_T,
	class Cell_T,
	class MHD_Flux_T,
	class Mass_Density_T,
	class Momentum_Density_T,
	class Total_Energy_Density_T,
	class Magnetic_Field_T
> void zero_fluxes(
	Grid_T& grid,
	const std::vector<uint64_t>& cells
) {
	for (const auto& cell_id: cells) {
		auto* cell_data = grid[cell_id];
		if (cell_data == NULL) {
			std::cerr <<  __FILE__ << "(" << __LINE__ << "): "
				<< " No data for cell " << cell_id
				<< std::endl;
			abort();
		}

		zero_fluxes<
			Cell_T,
			MHD_Flux_T,
			Mass_Density_T,
			Momentum_Density_T,
			Total_Energy_Density_T,
			Magnetic_Field_T
		>(*cell_data);
	}
}

}} // namespaces

#endif // ifndef PAMHD_MHD_SOLVE_HPP
