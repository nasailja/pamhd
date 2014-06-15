/*
Saves the MHD solution of PAMHD.

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

#ifndef PAMHD_MHD_SOLVE_HPP
#define PAMHD_MHD_SOLVE_HPP

#include "cmath"
#include "limits"
#include "vector"

#include "gensimcell.hpp"

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

		const double min_cell_length
			= std::min({
				cell_length[0],
				cell_length[1],
				cell_length[2]
			});

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

			if (neighbor_dim != 0) {
				std::cerr <<  __FILE__ << "(" << __LINE__ << ") Dimension"
					<< neighbor_dim << " not supported."
					<< std::endl;
				abort();
			}
			// TODO: rotate vector variables for other dimensions
			const typename MHD_T::data_type
				state_neg
					= [&](){
						if (neighbor_dir > 0) {
							return (*cell_data)[MHD_T()];
						} else {
							return (*neighbor_data)[MHD_T()];
						}
					}(),
				state_pos
					= [&](){
						if (neighbor_dir > 0) {
							return (*neighbor_data)[MHD_T()];
						} else {
							return (*cell_data)[MHD_T()];
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
			} else {
				std::cerr <<  __FILE__ << "(" << __LINE__ << ") Invalid solver: "
					<< solver << ", use --help to list available solvers"
					<< std::endl;
				abort();
			}
			max_dt = std::min(max_dt, min_cell_length / max_vel);

			// TODO: rotate vector variables of flux back

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

Zeros fluxes.
*/
template <
	class Grid_T,
	class Cell_T,
	class MHD_T,
	class MHD_Flux_T,
	class Mass_Density_T,
	class Momentum_Density_T,
	class Total_Energy_Density_T,
	class Magnetic_Field_T
> void apply_solution(
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

		apply_flux<
			Cell_T,
			MHD_T,
			MHD_Flux_T,
			Mass_Density_T,
			Momentum_Density_T,
			Total_Energy_Density_T,
			Magnetic_Field_T
		>(
			(*cell_data),
			inverse_volume
		);
	}
}

}} // namespaces

#endif // ifndef PAMHD_MHD_SOLVE_HPP
