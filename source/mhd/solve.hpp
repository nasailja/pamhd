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

*/

#ifndef PAMHD_MHD_SOLVE_HPP
#define PAMHD_MHD_SOLVE_HPP

#include "cmath"
#include "functional"
#include "limits"
#include "vector"


namespace pamhd {
namespace mhd {


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
	class MHD_T,
	class MHD_Flux_T
> using solver_t = std::function<
	std::pair<
		typename MHD_Flux_T::data_type,
		double
	>(
		typename MHD_T::data_type /* state_neg */,
		typename MHD_T::data_type /* state_pos */,
		const double /* shared_area */,
		const double /* dt */,
		const double /* adiabatic_index */,
		const double /* vacuum_permeability */
	)
>;


/*!
Advances MHD solution for one time step of length dt with given solver.

Returns the maximum allowed length of time step for the next step on this process.
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
	const solver_t<MHD_T, MHD_Flux_T>& solver,
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

	if (not std::isfinite(dt) or dt < 0) {
		throw std::domain_error(
			"Invalid time step: "
			+ boost::lexical_cast<std::string>(dt)
		);
	}

	const Mass_Density_T Rho{};
	const Momentum_Density_T M{};
	const Total_Energy_Density_T E{};
	const Magnetic_Field_T B{};

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
			try {
				std::tie(
					flux,
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
					<< (*cell_data)[MHD_T()][Mass_Density_T()] << ", "
					<< (*cell_data)[MHD_T()][Momentum_Density_T()] << ", "
					<< (*cell_data)[MHD_T()][Total_Energy_Density_T()] << ", "
					<< (*cell_data)[MHD_T()][Magnetic_Field_T()] << " and "
					<< (*neighbor_data)[MHD_T()][Mass_Density_T()] << ", "
					<< (*neighbor_data)[MHD_T()][Momentum_Density_T()] << ", "
					<< (*neighbor_data)[MHD_T()][Total_Energy_Density_T()] << ", "
					<< (*neighbor_data)[MHD_T()][Magnetic_Field_T()]
					<< " because: " << error.what()
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
	class MHD_Flux_T,
	class Mass_Density_T,
	class Momentum_Density_T,
	class Total_Energy_Density_T,
	class Magnetic_Field_T
> void apply_fluxes(
	Grid_T& grid,
	const std::vector<uint64_t>& cells,
	const double adiabatic_index,
	const double vacuum_permeability
) {
	static_assert(
		std::is_same<typename MHD_T::data_type, typename MHD_Flux_T::data_type>::value,
		"The data types of variables MHD_T and MHD_Flux_T must be equal"
	);

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

		const auto previous_state = (*cell_data)[MHD_T()];
		const auto& flux = (*cell_data)[MHD_Flux_T()];
		try {
			apply_fluxes<
				Cell_T,
				MHD_T,
				MHD_Flux_T,
				Mass_Density_T,
				Momentum_Density_T,
				Total_Energy_Density_T,
				Magnetic_Field_T
			>(
				*cell_data,
				inverse_volume,
				adiabatic_index,
				vacuum_permeability
			);
		} catch (const std::domain_error& error) {
			std::cerr <<  __FILE__ << "(" << __LINE__
				<< ") New MHD state for cell " << cell_id
				<< " at " << grid.geometry.get_center(cell_id)
				<< " would be unphysical because (" << error.what()
				<< ") with old state:\n"
				<< previous_state[Mass_Density_T()] << ", "
				<< previous_state[Momentum_Density_T()] << ", "
				<< previous_state[Total_Energy_Density_T()] << ", "
				<< previous_state[Magnetic_Field_T()]
				<< "\nand flux * " << inverse_volume << ":\n"
				<< flux[Mass_Density_T()] * inverse_volume << ", "
				<< flux[Momentum_Density_T()] * inverse_volume << ", "
				<< flux[Total_Energy_Density_T()] * inverse_volume << ", "
				<< flux[Magnetic_Field_T()] * inverse_volume
				<< "\ngiving:\n"
				<< (*cell_data)[MHD_T()][Mass_Density_T()] << ", "
				<< (*cell_data)[MHD_T()][Momentum_Density_T()] << ", "
				<< (*cell_data)[MHD_T()][Total_Energy_Density_T()] << ", "
				<< (*cell_data)[MHD_T()][Magnetic_Field_T()]
				<< "\nwith pressure: "
				<< get_pressure<
					typename MHD_T::data_type,
					Mass_Density_T,
					Momentum_Density_T,
					Total_Energy_Density_T,
					Magnetic_Field_T
				>(
					(*cell_data)[MHD_T()],
					adiabatic_index,
					vacuum_permeability
				)
				<< std::endl;
			abort();
		}
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
		auto* const cell_data = grid[cell_id];
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
