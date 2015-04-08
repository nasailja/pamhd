/*
Tests parallel particle solver of PAMHD in 2 dimensions.

Copyright 2015 Ilja Honkonen
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

#include "array"
#include "cmath"
#include "cstdlib"
#include "iostream"

#include "boost/numeric/odeint.hpp"
#include "dccrg.hpp"
#include "dccrg_cartesian_geometry.hpp"
#include "Eigen/Core" // must be included before gensimcell
#include "Eigen/Geometry"
#include "mpi.h" // must be included before gensimcell
#include "gensimcell.hpp"

#include "particle/solve_dccrg.hpp"
#include "particle/variables.hpp"

using namespace std;
using namespace pamhd::particle;


int main(int argc, char* argv[])
{
	using Cell_T = pamhd::particle::Cell;
	using Grid_T = dccrg::Dccrg<Cell_T, dccrg::Cartesian_Geometry>;

	/*
	Initialize MPI
	*/

	if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
		std::cerr << "Couldn't initialize MPI." << std::endl;
		abort();
	}

	MPI_Comm comm = MPI_COMM_WORLD;

	int rank = 0, comm_size = 0;
	if (MPI_Comm_rank(comm, &rank) != MPI_SUCCESS) {
		std::cerr << "Couldn't obtain MPI rank." << std::endl;
		abort();
	}
	if (MPI_Comm_size(comm, &comm_size) != MPI_SUCCESS) {
		std::cerr << "Couldn't obtain size of MPI communicator." << std::endl;
		abort();
	}

	// intialize Zoltan
	float zoltan_version;
	if (Zoltan_Initialize(argc, argv, &zoltan_version) != ZOLTAN_OK) {
		std::cerr << "Zoltan_Initialize failed." << std::endl;
		abort();
	}


	Grid_T grid_x, grid_y, grid_z;

	const unsigned int neighborhood_size = 1;
	const std::array<uint64_t, 3>
		number_of_cells_x{{ 1, 10, 10}},
		number_of_cells_y{{10,  1, 10}},
		number_of_cells_z{{10, 10,  1}};
	if (
		not grid_x.initialize(
			number_of_cells_x,
			comm,
			"RCB",
			neighborhood_size,
			0,
			false, false, false)
		or not grid_y.initialize(
			number_of_cells_y,
			comm,
			"RCB",
			neighborhood_size,
			0,
			false, false, false)
		or not grid_z.initialize(
			number_of_cells_z,
			comm,
			"RCB",
			neighborhood_size,
			0,
			false, false, false
		)
	) {
		std::cerr << __FILE__ << "(" << __LINE__ << "): "
			<< "Couldn't initialize one or more grids."
			<< std::endl;
		abort();
	}

	// set grid geometry
	dccrg::Cartesian_Geometry::Parameters geom_params;
	geom_params.start = {{0.0, 0.0, 0.0}};
	geom_params.level_0_cell_length = {{1.0, 1.0, 1.0}};

	if (
		not grid_x.set_geometry(geom_params)
		or not grid_y.set_geometry(geom_params)
		or not grid_z.set_geometry(geom_params)
	) {
		std::cerr << __FILE__ << "(" << __LINE__ << "): "
			<< "Couldn't set geometry of one or more grids."
			<< std::endl;
		abort();
	}

	// use same domain decomposition in all grids
	for (size_t cell_id = 1; cell_id <= 100; cell_id++) {
		grid_x.pin(cell_id, cell_id % grid_x.get_comm_size());
		grid_y.pin(cell_id, cell_id % grid_y.get_comm_size());
		grid_z.pin(cell_id, cell_id % grid_z.get_comm_size());
	}
	grid_x.balance_load(false);
	grid_y.balance_load(false);
	grid_z.balance_load(false);

	const auto
		inner_cell_ids = grid_x.get_local_cells_not_on_process_boundary(),
		outer_cell_ids = grid_x.get_local_cells_on_process_boundary(),
		remote_cell_ids = grid_x.get_remote_cells_on_process_boundary(),
		cell_ids = grid_x.get_cells();

	// initial condition
	for (const auto& cell_id: cell_ids) {
		auto
			*const cell_ptr_x = grid_x[cell_id],
			*const cell_ptr_y = grid_y[cell_id],
			*const cell_ptr_z = grid_z[cell_id];
		if (
			cell_ptr_x == nullptr
			or cell_ptr_y == nullptr
			or cell_ptr_z == nullptr
		) {
			std::cerr << __FILE__ << "(" << __LINE__ << ")" << std::endl;
			abort();
		}
		auto
			&cell_data_x = *cell_ptr_x,
			&cell_data_y = *cell_ptr_y,
			&cell_data_z = *cell_ptr_z;

		const auto
			cell_center_x = grid_x.geometry.get_center(cell_id),
			cell_center_y = grid_y.geometry.get_center(cell_id),
			cell_center_z = grid_z.geometry.get_center(cell_id);

		Particle_Internal particle_x, particle_y, particle_z;
		particle_x[Position()] = {
			cell_center_x[0],
			cell_center_x[1],
			cell_center_x[2]
		};
		particle_y[Position()] = {
			cell_center_y[0],
			cell_center_y[1],
			cell_center_y[2]
		};
		particle_z[Position()] = {
			cell_center_z[0],
			cell_center_z[1],
			cell_center_z[2]
		};
		particle_x[Velocity()] = {  0, 1.0, 1.0};
		particle_y[Velocity()] = {1.0,   0, 1.0};
		particle_z[Velocity()] = {1.0, 1.0,   0};

		particle_x[Mass()] =
		particle_y[Mass()] =
		particle_z[Mass()] =
		particle_x[Charge_Mass_Ratio()] =
		particle_y[Charge_Mass_Ratio()] =
		particle_z[Charge_Mass_Ratio()] = 0;

		cell_data_x[Particles_Internal()].push_back(particle_x);
		cell_data_y[Particles_Internal()].push_back(particle_y);
		cell_data_z[Particles_Internal()].push_back(particle_z);

		cell_data_x[Electric_Field()] =
		cell_data_x[Magnetic_Field()] =
		cell_data_y[Electric_Field()] =
		cell_data_y[Magnetic_Field()] =
		cell_data_z[Electric_Field()] =
		cell_data_z[Magnetic_Field()] = {0, 0, 0};

		cell_data_x[Nr_Particles_External()] = cell_data_x[Particles_External()].size();
		cell_data_y[Nr_Particles_External()] = cell_data_y[Particles_External()].size();
		cell_data_z[Nr_Particles_External()] = cell_data_z[Particles_External()].size();
	}
	// allocate copies of remote neighbor cells
	grid_x.update_copies_of_remote_neighbors();
	grid_y.update_copies_of_remote_neighbors();
	grid_z.update_copies_of_remote_neighbors();

	// short hand notation for calling solvers
	auto solve = [](
		const std::vector<uint64_t>& cell_ids,
		Grid_T& grid
	) {
		pamhd::particle::solve<
			Cell_T,
			pamhd::particle::Electric_Field,
			pamhd::particle::Magnetic_Field,
			pamhd::particle::Nr_Particles_External,
			pamhd::particle::Particles_Internal,
			pamhd::particle::Particles_External,
			pamhd::particle::Position,
			pamhd::particle::Velocity,
			pamhd::particle::Charge_Mass_Ratio,
			pamhd::particle::Mass,
			pamhd::particle::Destination_Cell,
			boost::numeric::odeint::runge_kutta_fehlberg78<state_t>
		>(1.0, cell_ids, grid);
	};

	auto resize_receiving = [](
		const std::vector<uint64_t>& cell_ids,
		Grid_T& grid
	) {
		pamhd::particle::resize_receiving_containers<
			Cell_T,
			pamhd::particle::Nr_Particles_External,
			pamhd::particle::Particles_External
		>(cell_ids, grid);
	};

	auto incorporate_external = [](
		const std::vector<uint64_t>& cell_ids,
		Grid_T& grid
	) {
		pamhd::particle::incorporate_external_particles<
			Cell_T,
			pamhd::particle::Nr_Particles_External,
			pamhd::particle::Particles_Internal,
			pamhd::particle::Particles_External,
			pamhd::particle::Destination_Cell
		>(cell_ids, grid);
	};

	auto remove_external = [](
		const std::vector<uint64_t>& cell_ids,
		Grid_T& grid
	) {
		pamhd::particle::remove_external_particles<
			Cell_T,
			pamhd::particle::Nr_Particles_External,
			pamhd::particle::Particles_External
		>(cell_ids, grid);
	};


	for (size_t step = 0; step < 10; step++) {
		solve(outer_cell_ids, grid_x);
		solve(outer_cell_ids, grid_y);
		solve(outer_cell_ids, grid_z);

		Cell_T::set_transfer_all(
			true,
			pamhd::particle::Electric_Field(),
			pamhd::particle::Magnetic_Field(),
			pamhd::particle::Nr_Particles_External()
		);
		grid_x.start_remote_neighbor_copy_updates();
		grid_y.start_remote_neighbor_copy_updates();
		grid_z.start_remote_neighbor_copy_updates();

		solve(inner_cell_ids, grid_x);
		solve(inner_cell_ids, grid_y);
		solve(inner_cell_ids, grid_z);

		grid_x.wait_remote_neighbor_copy_update_receives();
		grid_y.wait_remote_neighbor_copy_update_receives();
		grid_z.wait_remote_neighbor_copy_update_receives();
		resize_receiving(remote_cell_ids, grid_x);
		resize_receiving(remote_cell_ids, grid_y);
		resize_receiving(remote_cell_ids, grid_z);

		grid_x.wait_remote_neighbor_copy_update_sends();
		grid_y.wait_remote_neighbor_copy_update_sends();
		grid_z.wait_remote_neighbor_copy_update_sends();

		Cell_T::set_transfer_all(
			false,
			pamhd::particle::Electric_Field(),
			pamhd::particle::Magnetic_Field(),
			pamhd::particle::Nr_Particles_External()
		);
		Cell_T::set_transfer_all(
			true,
			pamhd::particle::Particles_External()
		);

		grid_x.start_remote_neighbor_copy_updates();
		grid_y.start_remote_neighbor_copy_updates();
		grid_z.start_remote_neighbor_copy_updates();

		incorporate_external(inner_cell_ids, grid_x);
		incorporate_external(inner_cell_ids, grid_y);
		incorporate_external(inner_cell_ids, grid_z);

		grid_x.wait_remote_neighbor_copy_update_receives();
		grid_y.wait_remote_neighbor_copy_update_receives();
		grid_z.wait_remote_neighbor_copy_update_receives();

		incorporate_external(outer_cell_ids, grid_x);
		incorporate_external(outer_cell_ids, grid_y);
		incorporate_external(outer_cell_ids, grid_z);

		remove_external(inner_cell_ids, grid_x);
		remove_external(inner_cell_ids, grid_y);
		remove_external(inner_cell_ids, grid_z);

		grid_x.wait_remote_neighbor_copy_update_sends();
		grid_y.wait_remote_neighbor_copy_update_sends();
		grid_z.wait_remote_neighbor_copy_update_sends();
		Cell_T::set_transfer_all(
			false,
			pamhd::particle::Particles_External()
		);

		remove_external(outer_cell_ids, grid_x);
		remove_external(outer_cell_ids, grid_y);
		remove_external(outer_cell_ids, grid_z);


		// check that solution is correct
		std::array<int, 3>
			total_particles_local{{0, 0, 0}},
			total_particles{{0, 0, 0}};

		for (const auto& cell_id: cell_ids) {
			auto
				*const cell_ptr_x = grid_x[cell_id],
				*const cell_ptr_y = grid_y[cell_id],
				*const cell_ptr_z = grid_z[cell_id];
			if (
				cell_ptr_x == nullptr
				or cell_ptr_y == nullptr
				or cell_ptr_z == nullptr
			) {
				std::cerr << __FILE__ << "(" << __LINE__ << ")" << std::endl;
				abort();
			}
			auto
				&cell_data_x = *cell_ptr_x,
				&cell_data_y = *cell_ptr_y,
				&cell_data_z = *cell_ptr_z;

			total_particles_local[0]
				+= cell_data_x[pamhd::particle::Particles_Internal()].size()
				+ cell_data_x[pamhd::particle::Particles_External()].size();
			total_particles_local[1]
				+= cell_data_y[pamhd::particle::Particles_Internal()].size()
				+ cell_data_y[pamhd::particle::Particles_External()].size();
			total_particles_local[2]
				+= cell_data_z[pamhd::particle::Particles_Internal()].size()
				+ cell_data_z[pamhd::particle::Particles_External()].size();

			if (cell_data_x[pamhd::particle::Particles_Internal()].size() > 1) {
				std::cerr << __FILE__ << "(" << __LINE__ << ") "
					<< "Incorrect number of internal particles in cell " << cell_id
					<< ": " << cell_data_x[pamhd::particle::Particles_Internal()].size()
					<< std::endl;
				abort();
			}
			if (cell_data_y[pamhd::particle::Particles_Internal()].size() > 1) {
				std::cerr << __FILE__ << "(" << __LINE__ << ") "
					<< "Incorrect number of internal particles in cell " << cell_id
					<< ": " << cell_data_y[pamhd::particle::Particles_Internal()].size()
					<< std::endl;
				abort();
			}
			if (cell_data_z[pamhd::particle::Particles_Internal()].size() > 1) {
				std::cerr << __FILE__ << "(" << __LINE__ << ") "
					<< "Incorrect number of internal particles in cell " << cell_id
					<< ": " << cell_data_z[pamhd::particle::Particles_Internal()].size()
					<< std::endl;
				abort();
			}
		}

		MPI_Allreduce(
			total_particles_local.data(),
			total_particles.data(),
			3,
			MPI_INT,
			MPI_SUM,
			comm
		);
		if (total_particles[0] != int(std::pow(9 - step, 2))) {
			if (grid_x.get_rank() == 0) {
				std::cerr << __FILE__ << "(" << __LINE__ << ") "
					<< "Incorrect total number of particles at end of step "
					<< step << ": " << total_particles[0]
					<< std::endl;
			}
			abort();
		}
		if (total_particles[1] != int(std::pow(9 - step, 2))) {
			if (grid_y.get_rank() == 0) {
				std::cerr << __FILE__ << "(" << __LINE__ << ") "
					<< "Incorrect total number of particles at end of step "
					<< step << ": " << total_particles[1]
					<< std::endl;
			}
			abort();
		}
		if (total_particles[2] != int(std::pow(9 - step, 2))) {
			if (grid_z.get_rank() == 0) {
				std::cerr << __FILE__ << "(" << __LINE__ << ") "
					<< "Incorrect total number of particles at end of step "
					<< step << ": " << total_particles[2]
					<< std::endl;
			}
			abort();
		}
	}

	MPI_Finalize();

	return EXIT_SUCCESS;
}
