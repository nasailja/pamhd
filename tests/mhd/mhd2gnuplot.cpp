/*
Program for plotting MHD output of PAMHD with gnuplot.

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

#include "boost/lexical_cast.hpp"
#include "cstdint"
#include "cstdlib"
#include "Eigen/Core" // must be included before gensimcell
#include "fstream"
#include "mpi.h"
#include "string"
#include "tuple"
#include "unordered_map"
#include "vector"

#include "dccrg_cartesian_geometry.hpp"
#include "dccrg_mapping.hpp"
#include "dccrg_topology.hpp"
#include "gensimcell.hpp"

#include "mhd/common.hpp"
#include "mhd/variables.hpp"

using namespace std;
using namespace pamhd::mhd;

int main(int argc, char* argv[])
{
	using Cell_T = MHD_Conservative;

	constexpr double
		adiabatic_index = 5.0 / 3.0,    
		vacuum_permeability = 4e-7 * M_PI;

	if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
		cerr << "Coudln't initialize MPI." << endl;
		abort();
	}

	MPI_Comm comm = MPI_COMM_WORLD;

	int rank = 0, comm_size = 0;
	MPI_Comm_rank(comm, &rank);
	MPI_Comm_size(comm, &comm_size);


	dccrg::Mapping mapping;
	dccrg::Grid_Topology topology;
	dccrg::Cartesian_Geometry geometry(mapping.length, mapping, topology);

	for (int i = 1; i < argc; i++) {

		if ((i - 1) % comm_size != rank) {
			continue;
		}

		const string argv_string(argv[i]);

		MPI_File file;
		if (
			MPI_File_open(
				MPI_COMM_SELF,
				const_cast<char*>(argv_string.c_str()),
				MPI_MODE_RDONLY,
				MPI_INFO_NULL,
				&file
			) != MPI_SUCCESS
		) {
			cerr << "Process " << rank
				<< " couldn't open file " << argv_string
				<< endl;
			continue;
		}

		// skip endianness check data
		MPI_Offset offset = sizeof(uint64_t);

		if (not mapping.read(file, offset)) {
			cerr << "Process " << rank
				<< " couldn't set cell id mapping from file " << argv_string
				<< endl;
			continue;
		}

		offset
			+= mapping.data_size()
			+ sizeof(unsigned int)
			+ topology.data_size();

		if (not geometry.read(file, offset)) {
			cerr << "Process " << rank
				<< " couldn't read geometry from file " << argv_string
				<< endl;
			continue;
		}
		offset += geometry.data_size();

		// read number of cells
		uint64_t total_cells = 0;
		MPI_File_read_at(
			file,
			offset,
			&total_cells,
			1,
			MPI_UINT64_T,
			MPI_STATUS_IGNORE
		);
		offset += sizeof(uint64_t);

		if (total_cells == 0) {
			MPI_File_close(&file);
			continue;
		}

		// read cell ids and data offsets
		vector<pair<uint64_t, uint64_t>> cells_offsets(total_cells);
		MPI_File_read_at(
			file,
			offset,
			cells_offsets.data(),
			2 * total_cells,
			MPI_UINT64_T,
			MPI_STATUS_IGNORE
		);

		// read cell data
		unordered_map<
			uint64_t,
			Cell_T
		> simulation_data;
		Cell_T::set_transfer_all(
			true,
			Mass_Density(),
			Momentum_Density(),
			Total_Energy_Density(),
			Magnetic_Field()
		);

		/*
		Sort cells for plotting, x, y or z
		coordinate increases with cell id
		*/
		std::vector<uint64_t> sorted_cells;
		sorted_cells.reserve(cells_offsets.size());

		for (const auto& item: cells_offsets) {
			const uint64_t
				cell_id = item.first,
				file_address = item.second;

			sorted_cells.push_back(cell_id);
			simulation_data[cell_id];

			/*
			dccrg writes cell data without padding so store the
			non-padded version of the cell's datatype into file_datatype
			*/
			void* memory_address = NULL;
			int memory_count = -1;
			MPI_Datatype
				memory_datatype = MPI_DATATYPE_NULL,
				file_datatype = MPI_DATATYPE_NULL;

			tie(
				memory_address,
				memory_count,
				memory_datatype
			) = simulation_data.at(cell_id).get_mpi_datatype();

			int sizeof_memory_datatype;
			MPI_Type_size(memory_datatype, &sizeof_memory_datatype);
			MPI_Type_contiguous(sizeof_memory_datatype, MPI_BYTE, &file_datatype);

			// interpret data from the file using the non-padded type
			MPI_Type_commit(&file_datatype);
			MPI_File_set_view(
				file,
				file_address,
				MPI_BYTE,
				file_datatype,
				const_cast<char*>("native"),
				MPI_INFO_NULL
			);

			MPI_Type_commit(&memory_datatype);
			MPI_File_read_at(
				file,
				0,
				memory_address,
				memory_count,
				memory_datatype,
				MPI_STATUS_IGNORE
			);

			MPI_Type_free(&memory_datatype);
			MPI_Type_free(&file_datatype);
		}
		if (sorted_cells.size() > 0) {
			std::sort(sorted_cells.begin(), sorted_cells.end());
		}

		MPI_File_close(&file);
		cells_offsets.clear();

		// get direction of the tube, first dimension with > 1 cell
		const auto grid_length = mapping.length.get();
		const size_t tube_dim = [&](){
			for (size_t i = 0; i < grid_length.size(); i++) {
				if (grid_length[i] > 1) {
					return i;
				}
			}
			std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
				<< "Grid must have > 1 cell in one dimension"
				<< std::endl;
			abort();
		}();

		const string
			gnuplot_file_name(argv_string + ".dat"),
			plot_file_name(argv_string + ".png"),
			B_plot_file_name("B_" + argv_string + ".png"),
			v_plot_file_name("v_" + argv_string + ".png");

		ofstream gnuplot_file(gnuplot_file_name);

		const double
			tube_start = geometry.get_start()[tube_dim],
			tube_end = geometry.get_end()[tube_dim];

		// mass density & pressure
		gnuplot_file
			<< "set term png enhanced\nset output '"
			<< plot_file_name
			<< "'\nset xlabel 'x (m)'\nset xrange ["
			<< boost::lexical_cast<std::string>(tube_start)
			<< " : " << boost::lexical_cast<std::string>(tube_end)
			<< "]\nset ylabel 'Mass density (kg m^{-3})'\n"
			   "set y2label 'Pressure (Pa)'\n"
			   "set key horizontal outside bottom\n"
			   "set ytics nomirror\n"
			   "set yrange [0 : *]\n"
			   "set y2range [0 : *]\n"
			   "set y2tics auto\n"
			   "plot "
			     "'-' using 1:2 axes x1y1 with line linewidth 2 title 'density', "
			     "'-' u 1:2 axes x1y2 w l lw 2 t 'pressure'\n";

		for (const auto& cell_id: sorted_cells) {
			const auto& mhd_data = simulation_data.at(cell_id);
			const double x = geometry.get_center(cell_id)[tube_dim];
			gnuplot_file << x << " " << mhd_data[Mass_Density()] << "\n";
		}
		gnuplot_file << "end\n";

		for (const auto& cell_id: sorted_cells) {
			const auto& mhd_data = simulation_data.at(cell_id);
			const double x = geometry.get_center(cell_id)[tube_dim];
			gnuplot_file
				<< x << " "
				<< get_pressure<
					Cell_T,
					Mass_Density,
					Momentum_Density,
					Total_Energy_Density,
					Magnetic_Field
				>(mhd_data, adiabatic_index, vacuum_permeability)
				<< "\n";
		}
		gnuplot_file << "end\n";

		// velocity
		gnuplot_file
			<< "set term png enhanced\nset output '"
			<< v_plot_file_name
			<< "'\nset ylabel 'Velocity (m s^{-1})'\n"
			   "set yrange [* : *]\n"
			   "unset y2label\n"
			   "unset y2tics\n"
			   "set ytics mirror\n"
			   "plot "
			     "'-' u 1:2 w l lw 2 t 'v_x', "
			     "'-' u 1:2 w l lw 2 t 'v_y', "
			     "'-' u 1:2 w l lw 2 t 'v_z'\n";

		for (const auto& cell_id: sorted_cells) {
			const auto& mhd_data = simulation_data.at(cell_id);
			const auto& m = mhd_data[Momentum_Density()];
			const auto& rho = mhd_data[Mass_Density()];
			const auto v(m / rho);
			const double x = geometry.get_center(cell_id)[tube_dim];
			gnuplot_file << x << " " << v[0] << "\n";
		}
		gnuplot_file << "end\n";

		for (const auto& cell_id: sorted_cells) {
			const auto& mhd_data = simulation_data.at(cell_id);
			const auto& m = mhd_data[Momentum_Density()];
			const auto& rho = mhd_data[Mass_Density()];
			const auto v(m / rho);
			const double x = geometry.get_center(cell_id)[tube_dim];
			gnuplot_file << x << " " << v[1] << "\n";
		}
		gnuplot_file << "end\n";

		for (const auto& cell_id: sorted_cells) {
			const auto& mhd_data = simulation_data.at(cell_id);
			const auto& m = mhd_data[Momentum_Density()];
			const auto& rho = mhd_data[Mass_Density()];
			const auto v(m / rho);
			const double x = geometry.get_center(cell_id)[tube_dim];
			gnuplot_file << x << " " << v[2] << "\n";
		}
		gnuplot_file << "end\n";

		// magnetic field
		gnuplot_file
			<< "set term png enhanced\nset output '"
			<< B_plot_file_name
			<< "'\nset ylabel 'Magnetic field (T)'\n"
			   "plot "
			     "'-' u 1:2 w l lw 2 t 'B_x', "
			     "'-' u 1:2 w l lw 2 t 'B_y', "
			     "'-' u 1:2 w l lw 2 t 'B_z'\n";

		for (const auto& cell_id: sorted_cells) {
			const auto& mhd_data = simulation_data.at(cell_id);
			const auto& B = mhd_data[Magnetic_Field()];
			const double x = geometry.get_center(cell_id)[tube_dim];
			gnuplot_file << x << " " << B[0] << "\n";
		}
		gnuplot_file << "end\n";

		for (const auto& cell_id: sorted_cells) {
			const auto& mhd_data = simulation_data.at(cell_id);
			const auto& B = mhd_data[Magnetic_Field()];
			const double x = geometry.get_center(cell_id)[tube_dim];
			gnuplot_file << x << " " << B[1] << "\n";
		}
		gnuplot_file << "end\n";

		for (const auto& cell_id: sorted_cells) {
			const auto& mhd_data = simulation_data.at(cell_id);
			const auto& B = mhd_data[Magnetic_Field()];
			const double x = geometry.get_center(cell_id)[tube_dim];
			gnuplot_file << x << " " << B[2] << "\n";
		}
		gnuplot_file << "end\n";

		gnuplot_file.close();

		system(("gnuplot " + gnuplot_file_name).c_str());
	}

	MPI_Finalize();

	return EXIT_SUCCESS;
}
