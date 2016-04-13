/*
Program for plotting bulk data of particle test of PAMHD with gnuplot.

Copyright 2015, 2016 Ilja Honkonen
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

#include "array"
#include "cstdint"
#include "cstdlib"
#include "functional"
#include "fstream"
#include "limits"
#include "string"
#include "tuple"
#include "unordered_map"
#include "vector"

#include "boost/lexical_cast.hpp"
#include "boost/optional.hpp"
#include "boost/program_options.hpp"
#include "dccrg_cartesian_geometry.hpp"
#include "dccrg_mapping.hpp"
#include "dccrg_topology.hpp"
#include "mpi.h" // must be included before gensimcell
#include "Eigen/Core" // must be included before gensimcell
#include "gensimcell.hpp"
#include "prettyprint.hpp"

#include "mhd/variables.hpp"
#include "particle/save.hpp"
#include "particle/variables.hpp"

using namespace std;
using namespace pamhd::particle;

/*
Reads simulation data from given file.

Fills out grid info and simulation data withing given volume.

On success returns vacuum permeability.
*/
boost::optional<std::array<double, 4>> read_data(
	dccrg::Mapping& cell_id_mapping,
	dccrg::Grid_Topology& topology,
	dccrg::Cartesian_Geometry& geometry,
	unordered_map<uint64_t, Cell>& simulation_data,
	const std::string& file_name,
	const int mpi_rank
) {
	MPI_File file;
	if (
		MPI_File_open(
			MPI_COMM_SELF,
			const_cast<char*>(file_name.c_str()),
			MPI_MODE_RDONLY,
			MPI_INFO_NULL,
			&file
		) != MPI_SUCCESS
	) {
		cerr << "Process " << mpi_rank
			<< " couldn't open file " << file_name
			<< endl;
		return boost::optional<std::array<double, 4>>();
	}

	MPI_Offset offset = 0;

	// read physical constants
	std::array<double, 4> metadata;
	MPI_File_read_at(
		file,
		offset,
		(void*) metadata.data(),
		metadata.size(),
		MPI_DOUBLE,
		MPI_STATUS_IGNORE
	);
	offset += sizeof(double) * metadata.size();

	// skip endianness check data
	offset += sizeof(uint64_t);

	if (not cell_id_mapping.read(file, offset)) {
		cerr << "Process " << mpi_rank
			<< " couldn't set cell id mapping from file " << file_name
			<< endl;
		return boost::optional<std::array<double, 4>>();
	}

	offset
		+= cell_id_mapping.data_size()
		+ sizeof(unsigned int)
		+ topology.data_size();

	if (not geometry.read(file, offset)) {
		cerr << "Process " << mpi_rank
			<< " couldn't read geometry from file " << file_name
			<< endl;
		return boost::optional<std::array<double, 4>>();
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
		return boost::optional<std::array<double, 4>>(metadata);
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
	for (const auto& item: cells_offsets) {
		uint64_t
			cell_id = item.first,
			file_address = item.second;

		simulation_data[cell_id];
		auto& cell_data = simulation_data.at(cell_id);

		void* memory_address = NULL;
		int memory_count = -1;
		MPI_Datatype
			memory_datatype = MPI_DATATYPE_NULL,
			file_datatype = MPI_DATATYPE_NULL;

		Cell::set_transfer_all(
			true,
			Electric_Field(),
			Magnetic_Field(),
			pamhd::mhd::Electric_Current_Density()
		);
		tie(
			memory_address,
			memory_count,
			memory_datatype
		) = cell_data.get_mpi_datatype();

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

	MPI_File_close(&file);

	return boost::optional<std::array<double, 4>>(metadata);
}


/*
Plots 1d data from given list in that order.
*/
int plot_1d(
	const dccrg::Cartesian_Geometry& geometry,
	const unordered_map<uint64_t, Cell>& simulation_data,
	const std::vector<uint64_t>& cells,
	const std::string& output_file_name_prefix
) {
	const string gnuplot_file_name(output_file_name_prefix + ".dat");
	ofstream gnuplot_file(gnuplot_file_name);

	const size_t tube_dim
		= [&](){
			for (size_t i = 0; i < geometry.length.get().size(); i++) {
				if (geometry.length.get()[i] > 1) {
					return i;
				}
			}
			return size_t(999);
		}();

	const double
		tube_start = geometry.get_start()[tube_dim],
		tube_end = geometry.get_end()[tube_dim];

	// magnetic field
	gnuplot_file
		<< "set term png enhanced\nset output '"
		<< output_file_name_prefix + "_B.png"
		<< "'\nset ylabel 'Magnetic field'\n"
		   "plot "
		     "'-' u 1:2 w l lw 2 t 'B_1', "
		     "'-' u 1:2 w l lw 2 t 'B_2', "
		     "'-' u 1:2 w l lw 2 t 'B_3'\n";

	for (const auto& cell_id: cells) {
		const auto& B = simulation_data.at(cell_id)[Magnetic_Field()];
		const double x = geometry.get_center(cell_id)[tube_dim];
		gnuplot_file << x << " " << B[0] << "\n";
	}
	gnuplot_file << "end\n";

	for (const auto& cell_id: cells) {
		const auto& B = simulation_data.at(cell_id)[Magnetic_Field()];
		const double x = geometry.get_center(cell_id)[tube_dim];
		gnuplot_file << x << " " << B[1] << "\n";
	}
	gnuplot_file << "end\n";

	for (const auto& cell_id: cells) {
		const auto& B = simulation_data.at(cell_id)[Magnetic_Field()];
		const double x = geometry.get_center(cell_id)[tube_dim];
		gnuplot_file << x << " " << B[2] << "\n";
	}
	gnuplot_file << "end\n";

	// electric field
	gnuplot_file
		<< "set term png enhanced\nset output '"
		<< output_file_name_prefix + "_E.png"
		<< "'\nset ylabel 'Electric field'\n"
		   "plot "
		     "'-' u 1:2 w l lw 2 t 'E_1', "
		     "'-' u 1:2 w l lw 2 t 'E_2', "
		     "'-' u 1:2 w l lw 2 t 'E_3'\n";

	for (const auto& cell_id: cells) {
		const auto& E = simulation_data.at(cell_id)[Electric_Field()];
		const double x = geometry.get_center(cell_id)[tube_dim];
		gnuplot_file << x << " " << E[0] << "\n";
	}
	gnuplot_file << "end\n";

	for (const auto& cell_id: cells) {
		const auto& E = simulation_data.at(cell_id)[Electric_Field()];
		const double x = geometry.get_center(cell_id)[tube_dim];
		gnuplot_file << x << " " << E[1] << "\n";
	}
	gnuplot_file << "end\n";

	for (const auto& cell_id: cells) {
		const auto& E = simulation_data.at(cell_id)[Electric_Field()];
		const double x = geometry.get_center(cell_id)[tube_dim];
		gnuplot_file << x << " " << E[2] << "\n";
	}
	gnuplot_file << "end\n";


	gnuplot_file.close();

	return system(("gnuplot '" + gnuplot_file_name + "'").c_str());
}


int main(int argc, char* argv[])
{
	if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
		cerr << "Coudln't initialize MPI." << endl;
		abort();
	}

	MPI_Comm comm = MPI_COMM_WORLD;

	int rank = 0, comm_size = 0;
	MPI_Comm_rank(comm, &rank);
	MPI_Comm_size(comm, &comm_size);

	// program options
	bool verbose = false;
	std::vector<std::string> input_files;
	std::string
		plot_command_prefix(
			"set term png enhanced size 800, 600\n"
			"set pal gray\n"
			"set format cb \"%.2e\""
		);

	boost::program_options::options_description
		options("Usage: program_name [options], where options are");
	options.add_options()
		("help", "Print this help message")
		("verbose", "Print run time information")
		("input-file",
			boost::program_options::value<>(&input_files)
				->default_value(input_files),
			"Results to plot, the --input-file part can be omitted")
		("common-2d",
			boost::program_options::value<std::string>(&plot_command_prefix)
				->default_value(plot_command_prefix),
			"String to prepend to gnuplot file");

	boost::program_options::positional_options_description positional_options;
	positional_options.add("input-file", -1);

	boost::program_options::variables_map var_map;
	try {
		boost::program_options::store(
			boost::program_options::command_line_parser(argc, argv)
				.options(options)
				.positional(positional_options)
				.run(),
			var_map
		);
	} catch (std::exception& e) {
		if (rank == 0) {
			std::cerr <<  __FILE__ << "(" << __LINE__ << "): "
				<< "Couldn't parse command line options: " << e.what()
				<< std::endl;
		}
		abort();
	}
	boost::program_options::notify(var_map);

	if (var_map.count("help") > 0) {
		if (rank == 0) {
			cout << options <<
				"To include newlines in arguments use e.g. in bash "
				"./mhd2gnuplot --opt $'set title \"T\"\\nset xrange...' "
				"which expands \\n etc."
				<< endl;
		}
		MPI_Finalize();
		return EXIT_SUCCESS;
	}

	if (var_map.count("verbose") > 0) {
		verbose = true;
	}

	for (size_t i = 0; i < input_files.size(); i++) {

		if (int(i) % comm_size != rank) {
			continue;
		}

		if (verbose) {
			std::cout << "MPI rank " << rank
				<< " processing file " << input_files[i]
				<< std::endl;
		}

		dccrg::Mapping cell_id_mapping;
		dccrg::Grid_Topology topology;
		dccrg::Cartesian_Geometry geometry(cell_id_mapping.length, cell_id_mapping, topology);
		unordered_map<uint64_t, Cell> simulation_data;

		boost::optional<std::array<double, 4>> metadata = read_data(
			cell_id_mapping,
			topology,
			geometry,
			simulation_data,
			input_files[i],
			rank
		);
		if (not metadata) {
			std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
				<< "Couldn't read simulation data from file " << input_files[i]
				<< std::endl;
			continue;
		}

		std::vector<uint64_t> cell_ids;
		cell_ids.reserve(simulation_data.size());
		for (const auto& item: simulation_data) {
			cell_ids.push_back(item.first);
		}
		sort(cell_ids.begin(), cell_ids.end());

		plot_1d(
			geometry,
			simulation_data,
			cell_ids,
			input_files[i]
		);
	}

	MPI_Finalize();

	return EXIT_SUCCESS;
}
