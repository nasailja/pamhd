/*
Program for converting MHD output of PAMHD to CDF format.

Copyright 2014, 2015, 2016 Ilja Honkonen
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
#include "string"
#include "unordered_map"
#include "vector"

#include "boost/filesystem.hpp"
#include "boost/optional.hpp"
#include "boost/program_options.hpp"
#include "cdf.h"
#include "dccrg_cartesian_geometry.hpp"
#include "dccrg_mapping.hpp"
#include "dccrg_topology.hpp"
#include "mpi.h" // must be included before gensimcell
#include "Eigen/Core" // must be included before gensimcell
#include "gensimcell.hpp"

#include "mhd/common.hpp"
#include "mhd/save.hpp"
#include "mhd/variables.hpp"


using namespace std;
using namespace pamhd::mhd;


/*
Reads simulation data from given file.

Fills out grid info and simulation data.

On success returns physical constants used by
simulation, on failure returns an uninitialized value.
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

	// read file version
	uint64_t file_version = 0;
	MPI_File_read_at( // TODO: add error checking
		file,
		offset,
		&file_version,
		1,
		MPI_UINT64_T,
		MPI_STATUS_IGNORE
	);
	offset += sizeof(uint64_t);
	if (file_version != 1) {
		cerr << "Process " << mpi_rank
			<< " Unsupported file version: " << file_version
			<< endl;
		return boost::optional<std::array<double, 4>>();
	}

	// read simulation parameters
	std::array<double, 4> simulation_parameters;
	MPI_File_read_at(
		file,
		offset,
		simulation_parameters.data(),
		4,
		MPI_DOUBLE,
		MPI_STATUS_IGNORE
	);
	offset += 4 * sizeof(double);

	// check endianness
	uint64_t endianness = 0;
	MPI_File_read_at(
		file,
		offset,
		&endianness,
		1,
		MPI_UINT64_T,
		MPI_STATUS_IGNORE
	);
	offset += sizeof(uint64_t);
	if (endianness != 0x1234567890abcdef) {
		cerr << "Process " << mpi_rank
			<< " Unsupported endianness: " << endianness
			<< ", should be " << 0x1234567890abcdef
			<< endl;
		return boost::optional<std::array<double, 4>>();
	}

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
		return boost::optional<std::array<double, 4>>(simulation_parameters);
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
	Cell::set_transfer_all(
		true,
		MHD_State_Conservative(),
		Electric_Current_Density(),
		Solver_Info(),
		MPI_Rank(),
		Resistivity(),
		Bg_Magnetic_Field_Pos_X(),
		Bg_Magnetic_Field_Pos_Y(),
		Bg_Magnetic_Field_Pos_Z()
	);
	for (const auto& item: cells_offsets) {
		const uint64_t
			cell_id = item.first,
			file_address = item.second;

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

	MPI_File_close(&file);

	return boost::optional<std::array<double, 3>>(simulation_parameters);
}


std::string get_cdf_error(CDFstatus& status)
{
	char error[CDF_STATUSTEXT_LEN + 1] = {'\0'};
	CDFgetStatusText(status, error);
	return std::string(error);
}


/*
Writes given data in cdf format to given file appended with .cdf.
*/
void convert(
	const dccrg::Cartesian_Geometry& geometry,
	const unordered_map<uint64_t, Cell>& simulation_data,
	const std::string& output_file_name_prefix,
	const double adiabatic_index,
	const double vacuum_permeability
) {
	const MHD_State_Conservative MHD{};

	// write each variable in identical order
	std::vector<uint64_t> cells;
	for (const auto& i: simulation_data) {
		cells.push_back(i.first);
	}

	CDFid cdf_id;
	CDFstatus status;

	std::vector<char> mutable_output_name(
		output_file_name_prefix.cbegin(),
		output_file_name_prefix.cend()
	);
	mutable_output_name.push_back('\0');
	status = CDFcreate(
		mutable_output_name.data(),
		0,
		nullptr,
		HOST_ENCODING,
		ROW_MAJOR,
		&cdf_id
	);
	if (status != CDF_OK) {
		std::cerr
			<< "CDF could not create output file "
			<< output_file_name_prefix << "(.cdf?): "
			<< get_cdf_error(status)
			<< std::endl;
		return;
	}

	/*
	Write global parameters as attributes
	*/
	long int attr_id;

	// adiabatic index
	status = CDFattrCreate(
		cdf_id,
		"adiabatic index",
		GLOBAL_SCOPE,
		(void*) &attr_id
	);
	if (status != CDF_OK) {
		std::cerr
			<< "Could not create attribute for adiabatic index in output file "
			<< output_file_name_prefix << "(.cdf?): "
			<< get_cdf_error(status)
			<< std::endl;
		abort();
	}

	status = CDFattrPut(
		cdf_id,
		attr_id,
		0,
		CDF_DOUBLE,
		1,
		(void*) &adiabatic_index
	);
	if (status != CDF_OK) {
		std::cerr
			<< "Could not write adiabatic index to output file "
			<< output_file_name_prefix << "(.cdf?): "
			<< get_cdf_error(status)
			<< std::endl;
		abort();
	}

	// vacuum permeability
	status = CDFattrCreate(
		cdf_id,
		"vacuum permeability",
		GLOBAL_SCOPE,
		&attr_id
	);
	if (status != CDF_OK) {
		std::cerr
			<< "Could not create attribute for vacuum permeability in output file "
			<< output_file_name_prefix << "(.cdf?): "
			<< get_cdf_error(status)
			<< std::endl;
		abort();
	}

	status = CDFattrPut(
		cdf_id,
		attr_id,
		0,
		CDF_DOUBLE,
		1,
		(void*) &vacuum_permeability
	);
	if (status != CDF_OK) {
		std::cerr
			<< "Could not write vacuum permeability to output file "
			<< output_file_name_prefix << "(.cdf?): "
			<< get_cdf_error(status)
			<< std::endl;
		abort();
	}

	/*
	Write cell data
	*/

	// data of currently written variable
	std::vector<double> data;
	long int var_id;
	std::array<long int, 1>
		vector_dim_sizes{{3}},
		variances{{VARY}},
		record_intervals{{1}},
		indices{{0}};

	// cell length
	status = CDFcreatezVar(
		cdf_id,
		"cell length",
		CDF_DOUBLE,
		1,
		1,
		vector_dim_sizes.data(),
		VARY,
		variances.data(),
		&var_id
	);
	if (status != CDF_OK) {
		std::cerr
			<< "Could not create cell length variable in output file "
			<< output_file_name_prefix << "(.cdf?): "
			<< get_cdf_error(status)
			<< std::endl;
		abort();
	}

	for (const auto& cell: cells) {
		const auto cell_length = geometry.get_length(cell);
		data.push_back(cell_length[0]);
		data.push_back(cell_length[1]);
		data.push_back(cell_length[2]);
	}
	status = CDFhyperPutzVarData(
		cdf_id,
		var_id,
		0,
		data.size() / 3,
		1,
		indices.data(),
		vector_dim_sizes.data(),
		record_intervals.data(),
		data.data()
	);
	data.clear();
	if (status != CDF_OK) {
		std::cerr
			<< "Could not write cell length variable to output file "
			<< output_file_name_prefix << "(.cdf?): "
			<< get_cdf_error(status)
			<< std::endl;
		abort();
	}


	// cell center
	status = CDFcreatezVar(
		cdf_id,
		"cell center",
		CDF_DOUBLE,
		1,
		1,
		vector_dim_sizes.data(),
		VARY,
		variances.data(),
		&var_id
	);
	if (status != CDF_OK) {
		std::cerr
			<< "Could not create cell center variable in output file "
			<< output_file_name_prefix << "(.cdf?): "
			<< get_cdf_error(status)
			<< std::endl;
		abort();
	}

	for (const auto& cell: cells) {
		const auto cell_center = geometry.get_center(cell);
		data.push_back(cell_center[0]);
		data.push_back(cell_center[1]);
		data.push_back(cell_center[2]);
	}
	status = CDFhyperPutzVarData(
		cdf_id,
		var_id,
		0,
		data.size() / 3,
		1,
		indices.data(),
		vector_dim_sizes.data(),
		record_intervals.data(),
		data.data()
	);
	data.clear();
	if (status != CDF_OK) {
		std::cerr
			<< "Could not write cell center variable to output file "
			<< output_file_name_prefix << "(.cdf?): "
			<< get_cdf_error(status)
			<< std::endl;
		abort();
	}


	// mass density
	status = CDFcreatezVar(
		cdf_id,
		"mass density",
		CDF_DOUBLE,
		1,
		0,
		nullptr,
		VARY,
		variances.data(),
		&var_id
	);
	if (status != CDF_OK) {
		std::cerr
			<< "Could not create mass density variable in output file "
			<< output_file_name_prefix << "(.cdf?): "
			<< get_cdf_error(status)
			<< std::endl;
		abort();
	}

	for (const auto& cell: cells) {
		data.push_back(simulation_data.at(cell)[MHD][Mass_Density()]);
	}
	status = CDFhyperPutzVarData(
		cdf_id,
		var_id,
		0,
		data.size(),
		1,
		indices.data(),
		nullptr,
		record_intervals.data(),
		data.data()
	);
	data.clear();
	if (status != CDF_OK) {
		std::cerr
			<< "Could not write mass density variable to output file "
			<< output_file_name_prefix << "(.cdf?): "
			<< get_cdf_error(status)
			<< std::endl;
		abort();
	}

	// momentum
	status = CDFcreatezVar(
		cdf_id,
		"momentum density",
		CDF_DOUBLE,
		1,
		1,
		vector_dim_sizes.data(),
		VARY,
		variances.data(),
		&var_id
	);
	if (status != CDF_OK) {
		std::cerr
			<< "Could not create momentum density variable in output file "
			<< output_file_name_prefix << "(.cdf?): "
			<< get_cdf_error(status)
			<< std::endl;
		abort();
	}

	for (const auto& cell: cells) {
		const auto& mom = simulation_data.at(cell)[MHD][Momentum_Density()];
		data.push_back(mom[0]);
		data.push_back(mom[1]);
		data.push_back(mom[2]);
	}
	status = CDFhyperPutzVarData(
		cdf_id,
		var_id,
		0,
		data.size() / 3,
		1,
		indices.data(),
		vector_dim_sizes.data(),
		record_intervals.data(),
		data.data()
	);
	data.clear();
	if (status != CDF_OK) {
		std::cerr
			<< "Could not write momentum density variable to output file "
			<< output_file_name_prefix << "(.cdf?): "
			<< get_cdf_error(status)
			<< std::endl;
		abort();
	}

	// total energy
	status = CDFcreatezVar(
		cdf_id,
		"total energy density",
		CDF_DOUBLE,
		1,
		0,
		nullptr,
		VARY,
		variances.data(),
		&var_id
	);
	if (status != CDF_OK) {
		std::cerr
			<< "Could not create total energy density variable in output file "
			<< output_file_name_prefix << "(.cdf?): "
			<< get_cdf_error(status)
			<< std::endl;
		abort();
	}

	for (const auto& cell: cells) {
		data.push_back(simulation_data.at(cell)[MHD][Total_Energy_Density()]);
	}
	status = CDFhyperPutzVarData(
		cdf_id,
		var_id,
		0,
		data.size(),
		1,
		indices.data(),
		nullptr,
		record_intervals.data(),
		data.data()
	);
	data.clear();
	if (status != CDF_OK) {
		std::cerr
			<< "Could not write total energy density variable to output file "
			<< output_file_name_prefix << "(.cdf?): "
			<< get_cdf_error(status)
			<< std::endl;
		abort();
	}

	// magnetic field
	status = CDFcreatezVar(
		cdf_id,
		"magnetic field",
		CDF_DOUBLE,
		1,
		1,
		vector_dim_sizes.data(),
		VARY,
		variances.data(),
		&var_id
	);
	if (status != CDF_OK) {
		std::cerr
			<< "Could not create magnetic field variable in output file "
			<< output_file_name_prefix << "(.cdf?): "
			<< get_cdf_error(status)
			<< std::endl;
		abort();
	}

	for (const auto& cell: cells) {
		const auto& mag = simulation_data.at(cell)[MHD][Magnetic_Field()];
		data.push_back(mag[0]);
		data.push_back(mag[1]);
		data.push_back(mag[2]);
	}
	status = CDFhyperPutzVarData(
		cdf_id,
		var_id,
		0,
		data.size() / 3,
		1,
		indices.data(),
		vector_dim_sizes.data(),
		record_intervals.data(),
		data.data()
	);
	data.clear();
	if (status != CDF_OK) {
		std::cerr
			<< "Could not write magnetic field variable to output file "
			<< output_file_name_prefix << "(.cdf?): "
			<< get_cdf_error(status)
			<< std::endl;
		abort();
	}

	status = CDFcloseCDF(cdf_id);
	if (status != CDF_OK) {
		std::cerr
			<< "Could not close output file "
			<< output_file_name_prefix << "(.cdf?)"
			<< std::endl;
		abort();
	}
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

	boost::program_options::options_description
		options("Usage: program_name [options], where options are");
	options.add_options()
		("help", "Print this help message")
		("verbose", "Print run time information")
		("input-file",
			boost::program_options::value<std::vector<std::string>>(&input_files),
			"Results to plot, the --input-file part can be omitted");

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
			cout
				<< "Converts given .dc files from MHD test program to CDF format\n"
				<< options
				<< endl;
		}
		MPI_Finalize();
		return EXIT_SUCCESS;
	}

	if (input_files.size() == 0) {
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

		boost::optional<std::array<double, 4> header
			= read_data(
				cell_id_mapping,
				topology,
				geometry,
				simulation_data,
				input_files[i],
				rank
			);
		if (not header) {
			std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
				<< "Couldn't read simulation data from file " << input_files[i]
				<< std::endl;
			continue;
		}

		// cdf doesn't support overwriting files
		const std::string output_prefix(
			input_files[i].substr(0, input_files[i].size() - 3)
		);
		boost::filesystem::remove(output_prefix + ".cdf");

		convert(
			geometry,
			simulation_data,
			output_prefix,
			(*header)[1],
			(*header)[3]
		);
	}

	MPI_Finalize();

	return EXIT_SUCCESS;
}
