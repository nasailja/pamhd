#! /usr/bin/env python3
'''
Extracts subvolumes of MHD output of PAMHD.

Copyright 2016 Ilja Honkonen
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


If called as the main program...
'''

'''
Reads simulation data from infile_name and writes it into outfile_name while excluding cells whose center is outside of given volume
'''
def extract(
	infile_name,
	outfile_name,
	min_x = float('-inf'),
	max_x = float('+inf'),
	min_y = float('-inf'),
	max_y = float('+inf'),
	min_z = float('-inf'),
	max_z = float('+inf')
):
	from os.path import exists, isfile
	from numpy import dtype, fromfile, int32, uint64, zeros

	if not exists(infile_name):
		exit('Path ' + infile_name + " doesn't exist")
	if not isfile(infile_name):
		exit('Path ' + infile_name + ' is not a file')

	with open(infile_name, 'rb') as infile:
		with open(outfile_name, 'wb') as outfile:

			file_version = fromfile(infile, dtype = 'uint64', count = 1)
			if file_version[0] > 1:
				exit('File ' + infile_name + ' is of unsupported version: ' + str(file_version[0]))
			sim_params = fromfile(infile, dtype = '4double', count = 1)
			endianness = fromfile(infile, dtype = 'uint64', count = 1)
			ref_lvl_0_cells = fromfile(infile, dtype = '3uint64', count = 1)
			max_ref_lvl = fromfile(infile, dtype = 'intc', count = 1)
			if max_ref_lvl[0] > 0:
				exit('Refinement level > 0 not supported')
			neighborhood_length = fromfile(infile, dtype = 'uintc', count = 1)
			periodicity = fromfile(infile, dtype = '3uint8', count = 1)
			geometry_id = fromfile(infile, dtype = 'intc', count = 1)
			if geometry_id[0] != int32(1):
				exit('Unsupported geometry: ' + str(geometry_id[0]))
			grid_start = fromfile(infile, dtype = '3double', count = 1)
			lvl_0_cell_length = fromfile(infile, dtype = '3double', count = 1)
			total_cells = fromfile(infile, dtype = 'uint64', count = 1)
			cell_ids_data_offsets = fromfile(infile, dtype = '2uint64', count = total_cells[0])

			# check which cells are inside requested volume
			max_x_i, max_y_i, max_z_i = 0, 0, 0
			min_x_i, min_y_i, min_z_i = ref_lvl_0_cells[0][0], ref_lvl_0_cells[0][1], ref_lvl_0_cells[0][2]
			out_cells = []
			out_cell_data = []
			cell_data_t = 'double, 3double, double, 3double, 3double, intc, intc, double, 3double, 3double, 3double'
			for item in cell_ids_data_offsets:
				cell_id = uint64(item[0])

				cell_id -= 1
				cell_index = (
					uint64(cell_id % ref_lvl_0_cells[0][0]),
					uint64(cell_id / ref_lvl_0_cells[0][0] % ref_lvl_0_cells[0][1]),
					uint64(cell_id / (ref_lvl_0_cells[0][0] * ref_lvl_0_cells[0][1]))
				)

				cell_center = (
					grid_start[0][0] + lvl_0_cell_length[0][0] * (0.5 + cell_index[0]),
					grid_start[0][1] + lvl_0_cell_length[0][1] * (0.5 + cell_index[1]),
					grid_start[0][2] + lvl_0_cell_length[0][2] * (0.5 + cell_index[2])
				)
				cell_id += 1

				if cell_center[0] < min_x or cell_center[0] > max_x or cell_center[1] < min_y or cell_center[1] > max_y or cell_center[2] < min_z or cell_center[2] > max_z:
					continue

				min_x_i = min(min_x_i, cell_index[0])
				min_y_i = min(min_y_i, cell_index[1])
				min_z_i = min(min_z_i, cell_index[2])
				max_x_i = max(max_x_i, cell_index[0])
				max_y_i = max(max_y_i, cell_index[1])
				max_z_i = max(max_z_i, cell_index[2])

				out_cells.append(cell_id)

				infile.seek(item[1], 0)
				out_cell_data.append(
					fromfile(
						infile,
						dtype = cell_data_t,
						count = 1
					)
				)

			file_version.tofile(outfile)
			sim_params.tofile(outfile)
			endianness.tofile(outfile)

			# write grid parameters corresponding to filtered grid
			new_ref_lvl_0_cells = ref_lvl_0_cells.copy()
			new_ref_lvl_0_cells[0][0] = max_x_i - min_x_i + 1
			new_ref_lvl_0_cells[0][1] = max_y_i - min_y_i + 1
			new_ref_lvl_0_cells[0][2] = max_z_i - min_z_i + 1
			new_ref_lvl_0_cells.tofile(outfile)

			max_ref_lvl.tofile(outfile)
			neighborhood_length.tofile(outfile)
			periodicity.tofile(outfile)
			geometry_id.tofile(outfile)

			new_grid_start = grid_start.copy()
			new_grid_start[0][0] += lvl_0_cell_length[0][0] * min_x_i
			new_grid_start[0][1] += lvl_0_cell_length[0][1] * min_y_i
			new_grid_start[0][2] += lvl_0_cell_length[0][2] * min_z_i
			new_grid_start.tofile(outfile)

			lvl_0_cell_length.tofile(outfile)

			new_total_cells = total_cells.copy()
			new_total_cells[0] = len(out_cells)
			new_total_cells.tofile(outfile)

			cell_data_start = outfile.tell()
			out_cell_ids_offsets = zeros(len(out_cells), dtype = '2uint64')
			for i in range(len(out_cells)):

				old_id = out_cells[i] - 1
				old_index = (
					uint64(old_id % ref_lvl_0_cells[0][0]),
					uint64(old_id / ref_lvl_0_cells[0][0] % ref_lvl_0_cells[0][1]),
					uint64(old_id / (ref_lvl_0_cells[0][0] * ref_lvl_0_cells[0][1]))
				)
				new_index = (
					old_index[0] - min_x_i,
					old_index[1] - min_y_i,
					old_index[2] - min_z_i
				)
				new_id = uint64(
					new_index[2] * new_ref_lvl_0_cells[0][0] * new_ref_lvl_0_cells[0][1]
					+ new_index[1] * new_ref_lvl_0_cells[0][0]
					+ new_index[0]
				)
				new_id += 1

				out_cell_ids_offsets[i][0] = new_id
				out_cell_ids_offsets[i][1] = cell_data_start + out_cell_ids_offsets.nbytes + dtype(cell_data_t).itemsize * i

			out_cell_ids_offsets.tofile(outfile)
			for d in out_cell_data:
				d.tofile(outfile)


if __name__ == '__main__':
	from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
	from sys import argv

	parser = ArgumentParser(
		formatter_class = ArgumentDefaultsHelpFormatter
	)
	parser.add_argument('files', metavar = 'F', nargs = '*', help = 'Names of input and output files to process (in1 out1 in2 out2 ...)')
	parser.add_argument('--min-x', type = float, default = float('-inf'), help = 'Include in output file(s) only cells with center x coordinate larger than MIN_X')
	parser.add_argument('--max-x', type = float, default = float('+inf'), help = 'Include in output file(s) only cells with center x coordinate smaller than MAX_X')
	parser.add_argument('--min-y', type = float, default = float('-inf'), help = 'Include in output file(s) only cells with center y coordinate larger than MIN_Y')
	parser.add_argument('--max-y', type = float, default = float('+inf'), help = 'Include in output file(s) only cells with center y coordinate smaller than MAX_Y')
	parser.add_argument('--min-z', type = float, default = float('-inf'), help = 'Include in output file(s) only cells with center z coordinate larger than MIN_Z')
	parser.add_argument('--max-z', type = float, default = float('+inf'), help = 'Include in output file(s) only cells with center z coordinate smaller than MAX_Z')
	args = parser.parse_args()

	if len(args.files) % 2 > 0:
		exit('Even number of file arguments required (in1 out1 in2 out2 ...)')

	if args.min_x >= args.max_x:
		exit('Minimum x coordinate must be less than maximum')
	if args.min_y >= args.max_y:
		exit('Minimum y coordinate must be less than maximum')
	if args.min_z >= args.max_z:
		exit('Minimum z coordinate must be less than maximum')

	for i in range(0, len(args.files), 2):
		extract(
			args.files[i],
			args.files[i + 1],
			min_x = args.min_x,
			max_x = args.max_x,
			min_y = args.min_y,
			max_y = args.max_y,
			min_z = args.min_z,
			max_z = args.max_z
		)
