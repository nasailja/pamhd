#! /usr/bin/env python3
'''
Function for reading output of PAMHD particle test program.

Copyright 2015 Ilja Honkonen
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


If called as the main program loads particle data from first given file into
dictionary particle_data and appends to it data from any subsequent files.
Requires NUMPY. Example listing the id of every particle in one file:
python -i particle2numpy.py particle_0.000e+00_s.dc
for key in particle_data:
	print(key)
exit()
'''

import os

try:
	import numpy
except:
	exit("Couldn't import numpy")


'''
Adds particle data produced by massless.exe of PAMHD
from given file into given dictionary.

Already existing particles' data is appended to from given file.
'''
def load(file_name, particle_data):
	import binascii

	if not os.path.isfile(file_name):
		raise Exception('Given file name (' + file_name + ')is not a file.')

	infile = open(file_name, 'rb')

	# read simulation header
	vacuum_permeability = numpy.fromfile(infile, dtype = numpy.float64, count = 1)[0]

	# check file endiannes
	endianness = numpy.fromfile(infile, dtype = 'uint64', count = 1)[0]
	if endianness != numpy.uint64(0x1234567890abcdef):
		raise Exception(
			'Wrong endiannes in given file, expected ' + str(hex(0x1234567890abcdef)) \
			+ ' but got ' + str(hex(endianness))
		)

	# number of refinement level 0 cells in each dimension
	ref_lvl_0_cells = numpy.fromfile(infile, dtype = '3uint64', count = 1)[0]
	#print(ref_lvl_0_cells)

	# maximum refinement level of grid cells
	max_ref_lvl = numpy.fromfile(infile, dtype = 'int32', count = 1)[0]
	if max_ref_lvl > numpy.uint32(0):
		raise Exception('Refinement level > 0 not supported')

	# length of every cells' neighborhood in cells of identical size
	neighborhood_length = numpy.fromfile(infile, dtype = 'uint32', count = 1)[0]
	#print(neighborhood_length)

	# whether grid is periodic each dimension (0 == no, 1 == yes)
	periodicity = numpy.fromfile(infile, dtype = '3uint8', count = 1)[0]
	#print(periodicity)

	geometry_id = numpy.fromfile(infile, dtype = 'int32', count = 1)[0]
	if geometry_id != numpy.int32(1):
		raise Exception('Unsupported geometry')

	# starting coordinate of grid
	grid_start = numpy.fromfile(infile, dtype = '3float64', count = 1)[0]
	#print(grid_start)

	# length of cells of refinement level 0
	lvl_0_cell_length = numpy.fromfile(infile, dtype = '3float64', count = 1)[0]
	#print(lvl_0_cell_length)

	# total number of cells in grid
	total_cells = numpy.fromfile(infile, dtype = 'uint64', count = 1)[0]
	#print(total_cells)

	# id of each cell and offset in bytes to data of each cell
	cell_ids_data_offsets = numpy.fromfile(infile, dtype = '2uint64', count = total_cells)
	#print(cell_ids_data_offsets)

	# read particle data
	for item in cell_ids_data_offsets:
		infile.seek(item[1], 0)
		bulk_data = numpy.fromfile(
			infile,
			dtype = '3float, 3float64, uint64',
			count = 1
		)[0]
		temp_particle_data = numpy.fromfile(
			infile,
			dtype = '3float64, 3float64, float64, float64, uint64',
			count = bulk_data[2]
		)
		for particle in temp_particle_data:
			position, velocity, mass, charge_mass_ratio, particle_id \
				= particle[0], particle[1], particle[2], particle[3], particle[4]

			if particle_id in particle_data:
				particle_data[particle_id].append(
					(position, velocity)
				)
			else:
				particle_data[particle_id] = [
					(position, velocity)
				]


if __name__ == '__main__':
	import sys

	particle_data = dict()
	for arg in sys.argv[1:]:
		load(arg, particle_data)
