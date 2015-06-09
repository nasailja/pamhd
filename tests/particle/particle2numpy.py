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


If called as the main program creates empty dictionaries bulk_data
and particle_data and fills them by calling load() on each file
given as an argument to this program.

Example printing the id of every particle in one file:
python -i particle2numpy.py particle_0.000e+00_s.dc
for key in particle_data:
	print(key)
exit()

Example printing the electric field in cell 3 in the first given file:
python -i particle2numpy.py particle_0.000e+00_s.dc
print(bulk_data[3][0][0])

Example printing the velocity of particle with id == 4 from the last given
file in which the particle existed:
python -i particle2numpy.py particle_0.000e+00_s.dc
print(particle_data[4][-1][1])
exit()
'''

import os

try:
	import numpy
except:
	exit("Couldn't import numpy")


'''
Returns bulk and particle data of PAMHD particle programs.

If bulk_data or particle_data are None will skip reading
the corresponding data otherwise bulk_data and particle_data
must be dictionaries to which already existing data is
appended from given file. If both are None will do nothing
and return None

Returns the following in a numpy array:
simulation time,
adiabatic index,
vacuum permeability,
ratio of particle temperature to energy (Boltzmann constant).

bulk_data will have cell ids as keys and each value will be
a list of tuples of numpy arrays of electric and magnetic field.
particle_data will have particle ids as keys and each value
will be a list of tuples of numpy arrays of particle position
and velocity.
'''
def load(file_name, bulk_data, particle_data):
	if bulk_data == None and particle_data == None:
		return None

	if not os.path.isfile(file_name):
		raise Exception('Given file name (' + file_name + ')is not a file.')

	infile = open(file_name, 'rb')

	# read simulation header, given by source/particle/save.hpp
	header = numpy.fromfile(infile, dtype = '4float64', count = 1)

	# from this point onward format given in by dccrg's save_grid_data()
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

	# until this point format decided by dccrg
	# read particle data
	for item in cell_ids_data_offsets:
		cell_id = item[0]

		if bulk_data == None:
			infile.seek(item[1] + 7 * 8, 0)
		else:
			infile.seek(item[1], 0)
			# given by source/particle/save.hpp
			temp_bulk_data = numpy.fromfile(
				infile,
				dtype = '3float, 3float64, uint64',
				count = 1
			)[0]
			if cell_id in bulk_data:
				bulk_data[cell_id].append((temp_bulk_data[0], temp_bulk_data[1]))
			else:
				bulk_data[cell_id] = [(temp_bulk_data[0], temp_bulk_data[1])]

		# given by Particle_Internal in source/particle/variables.hpp
		if particle_data == None:
			continue

		temp_particle_data = numpy.fromfile(
			infile,
			dtype = '3float64, 3float64, float64, float64, float64, uint64',
			count = temp_bulk_data[2]
		)
		for particle in temp_particle_data:
			position, velocity, mass, species_mass, charge_mass_ratio, particle_id \
				= particle[0], particle[1], particle[2], particle[3], particle[4], particle[5]

			if particle_id in particle_data:
				particle_data[particle_id].append((position, velocity))
			else:
				particle_data[particle_id] = [(position, velocity)]

	return header


if __name__ == '__main__':
	import sys

	bulk_data = dict()
	particle_data = dict()
	for arg in sys.argv[1:]:
		load(arg, bulk_data, particle_data)
