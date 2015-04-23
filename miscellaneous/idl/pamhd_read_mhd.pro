; GDL/IDL function for loading MHD test program data of PAMHD.
;
; Copyright 2015 Ilja Honkonen
; All rights reserved.
;
; Redistribution and use in source and binary forms, with or without modification,
; are permitted provided that the following conditions are met:
;
; * Redistributions of source code must retain the above copyright notice, this
;   list of conditions and the following disclaimer.
;
; * Redistributions in binary form must reproduce the above copyright notice, this
;   list of conditions and the following disclaimer in the documentation and/or
;   other materials provided with the distribution.
;
; * Neither the name of copyright holders nor the names of their contributors
;   may be used to endorse or promote products derived from this software
;   without specific prior written permission.
;
; THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
; ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
; WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
; DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
; ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
; (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
; ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
; (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
; SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
;
;
; Fills variable_names and variable_descriptions with
; comma-separated names and descriptions of variables
; stored in given file.
;
; Fills data_field with an array of arrays of MHD data,
; data_field[0:*, N] is simulation data of N+1th cell stored in file.
;
; Fills meta_data with adiabatic index, proton mass and vacuum premeability. 
;
; Returns 0 on success, > 0 on failure
;
function pamhd_read_mhd, filename, variable_names = variable_names, variable_descriptions = variable_descriptions, data_field = data_field, meta_data = meta_data

	variable_names = 'rx,ry,rz,dx,dy,dz,mas,momx,momy,momz,nrj,magx,magy,magz'
	variable_descriptions = 'x coord of cell center, '
	variable_descriptions += 'y coord of cell center, '
	variable_descriptions += 'z coord of cell center, '
	variable_descriptions += 'length of cell in x dimension, '
	variable_descriptions += 'length of cell in y dimension, '
	variable_descriptions += 'length of cell in z dimension, '
	variable_descriptions += 'mass density, '
	variable_descriptions += 'x component of momentum density, '
	variable_descriptions += 'y component of momentum density, '
	variable_descriptions += 'z component of momentum density, '
	variable_descriptions += 'total energy density, '
	variable_descriptions += 'x component of magnetic field, '
	variable_descriptions += 'y component of magnetic field, '
	variable_descriptions += 'z component of magnetic field'

	openr, in_file, filename, /get_lun

	items_read = ulong64(0)

	; whether fluxes were saved
	fluxes_header = bytarr(11)
	readu, in_file, fluxes_header, transfer_count = items_read
	if (items_read ne 11) then begin
		close, in_file
		free_lun, in_file
		print, "Couldn't read header", items_read
		return, 1
	endif

	fluxes_header = string(fluxes_header)

	have_fluxes = strmatch(fluxes_header, 'fluxes = n' + string(10B))
	if (have_fluxes ne 1) then begin
		close, in_file
		free_lun, in_file
		print, "File with MHD fluxes not supported."
		return, 2
	endif


	; metadata
	simulation_time = double(-1)
	adiabatic_index = double(-1)
	proton_mass = double(-1)
	vacuum_permeability = double(-1)
	readu, in_file, simulation_time, adiabatic_index, proton_mass, vacuum_permeability, transfer_count = items_read
	if (items_read ne 1) then begin
		close, in_file
		free_lun, in_file
		print, "Couldn't read metadata", items_read
		return, 3
	endif

	meta_data = make_array(4, /double)
	meta_data[0] = simulation_time
	meta_data[1] = adiabatic_index
	meta_data[2] = proton_mass
	meta_data[3] = vacuum_permeability


	; endianness check
	endianness = ulong64(1311768467294899695) ; 0x1234567890abcdef
	endianness_file = ulong64(0)
	readu, in_file, endianness_file, transfer_count = items_read
	if (items_read ne 1) then begin
		close, in_file
		free_lun, in_file
		print, "Couldn't read file endianness", items_read
		return, 4
	endif
	if (endianness ne endianness_file) then begin
		close, in_file
		free_lun, in_file
		print, "Unsupported endianness."
		return, 5
	endif


	; number of refinement level 0 cells
	grid_size_x = ulong64(0)
	grid_size_y = ulong64(0)
	grid_size_z = ulong64(0)
	readu, in_file, grid_size_x, grid_size_y, grid_size_z, transfer_count = items_read
	if (items_read ne 1) then begin
		close, in_file
		free_lun, in_file
		print, "Couldn't read grid size", items_read
		return, 6
	endif
	;print, "Grid size in ref. lvl. 0 cells", grid_size_x, grid_size_y, grid_size_z

	max_refinement_level = long(-1)
	readu, in_file, max_refinement_level, transfer_count = items_read
	if (items_read ne 1) then begin
		close, in_file
		free_lun, in_file
		print, "Couldn't read maximum refinement level of grid", items_read
		return, 7
	endif
	;print, "Maximum refinement level of grid cells", max_refinement_level
	if (max_refinement_level ne 0) then begin
		print, "Maximum refinement level > 0 not supported."
		return, 8
		close, in_file
		free_lun, in_file
	endif

	; size of cells' neighborhoods
	neighborhood_length = ulong(0)
	readu, in_file, neighborhood_length, transfer_count = items_read
	if (items_read ne 1) then begin
		close, in_file
		free_lun, in_file
		print, "Couldn't read length of neighborhood", items_read
		return, 9
	endif
	;print, "Length of cells' neighborhood: ", neighborhood_length

	; grid periodicity
	periodic_x = byte(0)
	periodic_y = byte(0)
	periodic_z = byte(0)
	readu, in_file, periodic_x, periodic_y, periodic_z, transfer_count = items_read
	if (items_read ne 1) then begin
		close, in_file
		free_lun, in_file
		print, "Couldn't read grid periodicity", items_read
		return, 10
	endif
	;print, "Grid periodicity info", periodic_x, periodic_y, periodic_z

	; grid geometry
	geometry_id = long(-1)
	readu, in_file, geometry_id, transfer_count = items_read
	if (items_read ne 1) then begin
		close, in_file
		free_lun, in_file
		print, "Couldn't read grid geometry id", items_read
		return, 11
	endif
	if (geometry_id ne 1) then begin
		close, in_file
		free_lun, in_file
		print, "Unsupported geometry: ", geometry_id
		return, 12
	endif

	grid_start_x = double(0)
	grid_start_y = double(0)
	grid_start_z = double(0)
	readu, in_file, grid_start_x, grid_start_y, grid_start_z, transfer_count = items_read
	if (items_read ne 1) then begin
		close, in_file
		free_lun, in_file
		print, "Couldn't read starting coordinate of grid", items_read
		return, 13
	endif
	;print, "Starting coordinate of grid", grid_start_x, grid_start_y, grid_start_z

	lvl_0_cell_length_x = double(0)
	lvl_0_cell_length_y = double(0)
	lvl_0_cell_length_z = double(0)
	readu, in_file, lvl_0_cell_length_x, lvl_0_cell_length_y, lvl_0_cell_length_z, transfer_count = items_read
	if (items_read ne 1) then begin
		close, in_file
		free_lun, in_file
		print, "Couldn't read length of ref. lvl. 0 cells", items_read
		return, 14
	endif
	;print, "Length of ref. lvl. 0 cells", lvl_0_cell_length_x, lvl_0_cell_length_y, lvl_0_cell_length_z


	; total number of cells
	total_cells = ulong64(0)
	readu, in_file, total_cells, transfer_count = items_read
	if (items_read ne 1) then begin
		close, in_file
		free_lun, in_file
		print, "Couldn't read total number of cells", items_read
		return, 15
	endif
	;print, "Total number of cells in grid", total_cells

	; read in cell ids and their datas' byte offsets
	cell_ids_data_offsets = make_array(2, total_cells, /ul64)
	readu, in_file, cell_ids_data_offsets, transfer_count = items_read
	if (items_read ne 2 * total_cells) then begin
		close, in_file
		free_lun, in_file
		print, "Couldn't read cell ids and data offsets", items_read
		return, 16
	endif


	; read cell data into final array
	data_field = make_array(14, total_cells, /double)
	cell_data = make_array(11, /double)
	for i = ulong64(0), total_cells - ulong64(1) do begin
		cell_id = cell_ids_data_offsets[0, i]
		data_offset = cell_ids_data_offsets[1, i]

		point_lun, in_file, data_offset
		readu, in_file, cell_data, transfer_count = items_read
		if (items_read ne 11) then begin
			close, in_file
			free_lun, in_file
			print, "Couldn't read cell data at byte offset", data_offset
			print, items_read, "items read, should be 11"
			return, 17
		endif

		data_field[6:13, i] = cell_data[0:7]

		; add derived data
		cell_id = cell_id - 1
		cell_x_index = cell_id mod grid_size_x
		cell_y_index = cell_id / grid_size_x mod grid_size_y
		cell_z_index = cell_id / (grid_size_x * grid_size_y)

		data_field[0, i] = grid_start_x + lvl_0_cell_length_x * (0.5 + cell_x_index)
		data_field[1, i] = grid_start_y + lvl_0_cell_length_y * (0.5 + cell_y_index)
		data_field[2, i] = grid_start_z + lvl_0_cell_length_z * (0.5 + cell_z_index)

		data_field[3, i] = lvl_0_cell_length_x
		data_field[4, i] = lvl_0_cell_length_y
		data_field[5, i] = lvl_0_cell_length_z
	endfor

	close, in_file
	free_lun, in_file

	return, 0
end
