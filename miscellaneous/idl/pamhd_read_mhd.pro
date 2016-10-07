; GDL/IDL function for loading MHD test program data of PAMHD.
;
; Copyright 2015, 2016 Ilja Honkonen
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
; Fills meta_data with simulation time, adiabatic index,
; proton mass and vacuum premeability.
;
; If given sort_data, returned data is sorted in ascending order
; by rz then items with identical rz are sorted by ry and items
; with identical rz and ry are sorted by rx.
;
; If given, volume must be an array of 6 doubles describing the volume
; from which to return data. Array items are interpreted in order
; min_x, max_x, min_y, ..., max_z. Simulation cells whose center
; is outside given ranges are excluded from returned data.
;
; Returns 0 on success, > 0 on failure
;
; Example when GDL/IDL started from this directory:
; .run pamhd_read_mhd
; vol = make_array(6,/double)
; vol[0] = 0
; vol[1] = 1
; vol[2] = -!values.d_infinity
; vol[3] = !values.d_infinity
; vol[4] = -!values.d_infinity
; vol[5] = !values.d_infinity
; pamhd_read_mhd('mhd_0.000e+00_s.dc', variable_names = names, variable_descriptions = descriptions, data_field = data, meta_data = meta, volume = vol)
;
function pamhd_read_mhd, filename, variable_names = variable_names, variable_descriptions = variable_descriptions, data_field = data_field, meta_data = meta_data, sort_data = sort_data, volume = volume

	variable_names = 'rx,ry,rz,dx,dy,dz,mas,momx,momy,momz,nrj,'
	variable_names += 'magx,magy,magz,curx,cury,curz,pressure,'
	variable_names += 'res,rank,type,pxbgmagx,pxbgmagy,pxbgmagz,'
	variable_names += 'pybgmagx,pybgmagy,pybgmagz,'
	variable_names += 'pzbgmagx,pzbgmagy,pzbgmagz'

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
	variable_descriptions += 'x component of perturbed magnetic field, '
	variable_descriptions += 'y component of perturbed magnetic field, '
	variable_descriptions += 'z component of perturbed magnetic field, '
	variable_descriptions += 'x component of electric current density, '
	variable_descriptions += 'y component of electric current density, '
	variable_descriptions += 'z component of electric current density, '
	variable_descriptions += 'thermal pressure, '
	variable_descriptions += 'electrical resistivity, '
	variable_descriptions += 'owner of cell (MPI rank), '
	variable_descriptions += 'cell type, '
	variable_descriptions += 'x component of background magnetic field on positive x side cell face, '
	variable_descriptions += 'y component of background magnetic field on pos x face, '
	variable_descriptions += 'z component of background magnetic field on pos x face, '
	variable_descriptions += 'background Bx on positive y side cell face, '
	variable_descriptions += 'background By on pos y face, '
	variable_descriptions += 'background Bz on pos y face, '
	variable_descriptions += 'background Bx on pos z face, '
	variable_descriptions += 'background By on pos z face, '
	variable_descriptions += 'background Bz on pos z face'

	openr, in_file, filename, /get_lun

	items_read = ulong64(0)

	; file version
	file_version = ulong64(0)
	readu, in_file, file_version, transfer_count = items_read
	if (items_read ne 1) then begin
		close, in_file
		free_lun, in_file
		print, "Couldn't read file version", items_read
		return, 1
	endif
	if (file_version ne 2) then begin
		close, in_file
		free_lun, in_file
		print, "Unsupported file version."
		return, 2
	endif

	; simulation parameters
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

	if (~keyword_set(volume)) then begin
		volume = make_array(6, /double)
		volume[0] = -!values.d_infinity
		volume[1] = +!values.d_infinity
		volume[2] = -!values.d_infinity
		volume[3] = +!values.d_infinity
		volume[4] = -!values.d_infinity
		volume[5] = +!values.d_infinity
	endif

	; read cell data into final array
	temp_data = make_array(30, total_cells, /double)
	cell_data1 = make_array(11, /double)
	cell_data2 = make_array(2, /long)
	cell_data3 = make_array(1, /double)
	cell_data4 = make_array(9, /double)
	loaded_cells = ulong64(0)
	B = make_array(3, /double) ; magnetic field
	M = make_array(3, /double) ; momentum density
	for i = ulong64(0), total_cells - ulong64(1) do begin
		cell_id = cell_ids_data_offsets[0, i]

		; cell center, defined by get_center() function in
		; dccrg_cartesian_geometry.hpp file of dccrg
		cell_id = cell_id - 1
		cell_x_index = cell_id mod grid_size_x
		cell_y_index = cell_id / grid_size_x mod grid_size_y
		cell_z_index = cell_id / (grid_size_x * grid_size_y)
		cell_id = cell_id + 1

		cell_center_x = grid_start_x + lvl_0_cell_length_x * (0.5 + cell_x_index)
		cell_center_y = grid_start_y + lvl_0_cell_length_y * (0.5 + cell_y_index)
		cell_center_z = grid_start_z + lvl_0_cell_length_z * (0.5 + cell_z_index)

		; don't load if cell outside requested volume
		if (cell_center_x lt volume[0]) then begin
			continue
		endif
		if (cell_center_x gt volume[1]) then begin
			continue
		endif
		if (cell_center_y lt volume[2]) then begin
			continue
		endif
		if (cell_center_y gt volume[3]) then begin
			continue
		endif
		if (cell_center_z lt volume[4]) then begin
			continue
		endif
		if (cell_center_z gt volume[5]) then begin
			continue
		endif

		temp_data[0, loaded_cells] = cell_center_x
		temp_data[1, loaded_cells] = cell_center_y
		temp_data[2, loaded_cells] = cell_center_z
		temp_data[3, loaded_cells] = lvl_0_cell_length_x
		temp_data[4, loaded_cells] = lvl_0_cell_length_y
		temp_data[5, loaded_cells] = lvl_0_cell_length_z

		data_offset = cell_ids_data_offsets[1, i]

		point_lun, in_file, data_offset
		readu, in_file, cell_data1, cell_data2, cell_data3, cell_data4, transfer_count = items_read
		if (items_read ne 9) then begin
			close, in_file
			free_lun, in_file
			print, "Couldn't read cell ", cell_id, " data at byte offset ", data_offset
			print, items_read, " item sets read, should be 1"
			return, 17
		endif

		temp_data[6, loaded_cells] = cell_data1
		B[0] = cell_data1[5]
		B[1] = cell_data1[6]
		B[2] = cell_data1[7]
		M[0] = cell_data1[1]
		M[1] = cell_data1[2]
		M[2] = cell_data1[3]
		if (cell_data1[0] le 0) then begin
			print, "Non-positive density in cell ", cell_id
		end
		pressure $
			= (adiabatic_index - 1) $
			* ( $
				cell_data1[4] $
				- (M[0]*M[0] + M[1]*M[1] + M[2]*M[2]) / 2 / cell_data1[0] $
				- (B[0]*B[0] + B[1]*B[1] + B[2]*B[2]) / 2 / vacuum_permeability $
			)
		temp_data[17, loaded_cells] = pressure
		temp_data[18, loaded_cells] = cell_data3[0]
		temp_data[19, loaded_cells] = double(cell_data2[1])
		temp_data[20, loaded_cells] = double(cell_data2[0])
		temp_data[21, loaded_cells] = cell_data4

		loaded_cells = loaded_cells + 1
	endfor

	; don't return too large array
	data_field = make_array(30, loaded_cells, /double)
	data_field = temp_data[*, 0:loaded_cells - 1]

	close, in_file
	free_lun, in_file

	if (keyword_set(sort_data)) then begin
		; sort by rz
		data_field = data_field[*, sort(data_field[2, *])]

		; sort by ry each set of rows with same value in rz
		subranges_y = uniq(data_field[2, *])
		range_y_start = 0
		foreach range_y_end, subranges_y do begin
			;print, 'processing subrange with Z value of', data_field[2, range_y_start]
			data_field[*, range_y_start : range_y_end] $
				= data_field[*, range_y_start + sort(data_field[1, range_y_start : range_y_end])]

			; sort by rx each set of rows with same value in ry
			subranges_x = range_y_start + uniq(data_field[1, range_y_start : range_y_end])
			range_x_start = range_y_start
			foreach range_x_end, subranges_x do begin
				;print, '   processing subrange with Y value of', data_field[1, range_x_start]
				data_field[*, range_x_start : range_x_end] $
					= data_field[*, range_x_start + sort(data_field[0, range_x_start : range_x_end])]

				range_x_start = range_x_end + 1
			endforeach

			range_y_start = range_y_end + 1
		endforeach

	endif

	return, 0
end
