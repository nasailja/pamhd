#! /usr/bin/env python

'''
CGI script for starting MHD simulations of PAMHD via a webserver.

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
'''

mpiexec = 'mpiexec'
processes = 4
pamhd_mhd = '../../tests/mhd/test.exe'
mhd2gnuplot = '../../tests/mhd/mhd2gnuplot.exe'
run_dir = '../../tests/mhd/run/'
log_dir = '../../tests/mhd/run/'


def print_form_data(form_data):
	for key in form_data:
		print key, '=', form_data.getfirst(key), '<br>'


from sys import stdout
import cgi
import cgitb
import glob
import os
import subprocess


print 'Content-type: text/html'
print
stdout.flush()

cgitb.enable(display = 0, logdir = log_dir)

print 'Preparing configuration file... '
stdout.flush()

if not os.path.exists(log_dir):
	os.makedirs(log_dir)

if not os.path.exists(run_dir):
	os.makedirs(run_dir)

form_data = cgi.FieldStorage()


'''
Prepare configuration file
'''
config_name = run_dir + '/config'
config = open(config_name, 'w')

config.write(
	'output-directory = ' + run_dir + '\n'
	+ 'load-balancer = RCB\n'
)

# simulation parameters
if not 'time_start' in form_data:
	print 'time_start missing from form data:<br>'
	print_form_data(form_data)
	exit()
config.write('time-start = ' + form_data.getfirst('time_start') + '\n')

if not 'time_length' in form_data:
	print 'time_length missing from form data:<br>'
	print_form_data(form_data)
	exit()
config.write('time-length = ' + form_data.getfirst('time_length') + '\n')

if not 'solver_mhd' in form_data:
	print 'solver_mhd missing from form data:<br>'
	print_form_data(form_data)
	exit()
config.write('solver-mhd = ' + form_data.getfirst('solver_mhd') + '\n')

if not 'save_mhd_n' in form_data:
	print 'save_mhd_n missing from form data:<br>'
	print_form_data(form_data)
	exit()
config.write('save-mhd-n = ' + form_data.getfirst('save_mhd_n') + '\n')

if not 'remove_div_b_n' in form_data:
	print 'remove_div_b_n missing from form data:<br>'
	print_form_data(form_data)
	exit()
config.write('remove-div-B-n = ' + form_data.getfirst('remove_div_b_n') + '\n')

if not 'adiabatic_index' in form_data:
	print 'adiabatic_index missing from form data:<br>'
	print_form_data(form_data)
config.write('adiabatic-index = ' + form_data.getfirst('adiabatic_index') + '\n')

if not 'vacuum_permeability' in form_data:
	print 'vacuum_permeability missing from form data:<br>'
	print_form_data(form_data)
config.write('vacuum-permeability = ' + form_data.getfirst('vacuum_permeability') + '\n')

if not 'proton_mass' in form_data:
	print 'proton_mass missing from form data:<br>'
	print_form_data(form_data)
config.write('proton-mass = ' + form_data.getfirst('proton_mass') + '\n')

if not 'time_step_factor' in form_data:
	print 'time_step_factor missing from form data:<br>'
	print_form_data(form_data)
config.write('time-step-factor = ' + form_data.getfirst('time_step_factor') + '\n')


config.write('[grid]\n')

if not 'periodic_x' in form_data:
	print 'periodic_x missing from form data:<br>'
	print_form_data(form_data)
	exit()
if not 'periodic_y' in form_data:
	print 'periodic_y missing from form data:<br>'
	print_form_data(form_data)
	exit()
if not 'periodic_z' in form_data:
	print 'periodic_z missing from form data:<br>'
	print_form_data(form_data)
	exit()
config.write(
	'periodic = {'
	+ form_data.getfirst('periodic_x') + ', '
	+ form_data.getfirst('periodic_y') + ', '
	+ form_data.getfirst('periodic_z') + '}\n'
)

if not 'cells_x' in form_data:
	print 'cells_x missing from form data:<br>'
	print_form_data(form_data)
	exit()
if not 'cells_y' in form_data:
	print 'cells_y missing from form data:<br>'
	print_form_data(form_data)
	exit()
if not 'cells_z' in form_data:
	print 'cells_z missing from form data:<br>'
	print_form_data(form_data)
	exit()
config.write(
	'nr-cells = {'
	+ form_data.getfirst('cells_x') + ', '
	+ form_data.getfirst('cells_y') + ', '
	+ form_data.getfirst('cells_z') + '}\n'
)

if not 'volume_x' in form_data:
	print 'volume_x missing from form data:<br>'
	print_form_data(form_data)
	exit()
if not 'volume_y' in form_data:
	print 'volume_y missing from form data:<br>'
	print_form_data(form_data)
	exit()
if not 'volume_z' in form_data:
	print 'volume_z missing from form data:<br>'
	print_form_data(form_data)
	exit()
config.write(
	'volume = {'
	+ form_data.getfirst('volume_x') + ', '
	+ form_data.getfirst('volume_y') + ', '
	+ form_data.getfirst('volume_z') + '}\n'
)

if not 'start_x' in form_data:
	print 'start_x missing from form data:<br>'
	print_form_data(form_data)
	exit()
if not 'start_y' in form_data:
	print 'start_y missing from form data:<br>'
	print_form_data(form_data)
	exit()
if not 'start_z' in form_data:
	print 'start_z missing from form data:<br>'
	print_form_data(form_data)
	exit()
config.write(
	'start = {'
	+ form_data.getfirst('start_x') + ', '
	+ form_data.getfirst('start_y') + ', '
	+ form_data.getfirst('start_z') + '}\n'
)


config.write('[initial]\n')

if not 'init_default_mass' in form_data:
	print 'init_default_mass missing from form data:<br>'
	print_form_data(form_data)
	exit()
config.write('default.number-density = ' + form_data.getfirst('init_default_mass') + '\n')

if not 'init_default_vx' in form_data:
	print 'init_default_vx missing from form data:<br>'
	print_form_data(form_data)
	exit()
if not 'init_default_vy' in form_data:
	print 'init_default_vy missing from form data:<br>'
	print_form_data(form_data)
	exit()
if not 'init_default_vz' in form_data:
	print 'init_default_vz missing from form data:<br>'
	print_form_data(form_data)
	exit()
config.write(
	'default.velocity = {'
	+ form_data.getfirst('init_default_vx') + ', '
	+ form_data.getfirst('init_default_vy') + ', '
	+ form_data.getfirst('init_default_vz') + '}\n'
)

if not 'init_default_pressure' in form_data:
	print 'init_default_pressure missing from form data:<br>'
	print_form_data(form_data)
	exit()
config.write('default.pressure = ' + form_data.getfirst('init_default_pressure') + '\n')

if not 'init_default_Bx' in form_data:
	print 'init_default_Bx missing from form data:<br>'
	print_form_data(form_data)
	exit()
if not 'init_default_By' in form_data:
	print 'init_default_By missing from form data:<br>'
	print_form_data(form_data)
	exit()
if not 'init_default_Bz' in form_data:
	print 'init_default_Bz missing from form data:<br>'
	print_form_data(form_data)
	exit()
config.write(
	'default.magnetic-field = {'
	+ form_data.getfirst('init_default_Bx') + ', '
	+ form_data.getfirst('init_default_By') + ', '
	+ form_data.getfirst('init_default_Bz') + '}\n'
)

# figure out number of initial boundary geometries
init_boxes = 0
init_spheres = 0
for key in form_data:
	if key.startswith('init_box_'):
		init_boxes = max(init_boxes, int(key.split('_')[2]))
	if key.startswith('init_sphere_'):
		init_spheres = max(init_spheres, int(key.split('_')[2]))

if init_boxes > 99:
	print 'too many initial condition boxes in form data:<br>'
	print_form_data(form_data)
	exit()
if init_spheres > 99:
	print 'too many initial condition spheres in form data:<br>'
	print_form_data(form_data)
	exit()

if init_boxes > 0:
	config.write('nr-boxes = ' + str(init_boxes) + '\n')
if init_spheres > 0:
	config.write('nr-spheres = ' + str(init_spheres) + '\n')

# write initial condition geometries
for init_box in range(1, init_boxes + 1):
	box_str = 'init_box_' + str(init_box) + '_'
	box_variables = ['startx', 'starty', 'startz', 'endx', 'endy', 'endz', 'mass', 'vx', 'vy', 'vz', 'pressure', 'Bx', 'By', 'Bz']
	for variable in box_variables:
		form_variable = box_str + variable
		if not form_variable in form_data:
			print form_variable + ' missing from form data:<br>'
			print_form_data(form_data)
			exit()
	box_nr_str = 'box' + str(init_box)
	config.write(
		box_nr_str + '.start = {'
			+ form_data.getfirst(box_str + 'startx') + ', '
			+ form_data.getfirst(box_str + 'starty') + ', '
			+ form_data.getfirst(box_str + 'startz') + '}\n'
		+ box_nr_str + '.end = {'
			+ form_data.getfirst(box_str + 'endx') + ', '
			+ form_data.getfirst(box_str + 'endy') + ', '
			+ form_data.getfirst(box_str + 'endz') + '}\n'
		+ box_nr_str + '.number-density = ' + form_data.getfirst(box_str + 'mass') + '\n'
		+ box_nr_str + '.velocity = {'
			+ form_data.getfirst(box_str + 'vx') + ', '
			+ form_data.getfirst(box_str + 'vy') + ', '
			+ form_data.getfirst(box_str + 'vz') + '}\n'
		+ box_nr_str + '.pressure = ' + form_data.getfirst(box_str + 'pressure') + '\n'
		+ box_nr_str + '.magnetic-field = {'
			+ form_data.getfirst(box_str + 'Bx') + ', '
			+ form_data.getfirst(box_str + 'By') + ', '
			+ form_data.getfirst(box_str + 'Bz') + '}\n'
	)

for init_sphere in range(1, init_spheres + 1):
	sphere_str = 'init_sphere_' + str(init_sphere) + '_'
	sphere_variables = ['centerx', 'centerx', 'centerx', 'radius', 'mass', 'vx', 'vy', 'vz', 'pressure', 'Bx', 'By', 'Bz']
	for variable in sphere_variables:
		form_variable = sphere_str + variable
		if not form_variable in form_data:
			print form_variable + ' missing from form data:<br>'
			print_form_data(form_data)
			exit()
	sphere_nr_str = 'sphere' + str(init_sphere)
	config.write(
		sphere_nr_str + '.center = {'
			+ form_data.getfirst(sphere_str + 'centerx') + ', '
			+ form_data.getfirst(sphere_str + 'centery') + ', '
			+ form_data.getfirst(sphere_str + 'centerz') + '}\n'
		+ sphere_nr_str + '.radius = ' + form_data.getfirst(sphere_str + 'radius') + '\n'
		+ sphere_nr_str + '.number-density = ' + form_data.getfirst(sphere_str + 'mass') + '\n'
		+ sphere_nr_str + '.velocity = {'
			+ form_data.getfirst(sphere_str + 'vx') + ', '
			+ form_data.getfirst(sphere_str + 'vy') + ', '
			+ form_data.getfirst(sphere_str + 'vz') + '}\n'
		+ sphere_nr_str + '.pressure = ' + form_data.getfirst(sphere_str + 'pressure') + '\n'
		+ sphere_nr_str + '.magnetic-field = {'
			+ form_data.getfirst(sphere_str + 'Bx') + ', '
			+ form_data.getfirst(sphere_str + 'By') + ', '
			+ form_data.getfirst(sphere_str + 'Bz') + '}\n'
	)


# figure out number of value boundary geometries
value_boxes = 0
value_spheres = 0
for key in form_data:
	if key.startswith('value_box_'):
		value_boxes = max(value_boxes, int(key.split('_')[2]))
	if key.startswith('value_sphere_'):
		value_spheres = max(value_spheres, int(key.split('_')[2]))

if value_boxes > 99:
	print 'too many value boundary boxes in form data:<br>'
	print_form_data(form_data)
	exit()
if value_spheres > 99:
	print 'too many value boundary spheres in form data:<br>'
	print_form_data(form_data)
	exit()

if value_boxes > 0 or value_spheres > 0:
	config.write('[value-boundaries]\n')

if value_boxes > 0:
	config.write('nr-boxes = ' + str(value_boxes) + '\n')
if value_spheres > 0:
	config.write('nr-spheres = ' + str(value_spheres) + '\n')

# write value boundary condition geometries
for value_box in range(1, value_boxes + 1):
	box_str = 'value_box_' + str(value_box) + '_'
	box_variables = ['time', 'startx', 'starty', 'startz', 'endx', 'endy', 'endz', 'mass', 'vx', 'vy', 'vz', 'pressure', 'Bx', 'By', 'Bz']
	for variable in box_variables:
		form_variable = box_str + variable
		if not form_variable in form_data:
			print form_variable + ' missing from form data:<br>'
			print_form_data(form_data)
			exit()
	box_nr_str = 'box' + str(value_box)
	config.write(
		box_nr_str + '.start = {'
			+ form_data.getfirst(box_str + 'startx') + ', '
			+ form_data.getfirst(box_str + 'starty') + ', '
			+ form_data.getfirst(box_str + 'startz') + '}\n'
		+ box_nr_str + '.end = {'
			+ form_data.getfirst(box_str + 'endx') + ', '
			+ form_data.getfirst(box_str + 'endy') + ', '
			+ form_data.getfirst(box_str + 'endz') + '}\n'
		+ box_nr_str + '.time = ' + form_data.getfirst(box_str + 'time') + '\n'
		+ box_nr_str + '.number-density = ' + form_data.getfirst(box_str + 'mass') + '\n'
		+ box_nr_str + '.velocity = {'
			+ form_data.getfirst(box_str + 'vx') + ', '
			+ form_data.getfirst(box_str + 'vy') + ', '
			+ form_data.getfirst(box_str + 'vz') + '}\n'
		+ box_nr_str + '.pressure = ' + form_data.getfirst(box_str + 'pressure') + '\n'
		+ box_nr_str + '.magnetic-field = {'
			+ form_data.getfirst(box_str + 'Bx') + ', '
			+ form_data.getfirst(box_str + 'By') + ', '
			+ form_data.getfirst(box_str + 'Bz') + '}\n'
	)

for value_sphere in range(1, value_spheres + 1):
	sphere_str = 'value_sphere_' + str(value_sphere) + '_'
	sphere_variables = ['time', 'centerx', 'centerx', 'centerx', 'radius', 'mass', 'vx', 'vy', 'vz', 'pressure', 'Bx', 'By', 'Bz']
	for variable in sphere_variables:
		form_variable = sphere_str + variable
		if not form_variable in form_data:
			print form_variable + ' missing from form data:<br>'
			print_form_data(form_data)
			exit()
	sphere_nr_str = 'sphere' + str(value_sphere)
	config.write(
		sphere_nr_str + '.center = {'
			+ form_data.getfirst(sphere_str + 'centerx') + ', '
			+ form_data.getfirst(sphere_str + 'centery') + ', '
			+ form_data.getfirst(sphere_str + 'centerz') + '}\n'
		+ sphere_nr_str + '.radius = ' + form_data.getfirst(sphere_str + 'radius') + '\n'
		+ sphere_nr_str + '.time = ' + form_data.getfirst(sphere_str + 'time') + '\n'
		+ sphere_nr_str + '.number-density = ' + form_data.getfirst(sphere_str + 'mass') + '\n'
		+ sphere_nr_str + '.velocity = {'
			+ form_data.getfirst(sphere_str + 'vx') + ', '
			+ form_data.getfirst(sphere_str + 'vy') + ', '
			+ form_data.getfirst(sphere_str + 'vz') + '}\n'
		+ sphere_nr_str + '.pressure = ' + form_data.getfirst(sphere_str + 'pressure') + '\n'
		+ sphere_nr_str + '.magnetic-field = {'
			+ form_data.getfirst(sphere_str + 'Bx') + ', '
			+ form_data.getfirst(sphere_str + 'By') + ', '
			+ form_data.getfirst(sphere_str + 'Bz') + '}\n'
	)


# figure out number of copy boundary geometries
copy_boxes = 0
copy_spheres = 0
for key in form_data:
	if key.startswith('copy_box_'):
		copy_boxes = max(copy_boxes, int(key.split('_')[2]))
	if key.startswith('copy_sphere_'):
		copy_spheres = max(copy_spheres, int(key.split('_')[2]))

if copy_boxes > 99:
	print 'too many copy boundary boxes in form data:<br>'
	print_form_data(form_data)
	exit()
if copy_spheres > 99:
	print 'too many copy boundary spheres in form data:<br>'
	print_form_data(form_data)
	exit()

if copy_boxes > 0 or copy_spheres > 0:
	config.write('[copy-boundaries]\n')

if copy_boxes > 0:
	config.write('nr-boxes = ' + str(copy_boxes) + '\n')
if copy_spheres > 0:
	config.write('nr-spheres = ' + str(copy_spheres) + '\n')

# write copy boundary condition geometries
for copy_box in range(1, copy_boxes + 1):
	box_str = 'copy_box_' + str(copy_box) + '_'
	box_variables = ['startx', 'starty', 'startz', 'endx', 'endy', 'endz']
	for variable in box_variables:
		form_variable = box_str + variable
		if not form_variable in form_data:
			print form_variable + ' missing from form data:<br>'
			print_form_data(form_data)
			exit()
	box_nr_str = 'box' + str(copy_box)
	config.write(
		box_nr_str + '.start = {'
			+ form_data.getfirst(box_str + 'startx') + ', '
			+ form_data.getfirst(box_str + 'starty') + ', '
			+ form_data.getfirst(box_str + 'startz') + '}\n'
		+ box_nr_str + '.end = {'
			+ form_data.getfirst(box_str + 'endx') + ', '
			+ form_data.getfirst(box_str + 'endy') + ', '
			+ form_data.getfirst(box_str + 'endz') + '}\n'
	)

for copy_sphere in range(1, copy_spheres + 1):
	sphere_str = 'copy_sphere_' + str(copy_sphere) + '_'
	sphere_variables = ['centerx', 'centerx', 'centerx', 'radius']
	for variable in sphere_variables:
		form_variable = sphere_str + variable
		if not form_variable in form_data:
			print form_variable + ' missing from form data:<br>'
			print_form_data(form_data)
			exit()
	sphere_nr_str = 'sphere' + str(copy_sphere)
	config.write(
		sphere_nr_str + '.center = {'
			+ form_data.getfirst(sphere_str + 'centerx') + ', '
			+ form_data.getfirst(sphere_str + 'centery') + ', '
			+ form_data.getfirst(sphere_str + 'centerz') + '}\n'
		+ sphere_nr_str + '.radius = ' + form_data.getfirst(sphere_str + 'radius') + '\n'
		+ sphere_nr_str + '.time = ' + form_data.getfirst(sphere_str + 'time') + '\n'
	)

config.close()

print 'done.<br>'
stdout.flush()

'''
Run simulation and plot results
'''

print 'Running simulation... '
stdout.flush()
try:
	subprocess.check_output([mpiexec, '-n', str(processes), pamhd_mhd, '--config', config_name], universal_newlines = True)
except subprocess.CalledProcessError as error:
	print 'failed, with output:<br>', error.output
	exit()

print 'done.<br> Running visualization... '
stdout.flush()
try:
	files = glob.glob(run_dir + '*.dc')
	subprocess.check_output([mpiexec, '-n', str(processes), mhd2gnuplot] + files, universal_newlines = True)
except subprocess.CalledProcessError as error:
	print 'failed, with output:<br>', error.output
	exit()
print 'done.<br>'
stdout.flush()
