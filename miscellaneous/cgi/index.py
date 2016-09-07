#! /usr/bin/env python

'''
CGI script for starting MHD simulations of PAMHD via a webserver.

Copyright 2015, 2016 Ilja Honkonen
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

# Edit the following to correspond to your configuration,
# example values are given as comments
mpiexec = None # 'mpiexec'
processes = None # 4
pamhd_mhd = None # '../../tests/mhd/test.exe'
mhd2gnuplot = None # '../../tests/mhd/mhd2gnuplot.exe'
run_dir = None # '../../tests/mhd/run/'
log_dir = None # '../../tests/mhd/run/'

# you shouldn't have to edit the rest of the file
# but if you do consider forking the code on github
# and submitting a pull request for the changes

from time import sleep
# prevent users from spamming this too often
sleep(0.2)

def print_form_data(form_data):
	for key in form_data:
		print key, '=', form_data.getfirst(key), '<br>'


from sys import stdout
import cgi
import cgitb
import glob
from json import dumps
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


config_data = dict()


'''
Prepare configuration file
'''
config_name = run_dir + '/config'
config = open(config_name, 'w')

config_data['output-directory'] = run_dir
config_data['load-balancer'] = 'RCB'


# simulation parameters
if not 'time-start' in form_data:
	print 'time-start missing from form data:<br>'
	print_form_data(form_data)
	exit()
config_data['time-start'] = float(form_data.getfirst('time-start'))

if not 'time-length' in form_data:
	print 'time-length missing from form data:<br>'
	print_form_data(form_data)
	exit()
config_data['time-length'] = float(form_data.getfirst('time-length'))

if not 'minimum-pressure' in form_data:
	print 'minimum-pressure missing from form data:<br>'
	print_form_data(form_data)
	exit()
config_data['minimum-pressure'] = float(form_data.getfirst('minimum-pressure'))

if not 'solver-mhd' in form_data:
	print 'solver-mhd missing from form data:<br>'
	print_form_data(form_data)
	exit()
config_data['solver-mhd'] = str(form_data.getfirst('solver-mhd'))

if not 'save-mhd-n' in form_data:
	print 'save-mhd-n missing from form data:<br>'
	print_form_data(form_data)
	exit()
config_data['save-mhd-n'] = float(form_data.getfirst('save-mhd-n'))

if not 'resistivity' in form_data:
	print 'resistivity missing from form data:<br>'
	print_form_data(form_data)
	exit()
config_data['resistivity'] = str(form_data.getfirst('resistivity'))

if not 'remove-div-b-n' in form_data:
	print 'remove-div-b-n missing from form data:<br>'
	print_form_data(form_data)
	exit()
config_data['remove-div-B-n'] = float(form_data.getfirst('remove-div-b-n'))

if not 'adiabatic-index' in form_data:
	print 'adiabatic-index missing from form data:<br>'
	print_form_data(form_data)
	exit()
config_data['adiabatic-index'] = float(form_data.getfirst('adiabatic-index'))

if not 'vacuum-permeability' in form_data:
	print 'vacuum-permeability missing from form data:<br>'
	print_form_data(form_data)
	exit()
config_data['vacuum-permeability'] = float(form_data.getfirst('vacuum-permeability'))

if not 'proton-mass' in form_data:
	print 'proton-mass missing from form data:<br>'
	print_form_data(form_data)
	exit()
config_data['proton-mass'] = float(form_data.getfirst('proton-mass'))

if not 'time-step-factor' in form_data:
	print 'time-step-factor missing from form data:<br>'
	print_form_data(form_data)
	exit()
config_data['time-step-factor'] = float(form_data.getfirst('time-step-factor'))

config_data['poisson-norm-stop'] = 1e-10
config_data['poisson-norm-increase-max'] = 10.0

# grid options
config_data['grid-options'] = dict()
if not 'periodic-x' in form_data or not 'periodic-y' in form_data or not 'periodic-z' in form_data:
	print 'periodicity missing from form data:<br>'
	print_form_data(form_data)
	exit()
config_data['grid-options']['periodic'] = '{' + str(form_data.getfirst('periodic-x')) + ', ' + str(form_data.getfirst('periodic-y')) + ', ' + str(form_data.getfirst('periodic-z')) + '}'

if not 'cells' in form_data:
	print 'cells missing from form data:<br>'
	print_form_data(form_data)
	exit()
config_data['grid-options']['cells'] = str(form_data.getfirst('cells'))

if not 'volume' in form_data:
	print 'grid volume missing from form data:<br>'
	print_form_data(form_data)
	exit()
config_data['grid-options']['volume'] = str(form_data.getfirst('volume'))

if not 'start' in form_data:
	print 'grid start missing from form data:<br>'
	print_form_data(form_data)
	exit()
config_data['grid-options']['start'] = str(form_data.getfirst('start'))

# geometries
geom_nr = 0
while True:
	box_id = 'box' + str(geom_nr)
	if \
		box_id + '-startx' in form_data \
		and box_id + '-starty' in form_data \
		and box_id + '-startz' in form_data \
		and box_id + '-endx' in form_data \
		and box_id + '-endy' in form_data \
		and box_id + '-endz' in form_data \
	:
		box = dict()
		box['start'] = [
			float(form_data.getfirst(box_id + '-startx')),
			float(form_data.getfirst(box_id + '-starty')),
			float(form_data.getfirst(box_id + '-startz'))
		]
		box['end'] = [
			float(form_data.getfirst(box_id + '-endx')),
			float(form_data.getfirst(box_id + '-endy')),
			float(form_data.getfirst(box_id + '-endz'))
		]
		if not 'geometries' in config_data:
			config_data['geometries'] = []
		config_data['geometries'].append(dict(box = box))
		geom_nr += 1
		continue

	sphere_id = 'sphere' + str(geom_nr)
	if \
		sphere_id + '-centerx' in form_data \
		and sphere_id + '-centery' in form_data \
		and sphere_id + '-centerz' in form_data \
		and sphere_id + '-radius' in form_data \
	:
		sphere = dict()
		sphere['center'] = [
			float(form_data.getfirst(sphere_id + '-centerx')),
			float(form_data.getfirst(sphere_id + '-centery')),
			float(form_data.getfirst(sphere_id + '-centerz'))
		]
		sphere['radius'] = float(form_data.getfirst(sphere_id + '-radius'))
		if not 'geometries' in config_data:
			config_data['geometries'] = []
		config_data['geometries'].append(dict(sphere = sphere))
		geom_nr += 1
		continue

	break

# background magnetic field
if 'background-magnetic-field' in form_data:
	bg_B = dict()

	if not 'background-B-value-x' in form_data or not 'background-B-value-y' in form_data or not 'background-B-value-z' in form_data:
		print 'background magnetic field value missing from form data:<br>'
		print_form_data(form_data)
		exit()
	bg_B['values'] = [
		float(form_data.getfirst('background-B-value-x')),
		float(form_data.getfirst('background-B-value-y')),
		float(form_data.getfirst('background-B-value-z'))
	]

	if 'dipole-min-distance' in form_data:
		bg_B['dipole-min-distance'] = float(form_data.getfirst('dipole-min-distance'))

	dipoles = []
	dipole_nr = 0
	while True:
		dipole_str = 'dipole_' + str(dipole_nr)
		if \
			dipole_str + '-momentum-x' in form_data \
			and dipole_str + '-momentum-y' in form_data \
			and dipole_str + '-momentum-z' in form_data \
			and dipole_str + '-position-x' in form_data \
			and dipole_str + '-position-y' in form_data \
			and dipole_str + '-position-z' in form_data \
		:
			dipole = dict()
			dipole['moment'] = [
				float(form_data.getfirst(dipole_str + '-momentum-x')),
				float(form_data.getfirst(dipole_str + '-momentum-y')),
				float(form_data.getfirst(dipole_str + '-momentum-z'))
			]
			dipole['position'] = [
				float(form_data.getfirst(dipole_str + '-position-x')),
				float(form_data.getfirst(dipole_str + '-position-y')),
				float(form_data.getfirst(dipole_str + '-position-z'))
			]
			dipoles.append(dipole)
			dipole_nr += 1
			continue
		break

	bg_B['dipoles'] = dipoles
	config_data['background-magnetic-field'] = bg_B


# returns a list of expressions as strings
def get_strings(string):
	string = string.replace(' ', '')
	values = []
	# vector exprs
	if '{' in string:
		for i in string.split('}'):
			for j in i.split('{'):
				if j == '':
					continue
				values.append('{' + j + '}')
	# scalar exprs
	else:
		[values.append(i) for i in string.split(',')]
	return values

# returns a list of floats
def get_numbers(string):
	return [float(a) for a in string.split(',')]


# number density
config_data['number-density'] = dict()
if not 'default-mass' in form_data:
	print 'default-mass missing from form data:<br>'
	print_form_data(form_data)
	exit()
config_data['number-density']['default'] = str(form_data.getfirst('default-mass'))

init_nr = 0
while True:
	init_str = 'initial-mass-' + str(init_nr)
	if \
		init_str + '-geometry-id' in form_data \
		and init_str + '-value' in form_data \
	:
		init = dict()
		init['geometry-id'] = int(form_data.getfirst(init_str + '-geometry-id'))
		init['value'] = str(form_data.getfirst(init_str + '-value'))
		if not 'initial-conditions' in config_data['number-density']:
			config_data['number-density']['initial-conditions'] = []
		config_data['number-density']['initial-conditions'].append(init)
		init_nr += 1
		continue
	break

value_nr = 0
while True:
	value_str = 'value-mass-' + str(value_nr) + '-'
	if \
		value_str + 'geometry-id' in form_data \
		and value_str + 'time-stamps-id' in form_data \
		and value_str + 'value' in form_data \
	:
		value = dict()
		value['geometry-id'] = int(form_data.getfirst(value_str + 'geometry-id'))
		value['time-stamps'] = get_numbers(form_data.getfirst(value_str + 'time-stamps-id'))
		value['values'] = get_strings(form_data.getfirst(value_str + 'value'))
		if not 'value-boundaries' in config_data['number-density']:
			config_data['number-density']['value-boundaries'] = []
		config_data['number-density']['value-boundaries'].append(value)
		value_nr += 1
		continue
	break

copy_nr = 0
while True:
	copy_str = 'copy-mass-' + str(copy_nr)
	if copy_str + '-geometry-id' in form_data:
		copy = dict()
		copy['geometry-id'] = int(form_data.getfirst(copy_str + '-geometry-id'))
		if not 'copy-boundaries' in config_data['number-density']:
			config_data['number-density']['copy-boundaries'] = []
		config_data['number-density']['copy-boundaries'].append(copy)
		copy_nr += 1
		continue
	break

# rest of simulation variables
for var in ['velocity', 'pressure', 'magnetic-field']:
	config_data[var] = dict()
	if not 'default-' + var in form_data:
		print 'default-' + var + ' missing from form data:<br>'
		print_form_data(form_data)
		exit()
	config_data[var]['default'] = str(form_data.getfirst('default-' + var))

	init_nr = 0
	while True:
		init_str = 'initial-' + var + '-' + str(init_nr)
		if \
			init_str + '-geometry-id' in form_data \
			and init_str + '-value' in form_data \
		:
			init = dict()
			init['geometry-id'] = int(form_data.getfirst(init_str + '-geometry-id'))
			init['value'] = str(form_data.getfirst(init_str + '-value'))
			if not 'initial-conditions' in config_data[var]:
				config_data[var]['initial-conditions'] = []
			config_data[var]['initial-conditions'].append(init)
			init_nr += 1
			continue
		break

	value_nr = 0
	while True:
		value_str = 'value-' + var + '-' + str(value_nr) + '-'
		if \
			value_str + 'geometry-id' in form_data \
			and value_str + 'time-stamps-id' in form_data \
			and value_str + 'value' in form_data \
		:
			value = dict()
			value['geometry-id'] = int(form_data.getfirst(value_str + 'geometry-id'))
			value['time-stamps'] = get_numbers(form_data.getfirst(value_str + 'time-stamps-id'))
			value['values'] = get_strings(form_data.getfirst(value_str + 'value'))
			if not 'value-boundaries' in config_data[var]:
				config_data[var]['value-boundaries'] = []
			config_data[var]['value-boundaries'].append(value)
			value_nr += 1
			continue
		break

	copy_nr = 0
	while True:
		copy_str = 'copy-' + var + '-' + str(copy_nr)
		if copy_str + '-geometry-id' in form_data:
			copy = dict()
			copy['geometry-id'] = int(form_data.getfirst(copy_str + '-geometry-id'))
			if not 'copy-boundaries' in config_data[var]:
				config_data[var]['copy-boundaries'] = []
			config_data[var]['copy-boundaries'].append(copy)
			copy_nr += 1
			continue
		break


s = dumps(config_data, sort_keys = True, indent = 4, separators = (',', ': '))
config.write(s)
config.close()

print 'done.<br>'
stdout.flush()

'''
Run simulation
'''

print 'Running simulation...'
stdout.flush()

proc = None
stdout_ = open(os.path.join(log_dir, 'stdout'), 'w')
stderr_ = open(os.path.join(log_dir, 'stderr'), 'w')
try:
	proc = subprocess.Popen( \
		[mpiexec, '-n', str(processes), pamhd_mhd, config_name], \
		universal_newlines = True, \
		stdout = stdout_, \
		stderr = stderr_ \
	)
except:
	print 'error starting simulation.<br>'
	stdout.flush()
	exit()

while proc.poll() == None:
	sleep(5)
	stdout.write('.')
	stdout.flush()
if proc.poll() != 0:
	print ' failed.<br>'
	stdout.flush()
	exit()

print 'done.<br>Running visualization... '
stdout.flush()
try:
	files = glob.glob(run_dir + '*.dc')
	subprocess.check_output([mpiexec, '-n', str(processes), mhd2gnuplot] + files, universal_newlines = True)
except subprocess.CalledProcessError as error:
	print 'failed.<br>'
	exit()
print 'done.<br>'
stdout.flush()
