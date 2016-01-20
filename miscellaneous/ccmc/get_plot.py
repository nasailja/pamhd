#! /usr/bin/env python3
'''
Downloads a plot obtained through CCMC visualization web interface.

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

from argparse import ArgumentParser
import shutil
from urllib import parse, request

default_headers = {
	'Accept': 'text/plain',
	'User-Agent': 'PAMHD get_plot.py'
}


'''
Returns directory id of given run id of format first_last_yyyymmdd_n.
'''
def get_run_directory(run_id, domain = 'LP'):
	url \
		= 'http://ccmc.gsfc.nasa.gov/results/viewrun.php?domain=' \
		+ domain \
		+ '&runnumber=' \
		+ run_id
	req = request.Request(url, data = None, headers = default_headers)

	dir_id = None
	with request.urlopen(req) as response:
		the_page = response.read()
		the_page = str(the_page)
		for i in the_page.split():
			if 'run_idl3d.cgi?dir=' in i:
				dir_id = i.split('run_idl3d.cgi?dir=')[1].split('"')[0]

	if not dir_id:
		print('Could not determine dir id for run', run_id, 'from page:')
		print(the_page)
		exit(1)
	return dir_id


'''
Returns a dictionary of available plot times and the default.

Dictionary keys are available times as presented to user.
Values are html codes corresponding to times presented to user.
Returned default is the key i.e. time presented to user.
'''
def get_plotting_times(the_page):
	opt_started = False
	times = dict()
	selected = None
	for line in the_page:
		line = line.lower()
		if '<select name="itime">' in line:
			opt_started = True
		if opt_started and '/select' in line:
			opt_started = False
			break
		if opt_started and 'option' in line and 'value' in line:
			# time presented to user
			t = line.split('select')[-1].split('value=')[1].split('>')[1].split('<')[0].strip()
			# html of above selection
			value = line.split('select')[-1].split('value=')[1].split('>')[0].strip()
			times[t] = value
			if 'selected' in line:
				selected = t

	return selected, times


'''
Returns a dictionary of available plot modes and the default.

Dictionary keys are available modes as presented to user.
Values are html codes corresponding to modes presented to user.
Returned default is the key i.e. mode presented to user.
'''
def get_plotting_modes(the_page):
	opt_started = False
	modes = dict()
	selected = None
	for line in the_page:
		line = line.lower()
		items = line.split('<')
		for item in items:
			if 'select name="plotmode"' in item:
				opt_started = True
			if opt_started and '/select' in item:
				opt_started = False
				break
			if opt_started and 'option' in item and 'value' in item:
				# name presented to user
				mode = item.split('value=')[1].split('>')[1].strip().replace(')', '').replace('(', '').replace(' ','')
				# html code of above selection
				value = item.split('value=')[1].split('>')[0]
				modes[mode] = value
				if 'selected' in item:
					selected = mode

	return selected, modes


'''
Returns a dictionary of available plot quantities and the default.

Dictionary keys are available quantities as presented to user.
Values are html codes corresponding to quantities presented to user.
Returned default is the key i.e. quantity presented to user.
'''
def get_plot_quantities(the_page):
	opt_started = False
	quantities = dict()
	selected = None
	for line in the_page:
		line = line.lower()
		if '<select name="quantity1">' in line:
			opt_started = True
		if opt_started and '/select' in line:
			opt_started = False
			break
		if opt_started and 'option' in line and 'value' in line:
			# name presented to user
			qty = line.split('value=')[1].split('>')[1].split('<')[0].strip()
			# html code of above selection
			value = line.split('value=')[1].split('>')[0].strip()
			quantities[qty] = value
			if 'selected' in line:
				selected = qty

	return selected, quantities


def get_plot_ranges(the_page):
	ranges = [] # -x, +x, -y, +y, -z, +z
	for line in the_page:
		if 'Range:' in line:
			ranges.append(float(line.split()[1]))
			ranges.append(float(line.split()[3]))
	cut_plane = None
	for line in the_page:
		if 'name="CutPlane"' in line and 'CHECKED' in line:
			if 'value=0' in line:
				cut_plane = '0'
			elif 'value=1' in line:
				cut_plane = '1'
			elif 'value=2' in line:
				cut_plane = '2'
	cut_plane_pos = None
	for line in the_page:
		if \
			(cut_plane == '0' and 'name="X1CUT"' in line) \
			or (cut_plane == '1' and 'name="X2CUT"' in line) \
			or (cut_plane == '2' and 'name="X3CUT"' in line) \
		:
			cut_plane_pos = line.split('value=')[1].split('>')[0]

	return ranges, cut_plane, cut_plane_pos


def get_info(dir_id):
	url = 'http://ccmc.gsfc.nasa.gov/cgi-bin/run_idl3d.cgi'

	values = {'dir': dir_id}

	data = parse.urlencode(values)
	data = data.encode('ascii')
	req = request.Request(url, data = data, headers = default_headers)
	with request.urlopen(req) as response:
		the_page = response.read()
		the_page = str(the_page)
	the_page = the_page.split('\\n')
	return the_page


'''
Returns address of figure from default CCMC plot.

Given arguments customize plotted figure to some extent.
'''
def get_default_plot(dir_id, request_values):
	url = 'http://ccmc.gsfc.nasa.gov/cgi-bin/run_idl3d.cgi'

	data = parse.urlencode(request_values)
	data = data.encode('ascii')
	req = request.Request(url, data = data, headers = default_headers)
	with request.urlopen(req) as response:
		the_page = response.read()
		the_page = str(the_page)
	the_page = the_page.split('\\n')
	return the_page


parser = ArgumentParser(description = 'Download a plot from CCMC web interface.')

parser.add_argument('--info', action = 'store_true', help = 'Print available plot options for given run without downloading plots')
parser.add_argument('--run-number', dest = 'run_id', metavar = 'R', required = True, help = 'Download plot for run number(s) R')
parser.add_argument('--variable', metavar = 'V', help = "Download plot for variable V, if not given CCMC's default variable is downloaded")
parser.add_argument('--time', metavar = 'T', help = "Download plot for simulation time T, if not given CCMC's default time is used")
parser.add_argument('--output-prefix', dest = 'prefix', metavar = 'O', help = 'Save downloaded image with prefix O, if not given run number is used')

args = parser.parse_args()


print('Requesting info for run', args.run_id)

run_dir = get_run_directory(args.run_id)
info_page = get_info(run_dir)
default_time, times = get_plotting_times(info_page)
selected_mode, modes = get_plotting_modes(info_page)
default_qty, quantities = get_plot_quantities(info_page)
ranges, cut_plane, cut_plane_pos = get_plot_ranges(info_page)


if args.info:
	print('\nTimes available for plotting:')
	for t in times:
		if t == default_time:
			print(t, '(default)', end = ', ')
		else:
			print(t, end = ', ')
	print('\n\nVariables available for plotting:')
	for v in quantities:
		if v == default_qty:
			print(v, '(default)', end = ', ')
		else:
			print(v, end = ', ')
	print()
	exit()

# these are required and aren't set by program options
values = {
	'subdir': '',
	'dirname': '',
	'Runname': 'mhd',
	'PlotMagni': '0',
	'Linethick': '0',
	'Charthick': '0',
	'aspectratio': '1',
	'Grid': '0',
	'Interpolate': '1',
	'AX': '30',
	'AZ': '30',
	'Output_all_parameters': '0',
	'Output_Pointdata_in_grid': '0',
	'Output_Pointdata': '0',
}



values['dir'] = run_dir

requested_time = None
if args.time:
	if args.time in times:
		requested_time = args.time
	else:
		print('Time', args.time, "doesn't exist in available plotting times:")
		print(times)
		exit(2)
else:
	requested_time = default_time
values['Itime'] = times[requested_time]

values['PlotMode'] = modes[selected_mode]

values['X1MIN'] = ranges[0]
values['X1MAX'] = ranges[1]
values['X2MIN'] = ranges[2]
values['X2MAX'] = ranges[3]
values['X3MIN'] = ranges[4]
values['X3MAX'] = ranges[5]

values['CutPlane'] = cut_plane

values['X1CUT'] = values['X2CUT'] = values['X3CUT'] = cut_plane_pos

requested_variable = None
if args.variable:
	if args.variable in quantities:
		requested_variable = args.variable
	else:
		print('Variable', args.variable, "doesn't exist in available plot variables:")
		print(quantities)
		exit(3)
else:
	requested_variable = default_qty
values['Quantity1'] = values['Quantity2'] = values['Quantity3'] = quantities[requested_variable]

print('Requesting plot for run', args.run_id)
plot_page = get_default_plot(run_dir, values)
for line in plot_page:
	if 'idl_images' in line and 'gif' in line:
		image_name = line.split('/idl_images/')[1].split('"')[0]

		local_name = ''
		if args.prefix:
			local_name = args.prefix
		else:
			local_name = args.run_id
		local_name += '_' + requested_variable + '_' + requested_time + '.gif'

		url = 'http://ccmc.gsfc.nasa.gov/idl_images/' + image_name
		print('Downloading file', url, 'as', local_name, flush = True)
		req = request.Request(url, headers = default_headers)
		with request.urlopen(req) as response, open(local_name, 'wb') as image_file:
		    shutil.copyfileobj(response, image_file)
