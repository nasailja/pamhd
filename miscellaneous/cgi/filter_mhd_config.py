#! /usr/bin/env python

'''
Performs some filtering on given MHD test program configuration file.

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

from argparse import ArgumentParser
from os import makedirs
from os.path import dirname, exists
from sys import exit

parser = ArgumentParser(description = 'Performs some filtering on given MHD test program configuration file.')

parser.add_argument('--infile', metavar = 'I', required = True, help = 'Filter configuration in file named I')
parser.add_argument('--outfile', metavar = 'O', required = True, help = 'Write filtered configuration to file named O, creating any required directories')

args = parser.parse_args()

infile = open(args.infile, 'r')

outdir = dirname(args.outfile)
if len(outdir) > 0 and not exists(outdir):
	makedirs(outdir)
outfile = open(args.outfile, 'w')

for line in infile:
	line = line.strip().split('#')[0]
	if len(line) == 0:
		continue
	if line.count('output-directory') > 0:
		continue
	if line.count('=') == 0:
		if not (line.count('[') > 0 and line.count(']') > 0):
			continue

	outfile.write(line + '\n')

outfile.close()
infile.close()
exit(0)
