/*
Tests performance of muParserX as used by PAMHD.

Copyright 2014 Ilja Honkonen
All rights reserved.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.
*/


#include "array"
#include "chrono"
#include "iostream"
#include "vector"

#include "mpParser.h"

using namespace std;
using namespace std::chrono;

int main()
{
	mup::ParserX parser(mup::pckCOMMON | mup::pckNON_COMPLEX | mup::pckMATRIX | mup::pckUNIT);
	mup::Value value1{0};
	mup::Variable variable1{&value1};
	parser.DefineVar("v1", variable1);
	parser.SetExpr("v1 + 1");

	const auto time_start = high_resolution_clock::now();

	for (int i = 0; i < 1000000; i++) {
		value1 = i;
		if (parser.Eval().GetInteger() != i + 1) {
			std::cerr <<  __FILE__ << "(" << __LINE__<< ") "
				<< "Wrong value for variable: " << parser.Eval().GetInteger()
				<< ", should be: " << i + 1
				<< endl;
			abort();
		}
	}

	const auto time_end = high_resolution_clock::now();

	const auto total = duration_cast<duration<double>>(time_end - time_start).count();
	if (total > 0.01) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< ") "
			<< "muParserX performance test took too long (s): " << total
			<< endl;
		abort();
	}

	return 0;
}
