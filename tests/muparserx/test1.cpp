/*
Test for muParserX.

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
#include "iostream"
#include "vector"

#include "mpParser.h"

using namespace std;

int main()
{
	mup::ParserX parser(mup::pckALL_NON_COMPLEX);
	mup::Value value1{0};
	mup::Variable variable1{&value1};
	parser.DefineVar("v1", variable1);
	parser.SetExpr("v1 + 1");

	double sum = 0;
	for (const auto i: {1, 3, 5, 7, 11}) {
		value1 = i;
		sum += parser.Eval().GetInteger();
	}
	if (sum != 32) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< ") "
			<< "Wrong value for sum: " << sum
			<< ", should be 32"
			<< endl;
	}


	mup::Value value2{3, 0};
	mup::Variable variable2{&value2};
	parser.DefineVar("v2", variable2);
	parser.SetExpr("sqrt(v2[0]*v2[0] + v2[1]*v2[1] + v2[2]*v2[2])");

	sum = 0;
	const std::vector<std::array<double, 3>> test_vectors{
		{1, 3, 5},
		{7, 11, 13}
	};
	for (const auto v: test_vectors) {
		for (size_t i = 0; i < v.size(); i++) {
			value2.At(i) = v[i];
		}
		sum += parser.Eval().GetFloat();
	}
	if (int(sum) != 24) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< ") "
			<< "Wrong value for int(sum): " << sum
			<< ", should be 24"
			<< endl;
	}

	return 0;
}
