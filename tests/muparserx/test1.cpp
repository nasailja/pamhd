/*
Test for muParserX.

Copyright 2014, 2015 Ilja Honkonen
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
	mup::ParserX parser(mup::pckCOMMON | mup::pckNON_COMPLEX | mup::pckMATRIX | mup::pckUNIT);
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
		abort();
	}


	mup::Value value2{3, 0};
	mup::Variable variable2{&value2};
	parser.DefineVar("v2", variable2);
	parser.SetExpr("{v2[0]*v2[0], v2[1]*v2[1], v2[2]*v2[2]}");

	// values of variables used in expression above
	const std::vector<std::array<double, 3>> test_vectors{
		{{1, 3, 5}},
		{{7, 11, 13}},
		{{-7, -11, -13}}
	};

	for (const auto test_vector: test_vectors) {
		for (size_t i = 0; i < test_vector.size(); i++) {
			value2.At(i) = test_vector[i];
		}

		const auto& evaluated
			= [&parser](){
				try {
					return parser.Eval().GetArray();
				} catch(mup::ParserError &error) {
					std::cerr <<  __FILE__ << "(" << __LINE__<< ") "
						<< "Could not evaluate expression "
						<< parser.GetExpr() << ": " << error.GetMsg()
						<< endl;
					abort();
				}
			}();

		if (evaluated.GetRows() != 1) {
			std::cerr <<  __FILE__ << "(" << __LINE__<< ") "
				<< "Invalid number of rows in expression " << parser.GetExpr()
				<< ": " << evaluated.GetRows() << ", should be 1"
				<< endl;
			abort();
		}
		if (evaluated.GetCols() != 3) {
			std::cerr <<  __FILE__ << "(" << __LINE__<< ") "
				<< "Invalid number of columns in expression " << parser.GetExpr()
				<< ": " << evaluated.GetCols() << ", should be 3"
				<< endl;
			abort();
		}

		for (int i = 0; i < evaluated.GetCols(); i++) {
			const double value = evaluated.At(i).GetFloat();
			if (test_vector[i]*test_vector[i] != value) {
				std::cerr <<  __FILE__ << "(" << __LINE__<< ") "
					<< "Invalid evaluated value at index "
					<< i << ": " << value
					<< ", should be: " << test_vector[i]*test_vector[i]
					<< endl;
				abort();
			}
		}
	}

	return 0;
}
