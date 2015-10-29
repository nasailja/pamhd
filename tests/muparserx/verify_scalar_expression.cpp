/*
Verifies that scalar math expression given on stdin is well formed.

Copyright 2015 Ilja Honkonen
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


#include "iostream"
#include "string"

#include "mpParser.h"

using namespace std;

/*!
Returns 0 if expression given on stdin is well formed, return 1 otherwise.

Expression can use scalar variables t and J and vector variables
r (e.g. r[0] or r[2]), cells and volume.
*/
int main()
{
	mup::ParserX parser;
	mup::Value
		time_val{0.0},
		pos_val{3, 0},
		cur_val{0.0},
		cells_val{3, 0},
		vol_val{3, 0};
	mup::Variable
		time{&time_val},
		pos{&pos_val},
		cur{&cur_val},
		cells{&cells_val},
		vol{&vol_val};
	parser.DefineVar("t", time);
	parser.DefineVar("r", pos);
	parser.DefineVar("J", cur);
	parser.DefineVar("cells", cells);
	parser.DefineVar("volume", vol);

	std::string expression;
	std::cin >> expression;
	parser.SetExpr(expression);

	try {
		const double evaluated = parser.Eval().GetFloat();
		(void)evaluated;
	} catch(mup::ParserError &error) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< ") "
			<< "Could not evaluate expression '"
			<< parser.GetExpr() << "': " << error.GetMsg()
			<< endl;
		return 1;
	}

	return 0;
}
