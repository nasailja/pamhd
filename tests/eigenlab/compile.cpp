/*
Tests compilation with of EigenLab.

Copyright 2014, 2015 Ilja Honkonen
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
*/


#include "iostream"

#include "EigenLab.h"
#include "Eigen/Core"

using namespace std;

int main()
{
	EigenLab::ParserXd parser;
	Eigen::MatrixXd variable1(1, 1);
	variable1 << 0;
	parser.var("v1").setShared(variable1);
	parser.eval("v1 = v1 + 1");

	Eigen::MatrixXd variable2(3, 1);
	variable2 << 1, 2, 3;
	parser.var("v2").setShared(variable2);
	parser.eval("v1 = sqrt(v2(0)*v2(0) + v2(1)*v2(1) + v2(2)*v2(2))");

	EigenLab::ParserXi parser2;
	Eigen::MatrixXi variable3(1, 1);
	variable3 << 0;
	parser2.var("v3").setShared(variable3);
	parser2.eval("v3 = v3 + 1");

	Eigen::MatrixXi variable4(3, 1);
	variable4 << 1, 2, 3;
	parser2.var("v4").setShared(variable4);
	parser2.eval("v3 = sqrt(v4(0)*v4(0) + v4(1)*v4(1) + v4(2)*v4(2))");

	return 0;
}
