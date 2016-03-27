/*
Tests EigenLab.

Copyright 2014, 2015, 2016 Ilja Honkonen
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


#include "array"
#include "iostream"
#include "vector"

#include "EigenLab.h"
#include "Eigen/Core"

using namespace std;

int main()
{
	EigenLab::ParserXd parser;

	Eigen::MatrixXd variable1(1, 1);
	parser.var("v1").setShared(variable1);

	double sum = 0;
	for (const auto i: {1, 3, 5, 7, 11}) {
		variable1(0, 0) = i;
		parser.eval("v1 = v1 + 1");
		sum += variable1(0, 0);
	}
	if (sum != 32) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< ") "
			<< "Wrong value for sum: " << sum
			<< ", should be 32"
			<< endl;
		abort();
	}


	Eigen::MatrixXd variable2(3, 1);
	variable2 << 0, 0, 0;
	parser.var("v2").setShared(variable2);

	const std::vector<std::array<double, 3>> test_vectors{
		{{1, 3, 5}},
		{{7, 11, 13}},
		{{-7, -11, -13}}
	};

	for (const auto test_vector: test_vectors) {
		for (size_t i = 0; i < test_vector.size(); i++) {
			variable2(i) = test_vector[i];
		}

		parser.eval("v2 = v2.^2");

		for (int i = 0; i < variable2.size(); i++) {
			if (test_vector[i]*test_vector[i] != variable2(i)) {
				std::cerr <<  __FILE__ << "(" << __LINE__<< ") "
					<< "Invalid evaluated value at index "
					<< i << ": " << variable2(i)
					<< ", should be: " << test_vector[i]*test_vector[i]
					<< endl << variable2
					<< endl;
				abort();
			}
		}
	}

	return 0;
}
