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

	return 0;
}
