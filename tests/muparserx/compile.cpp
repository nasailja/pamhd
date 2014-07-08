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
	parser.Eval().GetFloat();

	mup::Value value2{3, 0};
	mup::Variable variable2{&value2};
	parser.DefineVar("v2", variable2);
	parser.SetExpr("sqrt(v2[0]*v2[0] + v2[1]*v2[1] + v2[2]*v2[2])");
	parser.Eval().GetFloat();

	return 0;
}
