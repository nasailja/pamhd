/*
Tests a PAMHD boundary with many variables as data.

Copyright 2014 Ilja Honkonen
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
#include "cstdlib"
#include "iostream"
#include "string"

#include "boost/program_options.hpp"
#ifdef HAVE_EIGEN
#include "Eigen/Core"
#endif

#include "boundaries/boundary.hpp"


using namespace pamhd::boundaries;


struct Mass_Density {
	using data_type = double;

	static const std::string get_name()
	{
		return {"mass density"};
	}
	static const std::string get_option_name()
	{
		return {"mass-density"};
	}
	static const std::string get_option_help()
	{
		return {"Mass density (kg / m^3)"};
	}
};

struct Momentum_Density {
	#ifdef HAVE_EIGEN
	using data_type = Eigen::Vector3d;
	#else
	using data_type = std::array<double, 3>;
	#endif

	static const std::string get_name()
	{
		return {"momentum density"};
	}
	static const std::string get_option_name()
	{
		return {"momentum-density"};
	}
	static const std::string get_option_help()
	{
		return {"Momentum density (kg s / m^2)"};
	}
};

struct Magnetic_Field {
	#ifdef HAVE_EIGEN
	using data_type = Eigen::Vector3d;
	#else
	using data_type = std::array<double, 3>;
	#endif

	static const std::string get_name()
	{
		return {"magnetic field"};
	}
	static const std::string get_option_name()
	{
		return {"magnetic-field"};
	}
	static const std::string get_option_help()
	{
		return {"Magnetic field (T)"};
	}
};


#ifdef HAVE_EIGEN
#define MAKE_VEC(x, y, z) Eigen::Vector3d(x, y, z)
#else
#define MAKE_VEC(x, y, z) std::array<double, 3>{x, y, z}
#endif


using Box_Boundary = Boundary<
	Box<
		#ifdef HAVE_EIGEN
		Eigen::Vector3d
		#else
		std::array<double, 3>
		#endif
	>,
	int,
	Mass_Density,
	Momentum_Density,
	Magnetic_Field
>;

using Sphere_Boundary = Boundary<
	Sphere<
		#ifdef HAVE_EIGEN
		Eigen::Vector3d,
		#else
		std::array<double, 3>,
		#endif
		double
	>,
	int,
	Mass_Density,
	Momentum_Density,
	Magnetic_Field
>;


int main(int argc, char* argv[])
{
	#ifdef HAVE_EIGEN
	Eigen::Vector3d result;
	#else
	std::array<double, 3> result;
	#endif

	Box_Boundary box;
	Sphere_Boundary sphere;

	// set defaults
	box.geometry.set_geometry({-1, -2, -3}, {1, 2, 3});
	box.boundary_data.set_expression(
		Mass_Density(),
		"r[0] + r[1] + r[2] + t"
	);
	box.boundary_data.set_expression(
		Momentum_Density(),
		"{r[0] + t, r[1] + t, r[2] + t}"
	);
	box.boundary_data.set_expression(
		Magnetic_Field(),
		"{r[0] + t + 1, r[1] + t + 1, r[2] + t + 1}"
	);

	sphere.geometry.set_geometry({1, 2, 3}, 1);
	sphere.boundary_data.set_expression(
		Mass_Density(),
		"r[0] + r[1] + r[2] + t + 1"
	);
	sphere.boundary_data.set_expression(
		Momentum_Density(),
		"{r[0] + t + 1, r[1] + t + 1, r[2] + t + 1}"
	);
	sphere.boundary_data.set_expression(
		Magnetic_Field(),
		"{r[0] + t + 2, r[1] + t + 2, r[2] + t + 2}"
	);


	boost::program_options::options_description options(
		"Usage: program_name [options], where options are:"
	);

	options.add_options()
		("help", "Print options and their descriptions");
	box.add_options("box.", options);
	sphere.add_options("sphere.", options);

	boost::program_options::options_description boundary_options(
		"Options for initial and boundary conditions"
	);

	boost::program_options::variables_map option_variables;
	boost::program_options::store(
		boost::program_options::parse_command_line(
			argc,
			argv,
			options
		),
		option_variables
	);
	boost::program_options::notify(option_variables);

	if (option_variables.count("help") > 0) {
		std::cout << options << std::endl;
		return EXIT_SUCCESS;
	}

	/*
	Test box boundary
	*/
	if (not box.add_cell(-2, MAKE_VEC(0, 0, 0), MAKE_VEC(0.5, 1.5, 2.5))) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
			<< "Cell not added to box."
			<< std::endl;
		return EXIT_FAILURE;
	}
	if (box.add_cell(-1, MAKE_VEC(-2, -3, -4), MAKE_VEC(-1.1, -2.1, -3.1))) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
			<< "Cell added to box."
			<< std::endl;
		return EXIT_FAILURE;
	}
	if (not box.add_cell(0, MAKE_VEC(0, 0, 0), MAKE_VEC(1.5, 2.5, 3.5))) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
			<< "Cell not added to box."
			<< std::endl;
		return EXIT_FAILURE;
	}

	if (box.get_cells().size() != 2) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
			<< "Invalid number of cells in boundary."
			<< std::endl;
		return EXIT_FAILURE;
	}
	if (box.get_cells()[0] != -2) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
			<< "Invalid first cell in boundary."
			<< std::endl;
		return EXIT_FAILURE;
	}
	if (box.get_cells()[1] != 0) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
			<< "Invalid second cell in boundary."
			<< std::endl;
		return EXIT_FAILURE;
	}

	double mass = box.boundary_data.get_data(Mass_Density(), {1, 2, 3}, 4);
	if (mass != 10) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< ")" << std::endl;
		return EXIT_FAILURE;
	}

	result = box.boundary_data.get_data(Momentum_Density(), {1, 2, 3}, 4);
	if (result[0] != 5) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< ")" << std::endl;
		return EXIT_FAILURE;
	}
	if (result[1] != 6) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< ")" << std::endl;
		return EXIT_FAILURE;
	}
	if (result[2] != 7) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< ")" << std::endl;
		return EXIT_FAILURE;
	}

	result = box.boundary_data.get_data(Magnetic_Field(), {2, 4, 6}, 8);
	if (result[0] != 11) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< ")" << std::endl;
		return EXIT_FAILURE;
	}
	if (result[1] != 13) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< ")" << std::endl;
		return EXIT_FAILURE;
	}
	if (result[2] != 15) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< ")" << std::endl;
		return EXIT_FAILURE;
	}

	/*
	Test sphere boundary
	*/
	if (not sphere.add_cell(-2, MAKE_VEC(0, 0, 0), MAKE_VEC(0.5, 1.5, 2.5))) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
			<< "Cell not added to sphere boundary."
			<< std::endl;
		return EXIT_FAILURE;
	}
	if (sphere.add_cell(-1, MAKE_VEC(0, 0, 0), MAKE_VEC(0.1, 0.1, 0.1))) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
			<< "Cell added to sphere boundary."
			<< std::endl;
		return EXIT_FAILURE;
	}
	if (not sphere.add_cell(0, MAKE_VEC(0, 0, 0), MAKE_VEC(1.5, 2.5, 3.5))) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
			<< "Cell not added to sphere boundary."
			<< std::endl;
		return EXIT_FAILURE;
	}

	if (sphere.get_cells().size() != 2) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
			<< "Invalid number of cells in boundary."
			<< std::endl;
		return EXIT_FAILURE;
	}
	if (sphere.get_cells()[0] != -2) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
			<< "Invalid first cell in boundary."
			<< std::endl;
		return EXIT_FAILURE;
	}
	if (sphere.get_cells()[1] != 0) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
			<< "Invalid second cell in boundary."
			<< std::endl;
		return EXIT_FAILURE;
	}

	mass = sphere.boundary_data.get_data(Mass_Density(), {1, 2, 3}, 4);
	if (mass != 11) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< ")" << std::endl;
		return EXIT_FAILURE;
	}

	result = sphere.boundary_data.get_data(Momentum_Density(), {1, 2, 3}, 4);
	if (result[0] != 6) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< ")" << std::endl;
		return EXIT_FAILURE;
	}
	if (result[1] != 7) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< ")" << std::endl;
		return EXIT_FAILURE;
	}
	if (result[2] != 8) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< ")" << std::endl;
		return EXIT_FAILURE;
	}

	result = sphere.boundary_data.get_data(Magnetic_Field(), {2, 4, 6}, 8);
	if (result[0] != 12) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< ")" << std::endl;
		return EXIT_FAILURE;
	}
	if (result[1] != 14) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< ")" << std::endl;
		return EXIT_FAILURE;
	}
	if (result[2] != 16) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< ")" << std::endl;
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}
