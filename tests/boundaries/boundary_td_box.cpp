/*
Tests time-dependent boundary condition class of PAMHD with box geometry.

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
#include "set"
#include "string"
#include "vector"

#include "boost/optional.hpp"
#include "boost/program_options.hpp"
#ifdef HAVE_EIGEN
#include "Eigen/Core" // must be included before boundary_time_dependent.hpp
#endif

#include "boundaries/boundary_time_dependent.hpp"
#include "boundaries/boxes.hpp"


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


#ifdef HAVE_EIGEN
#define MAKE_VEC(x, y, z) Eigen::Vector3d(x, y, z)
#define MAKE_BOX(a, b, c, d, e, f) MAKE_VEC(a, b, c), MAKE_VEC(d, e, f)
#else
#define MAKE_VEC(x, y, z) std::array<double, 3>{{x, y, z}}
#define MAKE_BOX(a, b, c, d, e, f) MAKE_VEC(a, b, c), MAKE_VEC(d, e, f)
#endif


int main(int argc, char* argv[])
{
	Boundary_Time_Dependent<
		Boxes<
			#ifdef HAVE_EIGEN
			Eigen::Vector3d
			#else
			std::array<double, 3>
			#endif
		>,
		int,
		double,
		Mass_Density,
		Momentum_Density
	> bdy1, bdy2;

	// set default values
	bdy1.set_number_of_instances(1);
	bdy1.set_time_stamp(0, 0);
	bdy1.geometries.set_geometry(0, MAKE_BOX(0, 0, 0, 1, 1, 1));
	bdy1.boundary_data.set_expression(
		Mass_Density(),
		0,
		"r[0] + r[1] + r[2] + t"
	);
	bdy1.boundary_data.set_expression(
		Momentum_Density(),
		0,
		"{r[0] + t, r[1] + t, r[2] + t}"
	);

	bdy2.set_number_of_instances(2);
	bdy2.set_time_stamp(0, 0);
	bdy2.set_time_stamp(1, 1);
	bdy2.geometries.set_geometry(0, MAKE_BOX(0, 0, 0, 1, 1, 1));
	bdy2.geometries.set_geometry(1, MAKE_BOX(1, 1, 1, 2, 2, 2));
	bdy2.boundary_data.set_expression(
		Mass_Density(),
		0,
		"r[0] + r[1] + r[2] + t + 1"
	);
	bdy2.boundary_data.set_expression(
		Momentum_Density(),
		0,
		"{r[0] + t + 1, r[1] + t + 1, r[2] + t + 1}"
	);
	bdy2.boundary_data.set_expression(
		Mass_Density(),
		1,
		"r[0] + r[1] + r[2] + t + 2"
	);
	bdy2.boundary_data.set_expression(
		Momentum_Density(),
		1,
		"{r[0] + t + 2, r[1] + t + 2, r[2] + t + 2}"
	);


	std::string boundary_condition_file_name("");

	boost::program_options::options_description
		general_options(
			"Usage: program_name [options], where options are"
		),
		boundary_options(
			"Options related to initial condition"
		);

	general_options.add_options()
		("help", "Print options and their descriptions")
		("boundary-help", "Print help for options of boundary condition")
		("boundary-file",
			boost::program_options::value<std::string>(&boundary_condition_file_name)
				->default_value(boundary_condition_file_name),
			"Read boundary conditions from file arg (do not read if empty string)");

	boost::program_options::variables_map option_variables;
	boost::program_options::store(
		boost::program_options::parse_command_line(
			argc,
			argv,
			general_options
		),
		option_variables
	);
	boost::program_options::notify(option_variables);

	if (option_variables.count("help") > 0) {
		std::cout << general_options << std::endl;
		return EXIT_SUCCESS;
	}

	bdy1.add_options("bdy1.", boundary_options);
	bdy2.add_options("bdy2.", boundary_options);

	if (option_variables.count("boundary-help") > 0) {
		std::cout << boundary_options << std::endl;
		return EXIT_SUCCESS;
	}

	if (boundary_condition_file_name != "") {
		boost::program_options::store(
			boost::program_options::parse_config_file<char>(
				boundary_condition_file_name.c_str(),
				boundary_options
			),
			option_variables
		);
		boost::program_options::notify(option_variables);
	}


	/*
	Test boundary 1
	*/

	// cell classification logic
	if (not bdy1.add_cell(-1, 0, MAKE_BOX(0, 0, 0, 1, 1, 1))) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< ")" << std::endl;
		abort();
	}
	if (not bdy1.add_cell(0, 1, MAKE_BOX(0, 0, 0, 1, 1, 1))) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< ")" << std::endl;
		abort();
	}
	if (not bdy1.add_cell(1, 2, MAKE_BOX(0, 0, 0, 1, 1, 1))) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< ")" << std::endl;
		abort();
	}
	if (bdy1.get_cells().size() != 3) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< ")" << std::endl;
		abort();
	}

	if (bdy1.add_cell(-1, 3, MAKE_BOX(-2, -2, -2, -1, -1, -1))) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< ")" << std::endl;
		abort();
	}
	if (bdy1.add_cell(0, 4, MAKE_BOX(2, 2, 2, 3, 3, 3))) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< ")" << std::endl;
		abort();
	}
	if (bdy1.add_cell(1, 5, MAKE_BOX(-2, 0, 0, -1, 1, 1))) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< ")" << std::endl;
		abort();
	}

	// expression logic
	if (bdy1.get_data(Mass_Density(), {{1, 2, 3}}, -1) != 5) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< ")" << std::endl;
		abort();
	}
	if (bdy1.get_data(Mass_Density(), {{1, 2, 3}}, 0) != 6) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< ")" << std::endl;
		abort();
	}
	if (bdy1.get_data(Mass_Density(), {{1, 2, 3}}, 1) != 7) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< ")" << std::endl;
		abort();
	}

	if (bdy1.get_data(Momentum_Density(), {{1, 2, 3}}, -1)[0] != 0) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< ")" << std::endl;
		abort();
	}
	if (bdy1.get_data(Momentum_Density(), {{1, 2, 3}}, 0)[1] != 2) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< ")" << std::endl;
		abort();
	}
	if (bdy1.get_data(Momentum_Density(), {{1, 2, 3}}, 1)[2] != 4) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< ")" << std::endl;
		abort();
	}


	/*
	Test boundary 2
	*/

	// cell classification logic
	if (not bdy2.add_cell(-1, 0, MAKE_BOX(0, 0, 0, 1, 1, 1))) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< ")" << std::endl;
		abort();
	}
	if (not bdy2.add_cell(0.5, 1, MAKE_BOX(0, 0, 0, 1, 1, 1))) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< ")" << std::endl;
		abort();
	}
	if (bdy2.add_cell(1.5, 2, MAKE_BOX(0, 0, 0, 1, 1, 1))) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< ")" << std::endl;
		abort();
	}
	if (bdy2.add_cell(2.5, 3, MAKE_BOX(0, 0, 0, 1, 1, 1))) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< ")" << std::endl;
		abort();
	}
	if (bdy2.get_cells().size() != 2) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< ")" << std::endl;
		abort();
	}
	bdy2.clear_cells();

	if (bdy2.add_cell(-1, 0, MAKE_BOX(1, 1, 1, 2, 2, 2))) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< ")" << std::endl;
		abort();
	}
	if (bdy2.add_cell(0.5, 1, MAKE_BOX(1, 1, 1, 2, 2, 2))) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< ")" << std::endl;
		abort();
	}
	if (not bdy2.add_cell(1.5, 2, MAKE_BOX(1, 1, 1, 2, 2, 2))) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< ")" << std::endl;
		abort();
	}
	if (not bdy2.add_cell(2.5, 3, MAKE_BOX(1, 1, 1, 2, 2, 2))) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< ")" << std::endl;
		abort();
	}
	if (bdy2.get_cells().size() != 2) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< ")" << std::endl;
		abort();
	}
	bdy2.clear_cells();

	// expression logic
	if (bdy2.get_data(Mass_Density(), {{1, 2, 3}}, -1) != 6) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< ")" << std::endl;
		abort();
	}
	if (bdy2.get_data(Mass_Density(), {{1, 2, 3}}, 0) != 7) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< ")" << std::endl;
		abort();
	}
	if (bdy2.get_data(Mass_Density(), {{1, 2, 3}}, 0.5) != 7.5) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< ")" << std::endl;
		abort();
	}
	if (bdy2.get_data(Mass_Density(), {{1, 2, 3}}, 1.5) != 9.5) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< ")" << std::endl;
		abort();
	}
	if (bdy2.get_data(Mass_Density(), {{1, 2, 3}}, 3) != 11) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< ")" << std::endl;
		abort();
	}

	if (bdy2.get_data(Momentum_Density(), {{1, 2, 3}}, -1)[0] != 1) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< ")" << std::endl;
		abort();
	}
	if (bdy2.get_data(Momentum_Density(), {{1, 2, 3}}, 0)[1] != 3) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< ")" << std::endl;
		abort();
	}
	if (bdy2.get_data(Momentum_Density(), {{1, 2, 3}}, 0.5)[2] != 4.5) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< ")" << std::endl;
		abort();
	}
	if (bdy2.get_data(Momentum_Density(), {{1, 2, 3}}, 1.5)[2] != 6.5) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< ")" << std::endl;
		abort();
	}
	if (bdy2.get_data(Momentum_Density(), {{1, 2, 3}}, 3)[2] != 8) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< ")" << std::endl;
		abort();
	}


	return EXIT_SUCCESS;
}
