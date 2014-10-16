#! /bin/sh

# if given, override RUN with first command line argument
RUN=${1:-}

make

mkdir -p \
    results/shock_tube1_x \
    results/shock_tube1_-x \
    results/shock_tube1_y \
    results/shock_tube1_-y \
    results/shock_tube1_z \
    results/shock_tube1_-z \
    results/bow_shock_reference \
    results/kelvin-helmholtz_reference \
    results/magnetohydrodynamic/kelvin-helmholtz/reference_2d/roe_athena

echo
echo Running tests...

echo shock_tube1_x
$RUN ./test.exe \
    --config-file config_files/shock_tube1_x_general.cfg \
    --boundary-file config_files/shock_tube1_x_boundaries.cfg

echo shock_tube1_-x
$RUN ./test.exe \
    --config-file config_files/shock_tube1_-x_general.cfg \
    --boundary-file config_files/shock_tube1_-x_boundaries.cfg

echo shock_tube1_y
$RUN ./test.exe \
    --config-file config_files/shock_tube1_y_general.cfg \
    --boundary-file config_files/shock_tube1_y_boundaries.cfg

echo shock_tube1_-y
$RUN ./test.exe \
    --config-file config_files/shock_tube1_-y_general.cfg \
    --boundary-file config_files/shock_tube1_-y_boundaries.cfg

echo shock_tube1_z
$RUN ./test.exe \
    --config-file config_files/shock_tube1_z_general.cfg \
    --boundary-file config_files/shock_tube1_z_boundaries.cfg

echo shock_tube1_-z
$RUN ./test.exe \
    --config-file config_files/shock_tube1_-z_general.cfg \
    --boundary-file config_files/shock_tube1_-z_boundaries.cfg

echo kelvin-helmholtz_reference
$RUN ./test.exe \
    --config-file config_files/kelvin-helmholtz_reference_general.cfg \
    --boundary-file config_files/kelvin-helmholtz_reference_boundaries.cfg

echo bow_shock_reference
$RUN ./test.exe \
    --config-file config_files/bow_shock_reference_general.cfg \
    --boundary-file config_files/bow_shock_reference_boundaries.cfg

echo mhd kelvin-helmholtz reference 2d roe
$RUN ./test.exe \
    --config-file config_files/magnetohydrodynamic/kelvin-helmholtz/reference_2d/general.cfg \
    --boundary-file config_files/magnetohydrodynamic/kelvin-helmholtz/reference_2d/boundaries.cfg


echo
echo Plotting results with gnuplot

echo shock_tube1_x
$RUN ./mhd2gnuplot.exe results/shock_tube1_x/*dc

echo shock_tube1_-x
$RUN ./mhd2gnuplot.exe results/shock_tube1_-x/*dc

echo shock_tube1_y
$RUN ./mhd2gnuplot.exe results/shock_tube1_y/*dc

echo shock_tube1_-y
$RUN ./mhd2gnuplot.exe results/shock_tube1_-y/*dc

echo shock_tube1_z
$RUN ./mhd2gnuplot.exe results/shock_tube1_z/*dc

echo shock_tube1_-z
$RUN ./mhd2gnuplot.exe results/shock_tube1_-z/*dc

echo kelvin-helmholtz_reference
$RUN ./mhd2gnuplot.exe results/kelvin-helmholtz_reference/*dc

echo bow_shock_reference
$RUN ./mhd2gnuplot.exe results/bow_shock_reference/*dc

echo mhd kelvin-helmholtz reference 2d roe
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/kelvin-helmholtz/reference_2d/roe_athena/*dc
