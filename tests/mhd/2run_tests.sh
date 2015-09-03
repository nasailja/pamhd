#! /bin/sh -e

# if given, override RUN with first command line argument
RUN=${1:-}

echo
echo Compiling programs
echo

make


echo
echo Running tests...
echo

# two-fluid magnetic field advection
echo magnetic_field_advection 1 1d +x hll_athena
$RUN ./two_test.exe \
    --solver-mhd hll_athena \
    --output-directory results/magnetohydrodynamic/magnetic_field_advection/1d/+x/2/1/hll_athena/ \
    --config-file config_files/magnetohydrodynamic/magnetic_field_advection/1d/+x/two_fluid1.cfg
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/magnetic_field_advection/1d/+x/2/1/hll_athena/mhd*dc 2> /dev/null
$RUN ./hd2gnuplot.exe results/magnetohydrodynamic/magnetic_field_advection/1d/+x/2/1/hll_athena/hd*dc 2> /dev/null

echo magnetic_field_advection 1 1d +x hlld_athena
$RUN ./two_test.exe \
    --solver-mhd hlld_athena \
    --output-directory results/magnetohydrodynamic/magnetic_field_advection/1d/+x/2/1/hlld_athena/ \
    --config-file config_files/magnetohydrodynamic/magnetic_field_advection/1d/+x/two_fluid1.cfg
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/magnetic_field_advection/1d/+x/2/1/hlld_athena/mhd*dc 2> /dev/null
$RUN ./hd2gnuplot.exe results/magnetohydrodynamic/magnetic_field_advection/1d/+x/2/1/hlld_athena/hd*dc 2> /dev/null

echo magnetic_field_advection 1 1d +x roe_athena
$RUN ./two_test.exe \
    --solver-mhd roe_athena \
    --output-directory results/magnetohydrodynamic/magnetic_field_advection/1d/+x/2/1/roe_athena/ \
    --config-file config_files/magnetohydrodynamic/magnetic_field_advection/1d/+x/two_fluid1.cfg
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/magnetic_field_advection/1d/+x/2/1/roe_athena/mhd*dc 2> /dev/null
$RUN ./hd2gnuplot.exe results/magnetohydrodynamic/magnetic_field_advection/1d/+x/2/1/roe_athena/hd*dc 2> /dev/null

echo magnetic_field_advection 1 1d -z hll_athena
$RUN ./two_test.exe \
    --solver-mhd hll_athena \
    --output-directory results/magnetohydrodynamic/magnetic_field_advection/1d/-z/2/1/hll_athena/ \
    --config-file config_files/magnetohydrodynamic/magnetic_field_advection/1d/-z/two_fluid1.cfg
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/magnetic_field_advection/1d/-z/2/1/hll_athena/mhd*dc 2> /dev/null
$RUN ./hd2gnuplot.exe results/magnetohydrodynamic/magnetic_field_advection/1d/-z/2/1/hll_athena/hd*dc 2> /dev/null

echo magnetic_field_advection 1 1d -z hlld_athena
$RUN ./two_test.exe \
    --solver-mhd hlld_athena \
    --output-directory results/magnetohydrodynamic/magnetic_field_advection/1d/-z/2/1/hlld_athena/ \
    --config-file config_files/magnetohydrodynamic/magnetic_field_advection/1d/-z/two_fluid1.cfg
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/magnetic_field_advection/1d/-z/2/1/hlld_athena/mhd*dc 2> /dev/null
$RUN ./hd2gnuplot.exe results/magnetohydrodynamic/magnetic_field_advection/1d/-z/2/1/hlld_athena/hd*dc 2> /dev/null

echo magnetic_field_advection 1 1d -z roe_athena
$RUN ./two_test.exe \
    --solver-mhd roe_athena \
    --output-directory results/magnetohydrodynamic/magnetic_field_advection/1d/-z/2/1/roe_athena/ \
    --config-file config_files/magnetohydrodynamic/magnetic_field_advection/1d/-z/two_fluid1.cfg
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/magnetic_field_advection/1d/-z/2/1/roe_athena/mhd*dc 2> /dev/null
$RUN ./hd2gnuplot.exe results/magnetohydrodynamic/magnetic_field_advection/1d/-z/2/1/roe_athena/hd*dc 2> /dev/null

echo magnetic_field_advection 2 1d +x hll_athena
$RUN ./two_test.exe \
    --solver-mhd hll_athena \
    --output-directory results/magnetohydrodynamic/magnetic_field_advection/1d/+x/2/2/hll_athena/ \
    --config-file config_files/magnetohydrodynamic/magnetic_field_advection/1d/+x/two_fluid2.cfg
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/magnetic_field_advection/1d/+x/2/2/hll_athena/mhd*dc 2> /dev/null
$RUN ./hd2gnuplot.exe results/magnetohydrodynamic/magnetic_field_advection/1d/+x/2/2/hll_athena/hd*dc 2> /dev/null

echo magnetic_field_advection 2 1d +x hlld_athena
$RUN ./two_test.exe \
    --solver-mhd hlld_athena \
    --output-directory results/magnetohydrodynamic/magnetic_field_advection/1d/+x/2/2/hlld_athena/ \
    --config-file config_files/magnetohydrodynamic/magnetic_field_advection/1d/+x/two_fluid2.cfg
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/magnetic_field_advection/1d/+x/2/2/hlld_athena/mhd*dc 2> /dev/null
$RUN ./hd2gnuplot.exe results/magnetohydrodynamic/magnetic_field_advection/1d/+x/2/2/hlld_athena/hd*dc 2> /dev/null

echo magnetic_field_advection 2 1d +x roe_athena
$RUN ./two_test.exe \
    --solver-mhd roe_athena \
    --output-directory results/magnetohydrodynamic/magnetic_field_advection/1d/+x/2/2/roe_athena/ \
    --config-file config_files/magnetohydrodynamic/magnetic_field_advection/1d/+x/two_fluid2.cfg
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/magnetic_field_advection/1d/+x/2/2/roe_athena/mhd*dc 2> /dev/null
$RUN ./hd2gnuplot.exe results/magnetohydrodynamic/magnetic_field_advection/1d/+x/2/2/roe_athena/hd*dc 2> /dev/null

# shock tube
echo shock_tube mhd 0 1d +x hll_athena
$RUN ./two_test.exe \
    --solver-mhd hll_athena \
    --output-directory results/magnetohydrodynamic/shock_tube/1d/+x/2/0/hll_athena/ \
    --config-file config_files/magnetohydrodynamic/shock_tube/1d/+x/two_fluid0.cfg
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/shock_tube/1d/+x/2/0/hll_athena/mhd*dc 2> /dev/null
$RUN ./hd2gnuplot.exe results/magnetohydrodynamic/shock_tube/1d/+x/2/0/hll_athena/hd*dc 2> /dev/null

echo shock_tube mhd 0 1d +x hlld_athena
$RUN ./two_test.exe \
    --solver-mhd hlld_athena \
    --output-directory results/magnetohydrodynamic/shock_tube/1d/+x/2/0/hlld_athena/ \
    --config-file config_files/magnetohydrodynamic/shock_tube/1d/+x/two_fluid0.cfg
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/shock_tube/1d/+x/2/0/hlld_athena/mhd*dc 2> /dev/null
$RUN ./hd2gnuplot.exe results/magnetohydrodynamic/shock_tube/1d/+x/2/0/hlld_athena/hd*dc 2> /dev/null

echo shock_tube mhd 0 1d +x roe_athena
$RUN ./two_test.exe \
    --solver-mhd roe_athena \
    --output-directory results/magnetohydrodynamic/shock_tube/1d/+x/2/0/roe_athena/ \
    --config-file config_files/magnetohydrodynamic/shock_tube/1d/+x/two_fluid0.cfg
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/shock_tube/1d/+x/2/0/roe_athena/mhd*dc 2> /dev/null
$RUN ./hd2gnuplot.exe results/magnetohydrodynamic/shock_tube/1d/+x/2/0/roe_athena/hd*dc 2> /dev/null


echo shock_tube mhd 1 1d +x hll_athena
$RUN ./two_test.exe \
    --solver-mhd hll_athena \
    --output-directory results/magnetohydrodynamic/shock_tube/1d/+x/2/1/hll_athena/ \
    --config-file config_files/magnetohydrodynamic/shock_tube/1d/+x/two_fluid1.cfg
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/shock_tube/1d/+x/2/1/hll_athena/mhd*dc 2> /dev/null
$RUN ./hd2gnuplot.exe results/magnetohydrodynamic/shock_tube/1d/+x/2/1/hll_athena/hd*dc 2> /dev/null

echo shock_tube mhd 1 1d +x hlld_athena
$RUN ./two_test.exe \
    --solver-mhd hlld_athena \
    --output-directory results/magnetohydrodynamic/shock_tube/1d/+x/2/1/hlld_athena/ \
    --config-file config_files/magnetohydrodynamic/shock_tube/1d/+x/two_fluid1.cfg
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/shock_tube/1d/+x/2/1/hlld_athena/mhd*dc 2> /dev/null
$RUN ./hd2gnuplot.exe results/magnetohydrodynamic/shock_tube/1d/+x/2/1/hlld_athena/hd*dc 2> /dev/null

echo shock_tube mhd 1 1d +x roe_athena
$RUN ./two_test.exe \
    --solver-mhd roe_athena \
    --output-directory results/magnetohydrodynamic/shock_tube/1d/+x/2/1/roe_athena/ \
    --config-file config_files/magnetohydrodynamic/shock_tube/1d/+x/two_fluid1.cfg
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/shock_tube/1d/+x/2/1/roe_athena/mhd*dc 2> /dev/null
$RUN ./hd2gnuplot.exe results/magnetohydrodynamic/shock_tube/1d/+x/2/1/roe_athena/hd*dc 2> /dev/null


echo shock_tube mhd 2 1d +x hll_athena
$RUN ./two_test.exe \
    --solver-mhd hll_athena \
    --output-directory results/magnetohydrodynamic/shock_tube/1d/+x/2/2/hll_athena/ \
    --config-file config_files/magnetohydrodynamic/shock_tube/1d/+x/two_fluid2.cfg
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/shock_tube/1d/+x/2/2/hll_athena/mhd*dc 2> /dev/null
$RUN ./hd2gnuplot.exe results/magnetohydrodynamic/shock_tube/1d/+x/2/2/hll_athena/hd*dc 2> /dev/null

echo shock_tube mhd 2 1d +x hlld_athena
$RUN ./two_test.exe \
    --solver-mhd hlld_athena \
    --output-directory results/magnetohydrodynamic/shock_tube/1d/+x/2/2/hlld_athena/ \
    --config-file config_files/magnetohydrodynamic/shock_tube/1d/+x/two_fluid2.cfg
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/shock_tube/1d/+x/2/2/hlld_athena/mhd*dc 2> /dev/null
$RUN ./hd2gnuplot.exe results/magnetohydrodynamic/shock_tube/1d/+x/2/2/hlld_athena/hd*dc 2> /dev/null

echo shock_tube mhd 2 1d +x roe_athena
$RUN ./two_test.exe \
    --solver-mhd roe_athena \
    --output-directory results/magnetohydrodynamic/shock_tube/1d/+x/2/2/roe_athena/ \
    --config-file config_files/magnetohydrodynamic/shock_tube/1d/+x/two_fluid2.cfg
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/shock_tube/1d/+x/2/2/roe_athena/mhd*dc 2> /dev/null
$RUN ./hd2gnuplot.exe results/magnetohydrodynamic/shock_tube/1d/+x/2/2/roe_athena/hd*dc 2> /dev/null


echo shock_tube mhd 3 1d +x hll_athena
$RUN ./two_test.exe \
    --solver-mhd hll_athena \
    --output-directory results/magnetohydrodynamic/shock_tube/1d/+x/2/3/hll_athena/ \
    --config-file config_files/magnetohydrodynamic/shock_tube/1d/+x/two_fluid3.cfg
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/shock_tube/1d/+x/2/3/hll_athena/mhd*dc 2> /dev/null
$RUN ./hd2gnuplot.exe results/magnetohydrodynamic/shock_tube/1d/+x/2/3/hll_athena/hd*dc 2> /dev/null

echo shock_tube mhd 3 1d +x hlld_athena
$RUN ./two_test.exe \
    --solver-mhd hlld_athena \
    --output-directory results/magnetohydrodynamic/shock_tube/1d/+x/2/3/hlld_athena/ \
    --config-file config_files/magnetohydrodynamic/shock_tube/1d/+x/two_fluid3.cfg
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/shock_tube/1d/+x/2/3/hlld_athena/mhd*dc 2> /dev/null
$RUN ./hd2gnuplot.exe results/magnetohydrodynamic/shock_tube/1d/+x/2/3/hlld_athena/hd*dc 2> /dev/null

echo shock_tube mhd 3 1d +x roe_athena
$RUN ./two_test.exe \
    --solver-mhd roe_athena \
    --output-directory results/magnetohydrodynamic/shock_tube/1d/+x/2/3/roe_athena/ \
    --config-file config_files/magnetohydrodynamic/shock_tube/1d/+x/two_fluid3.cfg
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/shock_tube/1d/+x/2/3/roe_athena/mhd*dc 2> /dev/null
$RUN ./hd2gnuplot.exe results/magnetohydrodynamic/shock_tube/1d/+x/2/3/roe_athena/hd*dc 2> /dev/null


echo shock_tube mhd 4 1d +x hll_athena
$RUN ./two_test.exe \
    --solver-mhd hll_athena \
    --output-directory results/magnetohydrodynamic/shock_tube/1d/+x/2/4/hll_athena/ \
    --config-file config_files/magnetohydrodynamic/shock_tube/1d/+x/two_fluid4.cfg
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/shock_tube/1d/+x/2/4/hll_athena/mhd*dc 2> /dev/null
$RUN ./hd2gnuplot.exe results/magnetohydrodynamic/shock_tube/1d/+x/2/4/hll_athena/hd*dc 2> /dev/null

echo shock_tube mhd 4 1d +x hlld_athena
$RUN ./two_test.exe \
    --solver-mhd hlld_athena \
    --output-directory results/magnetohydrodynamic/shock_tube/1d/+x/2/4/hlld_athena/ \
    --config-file config_files/magnetohydrodynamic/shock_tube/1d/+x/two_fluid4.cfg
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/shock_tube/1d/+x/2/4/hlld_athena/mhd*dc 2> /dev/null
$RUN ./hd2gnuplot.exe results/magnetohydrodynamic/shock_tube/1d/+x/2/4/hlld_athena/hd*dc 2> /dev/null

echo shock_tube mhd 4 1d +x roe_athena
$RUN ./two_test.exe \
    --solver-mhd roe_athena \
    --output-directory results/magnetohydrodynamic/shock_tube/1d/+x/2/4/roe_athena/ \
    --config-file config_files/magnetohydrodynamic/shock_tube/1d/+x/two_fluid4.cfg
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/shock_tube/1d/+x/2/4/roe_athena/mhd*dc 2> /dev/null
$RUN ./hd2gnuplot.exe results/magnetohydrodynamic/shock_tube/1d/+x/2/4/roe_athena/hd*dc 2> /dev/null


# Orszag-Tang
echo orszag-tang mhd 2d z1 hll_athena
$RUN ./two_test.exe \
    --solver-mhd hll_athena \
    --output-directory results/magnetohydrodynamic/orszag-tang/2d/z1/2/hll_athena/ \
    --config-file config_files/magnetohydrodynamic/orszag-tang/2d/z1/two_fluid1.cfg
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/orszag-tang/2d/z1/2/hll_athena/mhd*dc 2> /dev/null
$RUN ./hd2gnuplot.exe results/magnetohydrodynamic/orszag-tang/2d/z1/2/hll_athena/hd*dc 2> /dev/null

echo orszag-tang mhd 2d z1 roe_athena
$RUN ./two_test.exe \
    --solver-mhd roe_athena \
    --output-directory results/magnetohydrodynamic/orszag-tang/2d/z1/2/roe_athena/ \
    --config-file config_files/magnetohydrodynamic/orszag-tang/2d/z1/two_fluid1.cfg
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/orszag-tang/2d/z1/2/roe_athena/mhd*dc 2> /dev/null
$RUN ./hd2gnuplot.exe results/magnetohydrodynamic/orszag-tang/2d/z1/2/roe_athena/hd*dc 2> /dev/null
