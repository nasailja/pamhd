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

# magnetic field advection with resistivity
echo magnetic_field_advection 1d +x roe_athena
$RUN ./test.exe config_files/magnetohydrodynamic/magnetic_field_advection/1d/+x/config.json
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/magnetic_field_advection/1d/+x/roe_athena/*dc 2> /dev/null

echo magnetic_field_advection 1d -x roe_athena
$RUN ./test.exe config_files/magnetohydrodynamic/magnetic_field_advection/1d/-x/config.json
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/magnetic_field_advection/1d/-x/roe_athena/*dc 2> /dev/null

echo magnetic_field_advection 1d +y roe_athena
$RUN ./test.exe config_files/magnetohydrodynamic/magnetic_field_advection/1d/+y/config.json
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/magnetic_field_advection/1d/+y/roe_athena/*dc 2> /dev/null

echo magnetic_field_advection 1d -y roe_athena
$RUN ./test.exe config_files/magnetohydrodynamic/magnetic_field_advection/1d/-y/config.json
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/magnetic_field_advection/1d/-y/roe_athena/*dc 2> /dev/null

echo magnetic_field_advection 1d +z roe_athena
$RUN ./test.exe config_files/magnetohydrodynamic/magnetic_field_advection/1d/+z/config.json
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/magnetic_field_advection/1d/+z/roe_athena/*dc 2> /dev/null

echo magnetic_field_advection 1d -z roe_athena
$RUN ./test.exe config_files/magnetohydrodynamic/magnetic_field_advection/1d/-z/config.json
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/magnetic_field_advection/1d/-z/roe_athena/*dc 2> /dev/null

# rarefaction
echo rarefaction hd 1d +x hll_athena
$RUN ./test.exe config_files/hydrodynamic/rarefaction/1d/+x/hll_athena.json
$RUN ./mhd2gnuplot.exe results/hydrodynamic/rarefaction/1d/+x/hll_athena/*dc 2> /dev/null

echo rarefaction hd 1d +x hlld_athena
$RUN ./test.exe config_files/hydrodynamic/rarefaction/1d/+x/hlld_athena.json
$RUN ./mhd2gnuplot.exe results/hydrodynamic/rarefaction/1d/+x/hlld_athena/*dc 2> /dev/null

# TODO: fix
#echo rarefaction hd 1d +x roe_athena
#$RUN ./test.exe config_files/hydrodynamic/rarefaction/1d/+x/roe_athena.json
#$RUN ./mhd2gnuplot.exe results/hydrodynamic/rarefaction/1d/+x/roe_athena/*dc 2> /dev/null

# shock tube
echo shock_tube mhd 1d +x hll_athena
$RUN ./test.exe config_files/magnetohydrodynamic/shock_tube/1d/+x/hll_athena.json
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/shock_tube/1d/+x/hll_athena/*dc 2> /dev/null

echo shock_tube mhd 1d +x hlld_athena
$RUN ./test.exe config_files/magnetohydrodynamic/shock_tube/1d/+x/hlld_athena.json
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/shock_tube/1d/+x/hlld_athena/*dc 2> /dev/null

echo shock_tube mhd 1d +x roe_athena
$RUN ./test.exe config_files/magnetohydrodynamic/shock_tube/1d/+x/roe_athena.json
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/shock_tube/1d/+x/roe_athena/*dc 2> /dev/null


echo shock_tube mhd 1d -x hll_athena
$RUN ./test.exe config_files/magnetohydrodynamic/shock_tube/1d/-x/hll_athena.json
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/shock_tube/1d/-x/hll_athena/*dc 2> /dev/null

echo shock_tube mhd 1d -x hlld_athena
$RUN ./test.exe config_files/magnetohydrodynamic/shock_tube/1d/-x/hlld_athena.json
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/shock_tube/1d/-x/hlld_athena/*dc 2> /dev/null

echo shock_tube mhd 1d -x roe_athena
$RUN ./test.exe config_files/magnetohydrodynamic/shock_tube/1d/-x/roe_athena.json
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/shock_tube/1d/-x/roe_athena/*dc 2> /dev/null


echo shock_tube mhd 1d +y hll_athena
$RUN ./test.exe config_files/magnetohydrodynamic/shock_tube/1d/+y/hll_athena.json
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/shock_tube/1d/+y/hll_athena/*dc 2> /dev/null

echo shock_tube mhd 1d +y hlld_athena
$RUN ./test.exe config_files/magnetohydrodynamic/shock_tube/1d/+y/hlld_athena.json
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/shock_tube/1d/+y/hlld_athena/*dc 2> /dev/null

echo shock_tube mhd 1d +y roe_athena
$RUN ./test.exe config_files/magnetohydrodynamic/shock_tube/1d/+y/roe_athena.json
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/shock_tube/1d/+y/roe_athena/*dc 2> /dev/null


echo shock_tube mhd 1d -y hll_athena
$RUN ./test.exe config_files/magnetohydrodynamic/shock_tube/1d/-y/hll_athena.json
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/shock_tube/1d/-y/hll_athena/*dc 2> /dev/null

echo shock_tube mhd 1d -y hlld_athena
$RUN ./test.exe config_files/magnetohydrodynamic/shock_tube/1d/-y/hlld_athena.json
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/shock_tube/1d/-y/hlld_athena/*dc 2> /dev/null

echo shock_tube mhd 1d -y roe_athena
$RUN ./test.exe config_files/magnetohydrodynamic/shock_tube/1d/-y/roe_athena.json
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/shock_tube/1d/-y/roe_athena/*dc 2> /dev/null


echo shock_tube mhd 1d +z hll_athena
$RUN ./test.exe config_files/magnetohydrodynamic/shock_tube/1d/+z/hll_athena.json
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/shock_tube/1d/+z/hll_athena/*dc 2> /dev/null

echo shock_tube mhd 1d +z hlld_athena
$RUN ./test.exe config_files/magnetohydrodynamic/shock_tube/1d/+z/hlld_athena.json
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/shock_tube/1d/+z/hlld_athena/*dc 2> /dev/null

echo shock_tube mhd 1d +z roe_athena
$RUN ./test.exe config_files/magnetohydrodynamic/shock_tube/1d/+z/roe_athena.json
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/shock_tube/1d/+z/roe_athena/*dc 2> /dev/null


echo shock_tube mhd 1d -z hll_athena
$RUN ./test.exe config_files/magnetohydrodynamic/shock_tube/1d/-z/hll_athena.json
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/shock_tube/1d/-z/hll_athena/*dc 2> /dev/null

echo shock_tube mhd 1d -z hlld_athena
$RUN ./test.exe config_files/magnetohydrodynamic/shock_tube/1d/-z/hlld_athena.json
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/shock_tube/1d/-z/hlld_athena/*dc 2> /dev/null

echo shock_tube mhd 1d -z roe_athena
$RUN ./test.exe config_files/magnetohydrodynamic/shock_tube/1d/-z/roe_athena.json
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/shock_tube/1d/-z/roe_athena/*dc 2> /dev/null



# Fast magnetosonic wave, small amplitude
echo fast_magnetosonic_wave mhd 1d +x small roe_athena
$RUN ./test.exe config_files/magnetohydrodynamic/fast_magnetosonic_wave/1d/+x/small.json
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/fast_magnetosonic_wave/1d/+x/small/roe_athena/*dc 2> /dev/null


# Fast magnetosonic wave, medium amplitude
echo fast_magnetosonic_wave mhd 1d +x medium roe_athena
$RUN ./test.exe config_files/magnetohydrodynamic/fast_magnetosonic_wave/1d/+x/medium.json
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/fast_magnetosonic_wave/1d/+x/medium/roe_athena/*dc 2> /dev/null

# Fast magnetosonic wave, large amplitude
echo fast_magnetosonic_wave mhd 1d +x large roe_athena
$RUN ./test.exe config_files/magnetohydrodynamic/fast_magnetosonic_wave/1d/+x/large.json
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/fast_magnetosonic_wave/1d/+x/large/roe_athena/*dc 2> /dev/null

# Fast magnetosonic wave, huge amplitude
echo fast_magnetosonic_wave mhd 1d +x huge roe_athena
$RUN ./test.exe config_files/magnetohydrodynamic/fast_magnetosonic_wave/1d/+x/huge.json
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/fast_magnetosonic_wave/1d/+x/huge/roe_athena/*dc 2> /dev/null


# Kelvin-Helmholtz
echo kelvin-helmholtz hd 2d z1 hll_athena
$RUN ./test.exe config_files/hydrodynamic/kelvin-helmholtz/2d/z1/hll_athena.json
$RUN ./mhd2gnuplot.exe results/hydrodynamic/kelvin-helmholtz/2d/z1/hll_athena/*dc 2> /dev/null

echo kelvin-helmholtz hd 2d z1 hlld_athena
$RUN ./test.exe config_files/hydrodynamic/kelvin-helmholtz/2d/z1/hlld_athena.json
$RUN ./mhd2gnuplot.exe results/hydrodynamic/kelvin-helmholtz/2d/z1/hlld_athena/*dc 2> /dev/null

echo kelvin-helmholtz hd 2d z1 roe_athena
$RUN ./test.exe config_files/hydrodynamic/kelvin-helmholtz/2d/z1/roe_athena.json
$RUN ./mhd2gnuplot.exe results/hydrodynamic/kelvin-helmholtz/2d/z1/roe_athena/*dc 2> /dev/null


echo kelvin-helmholtz mhd 2d z1 hll_athena
$RUN ./test.exe config_files/magnetohydrodynamic/kelvin-helmholtz/2d/z1/hll_athena.json
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/kelvin-helmholtz/2d/z1/hll_athena/*dc 2> /dev/null

echo kelvin-helmholtz mhd 2d z1 hlld_athena
$RUN ./test.exe config_files/magnetohydrodynamic/kelvin-helmholtz/2d/z1/hlld_athena.json
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/kelvin-helmholtz/2d/z1/hlld_athena/*dc 2> /dev/null

echo kelvin-helmholtz mhd 2d z1 roe_athena
$RUN ./test.exe config_files/magnetohydrodynamic/kelvin-helmholtz/2d/z1/roe_athena.json
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/kelvin-helmholtz/2d/z1/roe_athena/*dc 2> /dev/null



# Bow shock
echo bow_shock hd 2d z1 hll_athena
$RUN ./test.exe config_files/hydrodynamic/bow_shock/2d/z1/hll_athena.json
$RUN ./mhd2gnuplot.exe results/hydrodynamic/bow_shock/2d/z1/hll_athena/*dc 2> /dev/null

echo bow_shock hd 2d z1 hlld_athena
$RUN ./test.exe config_files/hydrodynamic/bow_shock/2d/z1/hlld_athena.json
$RUN ./mhd2gnuplot.exe results/hydrodynamic/bow_shock/2d/z1/hlld_athena/*dc 2> /dev/null

#TODO: fix echo bow_shock hd 2d z1 roe_athena
#$RUN ./test.exe config_files/hydrodynamic/bow_shock/2d/z1/roe_athena.json
#$RUN ./mhd2gnuplot.exe results/hydrodynamic/bow_shock/2d/z1/roe_athena/*dc 2> /dev/null


echo bow_shock mhd 2d z1 hll_athena
$RUN ./test.exe config_files/magnetohydrodynamic/bow_shock/2d/z1/hll_athena.json
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/bow_shock/2d/z1/hll_athena/*dc 2> /dev/null

echo bow_shock mhd 2d z1 hlld_athena
$RUN ./test.exe config_files/magnetohydrodynamic/bow_shock/2d/z1/hlld_athena.json
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/bow_shock/2d/z1/hlld_athena/*dc 2> /dev/null

echo bow_shock mhd 2d z1 roe_athena
$RUN ./test.exe config_files/magnetohydrodynamic/bow_shock/2d/z1/roe_athena.json
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/bow_shock/2d/z1/roe_athena/*dc 2> /dev/null



# Blast wave
echo blast_wave hd 2d z1 hll_athena
$RUN ./test.exe config_files/hydrodynamic/blast_wave/2d/z1/hll_athena.json
$RUN ./mhd2gnuplot.exe results/hydrodynamic/blast_wave/2d/z1/hll_athena/*dc 2> /dev/null

echo blast_wave hd 2d z1 hlld_athena
$RUN ./test.exe config_files/hydrodynamic/blast_wave/2d/z1/hlld_athena.json
$RUN ./mhd2gnuplot.exe results/hydrodynamic/blast_wave/2d/z1/hlld_athena/*dc 2> /dev/null

echo blast_wave hd 2d z1 roe_athena
$RUN ./test.exe config_files/hydrodynamic/blast_wave/2d/z1/roe_athena.json
$RUN ./mhd2gnuplot.exe results/hydrodynamic/blast_wave/2d/z1/roe_athena/*dc 2> /dev/null


echo blast_wave mhd 2d z1 hll_athena
$RUN ./test.exe config_files/magnetohydrodynamic/blast_wave/2d/z1/hll_athena.json
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/blast_wave/2d/z1/hll_athena/*dc 2> /dev/null

echo blast_wave mhd 2d z1 hlld_athena
$RUN ./test.exe config_files/magnetohydrodynamic/blast_wave/2d/z1/hlld_athena.json
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/blast_wave/2d/z1/hlld_athena/*dc 2> /dev/null

echo blast_wave mhd 2d z1 roe_athena
$RUN ./test.exe config_files/magnetohydrodynamic/blast_wave/2d/z1/roe_athena.json
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/blast_wave/2d/z1/roe_athena/*dc 2> /dev/null


# Reconnection
echo reconnection mhd 2d y1 hll_athena
$RUN ./test.exe config_files/magnetohydrodynamic/reconnection/2d/y1/hll_athena.json
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/reconnection/2d/y1/hll_athena/*dc 2> /dev/null

echo reconnection mhd 2d y1 hlld_athena
$RUN ./test.exe config_files/magnetohydrodynamic/reconnection/2d/y1/hlld_athena.json
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/reconnection/2d/y1/hlld_athena/*dc 2> /dev/null

echo reconnection mhd 2d y1 roe_athena
$RUN ./test.exe config_files/magnetohydrodynamic/reconnection/2d/y1/roe_athena.json
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/reconnection/2d/y1/roe_athena/*dc 2> /dev/null


# Orszag-Tang
echo orszag-tang mhd 2d z1 hll_athena
$RUN ./test.exe config_files/magnetohydrodynamic/orszag-tang/2d/z1/hll_athena.json
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/orszag-tang/2d/z1/hll_athena/*dc 2> /dev/null

echo orszag-tang mhd 2d z1 roe_athena
$RUN ./test.exe config_files/magnetohydrodynamic/orszag-tang/2d/z1/roe_athena.json
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/orszag-tang/2d/z1/roe_athena/*dc 2> /dev/null
