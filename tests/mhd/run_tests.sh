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
$RUN ./test.exe \
    --output-directory results/magnetohydrodynamic/magnetic_field_advection/1d/+x/ \
    --config-file config_files/magnetohydrodynamic/magnetic_field_advection/1d/+x/config.cfg
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/magnetic_field_advection/1d/+x/*dc 2> /dev/null

echo magnetic_field_advection 1d -x roe_athena
$RUN ./test.exe \
    --output-directory results/magnetohydrodynamic/magnetic_field_advection/1d/-x/ \
    --config-file config_files/magnetohydrodynamic/magnetic_field_advection/1d/-x/config.cfg
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/magnetic_field_advection/1d/-x/*dc 2> /dev/null

echo magnetic_field_advection 1d +y roe_athena
$RUN ./test.exe \
    --output-directory results/magnetohydrodynamic/magnetic_field_advection/1d/+y/ \
    --config-file config_files/magnetohydrodynamic/magnetic_field_advection/1d/+y/config.cfg
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/magnetic_field_advection/1d/+y/*dc 2> /dev/null

echo magnetic_field_advection 1d -y roe_athena
$RUN ./test.exe \
    --output-directory results/magnetohydrodynamic/magnetic_field_advection/1d/-y/ \
    --config-file config_files/magnetohydrodynamic/magnetic_field_advection/1d/-y/config.cfg
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/magnetic_field_advection/1d/-y/*dc 2> /dev/null

echo magnetic_field_advection 1d +z roe_athena
$RUN ./test.exe \
    --output-directory results/magnetohydrodynamic/magnetic_field_advection/1d/+z/ \
    --config-file config_files/magnetohydrodynamic/magnetic_field_advection/1d/+z/config.cfg
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/magnetic_field_advection/1d/+z/*dc 2> /dev/null

echo magnetic_field_advection 1d -z roe_athena
$RUN ./test.exe \
    --output-directory results/magnetohydrodynamic/magnetic_field_advection/1d/-z/ \
    --config-file config_files/magnetohydrodynamic/magnetic_field_advection/1d/-z/config.cfg
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/magnetic_field_advection/1d/-z/*dc 2> /dev/null

# shock tube
echo shock_tube mhd 1d +x hll_athena
$RUN ./test.exe \
    --solver-mhd hll_athena \
    --output-directory results/magnetohydrodynamic/shock_tube/1d/+x/hll_athena/ \
    --config-file config_files/magnetohydrodynamic/shock_tube/1d/+x/config.cfg
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/shock_tube/1d/+x/hll_athena/*dc 2> /dev/null

echo shock_tube mhd 1d +x hlld_athena
$RUN ./test.exe \
    --solver-mhd hlld_athena \
    --output-directory results/magnetohydrodynamic/shock_tube/1d/+x/hlld_athena/ \
    --config-file config_files/magnetohydrodynamic/shock_tube/1d/+x/config.cfg
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/shock_tube/1d/+x/hlld_athena/*dc 2> /dev/null

echo shock_tube mhd 1d +x roe_athena
$RUN ./test.exe \
    --solver-mhd roe_athena \
    --output-directory results/magnetohydrodynamic/shock_tube/1d/+x/roe_athena/ \
    --config-file config_files/magnetohydrodynamic/shock_tube/1d/+x/config.cfg
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/shock_tube/1d/+x/roe_athena/*dc 2> /dev/null


echo shock_tube mhd 1d -x hll_athena
$RUN ./test.exe \
    --solver-mhd hll_athena \
    --output-directory results/magnetohydrodynamic/shock_tube/1d/-x/hll_athena/ \
    --config-file config_files/magnetohydrodynamic/shock_tube/1d/-x/config.cfg
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/shock_tube/1d/-x/hll_athena/*dc 2> /dev/null

echo shock_tube mhd 1d -x hlld_athena
$RUN ./test.exe \
    --solver-mhd hlld_athena \
    --output-directory results/magnetohydrodynamic/shock_tube/1d/-x/hlld_athena/ \
    --config-file config_files/magnetohydrodynamic/shock_tube/1d/-x/config.cfg
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/shock_tube/1d/-x/hlld_athena/*dc 2> /dev/null

echo shock_tube mhd 1d -x roe_athena
$RUN ./test.exe \
    --solver-mhd roe_athena \
    --output-directory results/magnetohydrodynamic/shock_tube/1d/-x/roe_athena/ \
    --config-file config_files/magnetohydrodynamic/shock_tube/1d/-x/config.cfg
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/shock_tube/1d/-x/roe_athena/*dc 2> /dev/null


echo shock_tube mhd 1d +y hll_athena
$RUN ./test.exe \
    --solver-mhd hll_athena \
    --output-directory results/magnetohydrodynamic/shock_tube/1d/+y/hll_athena/ \
    --config-file config_files/magnetohydrodynamic/shock_tube/1d/+y/config.cfg
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/shock_tube/1d/+y/hll_athena/*dc 2> /dev/null

echo shock_tube mhd 1d +y hlld_athena
$RUN ./test.exe \
    --solver-mhd hlld_athena \
    --output-directory results/magnetohydrodynamic/shock_tube/1d/+y/hlld_athena/ \
    --config-file config_files/magnetohydrodynamic/shock_tube/1d/+y/config.cfg
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/shock_tube/1d/+y/hlld_athena/*dc 2> /dev/null

echo shock_tube mhd 1d +y roe_athena
$RUN ./test.exe \
    --solver-mhd roe_athena \
    --output-directory results/magnetohydrodynamic/shock_tube/1d/+y/roe_athena/ \
    --config-file config_files/magnetohydrodynamic/shock_tube/1d/+y/config.cfg
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/shock_tube/1d/+y/roe_athena/*dc 2> /dev/null


echo shock_tube mhd 1d -y hll_athena
$RUN ./test.exe \
    --solver-mhd hll_athena \
    --output-directory results/magnetohydrodynamic/shock_tube/1d/-y/hll_athena/ \
    --config-file config_files/magnetohydrodynamic/shock_tube/1d/-y/config.cfg
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/shock_tube/1d/-y/hll_athena/*dc 2> /dev/null

echo shock_tube mhd 1d -y hlld_athena
$RUN ./test.exe \
    --solver-mhd hlld_athena \
    --output-directory results/magnetohydrodynamic/shock_tube/1d/-y/hlld_athena/ \
    --config-file config_files/magnetohydrodynamic/shock_tube/1d/-y/config.cfg
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/shock_tube/1d/-y/hlld_athena/*dc 2> /dev/null

echo shock_tube mhd 1d -y roe_athena
$RUN ./test.exe \
    --solver-mhd roe_athena \
    --output-directory results/magnetohydrodynamic/shock_tube/1d/-y/roe_athena/ \
    --config-file config_files/magnetohydrodynamic/shock_tube/1d/-y/config.cfg
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/shock_tube/1d/-y/roe_athena/*dc 2> /dev/null


echo shock_tube mhd 1d +z hll_athena
$RUN ./test.exe \
    --solver-mhd hll_athena \
    --output-directory results/magnetohydrodynamic/shock_tube/1d/+z/hll_athena/ \
    --config-file config_files/magnetohydrodynamic/shock_tube/1d/+z/config.cfg
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/shock_tube/1d/+z/hll_athena/*dc 2> /dev/null

echo shock_tube mhd 1d +z hlld_athena
$RUN ./test.exe \
    --solver-mhd hlld_athena \
    --output-directory results/magnetohydrodynamic/shock_tube/1d/+z/hlld_athena/ \
    --config-file config_files/magnetohydrodynamic/shock_tube/1d/+z/config.cfg
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/shock_tube/1d/+z/hlld_athena/*dc 2> /dev/null

echo shock_tube mhd 1d +z roe_athena
$RUN ./test.exe \
    --solver-mhd roe_athena \
    --output-directory results/magnetohydrodynamic/shock_tube/1d/+z/roe_athena/ \
    --config-file config_files/magnetohydrodynamic/shock_tube/1d/+z/config.cfg
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/shock_tube/1d/+z/roe_athena/*dc 2> /dev/null


echo shock_tube mhd 1d -z hll_athena
$RUN ./test.exe \
    --solver-mhd hll_athena \
    --output-directory results/magnetohydrodynamic/shock_tube/1d/-z/hll_athena/ \
    --config-file config_files/magnetohydrodynamic/shock_tube/1d/-z/config.cfg
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/shock_tube/1d/-z/hll_athena/*dc 2> /dev/null

echo shock_tube mhd 1d -z hlld_athena
$RUN ./test.exe \
    --solver-mhd hlld_athena \
    --output-directory results/magnetohydrodynamic/shock_tube/1d/-z/hlld_athena/ \
    --config-file config_files/magnetohydrodynamic/shock_tube/1d/-z/config.cfg
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/shock_tube/1d/-z/hlld_athena/*dc 2> /dev/null

echo shock_tube mhd 1d -z roe_athena
$RUN ./test.exe \
    --solver-mhd roe_athena \
    --output-directory results/magnetohydrodynamic/shock_tube/1d/-z/roe_athena/ \
    --config-file config_files/magnetohydrodynamic/shock_tube/1d/-z/config.cfg
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/shock_tube/1d/-z/roe_athena/*dc 2> /dev/null



# Fast magnetosonic wave, small amplitude
echo fast_magnetosonic_wave mhd 1d +x small hll_athena
$RUN ./test.exe \
    --solver-mhd hll_athena \
    --output-directory results/magnetohydrodynamic/fast_magnetosonic_wave/1d/+x/small/hll_athena/ \
    --config-file config_files/magnetohydrodynamic/fast_magnetosonic_wave/1d/+x/small.cfg
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/fast_magnetosonic_wave/1d/+x/small/hll_athena/*dc 2> /dev/null

echo fast_magnetosonic_wave mhd 1d +x small hlld_athena
$RUN ./test.exe \
    --solver-mhd hlld_athena \
    --output-directory results/magnetohydrodynamic/fast_magnetosonic_wave/1d/+x/small/hlld_athena/ \
    --config-file config_files/magnetohydrodynamic/fast_magnetosonic_wave/1d/+x/small.cfg
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/fast_magnetosonic_wave/1d/+x/small/hlld_athena/*dc 2> /dev/null

echo fast_magnetosonic_wave mhd 1d +x small roe_athena
$RUN ./test.exe \
    --solver-mhd roe_athena \
    --output-directory results/magnetohydrodynamic/fast_magnetosonic_wave/1d/+x/small/roe_athena/ \
    --config-file config_files/magnetohydrodynamic/fast_magnetosonic_wave/1d/+x/small.cfg
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/fast_magnetosonic_wave/1d/+x/small/roe_athena/*dc 2> /dev/null


# Fast magnetosonic wave, medium amplitude
echo fast_magnetosonic_wave mhd 1d +x medium hll_athena
$RUN ./test.exe \
    --solver-mhd hll_athena \
    --output-directory results/magnetohydrodynamic/fast_magnetosonic_wave/1d/+x/medium/hll_athena/ \
    --config-file config_files/magnetohydrodynamic/fast_magnetosonic_wave/1d/+x/medium.cfg
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/fast_magnetosonic_wave/1d/+x/medium/hll_athena/*dc 2> /dev/null

echo fast_magnetosonic_wave mhd 1d +x medium hlld_athena
$RUN ./test.exe \
    --solver-mhd hlld_athena \
    --output-directory results/magnetohydrodynamic/fast_magnetosonic_wave/1d/+x/medium/hlld_athena/ \
    --config-file config_files/magnetohydrodynamic/fast_magnetosonic_wave/1d/+x/medium.cfg
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/fast_magnetosonic_wave/1d/+x/medium/hlld_athena/*dc 2> /dev/null

echo fast_magnetosonic_wave mhd 1d +x medium roe_athena
$RUN ./test.exe \
    --solver-mhd roe_athena \
    --output-directory results/magnetohydrodynamic/fast_magnetosonic_wave/1d/+x/medium/roe_athena/ \
    --config-file config_files/magnetohydrodynamic/fast_magnetosonic_wave/1d/+x/medium.cfg
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/fast_magnetosonic_wave/1d/+x/medium/roe_athena/*dc 2> /dev/null


# Fast magnetosonic wave, large amplitude
echo fast_magnetosonic_wave mhd 1d +x large hll_athena
$RUN ./test.exe \
    --solver-mhd hll_athena \
    --output-directory results/magnetohydrodynamic/fast_magnetosonic_wave/1d/+x/large/hll_athena/ \
    --config-file config_files/magnetohydrodynamic/fast_magnetosonic_wave/1d/+x/large.cfg
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/fast_magnetosonic_wave/1d/+x/large/hll_athena/*dc 2> /dev/null

echo fast_magnetosonic_wave mhd 1d +x large hlld_athena
$RUN ./test.exe \
    --solver-mhd hlld_athena \
    --output-directory results/magnetohydrodynamic/fast_magnetosonic_wave/1d/+x/large/hlld_athena/ \
    --config-file config_files/magnetohydrodynamic/fast_magnetosonic_wave/1d/+x/large.cfg
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/fast_magnetosonic_wave/1d/+x/large/hlld_athena/*dc 2> /dev/null

echo fast_magnetosonic_wave mhd 1d +x large roe_athena
$RUN ./test.exe \
    --solver-mhd roe_athena \
    --output-directory results/magnetohydrodynamic/fast_magnetosonic_wave/1d/+x/large/roe_athena/ \
    --config-file config_files/magnetohydrodynamic/fast_magnetosonic_wave/1d/+x/large.cfg
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/fast_magnetosonic_wave/1d/+x/large/roe_athena/*dc 2> /dev/null


# Fast magnetosonic wave, huge amplitude
echo fast_magnetosonic_wave mhd 1d +x huge hll_athena
$RUN ./test.exe \
    --solver-mhd hll_athena \
    --output-directory results/magnetohydrodynamic/fast_magnetosonic_wave/1d/+x/huge/hll_athena/ \
    --config-file config_files/magnetohydrodynamic/fast_magnetosonic_wave/1d/+x/huge.cfg
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/fast_magnetosonic_wave/1d/+x/huge/hll_athena/*dc 2> /dev/null

echo fast_magnetosonic_wave mhd 1d +x huge hlld_athena
$RUN ./test.exe \
    --solver-mhd hlld_athena \
    --output-directory results/magnetohydrodynamic/fast_magnetosonic_wave/1d/+x/huge/hlld_athena/ \
    --config-file config_files/magnetohydrodynamic/fast_magnetosonic_wave/1d/+x/huge.cfg
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/fast_magnetosonic_wave/1d/+x/huge/hlld_athena/*dc 2> /dev/null

echo fast_magnetosonic_wave mhd 1d +x huge roe_athena
$RUN ./test.exe \
    --solver-mhd roe_athena \
    --output-directory results/magnetohydrodynamic/fast_magnetosonic_wave/1d/+x/huge/roe_athena/ \
    --config-file config_files/magnetohydrodynamic/fast_magnetosonic_wave/1d/+x/huge.cfg
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/fast_magnetosonic_wave/1d/+x/huge/roe_athena/*dc 2> /dev/null


# Kelvin-Helmholtz
echo kelvin-helmholtz hd 2d z1 hll_athena
$RUN ./test.exe \
    --time-length 1.75 \
    --save-mhd-n 0.25 \
    --solver-mhd hll_athena \
    --output-directory results/hydrodynamic/kelvin-helmholtz/2d/z1/hll_athena/ \
    --config-file config_files/hydrodynamic/kelvin-helmholtz/2d/z1/config.cfg
$RUN ./mhd2gnuplot.exe results/hydrodynamic/kelvin-helmholtz/2d/z1/hll_athena/*dc 2> /dev/null

echo kelvin-helmholtz hd 2d z1 hlld_athena
$RUN ./test.exe \
    --solver-mhd hlld_athena \
    --output-directory results/hydrodynamic/kelvin-helmholtz/2d/z1/hlld_athena/ \
    --config-file config_files/hydrodynamic/kelvin-helmholtz/2d/z1/config.cfg
$RUN ./mhd2gnuplot.exe results/hydrodynamic/kelvin-helmholtz/2d/z1/hlld_athena/*dc 2> /dev/null

echo kelvin-helmholtz hd 2d z1 roe_athena
$RUN ./test.exe \
    --solver-mhd roe_athena \
    --output-directory results/hydrodynamic/kelvin-helmholtz/2d/z1/roe_athena/ \
    --config-file config_files/hydrodynamic/kelvin-helmholtz/2d/z1/config.cfg
$RUN ./mhd2gnuplot.exe results/hydrodynamic/kelvin-helmholtz/2d/z1/roe_athena/*dc 2> /dev/null


echo kelvin-helmholtz mhd 2d z1 hll_athena
$RUN ./test.exe \
    --time-length 5 \
    --save-mhd-n 0.5 \
    --solver-mhd hll_athena \
    --output-directory results/magnetohydrodynamic/kelvin-helmholtz/2d/z1/hll_athena/ \
    --config-file config_files/magnetohydrodynamic/kelvin-helmholtz/2d/z1/config.cfg
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/kelvin-helmholtz/2d/z1/hll_athena/*dc 2> /dev/null

echo kelvin-helmholtz mhd 2d z1 hlld_athena
$RUN ./test.exe \
    --time-length 2 \
    --save-mhd-n 0.25 \
    --solver-mhd hlld_athena \
    --output-directory results/magnetohydrodynamic/kelvin-helmholtz/2d/z1/hlld_athena/ \
    --config-file config_files/magnetohydrodynamic/kelvin-helmholtz/2d/z1/config.cfg
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/kelvin-helmholtz/2d/z1/hlld_athena/*dc 2> /dev/null

echo kelvin-helmholtz mhd 2d z1 roe_athena
$RUN ./test.exe \
    --solver-mhd roe_athena \
    --output-directory results/magnetohydrodynamic/kelvin-helmholtz/2d/z1/roe_athena/ \
    --config-file config_files/magnetohydrodynamic/kelvin-helmholtz/2d/z1/config.cfg
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/kelvin-helmholtz/2d/z1/roe_athena/*dc 2> /dev/null



# Bow shock
echo bow_shock hd 2d z1 hll_athena
$RUN ./test.exe \
    --solver-mhd hll_athena \
    --output-directory results/hydrodynamic/bow_shock/2d/z1/hll_athena/ \
    --config-file config_files/hydrodynamic/bow_shock/2d/z1/config.cfg
$RUN ./mhd2gnuplot.exe results/hydrodynamic/bow_shock/2d/z1/hll_athena/*dc 2> /dev/null

echo bow_shock hd 2d z1 hlld_athena
$RUN ./test.exe \
    --solver-mhd hlld_athena \
    --output-directory results/hydrodynamic/bow_shock/2d/z1/hlld_athena/ \
    --config-file config_files/hydrodynamic/bow_shock/2d/z1/config.cfg
$RUN ./mhd2gnuplot.exe results/hydrodynamic/bow_shock/2d/z1/hlld_athena/*dc 2> /dev/null


echo bow_shock mhd 2d z1 hll_athena
$RUN ./test.exe \
    --solver-mhd hll_athena \
    --output-directory results/magnetohydrodynamic/bow_shock/2d/z1/hll_athena/ \
    --config-file config_files/magnetohydrodynamic/bow_shock/2d/z1/config.cfg
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/bow_shock/2d/z1/hll_athena/*dc 2> /dev/null



# Blast wave
echo blast_wave hd 2d z1 hll_athena
$RUN ./test.exe \
    --solver-mhd hll_athena \
    --output-directory results/hydrodynamic/blast_wave/2d/z1/hll_athena/ \
    --config-file config_files/hydrodynamic/blast_wave/2d/z1/config.cfg
$RUN ./mhd2gnuplot.exe results/hydrodynamic/blast_wave/2d/z1/hll_athena/*dc 2> /dev/null

echo blast_wave hd 2d z1 hlld_athena
$RUN ./test.exe \
    --solver-mhd hlld_athena \
    --output-directory results/hydrodynamic/blast_wave/2d/z1/hlld_athena/ \
    --config-file config_files/hydrodynamic/blast_wave/2d/z1/config.cfg
$RUN ./mhd2gnuplot.exe results/hydrodynamic/blast_wave/2d/z1/hlld_athena/*dc 2> /dev/null

echo blast_wave hd 2d z1 roe_athena
$RUN ./test.exe \
    --solver-mhd roe_athena \
    --output-directory results/hydrodynamic/blast_wave/2d/z1/roe_athena/ \
    --config-file config_files/hydrodynamic/blast_wave/2d/z1/config.cfg
$RUN ./mhd2gnuplot.exe results/hydrodynamic/blast_wave/2d/z1/roe_athena/*dc 2> /dev/null


echo blast_wave mhd 2d z1 hll_athena
$RUN ./test.exe \
    --solver-mhd hll_athena \
    --output-directory results/magnetohydrodynamic/blast_wave/2d/z1/hll_athena/ \
    --config-file config_files/magnetohydrodynamic/blast_wave/2d/z1/config.cfg
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/blast_wave/2d/z1/hll_athena/*dc 2> /dev/null

echo blast_wave mhd 2d z1 hlld_athena
$RUN ./test.exe \
    --solver-mhd hlld_athena \
    --output-directory results/magnetohydrodynamic/blast_wave/2d/z1/hlld_athena/ \
    --config-file config_files/magnetohydrodynamic/blast_wave/2d/z1/config.cfg
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/blast_wave/2d/z1/hlld_athena/*dc 2> /dev/null

echo blast_wave mhd 2d z1 roe_athena
$RUN ./test.exe \
    --solver-mhd roe_athena \
    --output-directory results/magnetohydrodynamic/blast_wave/2d/z1/roe_athena/ \
    --config-file config_files/magnetohydrodynamic/blast_wave/2d/z1/config.cfg
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/blast_wave/2d/z1/roe_athena/*dc 2> /dev/null



# Reconnection
echo reconnection mhd 2d y1 hll_athena
$RUN ./test.exe --config-file config_files/magnetohydrodynamic/reconnection/2d/y1/hll_athena.cfg
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/reconnection/2d/y1/hll_athena/*dc 2> /dev/null

echo reconnection mhd 2d y1 hlld_athena
$RUN ./test.exe --config-file config_files/magnetohydrodynamic/reconnection/2d/y1/hlld_athena.cfg
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/reconnection/2d/y1/hlld_athena/*dc 2> /dev/null

echo reconnection mhd 2d y1 roe_athena
$RUN ./test.exe --config-file config_files/magnetohydrodynamic/reconnection/2d/y1/roe_athena.cfg
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/reconnection/2d/y1/roe_athena/*dc 2> /dev/null



# Orszag-Tang
echo orszag-tang mhd 2d z1 hll_athena
$RUN ./test.exe \
    --solver-mhd hll_athena \
    --output-directory results/magnetohydrodynamic/orszag-tang/2d/z1/hll_athena/ \
    --config-file config_files/magnetohydrodynamic/orszag-tang/2d/z1/config.cfg
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/orszag-tang/2d/z1/hll_athena/*dc 2> /dev/null

echo orszag-tang mhd 2d z1 roe_athena
$RUN ./test.exe \
    --solver-mhd roe_athena \
    --output-directory results/magnetohydrodynamic/orszag-tang/2d/z1/roe_athena/ \
    --config-file config_files/magnetohydrodynamic/orszag-tang/2d/z1/config.cfg
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/orszag-tang/2d/z1/roe_athena/*dc 2> /dev/null
