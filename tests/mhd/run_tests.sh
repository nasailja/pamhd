#! /bin/sh

# if given, override RUN with first command line argument
RUN=${1:-}

make

# shock tube
mkdir -p \
    results/magnetohydrodynamic/shock_tube/1d/+x/hll_athena \
    results/magnetohydrodynamic/shock_tube/1d/+x/hlld_athena \
    results/magnetohydrodynamic/shock_tube/1d/+x/roe_athena \
    results/magnetohydrodynamic/shock_tube/1d/-x/hll_athena \
    results/magnetohydrodynamic/shock_tube/1d/-x/hlld_athena \
    results/magnetohydrodynamic/shock_tube/1d/-x/roe_athena \
    results/magnetohydrodynamic/shock_tube/1d/+y/hll_athena \
    results/magnetohydrodynamic/shock_tube/1d/+y/hlld_athena \
    results/magnetohydrodynamic/shock_tube/1d/+y/roe_athena \
    results/magnetohydrodynamic/shock_tube/1d/-y/hll_athena \
    results/magnetohydrodynamic/shock_tube/1d/-y/hlld_athena \
    results/magnetohydrodynamic/shock_tube/1d/-y/roe_athena \
    results/magnetohydrodynamic/shock_tube/1d/+z/hll_athena \
    results/magnetohydrodynamic/shock_tube/1d/+z/hlld_athena \
    results/magnetohydrodynamic/shock_tube/1d/+z/roe_athena \
    results/magnetohydrodynamic/shock_tube/1d/-z/hll_athena \
    results/magnetohydrodynamic/shock_tube/1d/-z/hlld_athena \
    results/magnetohydrodynamic/shock_tube/1d/-z/roe_athena

# Kelvin-Helmholtz
mkdir -p \
    results/hydrodynamic/kelvin-helmholtz/2d/z1/hll_athena \
    results/hydrodynamic/kelvin-helmholtz/2d/z1/hlld_athena \
    results/hydrodynamic/kelvin-helmholtz/2d/z1/roe_athena \
    results/magnetohydrodynamic/kelvin-helmholtz/2d/z1/hll_athena \
    results/magnetohydrodynamic/kelvin-helmholtz/2d/z1/hlld_athena \
    results/magnetohydrodynamic/kelvin-helmholtz/2d/z1/roe_athena

# Bow shock
mkdir -p \
    results/hydrodynamic/bow_shock/2d/z1/hll_athena \
    results/hydrodynamic/bow_shock/2d/z1/hlld_athena \
    results/hydrodynamic/bow_shock/2d/z1/roe_athena \
    results/magnetohydrodynamic/bow_shock/2d/z1/hll_athena \
    results/magnetohydrodynamic/bow_shock/2d/z1/hlld_athena

# Blast wave
mkdir -p \
    results/hydrodynamic/blast_wave/2d/z1/hll_athena \
    results/hydrodynamic/blast_wave/2d/z1/hlld_athena \
    results/hydrodynamic/blast_wave/2d/z1/roe_athena \
    results/magnetohydrodynamic/blast_wave/2d/z1/hll_athena \
    results/magnetohydrodynamic/blast_wave/2d/z1/hlld_athena \
    results/magnetohydrodynamic/blast_wave/2d/z1/roe_athena

# Reconnection
mkdir -p \
    results/magnetohydrodynamic/reconnection/2d/y1/hll_athena \
    results/magnetohydrodynamic/reconnection/2d/y1/hlld_athena \
    results/magnetohydrodynamic/reconnection/2d/y1/roe_athena


echo
echo Running tests...
echo

echo shock_tube mhd 1d +x hll_athena
$RUN ./test.exe \
    --config-file config_files/magnetohydrodynamic/shock_tube/1d/+x/hll_athena.cfg \
    --boundary-file config_files/magnetohydrodynamic/shock_tube/1d/+x/boundaries.cfg
echo shock_tube mhd 1d +x hlld_athena
$RUN ./test.exe \
    --config-file config_files/magnetohydrodynamic/shock_tube/1d/+x/hlld_athena.cfg \
    --boundary-file config_files/magnetohydrodynamic/shock_tube/1d/+x/boundaries.cfg
echo shock_tube mhd 1d +x roe_athena
$RUN ./test.exe \
    --config-file config_files/magnetohydrodynamic/shock_tube/1d/+x/roe_athena.cfg \
    --boundary-file config_files/magnetohydrodynamic/shock_tube/1d/+x/boundaries.cfg

echo shock_tube mhd 1d -x hll_athena
$RUN ./test.exe \
    --config-file config_files/magnetohydrodynamic/shock_tube/1d/-x/hll_athena.cfg \
    --boundary-file config_files/magnetohydrodynamic/shock_tube/1d/-x/boundaries.cfg
echo shock_tube mhd 1d -x hlld_athena
$RUN ./test.exe \
    --config-file config_files/magnetohydrodynamic/shock_tube/1d/-x/hlld_athena.cfg \
    --boundary-file config_files/magnetohydrodynamic/shock_tube/1d/-x/boundaries.cfg
echo shock_tube mhd 1d -x roe_athena
$RUN ./test.exe \
    --config-file config_files/magnetohydrodynamic/shock_tube/1d/-x/roe_athena.cfg \
    --boundary-file config_files/magnetohydrodynamic/shock_tube/1d/-x/boundaries.cfg


echo shock_tube mhd 1d +y hll_athena
$RUN ./test.exe \
    --config-file config_files/magnetohydrodynamic/shock_tube/1d/+y/hll_athena.cfg \
    --boundary-file config_files/magnetohydrodynamic/shock_tube/1d/+y/boundaries.cfg
echo shock_tube mhd 1d +y hlld_athena
$RUN ./test.exe \
    --config-file config_files/magnetohydrodynamic/shock_tube/1d/+y/hlld_athena.cfg \
    --boundary-file config_files/magnetohydrodynamic/shock_tube/1d/+y/boundaries.cfg
echo shock_tube mhd 1d +y roe_athena
$RUN ./test.exe \
    --config-file config_files/magnetohydrodynamic/shock_tube/1d/+y/roe_athena.cfg \
    --boundary-file config_files/magnetohydrodynamic/shock_tube/1d/+y/boundaries.cfg

echo shock_tube mhd 1d -y hll_athena
$RUN ./test.exe \
    --config-file config_files/magnetohydrodynamic/shock_tube/1d/-y/hll_athena.cfg \
    --boundary-file config_files/magnetohydrodynamic/shock_tube/1d/-y/boundaries.cfg
echo shock_tube mhd 1d -y hlld_athena
$RUN ./test.exe \
    --config-file config_files/magnetohydrodynamic/shock_tube/1d/-y/hlld_athena.cfg \
    --boundary-file config_files/magnetohydrodynamic/shock_tube/1d/-y/boundaries.cfg
echo shock_tube mhd 1d -y roe_athena
$RUN ./test.exe \
    --config-file config_files/magnetohydrodynamic/shock_tube/1d/-y/roe_athena.cfg \
    --boundary-file config_files/magnetohydrodynamic/shock_tube/1d/-y/boundaries.cfg


echo shock_tube mhd 1d +z hll_athena
$RUN ./test.exe \
    --config-file config_files/magnetohydrodynamic/shock_tube/1d/+z/hll_athena.cfg \
    --boundary-file config_files/magnetohydrodynamic/shock_tube/1d/+z/boundaries.cfg
echo shock_tube mhd 1d +z hlld_athena
$RUN ./test.exe \
    --config-file config_files/magnetohydrodynamic/shock_tube/1d/+z/hlld_athena.cfg \
    --boundary-file config_files/magnetohydrodynamic/shock_tube/1d/+z/boundaries.cfg
echo shock_tube mhd 1d +z roe_athena
$RUN ./test.exe \
    --config-file config_files/magnetohydrodynamic/shock_tube/1d/+z/roe_athena.cfg \
    --boundary-file config_files/magnetohydrodynamic/shock_tube/1d/+z/boundaries.cfg

echo shock_tube mhd 1d -z hll_athena
$RUN ./test.exe \
    --config-file config_files/magnetohydrodynamic/shock_tube/1d/-z/hll_athena.cfg \
    --boundary-file config_files/magnetohydrodynamic/shock_tube/1d/-z/boundaries.cfg
echo shock_tube mhd 1d -z hlld_athena
$RUN ./test.exe \
    --config-file config_files/magnetohydrodynamic/shock_tube/1d/-z/hlld_athena.cfg \
    --boundary-file config_files/magnetohydrodynamic/shock_tube/1d/-z/boundaries.cfg
echo shock_tube mhd 1d -z roe_athena
$RUN ./test.exe \
    --config-file config_files/magnetohydrodynamic/shock_tube/1d/-z/roe_athena.cfg \
    --boundary-file config_files/magnetohydrodynamic/shock_tube/1d/-z/boundaries.cfg



echo kelvin-helmholtz hd 2d z1 hll_athena
$RUN ./test.exe \
    --config-file config_files/hydrodynamic/kelvin-helmholtz/2d/z1/hll_athena.cfg \
    --boundary-file config_files/hydrodynamic/kelvin-helmholtz/2d/z1/boundaries.cfg
echo kelvin-helmholtz hd 2d z1 hlld_athena
$RUN ./test.exe \
    --config-file config_files/hydrodynamic/kelvin-helmholtz/2d/z1/hlld_athena.cfg \
    --boundary-file config_files/hydrodynamic/kelvin-helmholtz/2d/z1/boundaries.cfg
echo kelvin-helmholtz hd 2d z1 roe_athena
$RUN ./test.exe \
    --config-file config_files/hydrodynamic/kelvin-helmholtz/2d/z1/roe_athena.cfg \
    --boundary-file config_files/hydrodynamic/kelvin-helmholtz/2d/z1/boundaries.cfg


echo kelvin-helmholtz mhd 2d z1 hll_athena
$RUN ./test.exe \
    --config-file config_files/magnetohydrodynamic/kelvin-helmholtz/2d/z1/hll_athena.cfg \
    --boundary-file config_files/magnetohydrodynamic/kelvin-helmholtz/2d/z1/boundaries.cfg
echo kelvin-helmholtz mhd 2d z1 hlld_athena
$RUN ./test.exe \
    --config-file config_files/magnetohydrodynamic/kelvin-helmholtz/2d/z1/hlld_athena.cfg \
    --boundary-file config_files/magnetohydrodynamic/kelvin-helmholtz/2d/z1/boundaries.cfg
echo kelvin-helmholtz mhd 2d z1 roe_athena
$RUN ./test.exe \
    --config-file config_files/magnetohydrodynamic/kelvin-helmholtz/2d/z1/roe_athena.cfg \
    --boundary-file config_files/magnetohydrodynamic/kelvin-helmholtz/2d/z1/boundaries.cfg



echo bow_shock hd 2d z1 hll_athena
$RUN ./test.exe \
    --config-file config_files/hydrodynamic/bow_shock/2d/z1/hll_athena.cfg \
    --boundary-file config_files/hydrodynamic/bow_shock/2d/z1/boundaries.cfg
echo bow_shock hd 2d z1 hlld_athena
$RUN ./test.exe \
    --config-file config_files/hydrodynamic/bow_shock/2d/z1/hlld_athena.cfg \
    --boundary-file config_files/hydrodynamic/bow_shock/2d/z1/boundaries.cfg


echo bow_shock mhd 2d z1 hll_athena
$RUN ./test.exe \
    --config-file config_files/magnetohydrodynamic/bow_shock/2d/z1/hll_athena.cfg \
    --boundary-file config_files/magnetohydrodynamic/bow_shock/2d/z1/boundaries.cfg
echo bow_shock mhd 2d z1 hlld_athena
$RUN ./test.exe \
    --config-file config_files/magnetohydrodynamic/bow_shock/2d/z1/hlld_athena.cfg \
    --boundary-file config_files/magnetohydrodynamic/bow_shock/2d/z1/boundaries.cfg



echo blast_wave hd 2d z1 hll_athena
$RUN ./test.exe \
    --config-file config_files/hydrodynamic/blast_wave/2d/z1/hll_athena.cfg \
    --boundary-file config_files/hydrodynamic/blast_wave/2d/z1/boundaries.cfg
echo blast_wave hd 2d z1 hlld_athena
$RUN ./test.exe \
    --config-file config_files/hydrodynamic/blast_wave/2d/z1/hlld_athena.cfg \
    --boundary-file config_files/hydrodynamic/blast_wave/2d/z1/boundaries.cfg
echo blast_wave hd 2d z1 roe_athena
$RUN ./test.exe \
    --config-file config_files/hydrodynamic/blast_wave/2d/z1/roe_athena.cfg \
    --boundary-file config_files/hydrodynamic/blast_wave/2d/z1/boundaries.cfg


echo blast_wave mhd 2d z1 hll_athena
$RUN ./test.exe \
    --config-file config_files/magnetohydrodynamic/blast_wave/2d/z1/hll_athena.cfg \
    --boundary-file config_files/magnetohydrodynamic/blast_wave/2d/z1/boundaries.cfg
echo blast_wave mhd 2d z1 hlld_athena
$RUN ./test.exe \
    --config-file config_files/magnetohydrodynamic/blast_wave/2d/z1/hlld_athena.cfg \
    --boundary-file config_files/magnetohydrodynamic/blast_wave/2d/z1/boundaries.cfg



echo reconnection mhd 2d y1 hll_athena
$RUN ./test.exe \
    --config-file config_files/magnetohydrodynamic/reconnection/2d/y1/hll_athena.cfg \
    --boundary-file config_files/magnetohydrodynamic/reconnection/2d/y1/boundaries.cfg
echo reconnection mhd 2d y1 hlld_athena
$RUN ./test.exe \
    --config-file config_files/magnetohydrodynamic/reconnection/2d/y1/hlld_athena.cfg \
    --boundary-file config_files/magnetohydrodynamic/reconnection/2d/y1/boundaries.cfg
echo reconnection mhd 2d y1 roe_athena
$RUN ./test.exe \
    --config-file config_files/magnetohydrodynamic/reconnection/2d/y1/roe_athena.cfg \
    --boundary-file config_files/magnetohydrodynamic/reconnection/2d/y1/boundaries.cfg



echo orszag-tang mhd 2d z1 roe_athena
$RUN ./test.exe \
    --config-file config_files/magnetohydrodynamic/orszag-tang/2d/z1/roe_athena.cfg \
    --boundary-file config_files/magnetohydrodynamic/orszag-tang/2d/z1/boundaries.cfg



echo
echo Plotting results with gnuplot

echo shock_tube mhd 1d +x hll_athena
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/shock_tube/1d/+x/hll_athena/*dc
echo shock_tube mhd 1d +x hlld_athena
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/shock_tube/1d/+x/hlld_athena/*dc
echo shock_tube mhd 1d +x roe_athena
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/shock_tube/1d/+x/roe_athena/*dc

echo shock_tube mhd 1d -x hll_athena
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/shock_tube/1d/-x/hll_athena/*dc
echo shock_tube mhd 1d -x hlld_athena
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/shock_tube/1d/-x/hlld_athena/*dc
echo shock_tube mhd 1d -x roe_athena
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/shock_tube/1d/-x/roe_athena/*dc


echo shock_tube mhd 1d +y hll_athena
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/shock_tube/1d/+y/hll_athena/*dc
echo shock_tube mhd 1d +y hlld_athena
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/shock_tube/1d/+y/hlld_athena/*dc
echo shock_tube mhd 1d +y roe_athena
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/shock_tube/1d/+y/roe_athena/*dc

echo shock_tube mhd 1d -y hll_athena
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/shock_tube/1d/-y/hll_athena/*dc
echo shock_tube mhd 1d -y hlld_athena
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/shock_tube/1d/-y/hlld_athena/*dc
echo shock_tube mhd 1d -y roe_athena
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/shock_tube/1d/-y/roe_athena/*dc


echo shock_tube mhd 1d +z hll_athena
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/shock_tube/1d/+z/hll_athena/*dc
echo shock_tube mhd 1d +z hlld_athena
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/shock_tube/1d/+z/hlld_athena/*dc
echo shock_tube mhd 1d +z roe_athena
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/shock_tube/1d/+z/roe_athena/*dc

echo shock_tube mhd 1d -z hll_athena
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/shock_tube/1d/-z/hll_athena/*dc
echo shock_tube mhd 1d -z hlld_athena
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/shock_tube/1d/-z/hlld_athena/*dc
echo shock_tube mhd 1d -z roe_athena
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/shock_tube/1d/-z/roe_athena/*dc



echo kelvin-helmholtz hd 2d z1 hll_athena
$RUN ./mhd2gnuplot.exe results/hydrodynamic/kelvin-helmholtz/2d/z1/hll_athena/*dc
echo kelvin-helmholtz hd 2d z1 hlld_athena
$RUN ./mhd2gnuplot.exe results/hydrodynamic/kelvin-helmholtz/2d/z1/hlld_athena/*dc
echo kelvin-helmholtz hd 2d z1 roe_athena
$RUN ./mhd2gnuplot.exe results/hydrodynamic/kelvin-helmholtz/2d/z1/roe_athena/*dc


echo kelvin-helmholtz mhd 2d z1 hll_athena
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/kelvin-helmholtz/2d/z1/hll_athena/*dc
echo kelvin-helmholtz mhd 2d z1 hlld_athena
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/kelvin-helmholtz/2d/z1/hlld_athena/*dc
echo kelvin-helmholtz mhd 2d z1 roe_athena
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/kelvin-helmholtz/2d/z1/roe_athena/*dc



echo bow_shock hd 2d z1 hll_athena
$RUN ./mhd2gnuplot.exe results/hydrodynamic/bow_shock/2d/z1/hll_athena/*dc
echo bow_shock hd 2d z1 hlld_athena
$RUN ./mhd2gnuplot.exe results/hydrodynamic/bow_shock/2d/z1/hlld_athena/*dc


echo bow_shock mhd 2d z1 hll_athena
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/bow_shock/2d/z1/hll_athena/*dc
echo bow_shock mhd 2d z1 hlld_athena
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/bow_shock/2d/z1/hlld_athena/*dc



echo blast_wave hd 2d z1 hll_athena
$RUN ./mhd2gnuplot.exe results/hydrodynamic/blast_wave/2d/z1/hll_athena/*dc
echo blast_wave hd 2d z1 hlld_athena
$RUN ./mhd2gnuplot.exe results/hydrodynamic/blast_wave/2d/z1/hlld_athena/*dc
echo blast_wave hd 2d z1 roe_athena
$RUN ./mhd2gnuplot.exe results/hydrodynamic/blast_wave/2d/z1/roe_athena/*dc


echo blast_wave mhd 2d z1 hll_athena
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/blast_wave/2d/z1/hll_athena/*dc
echo blast_wave mhd 2d z1 hlld_athena
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/blast_wave/2d/z1/hlld_athena/*dc



echo reconnection mhd 2d y1 hll_athena
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/reconnection/2d/y1/hll_athena/*dc
echo reconnection mhd 2d y1 hlld_athena
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/reconnection/2d/y1/hlld_athena/*dc
echo reconnection mhd 2d y1 roe_athena
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/reconnection/2d/y1/roe_athena/*dc



echo orszag-tang mhd 2d z1 roe_athena
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/orszag-tang/2d/z1/roe_athena/*.dc
