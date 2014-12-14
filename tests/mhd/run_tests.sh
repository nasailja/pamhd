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

# shock tube
echo shock_tube mhd 1d +x hll_athena
$RUN ./test.exe \
    --solver-mhd hll_athena \
    --output-directory results/magnetohydrodynamic/shock_tube/1d/+x/hll_athena/ \
    --config-file config_files/magnetohydrodynamic/shock_tube/1d/+x/config.cfg
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/shock_tube/1d/+x/hll_athena/*dc

echo shock_tube mhd 1d +x hlld_athena
$RUN ./test.exe \
    --solver-mhd hlld_athena \
    --output-directory results/magnetohydrodynamic/shock_tube/1d/+x/hlld_athena/ \
    --config-file config_files/magnetohydrodynamic/shock_tube/1d/+x/config.cfg
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/shock_tube/1d/+x/hlld_athena/*dc

echo shock_tube mhd 1d +x roe_athena
$RUN ./test.exe \
    --solver-mhd roe_athena \
    --output-directory results/magnetohydrodynamic/shock_tube/1d/+x/roe_athena/ \
    --config-file config_files/magnetohydrodynamic/shock_tube/1d/+x/config.cfg
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/shock_tube/1d/+x/roe_athena/*dc


echo shock_tube mhd 1d -x hll_athena
$RUN ./test.exe \
    --solver-mhd hll_athena \
    --output-directory results/magnetohydrodynamic/shock_tube/1d/-x/hll_athena/ \
    --config-file config_files/magnetohydrodynamic/shock_tube/1d/-x/config.cfg
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/shock_tube/1d/-x/hll_athena/*dc

echo shock_tube mhd 1d -x hlld_athena
$RUN ./test.exe \
    --solver-mhd hlld_athena \
    --output-directory results/magnetohydrodynamic/shock_tube/1d/-x/hlld_athena/ \
    --config-file config_files/magnetohydrodynamic/shock_tube/1d/-x/config.cfg
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/shock_tube/1d/-x/hlld_athena/*dc

echo shock_tube mhd 1d -x roe_athena
$RUN ./test.exe \
    --solver-mhd roe_athena \
    --output-directory results/magnetohydrodynamic/shock_tube/1d/-x/roe_athena/ \
    --config-file config_files/magnetohydrodynamic/shock_tube/1d/-x/config.cfg
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/shock_tube/1d/-x/roe_athena/*dc


echo shock_tube mhd 1d +y hll_athena
$RUN ./test.exe \
    --solver-mhd hll_athena \
    --output-directory results/magnetohydrodynamic/shock_tube/1d/+y/hll_athena/ \
    --config-file config_files/magnetohydrodynamic/shock_tube/1d/+y/config.cfg
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/shock_tube/1d/+y/hll_athena/*dc

echo shock_tube mhd 1d +y hlld_athena
$RUN ./test.exe \
    --solver-mhd hlld_athena \
    --output-directory results/magnetohydrodynamic/shock_tube/1d/+y/hlld_athena/ \
    --config-file config_files/magnetohydrodynamic/shock_tube/1d/+y/config.cfg
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/shock_tube/1d/+y/hlld_athena/*dc

echo shock_tube mhd 1d +y roe_athena
$RUN ./test.exe \
    --solver-mhd roe_athena \
    --output-directory results/magnetohydrodynamic/shock_tube/1d/+y/roe_athena/ \
    --config-file config_files/magnetohydrodynamic/shock_tube/1d/+y/config.cfg
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/shock_tube/1d/+y/roe_athena/*dc


echo shock_tube mhd 1d -y hll_athena
$RUN ./test.exe \
    --solver-mhd hll_athena \
    --output-directory results/magnetohydrodynamic/shock_tube/1d/-y/hll_athena/ \
    --config-file config_files/magnetohydrodynamic/shock_tube/1d/-y/config.cfg
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/shock_tube/1d/-y/hll_athena/*dc

echo shock_tube mhd 1d -y hlld_athena
$RUN ./test.exe \
    --solver-mhd hlld_athena \
    --output-directory results/magnetohydrodynamic/shock_tube/1d/-y/hlld_athena/ \
    --config-file config_files/magnetohydrodynamic/shock_tube/1d/-y/config.cfg
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/shock_tube/1d/-y/hlld_athena/*dc

echo shock_tube mhd 1d -y roe_athena
$RUN ./test.exe \
    --solver-mhd roe_athena \
    --output-directory results/magnetohydrodynamic/shock_tube/1d/-y/roe_athena/ \
    --config-file config_files/magnetohydrodynamic/shock_tube/1d/-y/config.cfg
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/shock_tube/1d/-y/roe_athena/*dc


echo shock_tube mhd 1d +z hll_athena
$RUN ./test.exe \
    --solver-mhd hll_athena \
    --output-directory results/magnetohydrodynamic/shock_tube/1d/+z/hll_athena/ \
    --config-file config_files/magnetohydrodynamic/shock_tube/1d/+z/config.cfg
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/shock_tube/1d/+z/hll_athena/*dc

echo shock_tube mhd 1d +z hlld_athena
$RUN ./test.exe \
    --solver-mhd hlld_athena \
    --output-directory results/magnetohydrodynamic/shock_tube/1d/+z/hlld_athena/ \
    --config-file config_files/magnetohydrodynamic/shock_tube/1d/+z/config.cfg
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/shock_tube/1d/+z/hlld_athena/*dc

echo shock_tube mhd 1d +z roe_athena
$RUN ./test.exe \
    --solver-mhd roe_athena \
    --output-directory results/magnetohydrodynamic/shock_tube/1d/+z/roe_athena/ \
    --config-file config_files/magnetohydrodynamic/shock_tube/1d/+z/config.cfg
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/shock_tube/1d/+z/roe_athena/*dc


echo shock_tube mhd 1d -z hll_athena
$RUN ./test.exe \
    --solver-mhd hll_athena \
    --output-directory results/magnetohydrodynamic/shock_tube/1d/-z/hll_athena/ \
    --config-file config_files/magnetohydrodynamic/shock_tube/1d/-z/config.cfg
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/shock_tube/1d/-z/hll_athena/*dc

echo shock_tube mhd 1d -z hlld_athena
$RUN ./test.exe \
    --solver-mhd hlld_athena \
    --output-directory results/magnetohydrodynamic/shock_tube/1d/-z/hlld_athena/ \
    --config-file config_files/magnetohydrodynamic/shock_tube/1d/-z/config.cfg
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/shock_tube/1d/-z/hlld_athena/*dc

echo shock_tube mhd 1d -z roe_athena
$RUN ./test.exe \
    --solver-mhd roe_athena \
    --output-directory results/magnetohydrodynamic/shock_tube/1d/-z/roe_athena/ \
    --config-file config_files/magnetohydrodynamic/shock_tube/1d/-z/config.cfg
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/shock_tube/1d/-z/roe_athena/*dc



# Kelvin-Helmholtz
echo kelvin-helmholtz hd 2d z1 hll_athena
$RUN ./test.exe \
    --time-length 1.75 \
    --save-mhd-n 0.25 \
    --solver-mhd hll_athena \
    --output-directory results/hydrodynamic/kelvin-helmholtz/2d/z1/hll_athena/ \
    --config-file config_files/hydrodynamic/kelvin-helmholtz/2d/z1/config.cfg
$RUN ./mhd2gnuplot.exe results/hydrodynamic/kelvin-helmholtz/2d/z1/hll_athena/*dc

echo kelvin-helmholtz hd 2d z1 hlld_athena
$RUN ./test.exe \
    --solver-mhd hlld_athena \
    --output-directory results/hydrodynamic/kelvin-helmholtz/2d/z1/hlld_athena/ \
    --config-file config_files/hydrodynamic/kelvin-helmholtz/2d/z1/config.cfg
$RUN ./mhd2gnuplot.exe results/hydrodynamic/kelvin-helmholtz/2d/z1/hlld_athena/*dc

echo kelvin-helmholtz hd 2d z1 roe_athena
$RUN ./test.exe \
    --solver-mhd roe_athena \
    --output-directory results/hydrodynamic/kelvin-helmholtz/2d/z1/roe_athena/ \
    --config-file config_files/hydrodynamic/kelvin-helmholtz/2d/z1/config.cfg
$RUN ./mhd2gnuplot.exe results/hydrodynamic/kelvin-helmholtz/2d/z1/roe_athena/*dc


echo kelvin-helmholtz mhd 2d z1 hll_athena
$RUN ./test.exe \
    --time-length 5 \
    --save-mhd-n 0.5 \
    --solver-mhd hll_athena \
    --output-directory results/magnetohydrodynamic/kelvin-helmholtz/2d/z1/hll_athena/ \
    --config-file config_files/magnetohydrodynamic/kelvin-helmholtz/2d/z1/config.cfg
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/kelvin-helmholtz/2d/z1/hll_athena/*dc

echo kelvin-helmholtz mhd 2d z1 hlld_athena
$RUN ./test.exe \
    --time-length 2 \
    --save-mhd-n 0.25 \
    --solver-mhd hlld_athena \
    --output-directory results/magnetohydrodynamic/kelvin-helmholtz/2d/z1/hlld_athena/ \
    --config-file config_files/magnetohydrodynamic/kelvin-helmholtz/2d/z1/config.cfg
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/kelvin-helmholtz/2d/z1/hlld_athena/*dc

echo kelvin-helmholtz mhd 2d z1 roe_athena
$RUN ./test.exe \
    --solver-mhd roe_athena \
    --output-directory results/magnetohydrodynamic/kelvin-helmholtz/2d/z1/roe_athena/ \
    --config-file config_files/magnetohydrodynamic/kelvin-helmholtz/2d/z1/config.cfg
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/kelvin-helmholtz/2d/z1/roe_athena/*dc



# Bow shock
echo bow_shock hd 2d z1 hll_athena
$RUN ./test.exe \
    --solver-mhd hll_athena \
    --output-directory results/hydrodynamic/bow_shock/2d/z1/hll_athena/ \
    --config-file config_files/hydrodynamic/bow_shock/2d/z1/config.cfg
$RUN ./mhd2gnuplot.exe results/hydrodynamic/bow_shock/2d/z1/hll_athena/*dc

echo bow_shock hd 2d z1 hlld_athena
$RUN ./test.exe \
    --solver-mhd hlld_athena \
    --output-directory results/hydrodynamic/bow_shock/2d/z1/hlld_athena/ \
    --config-file config_files/hydrodynamic/bow_shock/2d/z1/config.cfg
$RUN ./mhd2gnuplot.exe results/hydrodynamic/bow_shock/2d/z1/hlld_athena/*dc


echo bow_shock mhd 2d z1 hll_athena
$RUN ./test.exe \
    --solver-mhd hll_athena \
    --output-directory results/magnetohydrodynamic/bow_shock/2d/z1/hll_athena/ \
    --config-file config_files/magnetohydrodynamic/bow_shock/2d/z1/config.cfg
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/bow_shock/2d/z1/hll_athena/*dc



# Blast wave
echo blast_wave hd 2d z1 hll_athena
$RUN ./test.exe \
    --solver-mhd hll_athena \
    --output-directory results/hydrodynamic/blast_wave/2d/z1/hll_athena/ \
    --config-file config_files/hydrodynamic/blast_wave/2d/z1/config.cfg
$RUN ./mhd2gnuplot.exe results/hydrodynamic/blast_wave/2d/z1/hll_athena/*dc

echo blast_wave hd 2d z1 hlld_athena
$RUN ./test.exe \
    --solver-mhd hlld_athena \
    --output-directory results/hydrodynamic/blast_wave/2d/z1/hlld_athena/ \
    --config-file config_files/hydrodynamic/blast_wave/2d/z1/config.cfg
$RUN ./mhd2gnuplot.exe results/hydrodynamic/blast_wave/2d/z1/hlld_athena/*dc

echo blast_wave hd 2d z1 roe_athena
$RUN ./test.exe \
    --solver-mhd roe_athena \
    --output-directory results/hydrodynamic/blast_wave/2d/z1/roe_athena/ \
    --config-file config_files/hydrodynamic/blast_wave/2d/z1/config.cfg
$RUN ./mhd2gnuplot.exe results/hydrodynamic/blast_wave/2d/z1/roe_athena/*dc


echo blast_wave mhd 2d z1 hll_athena
$RUN ./test.exe \
    --solver-mhd hll_athena \
    --output-directory results/magnetohydrodynamic/blast_wave/2d/z1/hll_athena/ \
    --config-file config_files/magnetohydrodynamic/blast_wave/2d/z1/config.cfg
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/blast_wave/2d/z1/hll_athena/*dc

echo blast_wave mhd 2d z1 hlld_athena
$RUN ./test.exe \
    --solver-mhd hlld_athena \
    --output-directory results/magnetohydrodynamic/blast_wave/2d/z1/hlld_athena/ \
    --config-file config_files/magnetohydrodynamic/blast_wave/2d/z1/config.cfg
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/blast_wave/2d/z1/hlld_athena/*dc

echo blast_wave mhd 2d z1 roe_athena
$RUN ./test.exe \
    --solver-mhd roe_athena \
    --output-directory results/magnetohydrodynamic/blast_wave/2d/z1/roe_athena/ \
    --config-file config_files/magnetohydrodynamic/blast_wave/2d/z1/config.cfg
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/blast_wave/2d/z1/roe_athena/*dc



# Reconnection
echo reconnection mhd 2d y1 hll_athena
$RUN ./test.exe --config-file config_files/magnetohydrodynamic/reconnection/2d/y1/hll_athena.cfg
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/reconnection/2d/y1/hll_athena/*dc

echo reconnection mhd 2d y1 hlld_athena
$RUN ./test.exe --config-file config_files/magnetohydrodynamic/reconnection/2d/y1/hlld_athena.cfg
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/reconnection/2d/y1/hlld_athena/*dc

echo reconnection mhd 2d y1 roe_athena
$RUN ./test.exe --config-file config_files/magnetohydrodynamic/reconnection/2d/y1/roe_athena.cfg
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/reconnection/2d/y1/roe_athena/*dc



# Orszag-Tang
echo orszag-tang mhd 2d z1 hll_athena
$RUN ./test.exe \
    --solver-mhd hll_athena \
    --output-directory results/magnetohydrodynamic/orszag-tang/2d/z1/hll_athena/ \
    --config-file config_files/magnetohydrodynamic/orszag-tang/2d/z1/config.cfg
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/orszag-tang/2d/z1/hll_athena/*dc

echo orszag-tang mhd 2d z1 roe_athena
$RUN ./test.exe \
    --solver-mhd roe_athena \
    --output-directory results/magnetohydrodynamic/orszag-tang/2d/z1/roe_athena/ \
    --config-file config_files/magnetohydrodynamic/orszag-tang/2d/z1/config.cfg
$RUN ./mhd2gnuplot.exe results/magnetohydrodynamic/orszag-tang/2d/z1/roe_athena/*dc
