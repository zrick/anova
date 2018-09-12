#!/bin/bash -ex

F90=gfortran
NC_INC=/opt/rrzk/lib/netcdf/4.1.3/include
NC_LIB=/opt/rrzk/lib/netcdf/4.1.3/lib
DBG=
OPT=-O2

$F90 $OPT $DBG -c src/module_timer.f90
$F90 $OPT $DBG -c src/module_anova_constants.f90
$F90 $OPT $DBG -c src/module_nc_readwrite.f90 -I$NC_INC
$F90 $OPT $DBG -c src/module_anova.f90 
$F90 $OPT $DBG -c src/anova_main.f90
$F90 $OPT $DBG -o anova.x anova_main.o module_anova.o module_nc_readwrite.o module_timer.o -L$NC_LIB -lnetcdff
rm *.o *.mod
