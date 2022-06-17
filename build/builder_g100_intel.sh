#! /bin/bash

module load autoload
module load intel/oneapi-2021 --binary
module load intelmpi/oneapi-2021--binary
module load netcdf/4.7.4--oneapi--2021.2.0-ifort
module load netcdff/4.5.3--oneapi--2021.2.0-ifort

ln -s ../macro/make.cal1_ifort Makefile.macro
make
make examples

