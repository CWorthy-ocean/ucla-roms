# ==============================================================
#
# SAMPLE SOURCE CODE - SUBJECT TO THE TERMS OF SAMPLE CODE LICENSE AGREEMENT,
# http://software.intel.com/en-us/articles/intel-sample-source-code-license-agreement/
#
# Copyright 2017 Intel Corporation
#
# THIS FILE IS PROVIDED "AS IS" WITH NO WARRANTIES, EXPRESS OR IMPLIED, INCLUDING BUT
# NOT LIMITED TO ANY IMPLIED WARRANTY OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
# PURPOSE, NON-INFRINGEMENT OF INTELLECTUAL PROPERTY RIGHTS.
#
# =============================================================
#
#******************************************************************************
#  Content:
#      NQueens Fortran sample
#******************************************************************************
#
# To compile with the Intel(R) Fortran Compiler for Linux*:
#
# source <path_to_compiler_bin>/compilervars.sh or compilervars.csh
#
# > make
#

SHELL = /bin/sh


FC = ifort
#FFLAG = -xHost -g -debug inline-debug-info -parallel-source-info=2 -qopt-report=3
#FFLAG = -O3 -g -debug inline-debug-info -parallel-source-info=2 -qopt-report=3
FFLAG = -O3 -qopt-report=3
#FFLAG = -O1 -g -debug inline-debug-info -parallel-source-info=2 -qopt-report=3
#FFLAG = -g -traceback -check all
#OMP_LIB = qopenmp

SRCDIR = .

PROGS = imp_loop exp_loop_do

all:: $(PROGS)

imp_loop: $(SRCDIR)/imp_loop.f90
	$(FC) $(FFLAG) $< -o $@
#       $(FC) $(FFLAG) -$(OMP_LIB) $< -o $@


exp_loop_do: $(SRCDIR)/exp_loop_do.f90
	$(FC) $(FFLAG) $< -o $@


clean:
	- rm -f core
	- rm -f $(PROGS) *.o

