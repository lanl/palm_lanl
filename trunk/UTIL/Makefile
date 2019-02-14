#--------------------------------------------------------------------------------#
# This file is part of PALM.
#
# PALM is free software: you can redistribute it and/or modify it under the terms
# of the GNU General Public License as published by the Free Software Foundation,
# either version 3 of the License, or (at your option) any later version.
#
# PALM is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# PALM. If not, see <http://www.gnu.org/licenses/>.
#
# Copyright 1997-2014  Leibniz Universitaet Hannover
#--------------------------------------------------------------------------------#
#
# Current revisions:
# -----------------
#
# Former revisions:
# -----------------
# $Id: Makefile 2795 2018-02-07 14:48:48Z hellstea $
#
# 1046 2012-11-09 14:38:45Z maronga
# code put under GPL (PALM 3.9)
#
# utility programs are created separately for each block in the configuration
# file,
# compare_palm_logs added
# added PROG3 in 'all:' statement (line 39, bugfix)
#
# Initial revision somewhere in 2003/04
#
# Description:
# ------------
# Makefile for generating the utility programs needed by mrun and palm
#--------------------------------------------------------------------------------#

PROG1 =  ../SCRIPTS/combine_plot_fields$(BLOCK).x
PROG2 =  ../SCRIPTS/compare_palm_logs$(BLOCK).x

OBJS1 =  combine_plot_fields.o
OBJS2 =  compare_palm_logs.o

CC = cc
CFLAGS = -O

F90 = mpif90
F90_SER = mpif90
COPT = -cpp -DMPI_REAL=MPI_DOUBLE_PRECISION -DMPI_2REAL=MPI_2DOUBLE_PRECISION -D__netcdf
F90FLAGS = 
LDFLAGS  = 

#PKG_CONFIG      := PKG_CONFIG_ALLOW_SYSTEM_CFLAGS=1 pkg-config
#NETCDF_F90FLAGS != $(PKG_CONFIG) --cflags netcdf-fortran
#NETCDF_LDFLAGS  != $(PKG_CONFIG) --libs netcdf-fortran
NETCDF_F90FLAGS = $(shell nc-config --fflags)
NETCDF_LDFLAGS  = $(shell nc-config --flibs) 
F90FLAGS        += $(NETCDF_F90FLAGS)
LDFLAGS         += $(NETCDF_LDFLAGS)

.SUFFIXES: $(SUFFIXES) .f90


all: $(PROG1) $(PROG2) $(PROG3)

$(PROG1): $(OBJS1)
	$(F90_SER) -o $(PROG1) $(OBJS1) $(LDFLAGS)

$(PROG2): $(OBJS2)
	$(F90_SER) -o $(PROG2) $(OBJS2) $(LDFLAGS)

.f90.o:
	$(F90_SER) $(F90FLAGS) $(COPT) -c $<


