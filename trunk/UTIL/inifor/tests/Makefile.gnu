#------------------------------------------------------------------------------#
# This file is part of the PALM model system.
#
# PALM is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# PALM is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# PALM. If not, see <http://www.gnu.org/licenses/>.
#
# Copyright 2017-2018 Leibniz Universitaet Hannover
# Copyright 2017-2018 Deutscher Wetterdienst Offenbach
#------------------------------------------------------------------------------#
#
# Current revisions:
# -----------------
# 
# 
# Former revisions:
# -----------------
# $Id: Makefile.gnu 2718 2018-01-02 08:49:38Z maronga $
# Initial revision
#
# 
#
# Authors:
# --------
# @author Eckhard Kadasch
#
# Description:
# ------------
# This file serves as a templete makefile for compiling INIFOR tests with
# gfortran.
#------------------------------------------------------------------------------!
PROJECT = inifor
PROJECT_PATH = ..
BIN_PATH = $(PROJECT_PATH)/bin
SRC_PATH = $(PROJECT_PATH)/src
TESTS    = $(patsubst %.f90,%,$(wildcard test-*.f90))

MODULES = util.mod $(SRC_PATH)/defs.mod $(SRC_PATH)/control.mod $(SRC_PATH)/util.mod $(SRC_PATH)/types.mod $(SRC_PATH)/transform.mod $(SRC_PATH)/io.mod $(SRC_PATH)/grid.mod
SOURCES = $(MODULES:%.mod=%.f90)
OBJECTS = $(SOURCES:%.f90=%.o)

FC      = gfortran
FFLAGS  = -g -Wall -fdefault-real-8 \
		  -ffpe-trap=invalid,zero,underflow,overflow
INCLUDE = -I/home/ekadasch/local/include -I$(SRC_PATH)
LIBRARY = -L/home/ekadasch/local/lib64 -lnetcdff -lutil
		 
.PHONY: clean test $(PROJECT)

test: information $(PROJECT) $(TESTS)

information:
	@echo $(TESTS)

$(PROJECT):
	$(MAKE) -C $(PROJECT_PATH) $(PROJECT)

test-grid: $(OBJECTS) $(MODULES) test-grid.o test-grid.mod util.mod util.o
	mkdir -p $(BIN_PATH)
	$(FC) $(FFLAGS) $(OBJECTS) $@.o -o $(BIN_PATH)/$@ -I$(SRC_PATH) $(LIBRARY)
	$(BIN_PATH)/$@

# $(TESTS) expands to files of this pattern:
test-%: $(OBJECTS) $(MODULES) test-%.o test-%.mod
	mkdir -p $(BIN_PATH)
	$(FC) $(FFLAGS) $(OBJECTS) $@.o -o $(BIN_PATH)/$@ -I$(SRC_PATH) $(LIBRARY)
	$(BIN_PATH)/$@

%.o: %.mod

util.o: util.f90
	$(FC) $(FFLAGS) -c $< -o $@

util.mod: util.f90
	$(FC) $(FFLAGS) -c $< -o $(@:%.mod=%.o)

test-%.o: test-%.f90
	$(FC) $(FFLAGS) -c $< -o $@ $(INCLUDE)

test-%.mod: test-%.f90
	$(FC) $(FFLAGS) -c $< -o $(@:%.mod=%.o) $(INCLUDE)

clean:
	rm -rf *.mod *.o
	rm -rf $(BIN_PATH)/test-*
