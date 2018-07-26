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
# $Id: Makefile 2718 2018-01-02 08:49:38Z maronga $
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
# This makefile may be used to exclusively compile INIFOR.
#------------------------------------------------------------------------------#
PROJECT = inifor
PROJECT_PATH = .
BIN_PATH  = $(PROJECT_PATH)/../../SCRIPTS
SRC_PATH  = $(PROJECT_PATH)/src
TEST_PATH = $(PROJECT_PATH)/tests

MODULES = $(SRC_PATH)/defs.mod $(SRC_PATH)/util.mod $(SRC_PATH)/control.mod $(SRC_PATH)/types.mod $(SRC_PATH)/transform.mod $(SRC_PATH)/io.mod $(SRC_PATH)/grid.mod
SOURCES = $(MODULES:%.mod=%.f90) $(SRC_PATH)/$(PROJECT).f90
OBJECTS = $(SOURCES:%.f90=%.o)

FC      = gfortran
WRNGS   = -Wall -Wextra -Wsurprising -Waliasing -Wcharacter-truncation \
		  -Wline-truncation -fbacktrace -fcheck=all -pedantic \
		  -ffpe-trap=invalid,zero,underflow,overflow
FFLAGS  = -fdefault-real-8 -Og -g $(WRNGS)
INCLUDE = -I/home/ekadasch/local/include
LIBRARY = -L/home/ekadasch/local/lib64 -lnetcdff

.PHONY: all clean doc run tags test test-verbose $(PROJECT)

$(PROJECT): $(OBJECTS) $(MODULES) tags
	@echo ""
	mkdir -p $(BIN_PATH)
	$(FC) $(FFLAGS) $(OBJECTS) -o $(BIN_PATH)/$(PROJECT) -I$(SRC_PATH) $(LIBRARY)

$(SRC_PATH)/%.o: $(SRC_PATH)/%.f90
	@echo ""
	$(FC) $(FFLAGS) -c $< -o $@ -J$(SRC_PATH) $(INCLUDE)

$(SRC_PATH)/%.mod: $(SRC_PATH)/%.f90
	@echo ""
	$(FC) $(FFLAGS) -c $< -o $(@:%.mod=%.o) -J$(SRC_PATH) $(INCLUDE)

all: clean $(PROJECT) test doc

clean:
	@echo ""
	rm -rf $(OBJECTS) $(MODULES) $(BIN_PATH)/$(PROJECT)
	rm -rf *.mod *.o
	rm -rf ./doc/latex ./doc/html
	rm -f $(SRC_PATH)/tags
	rm -f *.files
	rm -f palm-hsurf.nc
	$(MAKE) -C $(TEST_PATH) clean
	$(MAKE) -C $(SRC_PATH) clean

doc:
	@echo ""
	doxygen ./doc/doxygen.config

run: $(PROJECT)
	@echo ""
	$(BIN_PATH)/$(PROJECT)

tags:
	@echo ""
	$(MAKE) -C $(SRC_PATH) tags

test:
	@echo "Running tests"
	$(MAKE) -C $(TEST_PATH) test | grep Test
	@echo "...done testing"

test-verbose:
	@echo "Running tests"
	$(MAKE) -C $(TEST_PATH) test
	@echo "...done testing"
