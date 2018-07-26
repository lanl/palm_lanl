#--------------------------------------------------------------------------------#
# This file is part of the PALM model system.
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
# Copyright 1997-2018  Leibniz Universitaet Hannover
#--------------------------------------------------------------------------------#
#
# WARNING: don't write filenames with extension .f90 in this header!!!!
#
# Current revisions:
# ------------------
# 
# 
# Former revisions:
# -----------------
# $Id: Makefile 2718 2018-01-02 08:49:38Z maronga $
#
# 1585 2015-04-30 07:05:52Z maronga
# Initial revision
# 
# Description:
# ------------
# Makefile for the RRTMG library used in PALM
#------------------------------------------------------------------------------!

PROG =  librrtmg

SOURCES = rrsw_aer.f90 rrsw_kg18.f90 rrsw_kg23.f90 rrsw_kg28.f90 rrsw_vsn.f90  \
          rrtmg_sw_reftra.f90 rrtmg_sw_vrtqdr.f90             \
          rrsw_cld.f90 rrsw_kg19.f90 rrsw_kg24.f90    \
          rrsw_kg29.f90 rrsw_wvn.f90 rrtmg_sw_setcoef.f90                      \
          rrsw_con.f90 rrsw_kg20.f90 rrsw_kg25.f90     \
          rrsw_ncpar.f90 rrtmg_sw_rad.nomcica.f90                 \
          parkind.f90 rrsw_kg16.f90 rrsw_kg21.f90          \
          rrsw_kg26.f90 rrsw_ref.f90 rrtmg_sw_cldprop.f90                      \
          rrtmg_sw_spcvrt.f90 parrrsw.f90             \
          rrsw_kg17.f90 rrsw_kg22.f90 rrsw_kg27.f90 rrsw_tbl.f90              \
          rrtmg_sw_init.f90 rrtmg_sw_read_nc.f90 rrtmg_sw_taumol.f90 \
        rrlw_cld.f90  rrlw_kg01.f90  rrlw_kg03.f90  rrlw_kg05.f90  rrlw_kg07.f90  \
        rrlw_kg09.f90  rrlw_kg11.f90  rrlw_kg13.f90  rrlw_kg15.f90  rrlw_ncpar.f90   \
	rrlw_tbl.f90  rrlw_wvn.f90 rrlw_con.f90  rrlw_kg02.f90  rrlw_kg04.f90  rrlw_kg06.f90   \
        rrlw_kg08.f90  rrlw_kg10.f90  rrlw_kg12.f90  rrlw_kg14.f90  rrlw_kg16.f90  rrlw_ref.f90   \
        rrlw_vsn.f90 rrtmg_lw_cldprmc.f90  rrtmg_lw_init.f90 \
        rrtmg_lw_rad.nomcica.f90  rrtmg_lw_rtrn.f90    rrtmg_lw_rtrnmr.f90  \
        rrtmg_lw_taumol.f90 rrtmg_lw_cldprop.f90   parrrtm.f90 \
        rrtmg_lw_read_nc.f90 rrtmg_lw_rtrnmc.f90  rrtmg_lw_setcoef.f90


OBJS=$(SOURCES:.f90=.o)

CC = cc
CFLAGS = -O

F90 = 
COPT = 
F90FLAGS = 

#to_be_replaced_by_include

.SUFFIXES:
.SUFFIXES: .o .f90


all: $(PROG)

# create static library
$(PROG): $(OBJS)
	$(F90) $(F90FLAGS) $(COPT) -o $(PROG).so -shared $(OBJS)

clean:
	rm -f $(PROG) $(OBJS) *.mod *.i *.lst *.so *.o  *.a

.f90.o:
	$(F90) $(F90FLAGS) $(COPT) -fpic -c $<



parrrtm.o: parkind.o
rrtmg_lw_init.o: parkind.o rrlw_wvn.o rrtmg_lw_setcoef.o rrtmg_lw_rtrnmr.o rrtmg_lw_taumol.o
rrtmg_lw_rad.nomcica.o: parkind.o rrtmg_lw_cldprop.o rrtmg_lw_rtrn.o
rrtmg_lw_read_nc.o: parkind.o rrlw_kg01.o rrlw_kg02.o rrlw_kg03.o rrlw_kg04.o rrlw_kg05.o rrlw_kg06.o rrlw_kg07.o rrlw_kg08.o rrlw_kg09.o rrlw_kg10.o rrlw_kg11.o rrlw_kg12.o rrlw_kg13.o rrlw_kg14.o rrlw_kg15.o rrlw_kg16.o
rrtmg_lw_rtrnmr.o: parkind.o
rrtmg_lw_taumol.o: parkind.o
rrlw_kg01.o: parkind.o parrrtm.o
rrlw_kg02.o: parkind.o parrrtm.o
rrlw_kg03.o: parkind.o parrrtm.o
rrlw_kg04.o: parkind.o parrrtm.o
rrlw_kg05.o: parkind.o parrrtm.o
rrlw_kg06.o: parkind.o parrrtm.o
rrlw_kg07.o: parkind.o parrrtm.o
rrlw_kg08.o: parkind.o parrrtm.o
rrlw_kg09.o: parkind.o parrrtm.o
rrlw_kg10.o: parkind.o parrrtm.o
rrlw_kg11.o: parkind.o parrrtm.o
rrlw_kg12.o: parkind.o parrrtm.o
rrlw_kg13.o: parkind.o parrrtm.o
rrlw_kg14.o: parkind.o parrrtm.o
rrlw_kg15.o: parkind.o parrrtm.o
rrlw_kg16.o: parkind.o parrrtm.o

parrrsw.o: parkind.o
rrtmg_sw_cldprop.o: parkind.o rrsw_cld.o rrsw_wvn.o rrsw_vsn.o
rrtmg_sw_init.o: parkind.o rrsw_wvn.o rrtmg_sw_setcoef.o
rrtmg_sw_rad.nomcica.o: parkind.o rrtmg_sw_spcvrt.o rrtmg_sw_cldprop.o
rrtmg_sw_read_nc.o: parkind.o rrsw_kg16.o rrsw_kg17.o rrsw_kg18.o rrsw_kg19.o rrsw_kg20.o rrsw_kg21.o rrsw_kg22.o rrsw_kg23.o rrsw_kg24.o rrsw_kg25.o rrsw_kg26.o rrsw_kg27.o rrsw_kg28.o rrsw_kg29.o rrsw_ncpar.o
rrtmg_sw_reftra.o: parkind.o rrsw_tbl.o rrsw_vsn.o
rrtmg_sw_spcvrt.o: parkind.o parrrsw.o rrsw_tbl.o rrsw_vsn.o rrsw_wvn.o rrtmg_sw_reftra.o rrtmg_sw_taumol.o rrtmg_sw_vrtqdr.o
rrtmg_sw_setcoef.o: parkind.o rrsw_ref.o
rrtmg_sw_taumol.o: parkind.o rrsw_con.o rrsw_vsn.o rrsw_wvn.o rrsw_kg16.o rrsw_kg17.o rrsw_kg18.o rrsw_kg19.o rrsw_kg20.o rrsw_kg21.o rrsw_kg22.o rrsw_kg23.o rrsw_kg24.o rrsw_kg25.o rrsw_kg26.o rrsw_kg27.o rrsw_kg28.o rrsw_kg29.o
rrtmg_sw_vrtqdr.o: parkind.o
rrsw_aer.o: parkind.o parrrsw.o 
rrsw_cld.o: parkind.o
rrsw_con.o: parkind.o
rrsw_kg16.o: parkind.o parrrsw.o
rrsw_kg17.o: parkind.o parrrsw.o
rrsw_kg18.o: parkind.o parrrsw.o
rrsw_kg19.o: parkind.o parrrsw.o
rrsw_kg20.o: parkind.o parrrsw.o
rrsw_kg21.o: parkind.o parrrsw.o
rrsw_kg22.o: parkind.o parrrsw.o
rrsw_kg23.o: parkind.o parrrsw.o
rrsw_kg24.o: parkind.o parrrsw.o
rrsw_kg25.o: parkind.o parrrsw.o
rrsw_kg26.o: parkind.o parrrsw.o
rrsw_kg27.o: parkind.o parrrsw.o
rrsw_kg28.o: parkind.o parrrsw.o
rrsw_kg29.o: parkind.o parrrsw.o
rrsw_ncpar.o: parkind.o
rrsw_ref.o: parkind.o
rrsw_tbl.o: parkind.o
rrsw_vsn.o: parkind.o parrrsw.o
rrsw_wvn.o: parkind.o parrrsw.o 