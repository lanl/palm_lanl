#!/bin/ksh
#
# run_kpp4palm - script for running gasphase preprocessor

#------------------------------------------------------------------------------#
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
# Copyright 2017-2018  Klaus Ketelsen and MPI-CH (April 2007)
# Copyright 2017-2018  Karlsruhe Institute of Technology
# Copyright 2017-2018  Leibniz Universitaet Hannover
#------------------------------------------------------------------------------#
#
# Current revisions:
# ------------------
# 
# 
# Former revisions:
# -----------------
# $Id: run_kpp4palm.ksh 2718 2018-01-02 08:49:38Z maronga $ 
# Initial revision
# 
# 
#
#
# Other notes:
# ------------
# Small adaptations (added MECH as argument, re-introduced relative path,
#                    adaption to changing kp4 to kpp4palm)
# filename changed from run_kp4.ksh to run_kpp4palm.ksh
#
# Nov. 2016: Initial Version of KPP chemistry by Klaus Ketelsen
#

echo "\$1 = " $1 

if [[ $1 == "clean" ]]
then
  cd kpp
  make distclean
  cd ..

  cd kpp4palm
  make distclean                 
  cd ..
  exit
fi

export PATH=$PATH:`pwd`/kpp4palm/bin

# build kpp

cp mechanisms/UserRateLaws.f90 kpp/util
cd kpp
make
cd ..

# build kpp4palm.exe

cd kpp4palm
make install                    
cd ..

# run kpp4palm with default Setup
#    -d   DEFDIR=$OPTARG;;          # directory of definition files
#    -i   DE_INDEX="YES";;          # if set, deindexing
#    -f   DE_INDEX_FAST="YES";;     # if set, fast deindexing
#    -k   KEEP="NO";;               # keep Working directory
#    -o   OUTDIR=$OPTARG;;          # Output directory of Geneated Code
#    -p   PREFIX=$OPTARG;;          # Name Prefix
#    -s   KPP_SOLVER=$OPTARG;;      # Name Prefix

 
echo $PATH
# use smog mechnism as default
MECH=smog

while  getopts m:k  opt   # get options
do case $opt in
      m)   MECH=$OPTARG;;          # mechanism name as part of mechanisms/def_[mechanism_name]

      k)   KEEP="-k";;             # keep Working directory

      \?)  print ${0##*/} "unknown option:" $OPTARG
           print "USAGE: ${0##*/} [ -m mechanism_name] [ -k] "
           exit 1;;
   esac
done
echo $MECH

kpp4palm.ksh -d `pwd`/mechanisms/def_$MECH $KEEP
