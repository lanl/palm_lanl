#!/usr/bin/ksh

# kpp4palm - script for creating gasphase module

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
# $Id: kpp4palm.ksh 2718 2018-01-02 08:49:38Z maronga $ 
# Initial revision
# 
# 
#
#
# Other notes:
# ------------# 
# Re-introduced relative path for KPP_HOME
# Subroutine list adapted to lowercase subroutine names
# Added arr2, removed update_sun and k_3rd from subroutine list
# Renamed output file to chem_gasphase_mod
# Renamed this fikle from kp4/ksh to kpp4kpp.ksh
# changed location of def_mechanism directories to GASPHASE_PREPROC/mechanisms
#
# Nov. 2016: Initial Version of KPP chemistry convertor by Klaus Ketelsen
#
#

set -eu


########################### User SetUp ####################################

export KPP_HOME=`pwd`/kpp
export KPP=$KPP_HOME/bin/kpp

BASE=`pwd`/kpp4palm

########################## End User Setup ################################

WORK=tmp_kpp4palm

# Default

OUTDIR=`pwd`/../../../SOURCE
OUTFILE=chem_gasphase_mod
DEFDIR=`pwd`/mechanisms/def_smog
PREFIX=chem_gasphase_mod
MODE="scalar"
VLEN=1
KEEP="NO"
DE_INDEX="NO"
DE_INDEX_FAST="NO"

export KPP_SOLVER=Rosenbrock

# get Command line option

echo xxxxxxxxxx
while  getopts :d:ifkp:o:s:v:w:  c     # get options
do case $c in
      d)   DEFDIR=$OPTARG;;          # directory of definition files

      i)   DE_INDEX="YES";;          # if set, deindexing

      f)   DE_INDEX_FAST="YES";;     # if set, fast deindexing

      k)   KEEP="YES";;              # keep Working directory

      o)   OUTDIR=$OPTARG;;          # Output directory of Generated Code

      p)   PREFIX=$OPTARG;;          # Name Prefix

      s)   KPP_SOLVER=$OPTARG;;      # Name Prefix

      v)   MODE="vector"
           VLEN=$OPTARG;;            # Set to vector Mode

      w)   WORK=$OPTARG;;            # Working directory

      \?)  print ${0##*/} "unknown option:" $OPTARG
           print "USAGE: ${0##*/} [ -d dir -e -k -o dir -p name -s solver -v length -w dir ] "
           exit 1;;
   esac
done
shift OPTIND-1

echo $DEFDIR

DEF_PREFIX=${PREFIX}.kpp

# Create or clean working directory

MY_PWD=`pwd`
mkdir -p $WORK
rm -rf $WORK/*
cd $WORK

# kpp dependend, may be changed

KPP_FILE_LIST="Initialize Integrator LinearAlgebra Jacobian Function Rates Util"


KPP_SUBROUTINE_LIST="initialize"
KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST integrate fun"
KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST kppsolve kppdecomp wlamch wlamch_add"
KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST jac_sp k_arr "
KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST update_rconst arr2"
KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST initialize_kpp_ctrl error_output"

# if [[ $MODE = "vector" && $KPP_SOLVER = "ROS2" ]]
# then
#   cp $BASE/templates/${KPP_SOLVER}_vec.f90 ${KPP_SOLVER}.f90    # get vector Solver 
# else
# #  KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST FunTemplate JacTemplate Update_SUN "
#   KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST WCOPY WSCAL WAXPY"
#   if [[ $MODE = "vector" ]]
#   then
#     cp $BASE/templates/${KPP_SOLVER}_vec.f90 ${KPP_SOLVER}.f90  # get vector Solver 
#   else
#     KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST Rosenbrock  FunTemplate JacTemplate Update_SUN"
#   fi
# fi
 if [[ $MODE = "vector" ]]
 then
   # get vector Solver 
   cp $BASE/templates/${KPP_SOLVER}_vec.f90 ${KPP_SOLVER}.f90
fi

# Interface ignore list
KPP_INTERFACE_IGNORE="waxpy wcopy"

case $KPP_SOLVER in
    ROS2) ;;

    Rosenbrock)   
      KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST wcopy wscal waxpy"
    if [[ $MODE != "vector" ]]
    then
      KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST rosenbrock  funtemplate jactemplate"
    fi;;

    rosenbrock_mz)
      KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST WCOPY WSCAL WAXPY"
      KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST Rosenbrock  FunTemplate JacTemplate Update_SUN";;

    rosenbrock)
      KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST WCOPY WSCAL WAXPY"
      KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST rosenbrock  funtemplate jactemplate";;

    kpp_lsode)
      KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST WCOPY WSCAL WAXPY"
      KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST KppLsode DLSODE JAC_CHEM FUN_CHEM"
      KPP_INTERFACE_IGNORE="$KPP_INTERFACE_IGNORE JAC_CHEM KppDecomp KppSolve";;

    kpp_radau5)
      KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST WCOPY WSCAL WAXPY FUN_CHEM JAC_CHEM SET2ZERO"
      KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST RADAU5 Update_SUN"
      KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST KppSolveCmplx KppDecompCmplx";;

    kpp_sdirk)
       KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST WCOPY WSCAL WAXPY"
       KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST SDIRK JAC_CHEM SET2ZERO FUN_CHEM"
       KPP_INTERFACE_IGNORE="$KPP_INTERFACE_IGNORE Set2zero SET2ZERO FUN_CHEM";;

    kpp_seulex)
       KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST WCOPY WSCAL WAXPY"
       KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST ATMSEULEX"
       KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST SEULEX_ErrorMsg SEULEX_Integrator FUN_CHEM JAC_CHEM SEUL"
       KPP_INTERFACE_IGNORE="$KPP_INTERFACE_IGNORE SEULEX_Integrator SDIRK FUN_CHEM SEUL";;

   \?)  print "SORRY ONLY ROSENBROCK METHODS WORK AT THE MOMENT:" $KPP_SOLVER
        exit 1;;
esac
#mz-ak-20070509+

KPP_INCLUDE_LIST="Parameters Global JacobianSP Monitor"

#Get definition Files

cp $DEFDIR/*.eqn         .
cp $DEFDIR/*.spc         .
cp $DEFDIR/${PREFIX}.kpp     .

# Run kpp

$KPP $DEF_PREFIX

# Get templates for C++ program

cp $BASE/templates/module_header* .           # Use fixed Module_header
cp $BASE/templates/initialize_kpp_ctrl_template.f90 .  # CTRL kpp time stepping

# file with subroutine list for c++ program create_kpp_module

for i in $KPP_FILE_LIST
do
  echo ${PREFIX}_${i} >> file_list
done
echo initialize_kpp_ctrl_template >> file_list

# file with subroutine list for c++ program create_kpp_module

for i in $KPP_SUBROUTINE_LIST
do
  echo $i >> subroutine_list
done

# file with include list for c++ program create_kpp_module

for i in $KPP_INCLUDE_LIST
do
  echo ${PREFIX}_${i} >> include_list
done

touch interface_ignore_list
for i in $KPP_INTERFACE_IGNORE
do
  echo $i >> interface_ignore_list
done

$BASE/bin/kpp4palm.exe $PREFIX $MODE $VLEN $DE_INDEX $DE_INDEX_FAST


if [[ -e $OUTDIR/${OUTFILE}.f90 ]] 
then 
 mv $OUTDIR/${OUTFILE}.f90 $OUTDIR/${OUTFILE}.f90.sav
fi
cp -p kk_kpp.f90    $OUTDIR/${OUTFILE}.f90
#cp -p kk_kpp.f90    $MY_PWD/../SOURCE/${OUTFILE}.f90

echo " "
echo "Write kpp module -- > " $OUTDIR/${OUTFILE}.f90

if [[ $KEEP = "NO" ]]
then
  cd  $MY_PWD
  rm -rf $WORK
fi
exit

