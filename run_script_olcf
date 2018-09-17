#!/bin/bash

#PBS  -N PALM_LES
#PBS  -j oe
#PBS  -m ae
#PBS -A cli127
#PBS -l nodes=16
#PBS -q debug
#PBS -l walltime=01:00:00

#TODO: if threads, need a way to pass to PALM

exec_path="/path/to/case/directory"

cd $PBS_O_WORKDIR
export PATH="/ccs/home/vanroek/palm/current_version/trunk/SCRIPTS:${PATH}"
module load cray-netcdf-hdf5parallel/4.4.1.1 
source $MODULESHOME/init/bash # May already be included if using modules

#Calculate Walltime in seconds from PBS commands for palm run command
WTIME2=$(qstat -f $PBS_JOBID | sed -rn 's/.*Resource_List.walltime = (.*)/\1/p')
echo $WTIME2
hr=${WTIME2:0:2}
min=${WTIME2:3:2}
sec=${WTIME2:6:2}
WTIME=$((hr*3600+min*60+sec))
echo "WTIME (seconds) = $WTIME"

#Calculate number of PES and num tasks per node for palm as well
NP=$(wc -l $PBS_NODEFILE | awk '{print $1}')
echo "Total CPU count = $NP"
NNODES=`cat $PBS_NODEFILE | uniq -c | wc -l`
NN=$((NP/NNODES))
echo "n_tasks_per_node= $NN"

#Run PALM
palm_simple_run -b ifort.titan.nc4 -w $exec_path -c aprun -s test_oceanml -n $NN -p $NP -t 1 -T $WTIME 
