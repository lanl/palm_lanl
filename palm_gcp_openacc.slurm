#!/bin/bash
#SBATCH --job-name=palm_les    # Job name
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --ntasks-per-node=1
#SBATCH --partition=single-p100
#SBATCH --time=00:30:00               # Time limit hrs:min:sec
#SBATCH --output=palm_les_%j.log   # Standard output and error log

echo "HOST : " $(hostname)
date

module load pgi/18.10
module load openmpi/3.1.2/pgi/18.10
module load hdf5/1.10.3/parallel/openmpi/3.1.2/pgi/18.10
module load netcdf/4.6.1/parallel/openmpi/3.1.2/pgi/18.10
mpirun -np 1 ./palm

date
