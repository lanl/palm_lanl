#!/usr/bin/env python
"""
Name: diff_files.py
Author: Luke van Roekel

Computes difference in two files from PALM and returns any diffs above a tolerance

Example Call
	./diff_files.py 
	   -ref /usr/projects/climate/lvanroekel/palm_les_lanl/trunk/TESTS/JOBS/convection
	   -test /lustre/scratch4/turquoise/lvanroekel/palm/jobs/test
	   -tol 0
	   -time 2
"""

import numpy as np
import xarray

def check_files(basedir, testdir, tslice=-1, tol='0'):
	if basedir[-1] =='/':
		basepr = basedir+'DATA_1D_PR_NETCDF'
		base3d = basedir+'DATA_3D_NETCDF'
	else:
		basepr = basedir+'/DATA_1D_PR_NETCDF'
		base3d = basedir+'/DATA_3D_NETCDF'
	
	if testdir[-1] == '/':
		testpr = testdir+'DATA_1D_PR_NETCDF'
		test3d = testdir+'DATA_3D_NETCDF'
	else:
		testpr = testdir+'/DATA_1D_PR_NETCDF'
		test3d = testdir+'/DATA_3D_NETCDF'

	ds1pr = xarray.open_dataset(basepr).isel(time=tslice)
	ds2pr = xarray.open_dataset(testpr).isel(time=tslice)

	ds13d = xarray.open_dataset(base3d).isel(time=tslice)
	ds23d = xarray.open_dataset(test3d).isel(time=tslice)

	for val in ds1pr.dims.iterkeys():
		if ds1pr.dims[val] != ds2pr.dims[val]:
        		print 'ERROR, dimension mismatch in DATA_1D_PR for dimension = ',val,'check file inputs'
        
	for val in ds13d.dims.iterkeys():
    		if ds13d.dims[val] != ds23d.dims[val]:
        		print 'ERROR, dimension mismatch in DATA_3D for dimension = ',val,'check file inputs'

	passVal1d = True
	for key in ds1pr.iterkeys():
    		if key != 'time':
        		diff = abs(ds1pr[key].values - ds2pr[key].values)
		else: 
			diff = np.zeros(10)

        	if diff.max() > tol:
            		print 'failed diff in 1d profiles, field = ',key
           	 	passVal1d = False

	passVal3d = True
	for key in ds13d.iterkeys():
    		if key != 'time':
        		diff = abs(ds13d[key] - ds23d[key])
		else:
			diff = np.zeros(10)

        	if diff.max() > tol:
            		print 'failed diff in 3d profiles, field = ',key
            		passVal3d = False
            
	if passVal1d:
    		print '1D_PR files match within tolerance of ',tol
	if passVal3d:
    		print '3d files match within tolerance of ',tol

if __name__ == "__main__":
	import argparse
	parser = argparse.ArgumentParser(description=__doc__,
					 formatter_class=argparse.RawTextHelpFormatter)
	parser.add_argument("-ref", "--ref_run_directory", dest="reference_directory", 
			help="path to run directory of reference case", metavar="REF",
			required=True)
	parser.add_argument("-test", "--test_run_directory", dest="test_directory",
			help="path to test directory for comparison", metavar="TEST")
	parser.add_argument("-time", "--time_slice", dest="time_select",
			help="time slice to compare fields in ", metavar="TIME",
			required=False,default=-1)
	parser.add_argument("-tol", "--tolerance", dest="err_tol",
			help="error tolerance for comparison, defaults to zero", 
			metavar="TOLERANCE",default=0)
	args = parser.parse_args()

	check_files(basedir=args.reference_directory,testdir=args.test_directory, 
			tslice=args.time_select,tol=args.err_tol)
