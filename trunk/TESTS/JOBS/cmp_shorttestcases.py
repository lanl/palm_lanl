#!/usr/bin/env python2

# Author Carolyn Begeman

from diff_files import check_files

filedir1a = 'surface_momentumflux/'
filedir2a = 'surface_cooling/'
filedir1b = '' # enter path to new case with surface momentum fluxes
filedir2b = '' # enter path to new case with surface cooling

check_files(filedir1a, filedir1b, tslice=5, tol='0')
check_files(filedir2a, filedir2b, tslice=5, tol='0')