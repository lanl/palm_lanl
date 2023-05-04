#!/usr/bin/env python3
# -------------------------------------------------------------------- #
# This file is part of the PALM model system.
#
# PALM is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your
# option) any later version.
#
# PALM is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
# for more details.
#
# You should have received a copy of the GNU General Public License
# along with PALM. If not, see <http://www.gnu.org/licenses/>.
#
# Copyright 1997-2021 Leibniz Universitaet Hannover
# -------------------------------------------------------------------- #

import datetime
import math
import pytz

from netCDF4 import Dataset
import numpy as np


class StaticDriver:
    """This is an example script to generate static drivers for PALM.

    You can use it as a starting point for creating your setup specific
    driver.
    """

    def __init__(self):
        """Open the static driver as netCDF4 file. Here, you have to
        give the full path to the static driver that shall be created.
        Existing file with same name is deleted.
        """
        print('Opening file...')
        self.nc_file = Dataset('PIDS_STATIC', 'w', format='NETCDF4')

    def write_global_attributes(self):
        """Write global attributes to static driver."""
        print("Writing global attributes...")

        # Optional global attributes
        # --------------------------
        self.nc_file.title = 'PALM sea ice static driver'
        self.nc_file.author = 'cbegeman'
        self.nc_file.institution = 'Los Alamos National Laboratory'
        self.nc_file.comment = 'ice cover characteristic of sea ice lead'
        self.nc_file.creation_date = \
            pytz.utc.localize(datetime.datetime.utcnow()).strftime('%y-%m-%d %H:%M:%S %z')[:-2]
        self.nc_file.history = ''
        self.nc_file.keywords = ''
        self.nc_file.license = ''
        self.nc_file.palm_version = ''
        self.nc_file.references = ''
        self.nc_file.source = 'PALM trunk'
        self.nc_file.version = '1'

        # Mandatory global attributes
        # ---------------------------
        self.nc_file.Conventions = 'CF-1.7'
        self.nc_file.origin_lat = 52.50965  # (overwrite initialization_parameters)
        self.nc_file.origin_lon = 13.3139  # Used to initialize Coriolis parameter
        self.nc_file.origin_time = '2019-03-06 10:00:00 +00'
        self.nc_file.origin_x = 3455249.0
        self.nc_file.origin_y = 5424815.0
        self.nc_file.origin_z = 0.0
        self.nc_file.rotation_angle = 0.0

    def define_dimensions(self):
        """Set dimensions on which variables are defined."""
        print("Writing dimensions...")

        # Specify general grid parameters
        # These values must equal to those set in the initialization_parameters
        self.nx = 31
        self.ny = 31
        self.nz = 32
        dx = 2.5
        dy = 2.5
        dz = 2.5

        # Coordinates
        # -----------
        self.nc_file.createDimension('x', self.nx+1)
        self.x = self.nc_file.createVariable('x', 'f4', ('x',))
        self.x.long_name = 'distance to origin in x-direction'
        self.x.units = 'm'
        self.x.axis = 'X'
        lx = (self.nx+1)*dx
        self.x[:] = np.arange(0, lx, dx) + 0.5 * dx

        self.nc_file.createDimension('y', self.ny+1)
        self.y = self.nc_file.createVariable('y', 'f4', ('y',))
        self.y.long_name = 'distance to origin in y-direction'
        self.y.units = 'm'
        self.y.axis = 'Y'
        self.y[:] = np.arange(0, (self.ny+1)*dy, dy) + 0.5 * dy

        # NOTE if your simulation uses a stretched vertical grid, you
        # need to modify the z coordinates, e.g. z_array = (...)
        z_array = np.append(0, np.arange(dz/2, (self.nz)*dz, dz))
        self.nc_file.createDimension('z', self.nz+1)
        self.z = self.nc_file.createVariable('z', 'f4', ('z',))
        self.z.long_name = 'height above origin'
        self.z.units = 'm'
        self.z.axis = 'Z'
        self.z.positive = 'up'
        self.z[:] = z_array

        zlad_array = self.z[:6]
        self.nc_file.createDimension('zlad', len(zlad_array))
        self.zlad = self.nc_file.createVariable('zlad', 'f4', ('zlad',))
        self.zlad.long_name = 'height above ground'
        self.zlad.units = 'm'
        self.zlad.axis = 'Z'
        self.zlad.positive = 'up'
        self.zlad[:] = zlad_array

        self.nc_file.createDimension('ntop_surface_fraction', 1)
        self.ntop_surface_fraction = self.nc_file.createVariable(
            'ntop_surface_fraction', 'i4', ('ntop_surface_fraction',))
        self.ntop_surface_fraction[:] = np.arange(1)

    def define_variables(self):
        """Define variables for the static driver.

        Be aware that some variables depend on others. For a description
        of each variable, please have a look at the documentation at
        palm.muk.uni-hannover.de/trac/wiki/doc/app/iofiles/pids/static

        An example of how you modify the variables is given below:

        building_2d_array = np.ones((self.ny+1, self.nx+1)) * -9999.0
        south_wall, north_wall, left_wall, right_wall = 20, 25, 20, 25
        building_2d_array[
            south_wall:north_wall,
            left_wall:right_wall
            ] = 50
        nc_buildings_2d = self.nc_file.createVariable(
            'buildings_2d', 'f4', ('y', 'x'), fill_value=-9999.0)
        nc_buildings_2d.lod = 1
        nc_buildings_2d[:, :] = building_2d_array
        """
        print("Writing variables...")

        # Main surface clasification
        # --------------------------
        nc_top_surface_fraction = self.nc_file.createVariable(
            'top_surface_fraction', 'f4', ('ntop_surface_fraction', 'y', 'x'), fill_value=-9999.0)
        nc_top_surface_fraction.long_name = "top surface fraction of ice"
        nc_top_surface_fraction.units = "1"
        nc_top_surface_fraction[0, :, :] = nc_top_surface_fraction._FillValue  # ice fraction

        lx = self.x[-1]
        print(self.x[:])
        xarray, _ = np.meshgrid(self.x[:], self.y[:])
        ice_mask = np.logical_or(xarray < lx/3, xarray > 2*lx/3)
        nc_top_surface_fraction[0, :, :] = np.where(ice_mask, 1.0, 0.0)

    def finalize(self):
        """Close file."""
        print("Closing file...")

        self.nc_file.close()


if __name__ == '__main__':
    driver = StaticDriver()
    driver.write_global_attributes()
    driver.define_dimensions()
    driver.define_variables()
    driver.finalize()

