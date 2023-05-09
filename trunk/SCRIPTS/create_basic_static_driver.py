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
        self.nc_file = Dataset('example_static_file.nc', 'w', format='NETCDF4')

    def write_global_attributes(self):
        """Write global attributes to static driver."""
        print("Writing global attributes...")

        # Optional global attributes
        # --------------------------
        self.nc_file.title = 'Example PALM static driver'
        self.nc_file.author = 'PALM user'
        self.nc_file.institution = 'Institut of Meteorology and Climatology,' \
            'Leibniz University Hannover'
        self.nc_file.comment = 'Generic crossing example'
        self.nc_file.creation_date = \
            pytz.utc.localize(datetime.datetime.utcnow()).strftime('%y-%m-%d %H:%M:%S %z')[:-2]
        self.nc_file.history = ''
        self.nc_file.keywords = 'example, PALM-4U'
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
        self.nx = 19
        self.ny = 19
        self.nz = 20
        dx = 2
        dy = 2
        dz = 2

        # Create soil grid (only relevant if land surface module is used)
        dz_soil = np.array((0.01, 0.02, 0.04, 0.06, 0.14, 0.26, 0.54, 1.86))
        zsoil_fullLayers = np.zeros_like(dz_soil)
        zsoil_fullLayers = np.around(
            [np.sum(dz_soil[:zs]) for zs in np.arange(1, len(dz_soil)+1)], 2)
        zsoil_array = zsoil_fullLayers - dz_soil/2.

        # Coordinates
        # -----------
        self.nc_file.createDimension('x', self.nx+1)
        self.x = self.nc_file.createVariable('x', 'f4', ('x',))
        self.x.long_name = 'distance to origin in x-direction'
        self.x.units = 'm'
        self.x.axis = 'X'
        self.x[:] = np.arange(0, (self.nx+1)*dx, dx) + 0.5 * dx

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

        # self.nc_file.createDimension('zsoil', len(zsoil_array))
        # self.zsoil = self.nc_file.createVariable('zsoil', 'f4', ('zsoil',))
        # self.zsoil.long_name = 'depth in the soil'
        # self.zsoil.units = 'm'
        # self.zsoil.axis = 'Z'
        # self.zsoil.positive = 'down'
        # self.zsoil[:] = zsoil_array

        # Other dimensions
        # ----------------
        # self.nc_file.createDimension('nalbedo_pars', 3)
        # self.nalbedo_pars = self.nc_file.createVariable(
        #     'nalbedo_pars', 'i4', ('nalbedo_pars',))
        # self.nalbedo_pars[:] = np.arange(8)
        #
        # self.nc_file.createDimension('nbuilding_pars', 149)
        # self.nbuilding_pars = self.nc_file.createVariable(
        #     'nbuilding_pars', 'i4', ('nbuilding_pars',))
        # self.nbuilding_pars[:] = np.arange(149)
        #
        # self.nc_file.createDimension('npavement_pars', 4)
        # self.npavement_pars = self.nc_file.createVariable(
        #     'npavement_pars', 'i4', ('npavement_pars',))
        # self.npavement_pars[:] = np.arange(4)
        #
        # self.nc_file.createDimension('npavement_subsurface_pars', 2)
        # self.npavement_subsurface_pars = self.nc_file.createVariable(
        #     'npavement_subsurface_pars', 'i4', ('npavement_subsurface_pars',))
        # self.npavement_subsurface_pars[:] = np.arange(2)
        #
        # self.nc_file.createDimension('nsoil_pars', 8)
        # self.nsoil_pars = self.nc_file.createVariable(
        #     'nsoil_pars', 'i4', ('nsoil_pars',))
        # self.nsoil_pars[:] = np.arange(8)
        #
        self.nc_file.createDimension('nsurface_fraction', 3)
        self.nsurface_fraction = self.nc_file.createVariable(
            'nsurface_fraction', 'i4', ('nsurface_fraction',))
        self.nsurface_fraction[:] = np.arange(3)

        # self.nc_file.createDimension('nvegetation_pars', 12)
        # self.nvegetation_pars = self.nc_file.createVariable(
        #     'nvegetation_pars', 'i4', ('nvegetation_pars',))
        # self.nvegetation_pars[:] = np.arange(12)
        #
        # self.nc_file.createDimension('nwater_pars', 6)
        # self.nwater_pars = self.nc_file.createVariable(
        #     'nwater_pars', 'i4', ('nwater_pars',))
        # self.nwater_pars[:] = np.arange(7)

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

        # Topography set-up
        # -----------------
        nc_building_id = self.nc_file.createVariable(
            'building_id', 'i4', ('y', 'x'), fill_value=-9999)
        nc_building_id.long_name = "building id number"
        nc_building_id.units = "1"
        nc_building_id[:, :] = nc_building_id._FillValue

        nc_buildings_2d = self.nc_file.createVariable(
            'buildings_2d', 'f4', ('y', 'x'), fill_value=-9999.0)
        nc_buildings_2d.long_name = "building height"
        nc_buildings_2d.units = "m"
        nc_buildings_2d.lod = np.int32(1)
        nc_buildings_2d[:, :] = nc_buildings_2d._FillValue

        nc_buildings_3d = self.nc_file.createVariable(
            'buildings_3d', 'i1', ('z', 'y', 'x'), fill_value=-127)
        nc_buildings_3d.long_name = "building flag"
        nc_buildings_3d.units = "1"
        nc_buildings_3d.lod = np.int32(2)
        nc_buildings_3d[:, :, :] = nc_buildings_3d._FillValue

        nc_zt = self.nc_file.createVariable(
            'zt', 'f4', ('y', 'x'), fill_value=-9999.0)
        nc_zt.long_name = 'terrain height'
        nc_zt.units = 'm'
        nc_zt[:, :] = nc_zt._FillValue

        # nc_z0 = self.nc_file.createVariable(
        #     'z0', 'f4', ('y', 'x'), fill_value=-9999.0)
        # nc_z0.long_name = 'roughness length for momentum'
        # nc_z0.units = 'm'
        # nc_z0[:, :] = nc_z0._FillValue

        # Main surface clasification
        # --------------------------
        # nc_albedo_type = self.nc_file.createVariable(
        #     'albedo_type', 'f4', ('y', 'x'), fill_value=-9999.0)
        # nc_albedo_type.long_name = "albedo type"
        # nc_albedo_type.units = "1"
        # nc_albedo_type[:, :] = nc_albedo_type._FillValue
        #
        nc_building_type = self.nc_file.createVariable(
            'building_type', 'i1', ('y', 'x'), fill_value=-127)
        nc_building_type.long_name = "building type classification"
        nc_building_type.units = "1"
        nc_building_type[:, :] = nc_building_type._FillValue

        nc_pavement_type = self.nc_file.createVariable(
            'pavement_type', 'i1', ('y', 'x'), fill_value=-127)
        nc_pavement_type.long_name = "pavement type classification"
        nc_pavement_type.units = "1"
        nc_pavement_type[:, :] = nc_pavement_type._FillValue

        nc_soil_type = self.nc_file.createVariable(
            'soil_type', 'i1', ('y', 'x'), fill_value=-127)
        nc_soil_type.long_name = "soil type classification"
        nc_soil_type.units = "1"
        nc_soil_type.lod = np.int32(1)
        nc_soil_type[:, :] = nc_soil_type._FillValue

        # NOTE alternatively, define different soil types for different
        # soil layers
        # nc_soil_type = self.nc_file.createVariable(
        #     'soil_type', 'i1', ('zsoil', 'y', 'x'), fill_value=-127)
        # nc_soil_type.long_name = "soil type"
        # nc_soil_type.units = "1"
        # nc_soil_type.lod = np.int32(2)
        # nc_soil_type[:, :, :] = nc_soil_type._FillValue
        #
        # nc_street_crossing = self.nc_file.createVariable(
        #     'street_crossing', 'i1', ('y', 'x'),
        #     fill_value=-127)
        # nc_street_crossing.long_name = "street crossing"
        # nc_street_crossing.units = "1"
        # nc_street_crossing.valid_range = np.byte(1)
        # nc_street_crossing[:, :] = nc_street_crossing._FillValue
        #
        nc_street_type = self.nc_file.createVariable(
            'street_type', 'i1', ('y', 'x'),
            fill_value=-127)
        nc_street_type.long_name = "street type classification"
        nc_street_type.units = "1"
        nc_street_type[:, :] = nc_street_type._FillValue

        nc_surface_fraction = self.nc_file.createVariable(
            'surface_fraction', 'f4', ('nsurface_fraction', 'y', 'x'), fill_value=-9999.0)
        nc_surface_fraction.long_name = "surface fraction"
        nc_surface_fraction.units = "1"
        nc_surface_fraction[0, :, :] = nc_surface_fraction._FillValue  # vegetation fraction
        nc_surface_fraction[1, :, :] = nc_surface_fraction._FillValue  # pavement fraction
        nc_surface_fraction[2, :, :] = nc_surface_fraction._FillValue  # water fraction

        nc_vegetation_type = self.nc_file.createVariable(
            'vegetation_type', 'i1', ('y', 'x'), fill_value=-127)
        nc_vegetation_type.long_name = "vegetation type classification"
        nc_vegetation_type.units = "1"
        nc_vegetation_type[:, :] = nc_vegetation_type._FillValue

        nc_water_type = self.nc_file.createVariable(
            'water_type', 'i1', ('y', 'x'), fill_value=-127)
        nc_water_type.long_name = "water type classification"
        nc_water_type.units = "1"
        nc_water_type[:, :] = nc_water_type._FillValue

        # Pixel-based surface and sub-surface parameters
        # ----------------------------------------------
        # nc_albedo_pars = self.nc_file.createVariable(
        #     'albedo_pars', 'f4', ('nalbedo_pars', 'y', 'x'), fill_value=-9999.0)
        # nc_albedo_pars.long_name = 'albedo parameters'
        # nc_albedo_pars[:, :, :] = nc_albedo_pars._FillValue
        #
        # nc_building_pars = self.nc_file.createVariable(
        #     'building_pars', 'f4', ('nbuilding_pars', 'y', 'x'), fill_value=-9999.0)
        # nc_building_pars.long_name = 'building parameters'
        # nc_building_pars[:, :, :] = nc_building_pars._FillValue
        #
        # nc_pavement_pars = self.nc_file.createVariable(
        #     'pavement_pars', 'f4', ('npavement_pars', 'y', 'x'), fill_value=-9999.0)
        # nc_pavement_pars.long_name = "pavement parameters"
        # nc_pavement_pars[:, :, :] = nc_pavement_pars._FillValue
        #
        # nc_pavement_subsurface_pars = self.nc_file.createVariable(
        #     'pavement_subsurface_pars', 'i1',
        #     ('npavement_subsurface_pars','zsoil', 'y', 'x'), fill_value=-9999.0)
        # nc_pavement_subsurface_pars.long_name = "pavement subsurface parameters"
        # nc_pavement_subsurface_pars[:, :, :, :] = nc_pavement_subsurface_pars._FillValue
        #
        # nc_soil_pars = self.nc_file.createVariable(
        #     'soil_pars', 'f', ('nsoil_pars','zsoil', 'y', 'x'), fill_value=-9999.0)
        # nc_soil_pars.long_name = "soil parameters"
        # nc_soil_pars.lod = np.int32(2) # if set to 1, adjust dimensions of soil_pars
        # nc_soil_pars[:, :, :, :] = nc_soil_pars._FillValue
        #
        # nc_vegetation_pars = self.nc_file.createVariable(
        #     'vegetation_pars', 'f', ('nvegetation_pars', 'y', 'x'), fill_value=-9999.0)
        # nc_vegetation_pars.long_name = 'vegetation parameters'
        # nc_vegetation_pars[:, :, :] = nc_vegetation_pars._FillValue
        #
        # nc_water_pars = self.nc_file.createVariable(
        #     'water_pars', 'i1', ('nwater_pars', 'y', 'x'), fill_value=-9999.0)
        # nc_water_pars.long_name = "water parameters"
        # nc_water_pars[:, :, :] = nc_water_pars._FillValue

        # Vegetation parameters
        # ---------------------
        nc_lad = self.nc_file.createVariable(
            'lad', 'f4', ('zlad', 'y', 'x'), fill_value=-9999.0)
        nc_lad.long_name = "leaf area density"
        nc_lad.units = "m2 m-3"
        nc_lad[:, :, :] = nc_lad._FillValue

        # nc_bad = self.nc_file.createVariable(
        #     'bad', 'f4', ('zlad', 'y', 'x'), fill_value=-9999.0)
        # nc_bad.long_name = "basal area density"
        # nc_bad.units = "m2 m-3"
        # nc_bad[:, :, :] = nc_bad._FillValue
        #
        # nc_root_area_dens_r = self.nc_file.createVariable(
        #     'root_area_dens_s', 'f4', ('zsoil', 'y', 'x'), fill_value=-9999.0)
        # nc_root_area_dens_r.long_name = 'root area density of resolved vegetation'
        # nc_root_area_dens_r.units = '1'
        # nc_root_area_dens_r[:, :, :] = nc_root_area_dens_r._FillValue
        #
        # nc_root_area_dens_s = self.nc_file.createVariable(
        #     'root_area_dens_s', 'f4', ('zsoil', 'y', 'x'), fill_value=-9999.0)
        # nc_root_area_dens_s.long_name = 'root area density of parameterized vegetation'
        # nc_root_area_dens_s.units = '1'
        # nc_root_area_dens_s[:, :, :] = nc_root_area_dens_s._FillValue
        #
        # nc_tree_id = self.nc_file.createVariable(
        #     'tree_id', 'i4', ('y', 'x'), fill_value=-9999.0)
        # nc_tree_id.long_name = "tree id"
        # nc_tree_id.units = "1"
        # nc_tree_id[:, :, :] = nc_tree_id._FillValue

        # Set topography
        # --------------
        nc_building_id[:5, :5] = 1
        nc_building_id[-5:, :5] = 2
        nc_building_id[:5, -5:] = 3
        nc_building_id[-5:, -5:] = 4

        nc_buildings_2d[:, :] = np.where(
            nc_building_id[:, :] > nc_building_id._FillValue,
            nc_building_id[:, :] * 10.0,
            nc_buildings_2d[:, :])

        for k in range(self.nz+1):
            nc_buildings_3d[k, :, :] = np.where(nc_buildings_2d[:, :] > self.z[k], 1, 0)

        nc_zt[:, :] = 4.0

        # Set surface types
        # -----------------
        nc_building_type[:, :] = np.where(
            nc_building_id[:, :] > nc_building_id._FillValue,
            nc_building_id[:, :] * 1,
            nc_building_type[:, :])

        nc_pavement_type[9:11, :5] = 1
        nc_pavement_type[:, 7:13] = 2

        nc_vegetation_type[:, 5:7] = 3
        nc_vegetation_type[:, 13:15] = 3
        nc_vegetation_type[5:15, 15:] = 3
        nc_vegetation_type[7:9, 15:17] = 1
        nc_vegetation_type[11:13, 15:17] = 1

        nc_water_type[5:9, :5] = 2
        nc_water_type[11:15, :5] = 2

        nc_soil_type[:, :] = np.where(
            nc_vegetation_type[:, :] > nc_vegetation_type._FillValue,
            2,
            nc_soil_type[:, :])
        nc_soil_type[:, :] = np.where(
            nc_pavement_type[:, :] > nc_pavement_type._FillValue,
            2,
            nc_soil_type[:, :])

        nc_street_type[:, :] = np.where(nc_pavement_type[:, :] == 1, 11, nc_street_type[:, :])
        nc_street_type[:, :] = np.where(nc_pavement_type[:, :] == 2, 13, nc_street_type[:, :])

        nc_surface_fraction[0, :, :] = np.where(
            nc_building_id[:, :] > nc_building_id._FillValue,
            nc_surface_fraction[0, :, :],
            0)
        nc_surface_fraction[2, :, :] = nc_surface_fraction[1, :, :] = nc_surface_fraction[0, :, :]
        nc_surface_fraction[0, :, :] = np.where(
            nc_vegetation_type[:, :] > nc_vegetation_type._FillValue,
            1,
            nc_surface_fraction[0, :, :])
        nc_surface_fraction[1, :, :] = np.where(
            nc_pavement_type[:, :] > nc_pavement_type._FillValue,
            1,
            nc_surface_fraction[1, :, :])
        nc_surface_fraction[2, :, :] = np.where(
            nc_water_type[:, :] > nc_water_type._FillValue,
            1,
            nc_surface_fraction[2, :, :])

        # Set vegetation
        # --------------
        lad_profile = [0.0, 0.01070122, 0.1070122, 0.3130108, 0.3879193, 0.1712195]
        for k in range(len(self.zlad)):
            nc_lad[k, 7:9, 15:17] = lad_profile[k]
            nc_lad[k, 11:13, 15:17] = lad_profile[k]

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

