import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset

#--------------------------------
# LESProfile
#--------------------------------
class LESProfile(object):

    """LESProfile object"""

    def __init__(self, time=None, time_name='Time', time_units='s',
                       z=None, z_name='z', z_units='m',
                       data=None, data_name=None, data_units=None):
        """Initialize LESProfile

        :time: (1D numpy array/datetime object) time
        :time_name: (str, optional) name of time
        :time_units: (str, optional) units of time
        :z: (1D numpy array) vertical coordinate
        :z_name: (str) name of z
        :z_units: (str) units of z
        :data: (2D numpy array) data at each time and z
        :data_name: (str) name of variable
        :data_units: (str) units of variable

        """
        self.time = time
        self.time_name = time_name
        self.time_units = time_units
        self.z = z
        self.z_name = z_name
        self.z_units = z_units
        self.data = data
        self.data_name = data_name
        self.data_units = data_units
        try:
            self.data_mean = np.mean(data, axis=0)
        except TypeError:
            self.data_mean = None

    def plot(self, axis=None, xlim=None, ylim=None,
                   xlabel=None, ylabel=None, title=None,
                   ptype='contourf', **kwargs):
        """Plot the Hovmoller diagram (time - z)

        :axis: (matplotlib.axes, optional) axis to plot figure on
        :xlim: ([float, float], optional) upper and lower limits of the x-axis
        :ylim: ([float, float], optional) upper and lower limits of the y-axis
        :xlabel: (str, optional) x-label, 'Time' by default, 'off' to turn it off
        :ylabel: (str, optional) y-label, 'Depth (m)' by default, 'off' to turn it off
        :title: (str, optional) title
        :ptype: (str, optional) plot type, valid values: contourf (default), pcolor
        :**kwargs: (keyword arguments) to be passed to matplotlib.pyplot.contourf() or
                                       matplotlib.pyplot.pcolor() depending on ptype
        :returns: (matplotlib figure object) figure

        """
        # use curret axis if not specified
        if axis is None:
            axis = plt.gca()
        # plot type
        if ptype == 'contourf':
            fig = axis.contourf(self.time, self.z, np.transpose(self.data), **kwargs)
        elif ptype == 'pcolor':
            fig = axis.pcolor(self.time, self.z, np.transpose(self.data), **kwargs)
        else:
            raise ValueError('Plot type (ptype) should be \'contourf\' or \'pcolor\', got {}.'.format(ptype))
        # x- and y-label, turn off by passing in 'off'
        if xlabel is None:
            axis.set_xlabel(self.time_name+' ('+self.time_units+')')
        else:
            if xlabel != 'off':
                axis.set_xlabel(xlabel)
        if ylabel is None:
            axis.set_ylabel(self.z_name+' ('+self.z_units+')')
        else:
            if ylabel != 'off':
                axis.set_ylabel(ylabel)
        # x- and y-limits
        if xlim is not None:
            axis.set_xlim(xlim)
        if ylim is not None:
            axis.set_ylim(ylim)
        # return figure
        return fig

    def plot_mean(self, axis=None, norm=1.0, znorm=1.0, xlim=None, ylim=None,
                        xlabel=None, ylabel=None, title=None, **kwargs):
        """Plot the mean profile

        :axis: (matplotlib.axes, optional) axis to plot figure on
        :norm: (float) normalizing factor
        :znorm: (float) normalizing factor for vertical coordinate
        :xlim: ([float, float], optional) upper and lower limits of the x-axis
        :ylim: ([float, float], optional) upper and lower limits of the y-axis
        :xlabel: (str, optional) x-label, 'self.name' by default, 'off' to turn it off
        :ylabel: (str, optional) y-label, 'Depth (m)' by default, 'off' to turn it off
        :title: (str, optional) title
        :**kwargs: (keyword arguments) to be passed to matplotlib.pyplot.plot()
        :returns: (matplotlib figure object) figure

        """
        # use curret axis if not specified
        if axis is None:
            axis = plt.gca()
        # plot figure
        fig = axis.plot(self.data_mean*norm, self.z*znorm, **kwargs)
        # x- and y-label, turn off by passing in 'off'
        if xlabel is None:
            axis.set_xlabel(self.data_name+' ('+self.data_units+')')
        else:
            if xlabel != 'off':
                axis.set_xlabel(xlabel)
        if ylabel is None:
            axis.set_ylabel(self.z_name+' ('+self.z_units+')')
        else:
            if ylabel != 'off':
                axis.set_ylabel(ylabel)
        # x- and y-limits
        if xlim is not None:
            axis.set_xlim(xlim)
        if ylim is not None:
            axis.set_ylim(ylim)
        # return figure
        return fig


#--------------------------------
# LESTimeseries
#--------------------------------

class LESTimeseries(object):

    """LESTimeseries object. """

    def __init__(self, time=None, time_name='Time', time_units='s',
                       data=None, data_name=None, data_units=None):
        """Initialize LESTimeseries.

        :time: (1D numpy array/datetime object) time
        :time_name: (str, optional) name of time
        :time_units: (str, optional) units of time
        :data: (1D numpy array) data at each location
        :data_name: (str) name of variable
        :data_units: (str) units of variable

        """

        self.time = time
        self.time_name = time_name
        self.time_units = time_units
        self.data = data
        self.data_name = data_name
        self.data_units = data_units
        try:
            self.data_mean = np.mean(data, axis=0)
        except TypeError:
            self.data_mean = None

    def plot(self, axis=None, xlim=None, ylim=None,
                   xlabel=None, ylabel=None, title=None, **kwargs):
        """Plot timeseries

        :axis: (matplotlib.axes, optional) axis to plot figure on
        :xlim: ([float, float], optional) upper and lower limits of the x-axis
        :ylim: ([float, float], optional) upper and lower limits of the y-axis
        :xlabel: (str, optional) x-label, 'Time' by default, 'off' to turn it off
        :ylabel: (str, optional) y-label, 'self.name' by default, 'off' to turn it off
        :title: (str, optional) title
        :**kwargs: (keyword arguments) to be passed to matplotlib.pyplot.plot()
        :returns: (matplotlib figure object) figure

        """
        # use curret axis if not specified
        if axis is None:
            axis = plt.gca()
        # plot figure
        fig = axis.plot(self.time, self.data, **kwargs)
        # x- and y-label, turn off by passing in 'off'
        if xlabel is None:
            axis.set_xlabel(self.time_name+' ('+self.time_units+')')
        else:
            if xlabel != 'off':
                axis.set_xlabel(xlabel)
        if ylabel is None:
            axis.set_ylabel(self.data_name+' ('+self.data_units+')')
        else:
            if ylabel != 'off':
                axis.set_ylabel(ylabel)
        # x- and y-limits
        if xlim is not None:
            axis.set_xlim(xlim)
        if ylim is not None:
            axis.set_ylim(ylim)
        # return figure
        return fig


#--------------------------------
# LESSlice
#--------------------------------

class LESSlice(object):

    """LESSlice object. """

    def __init__(self, xx=None, xx_name=None, xx_units='m',
                       yy=None, yy_name=None, yy_units='m',
                       data=None, data_name=None, data_units=None):
        """Initialize LESSlice.

        :xx: (1D numpy array) x-coordinate
        :xx_name: (str) name of x-coordinate
        :xx_units: (str) units of x-coordinate
        :yy: (1D numpy array) y-coordinate
        :yy_name: (str) name of y-coordinate
        :yy_units: (str) units of y-coordinate
        :data: (2D numpy array) data at each xx and yy
        :data_name: (str) name of variable
        :data_units: (str) units of variable

        """
        self.xx = xx
        self.xx_name = xx_name
        self.xx_units = xx_units
        self.yy = yy
        self.yy_name = yy_name
        self.yy_units = yy_units
        self.data = data
        self.data_name = data_name
        self.data_units = data_units

    def plot(self, axis=None, ptype='contourf', clim=[None, None], **kwargs):
        """Plot slice.

        :axis: (matplotlib.axes, optional) axis to plot figure on
        :ptype: (str) plot type, contourf, pcolor
        :clim: ([float, float]) limits of colorbar
        :**kwargs: keyword arguments for plt.contourf() or plt.pcolor()
        :returns: (matplotlib figure object) figure

        """
        # use curret axis if not specified
        if axis is None:
            axis = plt.gca()
        if ptype == 'pcolor':
            fig = axis.pcolor(self.xx, self.yy, self.data, vmin=clim[0], vmax=clim[1], **kwargs)
        else:
            fig = axis.contourf(self.xx, self.yy, self.data, vmin=clim[0], vmax=clim[1], **kwargs)

        axis.set_xlabel(self.xx_name+' ('+self.xx_units+')')
        axis.set_ylabel(self.yy_name+' ('+self.yy_units+')')

        if self.xx_name == 'x' and self.yy_name == 'y':
            axis.set_aspect('equal', 'box')

        cb = plt.colorbar(fig)
        cb.set_clim(clim[0], clim[1])
        cb.formatter.set_powerlimits((-2, 2))
        cb.update_ticks()
        return fig


#--------------------------------
# LESnetCDF
#--------------------------------

class LESnetCDF(object):

    """LES output netCDF data. """

    def __init__(self, path):
        """Initialize data. """
        self.path = path
        self.dataset = Dataset(self.path, 'r')
        self.list_dimensions = list(self.dataset.dimensions.keys())
        list_variables = list(self.dataset.variables.keys())
        self.list_variables = [x for x in list_variables if x not in self.list_dimensions]


#--------------------------------
# PALMData1DPR
#--------------------------------

class PALMData1DPR(LESnetCDF):

    """PALM 1D Profile data. """

    # update self.list_variables
    def __init__(self, path):
        LESnetCDF.__init__(self, path)
        list_variables_ignore = ['NORM_wpt0', 'NORM_ws2', 'NORM_tsw2', 'NORM_ws3', 'NORM_ws2tsw', 'NORM_wstsw2', 'NORM_z_i']
        self.list_variables = [x for x in self.list_variables if x not in list_variables_ignore]

    def read_profile(self, varname, tidx_start=None, tidx_end=None):
        """ Read profile for the variable [varname].

        :varname: (str) variable name
        :tidx_start: (int) starting index
        :tidx_end: (int) ending index
        :returns: (LESProfile object) profile of variable

        """
        list_fluxes = ['wu', 'wv', 'wpt', 'wvpt', 'wq', 'wqv', 'ws', 'wsa']
        time = self.dataset.variables['time'][tidx_start:tidx_end]
        if varname in self.list_variables:
            var = self.dataset.variables[varname][tidx_start:tidx_end, :]
            varunits = self.dataset.variables[varname].units
            zname = 'z'+varname
            zcoord = self.dataset.variables[zname][:]
            zunits = self.dataset.variables[zname].units
        elif varname in list_fluxes:
            varname_sgs = varname[0]+'"'+varname[1:]+'"'
            varname_res = varname[0]+'*'+varname[1:]+'*'
            if varname_res in self.list_variables and varname_sgs in self.list_variables:
                var = self.dataset.variables[varname_res][tidx_start:tidx_end, :] \
                    + self.dataset.variables[varname_sgs][tidx_start:tidx_end, :]
                varunits = self.dataset.variables[varname_res].units
                zname = 'z'+varname_res
                zcoord = self.dataset.variables[zname][:]
                zunits = self.dataset.variables[zname].units
            else:
                raise ValueError('Variable {} or {} not found.'.format(varname_res, varname_sgs))
        else:
            raise ValueError('Variable \'{}\' not found. Available variables: {}'.format(varname, self.list_variables))

        out = LESProfile(time=time, z=zcoord, z_name=zname, z_units=zunits,
                          data=var, data_name=varname, data_units=varunits)
        return out


#--------------------------------
# PALMData1DTS
#--------------------------------

class PALMData1DTS(LESnetCDF):

    """PALM 1d Time series data. """

    def read_timeseries(self, varname, axis=None, tidx_start=None, tidx_end=None):
        """Read timeseries for variable [varname].

        :varname: (str) variable name
        :tidx_start: (int) starting index
        :tidx_end: (int) ending index
        :returns: (LESTimeseries object) time series of variable

        """
        if varname in self.list_variables:
            var = self.dataset.variables[varname][tidx_start:tidx_end]
            time = self.dataset.variables['time'][tidx_start:tidx_end]
            varunits = self.dataset.variables[varname].units
        else:
            raise ValueError('Variable \'{}\' not found. Available variables: {}'.format(varname, self.list_variables))

        out = LESTimeseries(time=time, data=var, data_name=varname, data_units=varunits )
        return out


#--------------------------------
# PALMData3D
#--------------------------------

class PALMData3D(LESnetCDF):

    """PALM 3D data. """

    def read_slice_xy(self, varname, tidx=-1, zidx=1):
        """Read x-y slice for variable [varname]

        :varname: (str) variable name
        :tidx: (int) time index
        :zidx: (int) z index
        :return: (LESSlice object) xy-slice of variable

        """
        if varname in self.list_variables:
            var = self.dataset.variables[varname][tidx, zidx, :, :]
            varunits = self.dataset.variables[varname].units
        else:
            raise ValueError('Variable \'{}\' not found. Available variables: {}'.format(varname, self.list_variables))

        if varname == 'u':
            xx = self.dataset.variables['xu'][:]
        else:
            xx = self.dataset.variables['x'][:]

        if varname == 'v':
            yy = self.dataset.variables['yv'][:]
        else:
            yy = self.dataset.variables['y'][:]

        out = LESSlice(xx=xx, xx_name='x', yy=yy, yy_name='y', data=var, data_name=varname, data_units=varunits)
        return out

    def read_slice_xz(self, varname, tidx=-1, yidx=1):
        """Read x-z slice for variable [varname]

        :varname: (str) variable name
        :tidx: (int) time index
        :yidx: (int) y index
        :return: (LESSlice object) xz-slice of variable

        """
        if varname in self.list_variables:
            var = self.dataset.variables[varname][tidx, :, yidx, :]
            varunits = self.dataset.variables[varname].units
        else:
            raise ValueError('Variable \'{}\' not found. Available variables: {}'.format(varname, self.list_variables))

        if varname == 'u':
            xx = self.dataset.variables['xu'][:]
        else:
            xx = self.dataset.variables['x'][:]

        if varname == 'w':
            zz = self.dataset.variables['zw_3d'][:]
        else:
            zz = self.dataset.variables['zu_3d'][:]

        out = LESSlice(xx=xx, xx_name='x', yy=zz, yy_name='z', data=var, data_name=varname, data_units=varunits)
        return out

    def read_slice_yz(self, varname, tidx=-1, xidx=1):
        """Read y-z slice for variable [varname]

        :varname: (str) variable name
        :tidx: (int) time index
        :xidx: (int) x index
        :return: (LESSlice object) yz-slice of variable

        """
        if varname in self.list_variables:
            var = self.dataset.variables[varname][tidx, :, :, xidx]
            varunits = self.dataset.variables[varname].units
        else:
            raise ValueError('Variable \'{}\' not found. Available variables: {}'.format(varname, self.list_variables))

        if varname == 'v':
            yy = self.dataset.variables['yv'][:]
        else:
            yy = self.dataset.variables['y'][:]

        if varname == 'w':
            zz = self.dataset.variables['zw_3d'][:]
        else:
            zz = self.dataset.variables['zu_3d'][:]

        out = LESSlice(xx=yy, xx_name='y', yy=zz, yy_name='z', data=var, data_name=varname, data_units=varunits)
        return out


#--------------------------------
# NCARLESData1DPR
#--------------------------------

class NCARLESData1DPR(LESnetCDF):

    """NCARLES 1D Profile data. """

    def _read_variable(self, varname, tidx_start=None, tidx_end=None, iscl=0):
        """ Read variable [varname].

        :varname: (str) variable name
        :tidx_start: (int) starting index
        :tidx_end: (int) ending index
        :returns: (LESProfile object) profile of variable

        """

        list_scalar = ['wt', 'ut', 'vt', 'tps', 'txym', 'utle', 'vtle', 'wtle', 'utsb', 'vtsb', 'wtsb']
        if varname in list_scalar:
            var = self.dataset.variables[varname][tidx_start:tidx_end, iscl, :]
        else:
            var = self.dataset.variables[varname][tidx_start:tidx_end, :]
        return var

    def read_profile(self, varname, tidx_start=None, tidx_end=None):
        """ Read profile for the variable [varname].

        :varname: (str) variable name
        :tidx_start: (int) starting index
        :tidx_end: (int) ending index
        :returns: (LESProfile object) profile of variable

        """
        list_fluxes = ['uw', 'vw', 'wt', 'ut', 'vt']
        list_zuvar = ['engsbz', 'ups', 'vps', 'tps', 'uxym', 'vxym', 'txym',
                      'tcube', 'uvle', 'udpdx', 'vdpdx', 'wdpdx', 'udpdy',
                      'vdpdy', 'wdpdy', 'ttau11', 'ttau12',
                      'ttau13', 'ttau22', 'ttau23', 'ttau33',
                      'dsle11', 'dsle12', 'dsle13', 'dsle22',
                      'dsle23', 'dsle33', 'utle', 'vtle', 'ut', 'vt',
                      'stokes']
        list_zwvar = ['wtle', 'wtsb', 'englez', 'engz', 't_dsle',
                      't_rprod', 't_sprod', 't_stokes', 't_tran',
                      't_wq', 't_wp', 'uwle', 'uwsb', 'vwle', 'vwsb',
                      'wps', 'wxym', 'wcube', 'wfour', 'shrz',
                      'dudz', 'dvdz', 't_diss', 'uuwle', 'uvwle',
                      'uwwle', 'vvwle', 'vwwle', 'udpdz','vdpdz', 'wdpdz',
                      'uw', 'vw', 'wt']
        time = self.dataset.variables['time'][tidx_start:tidx_end]
        if varname in self.list_variables:
            var = self._read_variable(varname, tidx_start=tidx_start, tidx_end=tidx_end)
        elif varname in list_fluxes:
            varname_sgs = varname+'sb'
            varname_res = varname+'le'
            if varname_res in self.list_variables and varname_sgs in self.list_variables:
                var_sgs = self._read_variable(varname_sgs, tidx_start=tidx_start, tidx_end=tidx_end)
                var_res = self._read_variable(varname_res, tidx_start=tidx_start, tidx_end=tidx_end)
                var = var_sgs + var_res
            else:
                raise ValueError('Variable \'{}\' or \'{}\' not found.'.format(varname_res, varname_sgs))
        else:
            raise ValueError('Variable \'{}\' not found. Available variables: {}'.format(varname, self.list_variables))

        if varname in list_zwvar:
            zname = 'z_w'
        elif varname in list_zuvar:
            zname = 'z_u'

        varunits = None
        zcoord = self.dataset.variables[zname][:]
        zunits = 'm'
        out = LESProfile(time=time, z=zcoord, z_name=zname, z_units=zunits,
                          data=var, data_name=varname, data_units=varunits)
        return out


