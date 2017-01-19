import pandas as pd
import numpy as np
import netCDF4

import oggm.cfg as cfg
from oggm.core.preprocessing import inversion

class GlaThiDa:
    """A class to store the GlaThiDa data and whith methods to compare GlaThiDa with OGGM output"""

    def __init__(self):
        self.POINT_LON = None
        self.POINT_LAT = None
        self.THICKNESS = None
        self.RGI_Reg = None
        self.GlaThiDa_ID = None

    def read_pickle(self, path, RGI_Reg, GlaThiDa_ID):
        """Reads the pickle file from the RGI GlaThiDa linkage"""
        name = path + '/' + RGI_Reg + '/' + GlaThiDa_ID + '/' + GlaThiDa_ID + '.p'
        df = pd.read_pickle(name)
        self.POINT_LON = df.POINT_LON
        self.POINT_LAT = df.POINT_LAT
        self.GTD_THICKNESS = df.THICKNESS
        self.RGI_Reg = RGI_Reg
        self.GlaThiDa_ID = GlaThiDa_ID
        #self.RGI_ID =
        return self

    def transform(self, gdir):
        """Transforms the lat long points to a local grid via salem, maybe."""
        crs = gdir.grid.center_grid
        self.i, self.j = crs.transform(self.POINT_LON.values, self.POINT_LAT.values, nearest=True)
        return self

    def Delta_Thickness(self, gdir):
        """Findes the difference of the oggm model output and the GlaThiDa points"""
        self.OGGM_THICK = np.zeros(self.GTD_THICKNESS.shape)
        self.DELTA = np.zeros(self.GTD_THICKNESS.shape)

        # Stuff from Fabi:
        how = 'per_interpolation'
        with netCDF4.Dataset(gdir.get_filepath('gridded_data')) as nc:
            topo = nc.variables['topo'][:]
            mask = nc.variables['glacier_mask'][:]

        grids_file = gdir.get_filepath('gridded_data', div_id=0)
        with netCDF4.Dataset(grids_file) as nc:
            vn = 'thickness'
            if how is not None:
                vn += '_' + how
            thick = nc.variables[vn][:]

        thick = np.where(mask, thick, np.NaN)



        for k in range(self.j.shape[0]):
            self.OGGM_THICK[k] = thick[self.j[k], self.i[k]]
            self.DELTA[k] = thick[self.j[k], self.i[k]] - self.GTD_THICKNESS.iloc[k]

        return self

    def volume_and_bias(self, gdir, start=0.9, stop=1.11, step=0.05):
        """Runs loops to be able to later plot different volumes and biases as a funciton of A \n
        the three kwargs refer to numpy arrange, that is what values is cfg.A multiplied with"""

        glen_as = cfg.A*np.arange(start, stop, step)
        volumes = np.zeros(glen_as.shape)
        thick_diffs = np.zeros([self.GTD_THICKNESS.shape[0], glen_as.shape[0]])
        for i in range(glen_as.shape[0]):
            volumes[i], _ = inversion.invert_parabolic_bed(
                gdir, fs=0., glen_a=glen_as[i], write=True)

            inversion.distribute_thickness(gdir, how='per_interpolation',
                                           add_slope=False,
                                           add_nc_name=True)

            self = self.Delta_Thickness(gdir=gdir)
            thick_diffs[:, i] = self.DELTA

        self.glens_As = glen_as
        self.volumes = volumes
        self.biases = np.nansum(thick_diffs, axis=0)

        return self

#    def optimize_bias(self, gdir):
#        from scipy import optimize
#        optimize.newton()
#        a = np.max(self.biases[self.biases < 0])

