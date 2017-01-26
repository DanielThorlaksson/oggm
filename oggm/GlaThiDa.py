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
        self.profile_masks = None

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

        # === The follwing functions started their life in "inv_test.ipynb" ===
    def find_measurement_profiles(self, gdir, d_tol=0.20, n_tol=5, sort=True):
        # Here a loop will be created which will seperate the different measurement profile
        mask = np.zeros([1, np.size(self.POINT_LAT)], dtype=bool)
        masks = mask.copy()
          # Tolerance in km
        mark = 0  # marks the begining of the next profile
        for i in range(0, len(self.POINT_LAT) - 1):
            dis = Haversine_Distance(
                self.POINT_LON.iloc[i], self.POINT_LAT.iloc[i],
                self.POINT_LON.iloc[i + 1], self.POINT_LAT.iloc[i + 1])
            if (dis > d_tol):  # & (i-mark > 0):
                mask[0, (mark):(i + 1)] = True
                masks = np.append(masks, mask, axis=0)
                mask[:, :] = False
                mark = i + 1



        # Remove 'Profiles' with less than a tolerence of points
        npoints = np.where(np.sum(masks, axis=1) > n_tol)
        masks = masks[npoints[0], :]

        if sort:
            # Load the oggm interpolated thickness
            how = 'per_interpolation'
            grids_file = gdir.get_filepath('gridded_data', div_id=0)
            with netCDF4.Dataset(grids_file) as nc:
                vn = 'thickness'
                if how is not None:
                    vn += '_' + how
                topo = nc.variables['topo'][:]

            # This cell will sort the masks into elevation order, highest elevation is the first mask
            elev = np.zeros([masks.shape[0], 2])  # elevation, profile number
            for k in range(0, masks.shape[0]):
                i = int(np.round(np.mean(self.i[masks[k, :]])))  # Find the characteristc point on the local grid
                j = int(np.round(np.mean(self.j[masks[k, :]])))
                elev[k, 0] = topo[j, i]  # Save the elevation
                elev[k, 1] = k  # Save the mask number

            # Sort it!
            elev = np.flipud(elev[elev[:, 0].argsort()])
            # Rearrange the masks
            temp = np.zeros(masks.shape, dtype=bool)
            n = 0
            for i in elev[:, 1]:
                temp[n, :] = masks[int(i), :]
                n += 1

            masks = temp.copy()

        self.profile_masks = masks

        # === The follwing functions started their life in "f_ttt2rgi.py" ===
    # path = /home/daniel/Dropbox/ttt2rgi/

def Haversine_Distance(phi1, lambda1, phi2, lambda2):
        """The Function calculates the distance on a spherical earth and gies the
        results in km, input is in degrees. The equation taken from wiki site,
        Great-circle distance"""
        # First convert to radians
        phi1 = np.deg2rad(phi1)
        lambda1 = np.deg2rad(lambda1)
        phi2 = np.deg2rad(phi2)
        lambda2 = np.deg2rad(lambda2)
        R_earth = 6371.  # Mean radios of the earth in km
        d = 2 * R_earth * np.arcsin(np.sqrt(np.sin((phi2 - phi1) / 2) ** 2 +
                                            np.cos(phi1) * np.cos(phi2) *
                                            np.sin((lambda2 - lambda1) / 2) ** 2))
        return d;


    #

 #   def