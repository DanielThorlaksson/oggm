import pandas as pd
import numpy as np
import netCDF4
from scipy import optimize as optimization

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
        """Findes the difference of the oggm model output and the GlaThiDa points
        ====> becoming obsolete, moved to 'best bias' """
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
        the three kwargs refer to numpy arrange, that is what values is cfg.A multiplied with
        ===> becoming obsolete"""

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
        self.biases = np.nanmean(thick_diffs, axis=0)

        # Save the best A of the run
        here = np.where(np.abs(self.biases) == np.min(np.abs(self.biases)))
        self.best_bias = self.biases[here]
        self.best_A = self.glens_As[here]
        self.cfg_A_multiplier = self.best_A/cfg.A
        print('cfg A Multiplier: ', self.cfg_A_multiplier)

        return self

    def best_bias(self, gdir, beta_min=0.1, beta_max=10, n_bands=10):
        """WIP: Finds the best bias with optimization, where the bias is a virtual bias with elevation weights. """

        # Read the major flowline
        fl = gdir.read_pickle('inversion_flowlines', div_id=-1)
        # sparsify to 10 bands
        # on second thought, I think this should be defined exactly as the flowlines
        # on a third thought, the centerlines-flowlines defnitions is not quite that simple ...
        indexes = np.round(np.linspace(0, len(fl[0].surface_h)-1, num=(n_bands+1))).astype(int)
        elevations = fl[0].surface_h[indexes]


        # load the topo, GlaThiDa usually has elevation also saved. However to be "locally constant" the OGGM DEM on
        # the glacier grid will be used.
        grids_file = gdir.get_filepath('gridded_data', div_id=0)
        with netCDF4.Dataset(grids_file) as nc:
            topo = nc.variables['topo'][:]

        df = pd.DataFrame()
        df['i'] = self.i
        df['j'] = self.j
        df['GTD_THICKNESS'] = self.GTD_THICKNESS.values
        df['ele'] = np.NaN
        for point in df.itertuples():
            df.loc[point.Index, 'ele'] = topo[point.j, point.i]

        # The elevations are from highest down to lowest!
        df['ele_band'] = -1
        for k in range(1, n_bands):
            df.loc[(df.ele <= elevations[k-1]) & (df.ele > elevations[k]), 'ele_band'] = k

        # Points still with -1 are outside of glacier... play no role in bias
        df = df.loc[df.ele_band != -1]
        # Is this nececcary? it will alsways be close to zero
        self.best_ele_bias = 100 # This will be updated within the optimization loop

        # Here starts the job starts
        def to_optimize(x):

            V, a = inversion.invert_parabolic_bed(
                    gdir, fs=0., glen_a=cfg.A*x, write=True)

            inversion.distribute_thickness(gdir, how='per_interpolation',
                                           add_slope=False,
                                           add_nc_name=True)

            how = 'per_interpolation'
            with netCDF4.Dataset(gdir.get_filepath('gridded_data')) as nc:
                mask = nc.variables['glacier_mask'][:]

            grids_file = gdir.get_filepath('gridded_data', div_id=0)
            with netCDF4.Dataset(grids_file) as nc:
                vn = 'thickness'
                if how is not None:
                    vn += '_' + how
                thick = nc.variables[vn][:]

            thick = np.where(mask, thick, np.NaN)

            df['OGGM_THICK'] = -1000
            for point in df.itertuples():
                df.loc[point.Index, 'OGGM_THICK'] = thick[point.j, point.i]

            df['Delta_thick'] = df.OGGM_THICK - df.GTD_THICKNESS

            runner = df.drop_duplicates('ele_band')
            bias = 0
            for ele in runner.itertuples():
                # is it worth anything to store the individual biases for each range?
                # here the bias is calculated, each "elevation band" is equally important for the total
                try:
                    weight = np.nanmean(df.loc[df.ele_band == ele.ele_band, 'Delta_thick'])
                except:
                    print('Glacier: {}, \n ele band {}, has no non NaN point'.format(RGI_ID, str(ele)))
                    weight = 0

                bias += weight

            if np.abs(self.best_ele_bias) > np.abs(bias):
                self.best_ele_bias = bias

            return bias

        try:
            opti = optimization.brenth(to_optimize, a=beta_min, b=beta_max, xtol=0.01, disp=True)
        except:
            # Wild guess here
            opti = 0.1

        self.best_A = cfg.A * opti
        self.multiplier = opti

        return self

        # === The follwing functions started their life in "inv_test.ipynb" ===
    def find_measurement_profiles(self, gdir, d_tol=0.20, n_tol=5, sort=True,
                                  sparsify=True, spars_max=40, max_ele_diff=None):
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



        # Remove 'Profiles' with more or less than a tolerence of points
        npoints = np.where((np.sum(masks, axis=1) > n_tol)) #& (np.sum(masks, axis=1) < n_max_tol))
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

        if sparsify:
            n_max = spars_max
            npoints = np.sum(masks, axis=1)
            sparse_masks = masks.copy()
            for i in range(0, masks.shape[0]):
                pro_start = np.where(masks[i, :])[0][0]
                pro_end = pro_start + npoints[i] - 1
                indexes = np.round(np.linspace(pro_start, pro_end, num=n_max))
                indexes[0] = pro_start
                indexes[-1] = pro_end
                sparse_masks[i, :] = False
                for index in indexes:
                    sparse_masks[i, int(index)] = True

                #Todo: Change to some sort of gradiant?
                if max_ele_diff is not None:
                    if np.abs(np.max(topo[self.j[masks[i, :]], self.i[masks[i, :]]]) -
                                      np.min(topo[self.j[masks[i, :]], self.i[masks[i, :]]])) > max_ele_diff:
                        sparse_masks[i, :] = False

            if max_ele_diff is not None:
                sparse_masks = sparse_masks[np.any(sparse_masks, axis=1), :]

            self.full_masks = masks
            self.profile_masks = sparse_masks


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

