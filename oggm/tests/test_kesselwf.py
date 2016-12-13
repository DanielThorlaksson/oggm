from __future__ import absolute_import, division

import oggm.utils

import unittest
import os
import shutil

import numpy as np
import pandas as pd
import geopandas as gpd
import netCDF4

# Local imports
from oggm.core.preprocessing import inversion
import oggm.cfg as cfg
from oggm.utils import get_demo_file
from oggm.tests import requires_py3

# Globals
current_dir = os.path.dirname(os.path.abspath(__file__))

class TestKesselWFInvert(unittest.TestCase):

    def setUp(self):

        # test directory
        self.testdir = os.path.join(current_dir, 'tmp_kesselwf')
        # self.rm_dir() # Uncomment when doing a dry run

        # Init
        cfg.initialize()
        cfg.PARAMS['temp_use_local_gradient'] = False
        cfg.PARAMS['use_multiple_flowlines'] = False
        cfg.PATHS['dem_file'] = get_demo_file('srtm_oetztal.tif')
        cfg.PATHS['climate_file'] = get_demo_file('HISTALP_oetztal.nc')
        cfg.PATHS['cru_dir'] = '~'

    def tearDown(self):
        pass
        # self.rm_dir()

    def rm_dir(self):
        if os.path.exists(self.testdir):
            shutil.rmtree(self.testdir)

    @requires_py3
    def test_invert(self):

        # rgi_file = get_demo_file('rgi_oetztal.shp')
        # entity = gpd.GeoDataFrame.from_file(rgi_file)
        # entity = entity.loc[entity.RGIID == 'RGI40-11.00787'].iloc[0]

        # # Load the RGI50 file:
        entity = gpd.read_file(
            '/home/daniel/Dropbox/dev/data/rgi50/11_rgi50_CentralEurope.shp')
        entity = entity.loc[entity.RGIId == 'RGI50-11.00787'].iloc[0]

        gdir = oggm.GlacierDirectory(entity, base_dir=self.testdir)
        # All the follwing line are only needed for dryruns
        # gis.define_glacier_region(gdir, entity=entity)
        # gis.glacier_masks(gdir)
        # centerlines.compute_centerlines(gdir)
        # geometry.initialize_flowlines(gdir)
        # geometry.catchment_area(gdir)
        # geometry.catchment_width_geom(gdir)
        # geometry.catchment_width_correction(gdir)
        # climate.distribute_climate_data([gdir])
        # climate.mu_candidates(gdir, div_id=0)
        # hef_file = get_demo_file('mbdata_RGI40-11.00787.csv')
        # mbdf = pd.read_csv(hef_file).set_index('YEAR')
        # t_star, bias = climate.t_star_from_refmb(gdir, mbdf['ANNUAL_BALANCE'])
        # t_star = t_star[-1]
        # bias = bias[-1]
        # climate.local_mustar_apparent_mb(gdir, tstar=t_star, bias=bias)

        # OK.
        inversion.prepare_for_inversion(gdir)
        glen_a = cfg.A
        v, _ = inversion.invert_parabolic_bed(gdir, fs=0., glen_a=glen_a, write=True)

        inversion.distribute_thickness(gdir, how='per_interpolation',
                                       add_slope=False,
                                       add_nc_name=True)

        # GlaThiDa data manipulation:
        # Read the oggm data:
        crs = gdir.grid.center_grid
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

        # Read the GlaThiDa data:
        # TODO: This needs to be some sort of funciton
        df = pd.read_pickle('/home/daniel/Dropbox/dev/data/ttt_2_rgi/11/1970/1970.p')
        gtd = pd.DataFrame()
        gtd['POINT_LON'] = df.POINT_LON
        gtd['POINT_LAT'] = df.POINT_LAT
        gtd['GTD_THICKNESS'] = df.THICKNESS
        gtd['i'], gtd['j'] = crs.transform(df.POINT_LON.values, df.POINT_LAT.values, nearest=True)
        gtd['OGGM_THICKNESS'] = np.nan
        gtd['DELTA'] = np.nan

        for p in gtd.itertuples():
            gtd.loc[p.Index, 'OGGM_THICKNESS'] = thick[p.j, p.i]
            gtd.loc[p.Index, 'DELTA'] = thick[p.j, p.i] - p.GTD_THICKNESS

        gdir.write_pickle(gtd, 'GlaThiDa')

        # oggm plotting routine:
        from oggm import plot_kesselwf
        import matplotlib.pyplot as plt

        plot_kesselwf.plot_distributed_thickness(gdir, how=how, GTD=True)

        #plt.show()
        plt.savefig('/home/daniel/tempfigs/fig_same.png')

        plot_kesselwf.plot_distributed_thickness(gdir, how='per_interpolation', Delta_GTD=True,
                                            title_comment='\n Left: OGGM inversion thickness [m] \n Right: OGGM inversion - GlaThiDa[m]')
        #plt.show()
        plt.savefig('/home/daniel/tempfigs/fig_Delta.png')
