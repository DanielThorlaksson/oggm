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
from oggm.core.preprocessing import  gis, geometry, climate, inversion, centerlines
import oggm.cfg as cfg
from oggm.utils import get_demo_file
from oggm.tests import requires_py3

# Globals
current_dir = os.path.dirname(os.path.abspath(__file__))

class TestKesselWFInvert(unittest.TestCase):
    print('Number 1!')

    def setUp(self):

        # test directory
        self.testdir = os.path.join(current_dir, 'tmp_kesselwf')
        self.rm_dir() # Uncomment when doing a dry run

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

        # Init:

        GlaThiDa_ID = str(1970)
        GlaThiDa_Path = '/home/daniel/Dropbox/dev/data/ttt_2_rgi'
        RGI_Region  = str(11)
        RGI_ID = 'RGI50-11.00787'
        RGI_Region_shp = '/home/daniel/Dropbox/dev/data/rgi50/11_rgi50_CentralEurope.shp'

        # rgi_file = get_demo_file('rgi_oetztal.shp')
        # entity = gpd.GeoDataFrame.from_file(rgi_file)
        # entity = entity.loc[entity.RGIID == 'RGI40-11.00787'].iloc[0]

        # # Load the RGI50 file:
        entity = gpd.read_file(RGI_Region_shp)
        entity = entity.loc[entity.RGIId == RGI_ID].iloc[0]

        gdir = oggm.GlacierDirectory(entity, base_dir=self.testdir)
        # All the follwing line are only needed for dryruns
        gis.define_glacier_region(gdir, entity=entity)
        gis.glacier_masks(gdir)
        centerlines.compute_centerlines(gdir)
        geometry.initialize_flowlines(gdir)
        geometry.catchment_area(gdir)
        geometry.catchment_width_geom(gdir)
        geometry.catchment_width_correction(gdir)
        climate.distribute_climate_data([gdir])
        climate.mu_candidates(gdir, div_id=0)
        hef_file = get_demo_file('mbdata_RGI40-11.00787.csv')
        mbdf = pd.read_csv(hef_file).set_index('YEAR')
        t_star, bias = climate.t_star_from_refmb(gdir, mbdf['ANNUAL_BALANCE'])
        t_star = t_star[-1]
        bias = bias[-1]
        climate.local_mustar_apparent_mb(gdir, tstar=t_star, bias=bias)

        # OK.
        # GlaThiDa data manipulation:
        # Read the oggm data:

        # Read the GlaThiDa data:
        # TODO: This needs to be some sort of funciton

        from oggm import GlaThiDa

        gtd = GlaThiDa.GlaThiDa()
        gtd = gtd.read_pickle(path=GlaThiDa_Path,
                              RGI_Reg=RGI_Region, GlaThiDa_ID=GlaThiDa_ID)
        gtd = gtd.transform(gdir=gdir)

        inversion.prepare_for_inversion(gdir)

        # ======= The necessary steps to plot a single glacier! =========
        v, _ = inversion.invert_parabolic_bed(
            gdir, fs=0., glen_a=cfg.A, write=True)

        inversion.distribute_thickness(gdir, how='per_interpolation',
                                       add_slope=False,
                                       add_nc_name=True)

        gtd = gtd.Delta_Thickness(gdir=gdir)
        # Here it ends

        # This line runs the inversions as a loop.
        # gtd.volume_and_bias(gdir=gdir)

        gdir.write_pickle(gtd, 'GlaThiDa')

        # oggm plotting routine:
        from oggm import plot_kesselwf
        import matplotlib.pyplot as plt

        plot_kesselwf.plot_distributed_thickness(gdir, how='per_interpolation', GTD=True)

        #plt.show()
        plt.savefig('/home/daniel/tempfigs/fig_same.png')

        plot_kesselwf.plot_distributed_thickness(gdir, how='per_interpolation', Delta_GTD=True,
                                            title_comment='\n Left: OGGM inversion thickness [m] \n Right: OGGM inversion - GlaThiDa[m]')
        #plt.show()
        plt.savefig('/home/daniel/tempfigs/fig_Delta.png')


class TestKesselWFInversions(unittest.TestCase):

    def setUp(self):
        # test directory
        self.testdir = os.path.join(current_dir, 'tmp_kesselwf')
        # self.rm_dir() # Uncomment when doing a dry run

        # Init
        cfg.initialize()
        cfg.PARAMS['temp_use_local_gradient'] = False
        cfg.PARAMS['use_multiple_flowlines'] = True
        cfg.PATHS['dem_file'] = '/home/daniel/Dropbox/dev/data/Climate_and_DEM/srtm_90m_v4_alps.tif'
        cfg.PATHS['climate_file'] = '/home/daniel/Dropbox/dev/data/Climate_and_DEM/histalp_merged_hydro_yrs.nc'
        cfg.PATHS['cru_dir'] = '~'

    def tearDown(self):
        pass
        # self.rm_dir()

    def rm_dir(self):
        if os.path.exists(self.testdir):
            shutil.rmtree(self.testdir)

    @requires_py3
    def test_invert(self):
        print('number two!')
        # rgi_file = get_demo_file('rgi_oetztal.shp')
        # entity = gpd.GeoDataFrame.from_file(rgi_file)
        # entity = entity.loc[entity.RGIID == 'RGI40-11.00787'].iloc[0]
        dryrun = False
        gtd_loops = False
        best_bias = 5.57
        # Init: === Preparing to make this into a function which takes these varibles to plot any GlaThiDa glacier

        GlaThiDa_ID = str(497)
        GlaThiDa_Path = '/home/daniel/Dropbox/dev/data/ttt_2_rgi'
        RGI_Region  = str(11)
        RGI_ID = 'RGI50-11.02789'
        RGI_Region_shp = '/home/daniel/Dropbox/dev/data/rgi50/11_rgi50_CentralEurope.shp'

        path = '/home/daniel/tempfigs/'
        path = path + 'FINDELEN/'

        # # Load the RGI50 file:
        entity = gpd.read_file(RGI_Region_shp)
        entity = entity.loc[entity.RGIId == RGI_ID].iloc[0]

        gdir = oggm.GlacierDirectory(entity, base_dir=self.testdir)
        # The follwing line are only needed for dryruns
        if dryrun:
            gis.define_glacier_region(gdir, entity=entity)
            gis.glacier_masks(gdir)
            centerlines.compute_centerlines(gdir)
            geometry.initialize_flowlines(gdir)
            geometry.catchment_area(gdir)
            geometry.catchment_width_geom(gdir)
            geometry.catchment_width_correction(gdir)
            climate.distribute_climate_data([gdir])
            climate.mu_candidates(gdir, div_id=0)
            hef_file = get_demo_file('mbdata_RGI40-11.00787.csv')
            mbdf = pd.read_csv(hef_file).set_index('YEAR')
            t_star, bias = climate.t_star_from_refmb(gdir, mbdf['ANNUAL_BALANCE'])
            t_star = t_star[-1]
            bias = bias[-1]
            climate.local_mustar_apparent_mb(gdir, tstar=t_star, bias=bias)

        # OK.
        # GlaThiDa data manipulation:
        # Read the oggm data:

        # Read the GlaThiDa data:
        # TODO: This needs to be some sort of funciton

        from oggm import GlaThiDa

        gtd = GlaThiDa.GlaThiDa()
        gtd = gtd.read_pickle(path=GlaThiDa_Path, RGI_Reg=RGI_Region, GlaThiDa_ID=GlaThiDa_ID)
        gtd = gtd.transform(gdir=gdir)

        inversion.prepare_for_inversion(gdir)

        # ======= The necessary steps to plot a single glacier! =========
        v, _ = inversion.invert_parabolic_bed(
            gdir, fs=0., glen_a=best_bias*cfg.A, write=True)

        inversion.distribute_thickness(gdir, how='per_interpolation',
                                       add_slope=False,
                                       add_nc_name=True,
                                       smooth=True)

        # # This line runs the inversions as a loop.
        if dryrun | gtd_loops:
            gtd.volume_and_bias(gdir=gdir, start=5.56, stop=5.58, step=0.01)

        gtd = gtd.Delta_Thickness(gdir=gdir)

        gtd.find_measurement_profiles(gdir, sort=True, d_tol=0.5, n_tol=3, max_ele_diff=200)

        gdir.write_pickle(gtd, 'GlaThiDa')



        # oggm plotting routine:
        from oggm import plot_kesselwf
        import matplotlib.pyplot as plt

        if dryrun | gtd_loops:
            plot_kesselwf.plot_As_vs_Volume(gdir)
            name = path + 'volume_vs_A.png'
            plt.savefig(name)
            plt.clf()

            plot_kesselwf.plot_As_vs_bias(gdir)
            name = path + 'bias_vs_A.png'
            plt.savefig(name)
            plt.clf()

        if not dryrun & gtd_loops:

            plot_kesselwf.plot_distributed_thickness(gdir, how='per_interpolation', GTD=True)
            name = path + 'fig_same.png'
            plt.savefig(name)
            plt.clf()

            plot_kesselwf.plot_distributed_thickness(gdir, how='per_interpolation', Delta_GTD=True,
                                                title_comment='\n Left: OGGM inversion thickness [m] \n Right: OGGM inversion - GlaThiDa[m]')
            name = path + 'fig_Delta.png'
            plt.savefig(name)
            plt.clf()

            plot_kesselwf.plot_bed_cross_sections(gdir)
            name = path + 'cross_sections.png'
            plt.savefig(name)
            plt.clf()

            plot_kesselwf.plot_profiles(gdir)
            name = path + 'profiles.png'
            plt.savefig(name, bbox_inches='tight')
            plt.clf()

            plot_kesselwf.plot_points_viridis(gdir)
            name = path + 'points.png'
            plt.savefig(name)
            plt.clf()

            plot_kesselwf.plot_catchment_width(gdir, corrected=True)
            name = path + 'catchment.png'
            plt.savefig(name)
            plt.clf()

class TestTASCHACHFERNERnvert(unittest.TestCase):
    print('Number three!')

    def setUp(self):

        # test directory
        self.testdir = os.path.join(current_dir, 'tmp_TASCHACHFERNER')
        # self.rm_dir() # Uncomment when doing a dry run

        # Init
        cfg.initialize()
        cfg.PARAMS['temp_use_local_gradient'] = False
        cfg.PARAMS['use_multiple_flowlines'] = True
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
        entity = entity.loc[entity.RGIId == 'RGI50-11.00687'].iloc[0]

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
        # GlaThiDa data manipulation:
        # Read the oggm data:

        # Read the GlaThiDa data:
        # TODO: This needs to be some sort of funciton

        from oggm import GlaThiDa

        gtd = GlaThiDa.GlaThiDa()
        gtd = gtd.read_pickle(path='/home/daniel/Dropbox/dev/data/ttt_2_rgi',
                              RGI_Reg=str(11), GlaThiDa_ID=str(1999))
        gtd = gtd.transform(gdir=gdir)

        inversion.prepare_for_inversion(gdir)

        # ======= The necessary steps to plot a single glacier! =========
        v, _ = inversion.invert_parabolic_bed(
            gdir, fs=0., glen_a=cfg.A, write=True)

        inversion.distribute_thickness(gdir, how='per_interpolation',
                                       add_slope=False,
                                       add_nc_name=True)

        gtd = gtd.Delta_Thickness(gdir=gdir)
        # Here it ends

        # This line runs the inversions as a loop.
        # gtd.volume_and_bias(gdir=gdir)

        gdir.write_pickle(gtd, 'GlaThiDa')

        # oggm plotting routine:
        from oggm import plot_kesselwf
        import matplotlib.pyplot as plt

        plot_kesselwf.plot_distributed_thickness(gdir, how='per_interpolation', GTD=True)

        #plt.show()
        plt.savefig('/home/daniel/tempfigs/fig_same2.png')

        plot_kesselwf.plot_distributed_thickness(gdir, how='per_interpolation', Delta_GTD=True,
                                            title_comment='\n Left: OGGM inversion thickness [m] \n Right: OGGM inversion - GlaThiDa[m]')
        #plt.show()
        plt.savefig('/home/daniel/tempfigs/fig_Delta2.png')

