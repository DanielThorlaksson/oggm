# This script started its life in "test_kesselwf.py" but here it will be turned into one single script
# and not in the test mode

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
from oggm import GlaThiDa

# Globals
# current_dir = os.path.dirname(os.path.abspath(__file__))
current_dir = '/home/daniel/Dropbox/dev/OGGM_vs_GlaThiDa'

def GlaTHiDa_with_TTT():

    # read the GlaThiDa_RGI linkage from the alps
    all_points = pd.read_csv('~/Dropbox/dev/data/GlaThiDa_2016/TTT_RGI_temporal_Alps.csv')
    glaciers = all_points.drop_duplicates('GlaThiDa_ID').copy()
    print('Raw number of glaciers: ', glaciers.shape[0])
    # I am not quite sure how mulitple RGI shapes work, for now I throw them out
    # which is not good, I am loosing more than half with that, need to find a work around

    glaciers['multiple_RGI'] = False
    # BEFORE DOING THIS: YOU NEED TO THROW OUT THE POINTS WHICH ARE OUTSIDE OF GLACIERS!!!
    for glacier in glaciers.itertuples():
        tmp = all_points.loc[all_points.GlaThiDa_ID == glacier.GlaThiDa_ID]
        tmp = tmp.drop_duplicates('rgi_id')
        if tmp.shape[0] > 1:
            glaciers.loc[glacier.Index, 'multiple_RGI'] = True
    glaciers = glaciers[~glaciers.multiple_RGI].copy()
    print('N single RGI shapes:', glaciers.shape[0])

    glaciers['Delta_Time'] = glaciers.SURVEY_DATE - glaciers.rgi_survey_date
    glaciers = glaciers.loc[np.abs(glaciers.Delta_Time) < 10]
    print('N temporally filtered glaciers:', glaciers.shape[0])

    RGI_Region_shp = '/home/daniel/Dropbox/dev/data/rgi50/11_rgi50_CentralEurope.shp'
    RGI_Region = str(11)
    fig_path = '/home/daniel/Dropbox/dev/OGGM_vs_GlaThiDa/figs/'
    fig_ext = '.png'
    GlaThiDa_Path = '/home/daniel/Dropbox/dev/data/ttt_2_rgi'
    run_dir = os.path.join(current_dir, 'GlaThiDa_run')
    entities = gpd.read_file(RGI_Region_shp)

    #glaciers = glaciers.iloc[1:3, :]

    for glacier in glaciers.itertuples():

        GlaThiDa_ID = str(glacier.GlaThiDa_ID)
        RGI_ID = str(glacier.rgi_id)


        if os.path.exists(run_dir):
            shutil.rmtree(run_dir)

        # Init, I th
        cfg.initialize()
        cfg.PARAMS['temp_use_local_gradient'] = False
        cfg.PARAMS['use_multiple_flowlines'] = True
        # Since I will stay in the alps then the DEM and histalp climate files are are constant
        cfg.PATHS['dem_file'] = '/home/daniel/Dropbox/dev/data/Climate_and_DEM/srtm_90m_v4_alps.tif'
        cfg.PATHS['climate_file'] = '/home/daniel/Dropbox/dev/data/Climate_and_DEM/histalp_merged_hydro_yrs.nc'
        cfg.PATHS['cru_dir'] = '~'


        # # Load the RGI50 file:
        entity = entities.loc[entities.RGIId == RGI_ID].iloc[0]

        gdir = oggm.GlacierDirectory(entity, base_dir=run_dir)
        # The follwing line are only needed for dryruns
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

        gtd = GlaThiDa.GlaThiDa()
        gtd = gtd.read_pickle(path=GlaThiDa_Path, RGI_Reg=RGI_Region, GlaThiDa_ID=GlaThiDa_ID)
        gtd = gtd.transform(gdir=gdir)

        inversion.prepare_for_inversion(gdir)

        # ======= The necessary steps to plot a single glacier! =========

        gtd.best_bias(gdir, beta_min=0.001, beta_max=1000)

        if gtd.multiplier < 0.1:
            multiplier = 0.1
        elif gtd.multiplier > 10:
            multiplier = 10
        else:
            multiplier = gtd.multiplier

        # Make sure that the last written run is run with the best bias
        v, a = inversion.invert_parabolic_bed(
            gdir, fs=0., glen_a=multiplier*cfg.A, write=True)

        inversion.distribute_thickness(gdir, how='per_interpolation',
                                       add_slope=False,
                                       add_nc_name=True,
                                       smooth=True)

        gtd.Delta_Thickness(gdir)
        gdir.write_pickle(gtd, 'GlaThiDa')

        to_write = pd.DataFrame(index=[1], columns=['GlaThiDa_ID'])
        to_write['GlaThiDa_ID'] = gtd.GlaThiDa_ID
        to_write['RGI_ID'] = RGI_ID
        to_write['best_A'] = gtd.best_A
        to_write['multiplier'] = gtd.multiplier
        to_write['Volume'] = v
        to_write['Area'] = gdir.rgi_area_km2 * 1e6
        name = current_dir + '/GlaThiDa_OGGM_data.csv'

        if os.path.isfile(name):
            with open(name, 'a') as f:
                to_write.to_csv(f, index=False, header=False)
        else:
            to_write.to_csv(name, index=False, header=True)

        # oggm plotting routine:
        from oggm import plot_kesselwf
        import matplotlib.pyplot as plt

        if not os.path.exists(fig_path):
            os.makedirs(fig_path)

        plot_kesselwf.plot_distributed_thickness(gdir, how='per_interpolation', GTD=True)
        name = '{}{}_same{}'.format(fig_path, GlaThiDa_ID, fig_ext)
        #name = fig_path + 'fig_same.png'
        plt.savefig(name)
        plt.clf()

        plot_kesselwf.plot_distributed_thickness(gdir, how='per_interpolation', Delta_GTD=True) #,
                                                    # title_comment='\n Left: OGGM inversion thickness [m] \n Right: OGGM inversion - GlaThiDa[m]')
        name = '{}{}_Delta{}'.format(fig_path, GlaThiDa_ID, fig_ext)
        plt.savefig(name)
        plt.clf()
        #
        # Only needed when plotting cross sections
        # gtd.find_measurement_profiles(gdir, sort=True, d_tol=0.5, n_tol=5, max_ele_diff=200,
        #                               sparsify=True)
        # gdir.write_pickle(gtd, 'GlaThiDa')


        #             plot_kesselwf.plot_bed_cross_sections(gdir)
        #             name = fig_path + 'cross_sections.pdf'
        #             plt.savefig(name)
        #             plt.clf()
        #
        #             plot_kesselwf.plot_profiles(gdir)
        #             name = fig_path + 'profiles.png'
        #             plt.savefig(name, bbox_inches='tight')
        #             plt.clf()
        #
        #             # plot_kesselwf.plot_points_viridis(gdir)
        #             # name = fig_path + 'points.png'
        #             # plt.savefig(name)
        #             # plt.clf()
        #             # #
        plot_kesselwf.plot_catchment_width(gdir, corrected=True, GlaThiDa_profiles=False, sparse=False)
        name = '{}{}_catchment{}'.format(fig_path, GlaThiDa_ID, fig_ext)
        plt.savefig(name)
        plt.clf()

if __name__ == "__main__":
    GlaTHiDa_with_TTT()
