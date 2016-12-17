"""Useful plotting functions"""
from __future__ import division

import logging
import functools

import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter

import numpy as np
import netCDF4
import salem

from oggm.utils import entity_task

# Local imports
import oggm.cfg as cfg

# Module logger
# from sandbox.run_alaska import gdir

log = logging.getLogger(__name__)

def _plot_map(plotfunc):
    """
    Decorator for common salem.Map plotting logic
    """
    commondoc = """
    Parameters
    ----------
    ax : matplotlib axes object, optional
        If None, uses own axis
    add_colorbar : Boolean, optional, default=True
        Adds colorbar to axis
    horizontal_colorbar : Boolean, optional, default=False
        Horizontal colorbar instead
    title : str, optional
        If left to None, the plot decides wether it writes a title or not. Set
        to '' for no title.
    title_comment : str, optional
        add something to the default title. Set to none to remove default
    lonlat_contours_kwargs: dict, optional
        pass kwargs to salem.Map.set_lonlat_contours
    """

    # Build on the original docstring
    plotfunc.__doc__ = '\n'.join((plotfunc.__doc__, commondoc))

    @functools.wraps(plotfunc)
    def newplotfunc(gdir, ax=None, add_colorbar=True, title=None,
                    title_comment=None, horizontal_colorbar=False,
                    lonlat_contours_kwargs=None,
                    **kwargs):

        dofig = False
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            dofig = True

        mp = salem.Map(gdir.grid, countries=False, nx=gdir.grid.nx)
        if lonlat_contours_kwargs is not None:
            mp.set_lonlat_contours(**lonlat_contours_kwargs)

        out = plotfunc(gdir, ax=ax, salemmap=mp, **kwargs)

        if add_colorbar and 'cbar_label' in out:
            cbprim = out.get('cbar_primitive', mp)
            if horizontal_colorbar:
                cb = cbprim.append_colorbar(ax, "bottom", size="5%", pad=0.4)
            else:
                cb = cbprim.append_colorbar(ax, "right", size="5%", pad=0.2)
            cb.set_label(out['cbar_label'])

        if title is None:
            if 'title' not in out:
                # Make a defaut one
                title = gdir.rgi_id
                if gdir.name is not None and gdir.name != '':
                    title += ': ' + gdir.name
                out['title'] = title

            if title_comment is None:
                title_comment = out.get('title_comment', '')

            out['title'] += title_comment
            ax.set_title(out['title'])
        else:
            ax.set_title(title)

        if dofig:
            plt.tight_layout()

    return newplotfunc

@entity_task(log)
@_plot_map
def plot_distributed_thickness(gdir, ax=None, salemmap=None, how=None, GTD=False, Delta_GTD=False):
    """Plots the result of the inversion out of a glacier directory.

    Method: 'alt' or 'interp', GTD = GlaTHiDa points, Delta_GTD is differences of oggm and GTD
    """

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

    salemmap.set_topography(topo)

    # TODO: center grid or corner grid???
    crs = gdir.grid.center_grid


    for i in gdir.divide_ids:
        geom = gdir.read_pickle('geometries', div_id=i)

        # Plot boundaries
        poly_pix = geom['polygon_pix']
        salemmap.set_geometry(poly_pix, crs=crs, fc='none', zorder=2, linewidth=.2)
        for l in poly_pix.interiors:
            salemmap.set_geometry(l, crs=crs, color='black', linewidth=0.5)

    salemmap.set_cmap(plt.get_cmap('viridis'))
    salemmap.set_plot_params(nlevels=256)
    salemmap.set_data(thick)


    # GlaThiDa

    if GTD:
        from salem import DataLevels
        gtd = gdir.read_pickle('GlaThiDa')
        x, y = salemmap.grid.transform(gtd.POINT_LON.values, gtd.POINT_LAT.values)
        dl = DataLevels(
            gtd.GTD_THICKNESS,
            nlevels=256,
            cmap=plt.get_cmap('viridis'),
            vmin=salemmap.vmin,
            vmax=salemmap.vmax,
        )
        ax.scatter(x, y, s=30, color=dl.to_rgb(), edgecolors='k', linewidths=1)

    if Delta_GTD:
        from salem import DataLevels
        gtd = gdir.read_pickle('GlaThiDa')
        x, y = salemmap.grid.transform(gtd.POINT_LON.values, gtd.POINT_LAT.values)
        dl = DataLevels(
            gtd.DELTA,
            extend='both',
            cmap=plt.get_cmap('coolwarm'),
            levels=np.arange(-100, 101, 20))

        ax.scatter(x, y, s=30, color=dl.to_rgb(), edgecolors='k', linewidths=1)
        dl.append_colorbar(ax, label=' ', position='right', pad=1.5)

    salemmap.plot(ax)

    return dict(cbar_label=' ')

def plot_As_vs_Volume(gdir, title=None):
    """Plots the different volumes of the glacier as a function of Glens A's"""
    gtd = gdir.read_pickle('GlaThiDa')

    if title is None:
        title = gdir.rgi_id
        if gdir.name is not None and gdir.name != '':
            title += ': ' + gdir.name


    plt.plot(gtd.glens_As, gtd.volumes/(10**9), '.b')
    plt.title(title, loc='left')
    plt.xlabel('Glens Creep Parameter [s$^{-1}$ Pa$^{-3}$]')
    plt.ylabel('Volume [km$^3$]')

    return

def plot_As_vs_bias(gdir, title=None):
    """Plots the different volumes of the glacier as a function of Glens A's"""
    gtd = gdir.read_pickle('GlaThiDa')

    if title is None:
        title = gdir.rgi_id
        if gdir.name is not None and gdir.name != '':
            title += ': ' + gdir.name


    plt.plot(gtd.glens_As, gtd.biases, '.b')
    plt.plot([0, np.max(gtd.glens_As)], [0, 0], '--k')
    #plt.ylim([-2000, 2000])
    plt.title(title, loc='left')
    plt.xlabel('Glens Creep Parameter [s$^{-1}$ Pa$^{-3}$]')
    plt.ylabel('Sum of all points on glacier, Bias [m]')

    return