"""Useful plotting functions"""
from __future__ import division

import logging
import functools

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.patches as mpatches


from matplotlib.ticker import NullFormatter

import numpy as np
import netCDF4
import salem
import geopandas as gpd
import shapely.geometry as shpg


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
def plot_catchment_width(gdir, ax=None, salemmap=None, corrected=False, GlaThiDa_profiles=False, sparse=True):
    """Plots the catchment widths out of a glacier directory.

    """

    with netCDF4.Dataset(gdir.get_filepath('gridded_data')) as nc:
        topo = nc.variables['topo'][:]

    salemmap.set_topography(topo)

    # TODO: center grid or corner grid???
    crs = gdir.grid.center_grid
    for i in gdir.divide_ids:
        geom = gdir.read_pickle('geometries', div_id=i)

        # Plot boundaries
        poly_pix = geom['polygon_pix']
        salemmap.set_geometry(poly_pix, crs=crs, fc='none', zorder=2,
                             linewidth=.2)
        for l in poly_pix.interiors:
            salemmap.set_geometry(l, crs=crs, color='black', linewidth=0.5)

        # plot Centerlines
        cls = gdir.read_pickle('inversion_flowlines', div_id=i)[::-1]
        color = gpd.plotting.gencolor(len(cls) + 1, colormap='Set1')
        for l, c in zip(cls, color):
            salemmap.set_geometry(l.line, crs=crs, color=c,
                                 linewidth=2.5, zorder=50)
            if corrected:
                for wi, cur, (n1, n2) in zip(l.widths, l.line.coords,
                                             l.normals):
                    l = shpg.LineString([shpg.Point(cur + wi / 2. * n1),
                                         shpg.Point(cur + wi / 2. * n2)])

                    salemmap.set_geometry(l, crs=crs, color=c,
                                         linewidth=0.6, zorder=50)
            else:
                for wl, wi in zip(l.geometrical_widths, l.widths):
                    col = c if np.isfinite(wi) else 'grey'
                    for w in wl:
                        salemmap.set_geometry(w, crs=crs, color=col,
                                             linewidth=0.6, zorder=50)

        if GlaThiDa_profiles:

            gtd = gdir.read_pickle('GlaThiDa')
            masks = gtd.profile_masks

            colors = iter(cm.rainbow(np.linspace(0, 1, masks.shape[0])))

            for i in range(masks.shape[0] + 1):
                if i < masks.shape[0]:
                    lab = 'n: {} \nN = {}'.format(str(i), np.sum(masks[i, :]))
                    ax.scatter(gtd.i[masks[i, :]], gtd.j[masks[i, :]],
                                color=next(colors), label=lab)
                if (i == masks.shape[0]) & (not sparse):
                    mask = np.sum(masks, axis=0)
                    mask = mask.astype(bool)
                    lab = 'undef. \nN = {}'.format(np.sum(~mask))
                    ax.scatter(gtd.j[~mask], gtd.i[~mask],
                                marker='x', color='black', label=lab)

                ax.legend(bbox_to_anchor=(1., 1.), loc=2)

    salemmap.plot(ax)

    return {}

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
        print('Bias for this run is:', np.nansum(gtd.DELTA))

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
    # here = np.where(np.abs(gtd.biases) == np.min(np.abs(gtd.biases)))
    plt.plot(gtd.best_A, gtd.best_bias, 'xr', ms=10)
    # plt.suptitle(title, loc='left')
    title = title + '\nBest bias = {} [m] when A = cfg,A * {} '.format(
        str(np.round(gtd.best_bias[0], decimals=1)), str(gtd.cfg_A_multiplier[0]))
    title = title + '[s$^{-1}$ Pa$^{-3}$]'
    plt.title(title, loc='left')
    plt.xlabel('Glens Creep Parameter [s$^{-1}$ Pa$^{-3}$]')
    plt.ylabel('Sum of all points on glacier, Bias [m]')

    return

# @_plot_map
def plot_bed_cross_sections(gdir):
    from oggm import GlaThiDa

    gtd = gdir.read_pickle('GlaThiDa')

    how = 'per_interpolation'
    grids_file = gdir.get_filepath('gridded_data', div_id=0)
    with netCDF4.Dataset(grids_file) as nc:
        vn = 'thickness'
        if how is not None:
            vn += '_' + how
        thick = nc.variables[vn][:]
        topo = nc.variables['topo'][:]
        mask = nc.variables['glacier_mask'][:]

    thick = np.where(mask, thick, np.NaN)

    
    # Plot the bed shape of the profiles
    fig = plt.figure(figsize=(7, 10))
    # ax1 = fig.add_subplot(2,1,1)
    # Cosmetics
    ymax = np.max(gtd.GTD_THICKNESS)
    if (ymax < np.max(thick)):
        ymax = np.max(thick)
    a_ymax = ymax + 0.2 * ymax  # Abs. difference of axis
    # find the data and xmax
    X = np.zeros(gtd.profile_masks.shape)
    xmax = 0
    delta_e = 0
    for i in range(0, gtd.profile_masks.shape[0]):
        lat = gtd.POINT_LAT[gtd.profile_masks[i, :]]
        lon = gtd.POINT_LON[gtd.profile_masks[i, :]]
        x = np.zeros(lat.shape)
        dis = 0
        for j in range(1, lat.shape[0]):
            dis = dis + GlaThiDa.Haversine_Distance(
                lon.iloc[j - 1], lat.iloc[j - 1],
                lon.iloc[j], lat.iloc[j])
            x[j] = dis
        x = x - np.mean(x)
        if xmax < x[-1]:
            xmax = x[-1]
        X[i, gtd.profile_masks[i, :]] = x

        topo_ = topo[gtd.j[gtd.profile_masks[i, :]], gtd.i[gtd.profile_masks[i, :]]]
        maxele = np.max(topo_)  # max elevation of profile
        minele = np.min(topo_)
        delta = maxele - minele
        if delta_e < delta:
            delta_e = delta

    a_ymax = a_ymax + delta_e

    for i in range(0, gtd.profile_masks.shape[0]):

        y3 = topo[gtd.j[gtd.profile_masks[i, :]], gtd.i[gtd.profile_masks[i, :]]]
        y = gtd.GTD_THICKNESS[gtd.profile_masks[i, :]]
        y = y3 - y
        y2 = thick[gtd.j[gtd.profile_masks[i, :]], gtd.i[gtd.profile_masks[i, :]]]
        y2 = y3 - y2
        x = X[i, gtd.profile_masks[i, :]]
        ax = plt.subplot(gtd.profile_masks.shape[0], 1, i + 1)
        l1 = ax.plot(x, y, 'x-')
        l2 = ax.plot(x, y2, '--x')
        l3 = ax.plot(x, y3, '-k', linewidth=2.0)

        # ax.legend()
        ax.set_xlim([-1.05 * xmax, 1.05 * xmax])
        if i < gtd.profile_masks.shape[0] - 1:
            ax.tick_params(
                axis='x',  # changes apply to the x-axis
                which='both',  # both major and minor ticks are affected
                labelbottom='off')  # labels along the bottom edge are off
        if i == 0:
            ax.set_title('Glacier cross sections \n Y axes: m a.s.l. and are not equal', loc='left')

        ymax = np.max(y3)
        ymin = (ymax - 1.1 * a_ymax)
        ymax = (ymax + 0.1 * a_ymax)

        #ax.set_ylim([ymin, ymax])
        ax.locator_params(axis='y', nbins=4)

        ylabel = 'n: {}'.format(str(i)) # is [m] necessary?
        ax.set_ylabel(ylabel)

    surface = mpatches.Patch(color='black', label=['Glacier Surface'])
    oggm = mpatches.Patch(color='green', label=['OGGM Inversion'])
    GlThDa = mpatches.Patch(color='blue')
    # plt.legend(handles=[red_patch])
    fig.legend(handles=[surface, oggm, GlThDa], labels=['Surface', 'OGGM', 'GlaThiDa'], loc=3,
               mode='expand', ncol=3)
               # bbox_to_anchor=(0., 1.02, 1., .102), loc=2,
               # ncol=3, mode="expand", borderaxespad=0.)
    ax.set_xlabel('[km]')

    return

def plot_profiles(gdir):

    fig = plt.figure(figsize=(7, 5))

    gtd = gdir.read_pickle('GlaThiDa')
    sparse = False
    try:
        test = gtd.full_masks
        sparse = True
    except:
        pass

    masks = gtd.profile_masks

    colors = iter(cm.rainbow(np.linspace(0, 1, masks.shape[0])))

    for i in range(masks.shape[0] + 1):
        if i < masks.shape[0]:
            lab = 'n: {}, N = {}'.format(str(i), np.sum(masks[i, :]))
            plt.scatter(gtd.POINT_LON[masks[i, :]], gtd.POINT_LAT[masks[i, :]],
                        color=next(colors), label=lab)
        if (i == masks.shape[0]) & (not sparse):
            mask = np.sum(masks, axis=0)
            mask = mask.astype(bool)
            lab = 'und. N = {}'.format(np.sum(~mask))
            plt.scatter(gtd.POINT_LON[~mask], gtd.POINT_LAT[~mask],
                        marker='x', color='black', label=lab)

        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        plt.axis('equal')

    return

def plot_points_viridis(gdir):

    fig = plt.figure(figsize=(5, 5))

    gtd = gdir.read_pickle('GlaThiDa')

    colors = iter(cm.viridis(np.linspace(0, 1, len(gtd.POINT_LON))))
    for i in range(0, len(gtd.POINT_LON)):
        plt.scatter(gtd.POINT_LON.iloc[i], gtd.POINT_LAT.iloc[i], color=next(colors))

    plt.axis('equal')

    return

