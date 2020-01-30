from metpy.calc import wind_components
import numpy as np
from typing import Optional
from matplotlib.colors import Normalize
import matplotlib.cm as mplcm


def wind(ax, x, wdir, wspd, wgst, *, vector_length: Optional[float] = 1.0, calm_thresh: float = 1.0,
         cm=mplcm.hsv, cm_dir_shift: float = 180., vector_multiplier: float = 1.0,
         wspd_kwargs: Optional[dict] = None, wgst_kwargs: Optional[dict] = None, quiver_kwargs: Optional[dict] = None,
         calms_kwargs: Optional[dict] = None):
    """
    Plots a timeseries of wind vectors.
    :param ax:  The axis on which to plot.
    :param x:  The x-values at which to plot.
    :param wdir:  The array of wind directions in degrees true.  Must have same length as x.
    :param wspd:  The array of wind speeds.  Unit-agnostic.  Must have same length as x.
    :param wgst:  The array of wind gusts.  Should have same unit as wspd.  Must have same length as x.
    :param vector_length:  If a float, the length of every vector (in data units).  If None, vector_multiplier is used
    instead.  Default 1.0.
    :param vector_multiplier:  A multiplier applied to the length of every vector.  Ignored if vector_length is a float.
    Default 1.0.
    :param calm_thresh:  The threshold below which an observation is considered calm.  Calm observations will not have
    vectors plotted, even if the speed is nonzero.  Default 1.0 - 1 m/s is roughly the calm threshold for ASOS (3 kt,
    more precisely).
    :param cm:  The colormap to use for coloring the vectors.  Must be a proper MPL colormap and not, e.g., a list.
    Default matplotlib.cm.hsv (the "hsv" colormap).
    :param cm_dir_shift:  The angle, in degrees, through which the colormap will be turned.  If 0, the first color on
    the map will correspond to southerly winds.  Default 180, which completely flips the colormap.
    :param wspd_kwargs:  kwargs to be passed to the plotter for the wind speed line.  Default None, which uses hardcoded
    defaults.
    :param wgst_kwargs:  kwargs to be passed to the plotter for the wind gust line.  Default None, which uses hardcoded
    defaults.
    :param quiver_kwargs:  kwargs to be passed to the plotter for the wind vectors.  Default None, which uses hardcoded
    defaults.
    :param calms_kwargs:  kwargs to be passed to the plotter for the calm wind dots.  Default None, which uses hardcoded
    defaults.
    :return:  None.
    """
    # Flip wind direction for colormapping.
    wdir_flipped = wdir + cm_dir_shift
    wdir_flipped[wdir_flipped > 360] -= 360

    # Get components of the wind.  Normalize all vectors to have unit magnitude.
    if vector_length is None:
        u1, v1 = wind_components(wspd * vector_multiplier, wdir)
    else:
        u1, v1 = wind_components(np.array([vector_length] * len(wdir)), wdir)
    # Create masks for all points where the wind is calm and where it is not calm.
    calms = wspd < calm_thresh
    not_calms = wspd >= calm_thresh

    norm = Normalize(0, 360)

    if wspd_kwargs is None:
        wspd_kwargs = dict(lw=0.5, color="grey")
    if wgst_kwargs is None:
        wgst_kwargs = dict(lw=0.75, color="black", markersize=4, marker=".")
    if quiver_kwargs is None:
        quiver_kwargs = dict(width=0.0025)
    if calms_kwargs is None:
        calms_kwargs = dict(color="grey")

    ax.plot(x, wspd.to("knots"), **wspd_kwargs)
    ax.plot(x, wgst.to("knots"), **wgst_kwargs)
    ax.quiver(x[not_calms], wspd.to("knots")[not_calms], u1[not_calms], v1[not_calms],
              color=cm(norm(wdir_flipped[not_calms])), **quiver_kwargs)
    ax.scatter(x[calms], [0] * len(x[calms]), **calms_kwargs)

    ax.set_ylim(0, ax.get_ylim()[1] + 2)
