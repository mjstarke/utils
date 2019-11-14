from datetime import datetime
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.gridspec as gridspec
from metpy.plots import SkewT, Hodograph
from metpy.units import pandas_dataframe_to_unit_arrays, units
import metpy.calc as mpcalc
import numpy as np
from typing import Optional
# Packages involved in downloading:
import shutil
from urllib.request import urlopen
from contextlib import closing
from os.path import exists


bufkit_filename = "rap_max.buf"
BARB_DENSITY = 2


class Sounding:
    def __init__(self, text):
        lines = text.split("\n")
        hdr_lines = []
        snd_lines = []
        stn_lines = []
        for line in lines:
            if len(line.strip()) > 0:
                if "=" in line:
                    stn_lines.append(line)
                else:
                    if line[0] not in "-0123456789.":
                        hdr_lines.append(line)
                    else:
                        snd_lines.append(line)

        headers = " ".join(hdr_lines).split(" ")
        data = " ".join(snd_lines).split(" ")
        station = ("STID = " + (" ".join(stn_lines))).split(" ")
        station = [item for item in station if len(item.strip()) > 0]

        levels = []
        for i in range(0, len(data), len(headers)):
            raw_level = data[i:i+len(headers)]
            level = dict()
            for header, datum in zip(headers, raw_level):
                level[header] = float(datum)
            levels.append(level)

        self.levels = levels

        params = dict()
        for i in range(0, len(station), 3):
            param, _, value = station[i:i+3]
            try:
                params[param] = float(value)
            except ValueError:
                params[param] = value

        self.params = params

    def __getitem__(self, item):
        return np.array([level[item] for level in self.levels])

    def wind_components(self):
        """
        :return: The zonal and meridional components of the wind.
        """
        mag = np.array(self["SKNT"])
        drc = np.array(self["DRCT"])
        drc = 270 - drc
        drc *= np.pi / 180.

        return np.cos(drc) * mag, np.sin(drc) * mag

    @property
    def u(self):
        return self.wind_components()[0]

    @property
    def v(self):
        return self.wind_components()[1]

    def parcel_trace(self, index_from):
        return mpcalc.parcel_profile(self["PRES"][index_from:] * units.hPa,
                                     self["TMPC"][index_from] * units.degC,
                                     self["DWPC"][index_from] * units.degC).to('degC')

    def cape_cin(self, index_from):
        return mpcalc.cape_cin(self["PRES"][index_from:] * units.hPa,
                               self["TMPC"][index_from] * units.degC,
                               self["DWPC"][index_from] * units.degC,
                               self.parcel_trace(index_from))

    def lcl(self, index):
        return mpcalc.lcl(self["PRES"][index] * units.hPa,
                          self["TMPC"][index] * units.degC,
                          self["DWPC"][index] * units.degC)

    def bunkers_storm_motion(self):
        return mpcalc.bunkers_storm_motion(self["PRES"] * units.hPa,
                                           self.u * units.knot,
                                           self.v * units.knot,
                                           self["HGHT"] * units.meter)

    def bulk_shear(self, depth=6000*units.meter):
        return mpcalc.bulk_shear(self["PRES"] * units.hPa,
                                 self.u * units.knot,
                                 self.v * units.knot,
                                 self["HGHT"] * units.meter,
                                 depth)

    def storm_relative_helicity(self):
        sm_u, sm_v = self.bunkers_storm_motion()[2]
        return mpcalc.storm_relative_helicity(self.u * units.knot,
                                              self.v * units.knot,
                                              self["HGHT"] * units.meter,
                                              1000*units.meter,
                                              0*units.meter,
                                              sm_u,
                                              sm_v)


########################################################################################################################
# Download BUFKIT file.
def download_rap_bufkit(sid: str, dt: datetime, check_existing: bool = True) -> str:
    """
    Downloads a RAP BUFKIT sounding from Iowa State's archive.
    :param sid: The ID of the sounding site.  See https://www.meteor.iastate.edu/~ckarsten/bufkit/data/
    :param dt: The datetime of the sounding.
    :param check_existing: Whether or not to check that the file exists locally before downloading.
    :return: The name of the downloaded file.
    """
    remote = "http://mtarchive.geol.iastate.edu/" \
             "{0.year}/{0.month:0>2}/{0.day:0>2}/bufkit/{0.hour:0>2}/rap/rap_{1}.buf".format(dt, sid)
    local = "rap_{1}_{0.year}-{0.month:0>2}-{0.day:0>2}-{0.hour:0>2}.buf".format(dt, sid)

    if not (check_existing and exists(local)):
        with closing(urlopen(remote)) as r:
            with open(local, 'wb') as f:
                shutil.copyfileobj(r, f)

    return local


def plot_skewt(snd: Sounding, save_to: Optional[str] = None):
    ####################################################################################################################
    # Data extraction and masking

    # Extract data from sounding.
    p = snd["PRES"]
    T = snd["TMPC"]
    Td = snd["DWPC"]
    Tw = snd["TMWC"]
    Te = snd["THTE"]
    z = snd["HGHT"]
    u, v = snd.wind_components()

    # Create masks to filter what data is plotted.
    mask_dewpoint = Td > -9000.  # Plot only non-missing dewpoints.
    mask_wetbulb = Tw > -9000.  # Plot only non-missing dewpoints.
    mask_thetae = Te > 0.  # Plot only non-missing theta-es.
    mask_barbs = p > 100.  # Plot only winds below 100mb.


    ####################################################################################################################
    # Define intervals of height for coloring barbs and hodograph.
    z_interval_levels = [1500, 3000, 6000, 9000, 12000, 99999]
    z_interval_colors = ["red", "orange", "green", "blue", "purple", "black"]

    z_colors = []

    for item in z:
        for color, interval in zip(z_interval_colors, z_interval_levels):
            if item <= interval:
                z_colors.append(color)
                break


    ####################################################################################################################
    # Plotting skew-T

    # Use gridspec to more easily locate auxilliary axes.
    fig = plt.figure(figsize=(11, 11))
    ax_hodo = fig.add_axes([0.70, 0.675, 0.25, 0.25])
    ax_thte = fig.add_axes([0.70, 0.375, 0.25, 0.25])
    skew = SkewT(fig, rotation=45, rect=[0.05, 0.05, 0.60, 0.9])

    # Plot temperature, dewpoint, and wet-bulb.
    skew.plot(p, T, 'r')
    skew.plot(p[mask_dewpoint], Td[mask_dewpoint], 'g')
    skew.plot(p[mask_wetbulb], Tw[mask_wetbulb], color='#009999', linewidth=1)

    # Calculate and plot surface parcel trace.
    sfc_trace = snd.parcel_trace(0)
    sfc_trace_plot = skew.plot(p, sfc_trace, 'k', linewidth=2, zorder=-10)

    # Calculate and plot MU parcel trace.
    mu_level_index = np.argmax(Te[p > 750.])
    mu_trace = snd.parcel_trace(mu_level_index)
    mu_trace_plot = skew.plot(p[mu_level_index:], mu_trace, c='gray', linewidth=2, zorder=-9)

    # Plot each barb individually for control over color.  Unfortunately, the c arg of plot_barbs doesn't work for this
    # purpose.
    for p_, u_, v_, c_ in zip(p[mask_barbs][::BARB_DENSITY],
                              u[mask_barbs][::BARB_DENSITY],
                              v[mask_barbs][::BARB_DENSITY],
                              np.array(z_colors)[mask_barbs][::BARB_DENSITY]):
        skew.plot_barbs(p_, u_, v_, y_clip_radius=0.03, barbcolor=c_)


    ####################################################################################################################
    # Tweaking

    # Set some appropriate axes limits for x and y
    skew.ax.set_xlim(-30, 40)
    skew.ax.set_ylim(1020, 100)

    # Add legend for the parcel traces.
    skew.ax.legend(handles=[
        mlines.Line2D([], [], color='black', label='Surface parcel'),
        mlines.Line2D([], [], color='gray', label=r"Max $\theta_e$ below 750mb")
    ], loc="upper center")

    # Add adiabats and isohumes.
    skew.plot_dry_adiabats(t0=np.arange(233, 533, 10) * units.K,
                           alpha=0.25, color='orangered')
    skew.plot_moist_adiabats(t0=np.arange(233, 400, 5) * units.K,
                             alpha=0.25, color='tab:green')
    skew.plot_mixing_lines(p=np.arange(1000, 99, -20) * units.hPa,
                           linestyle='dotted', color='tab:blue')

    plt.title('RAP sounding at {}'.format(snd.params["STID"]), loc='left')
    plt.title('{:.0f}-hour forecast valid at {}'.format(snd.params["STIM"], snd.params["TIME"]), loc='right')


    ####################################################################################################################
    # Theta-E plot

    # Set up axis for theta-e plot.
    ax_thte.plot(Te[mask_thetae], p[mask_thetae])

    ax_thte.set_xlim(300, 360)
    ax_thte.set_ylim(1020, 100)
    ax_thte.set_yscale("log")
    ax_thte.set_yticks(np.arange(100, 1001, 100))
    ax_thte.set_yticklabels(np.arange(100, 1001, 100))
    ax_thte.grid(axis="both")
    plt.text(0.5, 0.9, "Theta-E", ha="center", va="center", transform=ax_thte.transAxes)


    ####################################################################################################################
    # Hodograph

    # Set up axis for hodograph.
    h = Hodograph(ax_hodo)
    h.add_grid(20)

    # Plot each segment individually for control over color, reversed so that the full hodograph is plotted first,
    # followed by all but the last segment, etc.  Unfortunately, the plot_colormapped() function doesn't work for this
    # purpose.
    for color, interval in zip(reversed(z_interval_colors), reversed(z_interval_levels)):
        mask = z < interval
        h.plot(u[mask], v[mask], c=color)

    ax_hodo.set_xticks([])
    ax_hodo.set_yticks([])
    for a in range(20, 81, 20):
        plt.text(-a*0.71, -a*0.71, a, ha="center", va="center")

    if save_to is None:
        plt.show()
    else:
        plt.savefig(save_to)
        plt.close()


########################################################################################################################
