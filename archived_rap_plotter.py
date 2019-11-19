from datetime import datetime
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.gridspec as gridspec
from metpy.plots import SkewT, Hodograph
from metpy.units import pandas_dataframe_to_unit_arrays, units
import metpy.calc as mpcalc
import metpy.interpolate as mpint
import numpy as np
import pint
from typing import Optional
# Packages involved in downloading:
import shutil
from urllib.request import urlopen
from contextlib import closing
from os.path import exists


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
            raw_level = data[i:i + len(headers)]
            level = dict()
            for header, datum in zip(headers, raw_level):
                level[header] = float(datum)
            levels.append(level)

        self.levels = levels

        params = dict()
        for i in range(0, len(station), 3):
            param, _, value = station[i:i + 3]
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

        return np.cos(drc) * mag * units.knot, np.sin(drc) * mag * units.knot

    @property
    def u(self):
        return self.wind_components()[0]

    @property
    def v(self):
        return self.wind_components()[1]

    @property
    def p(self):
        return self["PRES"] * units.hPa

    @property
    def T(self):
        return self["TMPC"] * units.degC

    @property
    def Td(self):
        return self["DWPC"] * units.degC

    @property
    def Tw(self):
        return self["TMWC"] * units.degC

    @property
    def thetaE(self):
        return self["THTE"] * units.kelvin

    @property
    def z(self):
        return self["HGHT"] * units.meter

    @property
    def omega(self):
        return self["OMEG"] * units("microbar/second")

    def parcel_trace(self, index_from):
        return mpcalc.parcel_profile(self.p[index_from:],
                                     self.T[index_from],
                                     self.Td[index_from])

    def cape_cin(self, index_from):
        return mpcalc.cape_cin(self.p[index_from:],
                               self.T[index_from:],
                               self.Td[index_from:],
                               self.parcel_trace(index_from))

    def lcl(self, index):
        return mpcalc.lcl(self.p[index],
                          self.T[index],
                          self.Td[index])

    def bunkers_storm_motion(self):
        return mpcalc.bunkers_storm_motion(self.p,
                                           self.u,
                                           self.v,
                                           self.z)

    def bulk_shear(self, depth=6 * units.kilometer):
        return mpcalc.bulk_shear(self.p,
                                 self.u,
                                 self.v,
                                 self.z,
                                 depth=depth)

    def storm_relative_helicity(self):
        sm_u, sm_v = self.bunkers_storm_motion()[2]
        return mpcalc.storm_relative_helicity(self.u,
                                              self.v,
                                              self.z,
                                              1000 * units.meter,
                                              0 * units.meter,
                                              sm_u,
                                              sm_v)

    def significant_tornado(self):
        u, v = self.bulk_shear()
        return mpcalc.significant_tornado(self.cape_cin(0)[0],
                                          mpint.log_interpolate_1d(
                                              self.lcl(0)[0],
                                              self.p,
                                              self.z
                                          ),
                                          self.storm_relative_helicity()[2],
                                          (u**2 + v**2)**0.5)


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


def plot_skewt(snd: Sounding, save_to: Optional[str] = None, p_top: int = 100):
    """
    Plots a skew-T from the given Sounding.
    :param snd: The Sounding.
    :param save_to: Where to save the figure.  Default None, which does not save the figure, and instead shows it.
    :param p_top: Pressure at the top of the skew-T.  If you change this, Metpy might change the rotation of the
    isotherms.  No known fix yet.
    :return: None.
    """
    ####################################################################################################################
    # Data extraction and masking

    # Extract data from sounding.
    p = snd.p
    T = snd.T
    Td = snd.Td
    Tw = snd.Tw
    Te = snd.thetaE
    z = snd.z
    cf = snd["CFRL"]
    omega = snd.omega
    u, v = snd.wind_components()

    e = mpcalc.saturation_vapor_pressure(T)
    rv = mpcalc.mixing_ratio(e, p)
    w = mpcalc.vertical_velocity(omega, p, T, rv).to("cm/s")

    # Create masks to filter what data is plotted.
    mask_dewpoint = Td > -9000. * units.degC  # Plot only non-missing dewpoints.
    mask_wetbulb = Tw > -9000. * units.degC  # Plot only non-missing dewpoints.
    mask_thetae = Te > 0. * units.K  # Plot only non-missing theta-es.
    mask_barbs = p > p_top * units.hPa  # Plot only winds below the top of the sounding.

    ####################################################################################################################
    # Define intervals of height for coloring barbs and hodograph.
    z_interval_levels = [1500, 3000, 6000, 9000, 12000, 99999]
    z_interval_colors = ["red", "orange", "green", "blue", "purple", "grey"]

    z_colors = []

    for item in z:
        for color, interval in zip(z_interval_colors, z_interval_levels):
            if item <= interval*units.meter:
                z_colors.append(color)
                break

    ####################################################################################################################
    # Plotting skew-T

    fig = plt.figure(figsize=(11, 11))
    ax_hodo = fig.add_axes([0.70, 0.675, 0.25, 0.25])
    ax_thte = fig.add_axes([0.70, 0.375, 0.25, 0.25])
    skew = SkewT(fig, rotation=45, rect=[0.05, 0.05, 0.60, 0.9])

    # Plot temperature, dewpoint, and wet-bulb.
    skew.plot(p, T, 'r')
    skew.plot(p[mask_dewpoint], Td[mask_dewpoint], 'g')
    skew.plot(p[mask_wetbulb], Tw[mask_wetbulb], color='#009999', linewidth=1)

    # Calculate and plot surface parcel trace.
    sfc_trace = snd.parcel_trace(0).to('degC')
    sfc_trace_plot = skew.plot(p, sfc_trace, c='orange', linewidth=2, zorder=-10)

    # Calculate and plot MU parcel trace.
    mu_level_index = np.argmax(Te[p > 750.*units.hPa])
    mu_trace = snd.parcel_trace(mu_level_index).to('degC')
    mu_trace_plot = skew.plot(p[mu_level_index:], mu_trace, c='gray', linewidth=2, zorder=-9)

    # Plot each barb individually for control over color.  Unfortunately, the c arg of plot_barbs doesn't work for this
    # purpose.
    for p_, u_, v_, c_ in zip(p[mask_barbs][::BARB_DENSITY],
                              u[mask_barbs][::BARB_DENSITY],
                              v[mask_barbs][::BARB_DENSITY],
                              np.array(z_colors)[mask_barbs][::BARB_DENSITY]):
        skew.plot_barbs(p_, u_, v_, y_clip_radius=0.03, barbcolor=c_)

    ####################################################################################################################
    # Cloud fraction and omega
    zero_line = 1/15
    cf_plot = (cf * zero_line) / 100
    w_plot = (w.magnitude / 20) + zero_line
    skew.ax.plot(np.zeros(cf_plot.shape) + 1/15, snd.p, transform=skew.ax.get_yaxis_transform(), color="grey")
    skew.ax.plot(cf_plot, snd.p, transform=skew.ax.get_yaxis_transform(), color="black")
    skew.ax.plot(w_plot, snd.p, transform=skew.ax.get_yaxis_transform(), color="purple")

    skew.ax.text(np.max(w_plot), snd.p[np.argmax(w_plot)], " {:.1f}".format(np.max(w.magnitude)),
                 color="purple", ha="left", va="center", transform=skew.ax.get_yaxis_transform())
    skew.ax.text(max(np.min(w_plot), 0), snd.p[np.argmin(w_plot)], " {:.1f}".format(np.min(w.magnitude)),
                 color="purple", ha="left", va="center", transform=skew.ax.get_yaxis_transform())
    # skew.ax.fill_betweenx(snd.p, cloud_fractions, np.zeros(cloud_fractions.shape))

    ####################################################################################################################
    # Tweaking

    skew.ax.set_xlim(-30, 40)
    skew.ax.set_ylim(1020, p_top)
    skew.ax.set_xlabel("")
    skew.ax.set_ylabel("")

    # Add legend for the parcel traces.
    skew.ax.legend(handles=[
        mlines.Line2D([], [], color='orange', label='Surface parcel'),
        mlines.Line2D([], [], color='gray', label=r"Max $\theta_e$ below 750mb"),
        mlines.Line2D([], [], color='black', label=r"Cloud fraction"),
        mlines.Line2D([], [], color='purple', label=r"Vertical velocity (cm/s)"),
    ], loc="upper center")

    # Add adiabats and isohumes.
    skew.plot_dry_adiabats(t0=np.arange(233, 533, 10) * units.K,
                           alpha=0.25, color='orangered')
    skew.plot_moist_adiabats(t0=np.arange(233, 400, 5) * units.K,
                             alpha=0.25, color='tab:green')
    # Reshape required as a quirk of metpy.
    skew.plot_mixing_lines(w=np.array([1, 2, 3, 4, 6, 8, 10, 12, 16, 20, 24, 28, 36]).reshape(-1, 1) / 1000.,
                           p=np.arange(1000, 99, -100) * units.hPa,
                           linestyle='dotted', color='tab:blue')

    plt.title('RAP sounding at {}'.format(snd.params["STID"]), loc='left')
    plt.title('{:.0f}-hour forecast valid at {}'.format(snd.params["STIM"], snd.params["TIME"]), loc='right')

    ####################################################################################################################
    # Theta-E plot

    # Set up axis for theta-e plot.
    ax_thte.plot(Te[mask_thetae], p[mask_thetae])

    ax_thte.set_xlim(300, 360)
    ax_thte.set_ylim(1020, p_top)
    ax_thte.set_yscale("log")
    ax_thte.set_yticks(np.arange(p_top, 1001, 100))
    ax_thte.set_yticklabels(np.arange(p_top, 1001, 100))
    ax_thte.set_xlabel("")
    ax_thte.grid(axis="both")
    plt.text(0.5, 0.9, "Theta-E (Kelvins)", ha="center", va="center", transform=ax_thte.transAxes)

    ####################################################################################################################
    # Hodograph

    # Set up axis for hodograph.
    h = Hodograph(ax_hodo, component_range=100)
    h.add_grid(20)

    # Plot each segment individually for control over color, reversed so that the full hodograph is plotted first,
    # followed by all but the last segment, etc.  Unfortunately, the plot_colormapped() function doesn't work for this
    # purpose.
    for color, interval in zip(reversed(z_interval_colors), reversed(z_interval_levels)):
        mask = z < interval*units.meter
        h.plot(u.magnitude[mask], v.magnitude[mask], c=color)

    for vector in snd.bunkers_storm_motion():
        h.plot(vector[0], vector[1], c="black", markersize=3, marker="o")

    ax_hodo.set_xticks([])
    ax_hodo.set_yticks([])
    ax_hodo.set_xlim(-60, 100)
    ax_hodo.set_ylim(-60, 100)
    plt.text(0.1, 0.9, "Velocity (knots)", ha="left", va="center", transform=ax_hodo.transAxes)
    for a in range(20, 61, 20):
        ax_hodo.text(-a * 0.71, -a * 0.71, a, ha="center", va="center")

########################################################################################################################
    parameter_names = [
        "SB STP",
        "0-1 SRH",
        "SB CAPE",
        "SB CIN"
    ]
    parameters = [
        snd.significant_tornado()[0],
        snd.storm_relative_helicity()[2],
        snd.cape_cin(0)[0],
        snd.cape_cin(0)[1]
    ]

    for name, value, i in zip(parameter_names, parameters, range(len(parameters))):
        s = "{:15} {:10.3f}".format(name, value.magnitude)
        fig.text(0.70, 0.32 - (0.02*i), s, ha="left", va="top", family="monospace", transform=fig.transFigure)

########################################################################################################################
    if save_to is None:
        plt.show()
    else:
        plt.savefig(save_to)
        plt.close()
