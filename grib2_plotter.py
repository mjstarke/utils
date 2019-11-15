import pygrib
import cartopy
import cartopy.crs as ccrs
from cartopy.io import shapereader
import matplotlib.pyplot as plt
import numpy as np


def add_states(axis, projection):
    """
    Adds state outlines to the axis.
    :param axis: The axis.
    :param projection: The projection of the axis.
    :return: None.
    """
    shape = shapereader.natural_earth('110m', 'cultural', 'admin_1_states_provinces_lakes_shp')
    geometries = shapereader.Reader(shape).geometries()
    axis.add_geometries(geometries, projection, facecolor='none', edgecolor='black')


print("Opening GRB...")
grbs = pygrib.open("hrrr.t20z.wrfprsf00.grib2")
lats, lons = grbs.read(1)[0].latlons()
dpt850, _, _ = grbs(shortName="dpt", level=850)[0].data()
t850, _, _ = grbs(shortName="t", level=850)[0].data()
grbs.close()

fig = plt.figure()
ax = fig.add_axes([0, 0, 1, 1], projection=ccrs.PlateCarree())
ax.set_extent([-90, -82, 30, 36])
ax.contourf(lons, lats, t850 - dpt850, transform=ccrs.PlateCarree(), levels=np.arange(0, 10, 0.2))
ax.coastlines()
add_states(ax, ccrs.PlateCarree())

plt.show()
