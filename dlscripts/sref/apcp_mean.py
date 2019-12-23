import glob
import subprocess as sp
import csv
import numpy as np
import netCDF4 as nc4
import pygrib as pg
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import datetime
import scipy
import os
import sys

from mpl_toolkits.basemap import Basemap
from matplotlib.patches import Polygon
from matplotlib.colors import LinearSegmentedColormap
from scipy.spatial import Delaunay
from scipy.interpolate import LinearNDInterpolator
from shutil import copyfile

forecasthoursub = str(sys.argv[1])
plt.figure(figsize=(16,9))

m = Basemap(projection='lcc',lat_0=5,lon_0=-100,llcrnrlon=-126,llcrnrlat=23,urcrnrlon=-63,urcrnrlat=50,resolution='h')
shp_info = m.readshapefile('/gpfs_backup/stormtrack/jtradfor/ensemble_data/reference/st99_d00','states',drawbounds=False)

ax = plt.gca()

for nshape,seg in enumerate(m.states):
	poly = Polygon(seg,facecolor='white',edgecolor='white',zorder=1,linewidth=1)
	poly2 = Polygon(seg,facecolor='none',edgecolor='black',zorder=3,linewidth=1)
	ax.add_patch(poly)
	ax.add_patch(poly2)

totalPrecips = np.load('/gpfs_backup/stormtrack/jtradfor/ensemble_data/rawdata/sref/%s_apcp.npy' % (forecasthoursub))
totalPrecips_mask = np.load('/gpfs_backup/stormtrack/jtradfor/ensemble_data/reference/sref_arw_mask.npy')

totalPrecips_copy = []
for totalPrecip in totalPrecips:
	totalPrecip[totalPrecips_mask] = np.nan
	totalPrecip[totalPrecip>1000000] = np.nan
	totalPrecips_copy.append(totalPrecip)

totalPrecip_mean = np.mean(totalPrecips_copy,axis=0)
totalPrecip_mean[totalPrecip_mean<0.25] = np.nan
totalPrecip_mean[totalPrecip_mean>20] = 20.0

im = m.imshow(totalPrecip_mean,zorder=2,interpolation='none',cmap='jet',vmin=0,vmax=20)
cbar = plt.colorbar(im,fraction=0.023,ticks=[0,2,4,6,8,10,12,14,16,18,20])
cbar.ax.yaxis.set_tick_params(color='w')
cbar.ax.set_yticklabels([0,2,4,6,8,10,12,14,16,18,20],color='w')
plt.box(False)
apcpfil = '/gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/uploads/outimages/sref/%s_apcp_mean.png' % (forecasthoursub)
plt.savefig(apcpfil,facecolor='#101010',bbox_inches='tight',dpi=500)
plt.close()

