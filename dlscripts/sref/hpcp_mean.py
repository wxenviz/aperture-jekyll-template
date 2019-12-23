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

totalprecips = np.load('/gpfs_backup/stormtrack/jtradfor/ensemble_data/rawdata/sref/%s_apcp.npy' % (forecasthoursub))
totalPrecips_mask = np.load('/gpfs_backup/stormtrack/jtradfor/ensemble_data/reference/sref_arw_mask.npy')

forecasthoursubback = forecasthoursub[:9] + str(int(forecasthoursub[9:11]) - 3).zfill(2) + '00'
totalprecipsback = np.load('/gpfs_backup/stormtrack/jtradfor/ensemble_data/rawdata/sref/%s_apcp.npy' % (forecasthoursubback))

totalprecips_copy = []
for totalprecip in totalprecips:
	totalprecip[totalPrecips_mask] = np.nan
	totalprecip[totalprecip>1000000] = np.nan
	totalprecips_copy.append(totalprecip)

totalprecipsback_copy = []
for totalprecip in totalprecipsback:
	totalprecip[totalPrecips_mask] = np.nan
	totalprecip[totalprecip>1000000] = np.nan
	totalprecipsback_copy.append(totalprecip)

hprecips = []
for i in range(0,9):
	hprecips.append(totalprecips_copy[i] - totalprecipsback_copy[i])

hprecip_mean = np.mean(hprecips,axis=0)
hprecip_mean[hprecip_mean<0.25] = np.nan
hprecip_mean[hprecip_mean>1.0] = 1.0

im = m.imshow(hprecip_mean,zorder=2,interpolation='none',cmap='jet',vmin=0,vmax=1.0)
cbar = plt.colorbar(im,fraction=0.023,ticks=[0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1.0])
cbar.ax.yaxis.set_tick_params(color='w')
cbar.ax.set_yticklabels([0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1.0],color='w')
plt.box(False)
hpcpfil = '/gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/uploads/outimages/sref/%s_hpcp_mean.png' % (forecasthoursub)
plt.savefig(hpcpfil,facecolor='#101010',bbox_inches='tight',dpi=500)
plt.close()



