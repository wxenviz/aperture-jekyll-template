import ftplib
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

snowtotals = np.load('/gpfs_backup/stormtrack/jtradfor/ensemble_data/rawdata/sref/%s_asnow.npy' % (forecasthoursub))
totalsnow_mask = np.load('/gpfs_backup/stormtrack/jtradfor/ensemble_data/reference/sref_arw_mask.npy')
snowtotals_copy = []

for snowtotal in snowtotals:
	snowtotal[totalsnow_mask] = np.nan
	snowtotal[snowtotal>1000000] = np.nan
	snowtotals_copy.append(snowtotal)

snowtotals_mean = np.mean(snowtotals_copy,axis=0)

snowtotals_mean[snowtotals_mean<0.25] = np.nan
snowtotals_mean[snowtotals_mean>18] = 18.0

im = m.imshow(snowtotals_mean,zorder=2,aspect='equal',interpolation='none',cmap='Blues',vmin=0.0,vmax=18.0)
cbar = plt.colorbar(im,fraction=0.023,ticks=[0,2,4,6,8,10,12,14,16,18])
cbar.ax.yaxis.set_tick_params(color='w')
cbar.ax.set_yticklabels([0,2,4,6,8,10,12,14,16,18],color='w')
plt.box(False)
meanfil = '/gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/uploads/outimages/sref/%s_aweasd_mean.png' % (forecasthoursub)
plt.savefig(meanfil,facecolor='#101010',bbox_inches='tight',dpi=500)
plt.close()

