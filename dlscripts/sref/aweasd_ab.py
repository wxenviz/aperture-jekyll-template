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

totalsnow = np.load('/gpfs_backup/stormtrack/jtradfor/ensemble_data/rawdata/sref/%s_asnow.npy' % (forecasthoursub))
totalsnow_mask = np.load('/gpfs_backup/stormtrack/jtradfor/ensemble_data/reference/sref_arw_mask.npy')

totalsnow_copy = []
for snow in totalsnow:
	snow[totalsnow_mask] = np.nan
	snow[snow>1000000] = np.nan
	totalsnow_copy.append(snow)

totalsnow_ab = np.zeros_like(totalsnow_copy[0])

for i in range(0,len(totalsnow_ab)):
	for j in range(0,len(totalsnow_ab[0])):
   		c = 0
		for k in range(0,4):
			for l in range(0,4):
				if i%4==k and j%4==l:
					totalsnow_ab[i,j] = totalsnow_copy[c][i,j]
				c+=1
	
totalsnow_ab[totalsnow_ab<0.25] = np.nan
totalsnow_ab[totalsnow_ab>18] = 18.0

im = m.imshow(totalsnow_ab,zorder=2,aspect='equal',interpolation='none',cmap='Blues',vmin=0.0,vmax=18.0)
cbar = plt.colorbar(im,fraction=0.023,ticks=[0,2,4,6,8,10,12,14,16,18])
cbar.ax.yaxis.set_tick_params(color='w')
cbar.ax.set_yticklabels([0,2,4,6,8,10,12,14,16,18],color='w')
plt.box(False)
sdfil = '/gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/uploads/outimages/sref/%s_aweasd_ab.png' % (forecasthoursub)
plt.savefig(sdfil,facecolor='#101010',bbox_inches='tight',dpi=800)
plt.close()


