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

totalPrecips = np.load('/gpfs_backup/stormtrack/jtradfor/ensemble_data/rawdata/sref/%s_apcp.npy' % (forecasthoursub))
totalPrecips_mask = np.load('/gpfs_backup/stormtrack/jtradfor/ensemble_data/reference/sref_arw_mask.npy')

forecasthoursubback = forecasthoursub[:9] + str(int(forecasthoursub[9:11]) - 3).zfill(2) + '00'
totalPrecipsback = np.load('/gpfs_backup/stormtrack/jtradfor/ensemble_data/rawdata/sref/%s_apcp.npy' % (forecasthoursubback))

totalPrecips_copy = []
for totalPrecip in totalPrecips:
	totalPrecip[totalPrecips_mask] = np.nan
	totalPrecip[totalPrecip>1000000] = np.nan
	totalPrecips_copy.append(totalPrecip)

totalPrecipsback_copy = []
for totalPrecip in totalPrecipsback:
	totalPrecip[totalPrecips_mask] = np.nan
	totalPrecip[totalPrecip>1000000] = np.nan
	totalPrecipsback_copy.append(totalPrecip)

totalPrecips_ab = np.zeros_like(totalPrecips_copy[0])

for i in range(0,len(totalPrecips_ab)):
	for j in range(0,len(totalPrecips_ab[0])):
    		c = 0
    		for k in range(0,4):
			for l in range(0,4):
	    			if i%4==k and j%4==l:
					totalPrecips_ab[i,j] = totalPrecips_copy[c][i,j] - totalPrecipsback_copy[c][i,j]
	    			c+=1
	
totalPrecips_ab[totalPrecips_ab<0.25] = np.nan
totalPrecips_ab[totalPrecips_ab>1.0] = 1.0

im = m.imshow(totalPrecips_ab,zorder=2,aspect='equal',interpolation='none',cmap='jet',vmin=0.0,vmax=1.0)
cbar = plt.colorbar(im,fraction=0.023,ticks=[0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1.0])
cbar.ax.yaxis.set_tick_params(color='w')
cbar.ax.set_yticklabels([0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1.0],color='w')
plt.box(False)
sdfil = '/gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/uploads/outimages/sref/%s_hpcp_ab.png' % (forecasthoursub)
plt.savefig(sdfil,facecolor='#101010',bbox_inches='tight',dpi=800)
plt.close()


