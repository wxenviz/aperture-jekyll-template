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

apcps = np.load('/gpfs_backup/stormtrack/jtradfor/ensemble_data/rawdata/sref/%s_apcp.npy' % (forecasthoursub))
totalPrecips_mask = np.load('/gpfs_backup/stormtrack/jtradfor/ensemble_data/reference/sref_arw_mask.npy')

apcp_copy = []
for apcp in apcps:
	apcp[totalPrecips_mask] = np.nan
	apcp_copy.append(apcp)

fig, axes = plt.subplots(nrows=4,ncols=4,figsize=(16,9))
plt.subplots_adjust(wspace=0.05,hspace=-.10)

for i in range(0,16):

	ax = axes.flat[i]
	m = Basemap(projection='lcc',lat_0=5,lon_0=-100,llcrnrlon=-126,llcrnrlat=23,urcrnrlon=-63,urcrnrlat=50,resolution='l',ax=ax)
	shp_info = m.readshapefile('/gpfs_backup/stormtrack/jtradfor/ensemble_data/reference/st99_d00','states',drawbounds=False)

	for nshape,seg in enumerate(m.states):
		poly = Polygon(seg,facecolor='white',edgecolor='white',zorder=1,linewidth=.5)
		poly2 = Polygon(seg,facecolor='none',edgecolor='black',zorder=3,linewidth=.5)
		ax.add_patch(poly)
		ax.add_patch(poly2)

	apcp_copy[i][apcp_copy[i]<0.25] = np.nan
	apcp_copy[i][apcp_copy[i]>1000000] = np.nan
	apcp_copy[i][apcp_copy[i]>20.] = 20.

	im = m.imshow(apcp_copy[i],zorder=2,aspect='equal',interpolation='none',cmap='jet',vmin=0.0,vmax=20.0)

cbar = fig.colorbar(im, ax=axes.ravel().tolist(),fraction=0.025,ticks=[0,2,4,6,8,10,12,14,16,18,20])
cbar.ax.yaxis.set_tick_params(color='w')
cbar.ax.set_yticklabels([0,2,4,6,8,10,12,14,16,18,20],color='w')
postagefil = '/gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/uploads/outimages/sref/%s_apcp_postage.png' % (forecasthoursub)
plt.savefig(postagefil,facecolor='#101010',bbox_inches='tight',dpi=500)
plt.close()

