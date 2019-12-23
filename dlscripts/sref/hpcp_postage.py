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
totalprecips = np.load('/gpfs_backup/stormtrack/jtradfor/ensemble_data/rawdata/sref/%s_apcp.npy' % (forecasthoursub))
forecasthoursubback = forecasthoursub[:9] + str(int(forecasthoursub[9:11]) - 3).zfill(2) + '00'
totalprecipsback = np.load('/gpfs_backup/stormtrack/jtradfor/ensemble_data/rawdata/sref/%s_apcp.npy' % (forecasthoursubback))

totalprecip_mask = np.load('/gpfs_backup/stormtrack/jtradfor/ensemble_data/reference/sref_arw_mask.npy')

totalprecips_copy = []
for totalprecip in totalprecips:
	totalprecip[totalprecip_mask] = np.nan
	totalprecip[totalprecip>1000000] = np.nan
	totalprecips_copy.append(totalprecip)

totalprecipsback_copy = []
for totalprecip in totalprecipsback:
	totalprecip[totalprecip_mask] = np.nan
	totalprecip[totalprecip>1000000] = np.nan
	totalprecipsback_copy.append(totalprecip)


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

	hprecip = totalprecips_copy[i] - totalprecipsback_copy[i]

	hprecip[hprecip<0.25] = np.nan
	hprecip[hprecip>1.0] = 1.0

	im = m.imshow(hprecip,zorder=2,interpolation='none',cmap='jet',vmin=0,vmax=1)

cbar = fig.colorbar(im,ax=axes.ravel().tolist(),fraction=0.025,ticks=[0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1.0])
cbar.ax.yaxis.set_tick_params(color='w')
cbar.ax.set_yticklabels([0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1.0],color='w')
hpcpfil = '/gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/uploads/outimages/sref/%s_hpcp_postage.png' % (forecasthoursub)
plt.savefig(hpcpfil,facecolor='#101010',bbox_inches='tight',dpi=500)
plt.close()




