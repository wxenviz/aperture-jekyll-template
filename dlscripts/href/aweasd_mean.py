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
import matplotlib.colors as colors

from mpl_toolkits.basemap import Basemap
from matplotlib.patches import Polygon
from matplotlib.colors import LinearSegmentedColormap
from scipy.spatial import Delaunay
from scipy.interpolate import LinearNDInterpolator
from shutil import copyfile

datesub = str(sys.argv[1])
models = ["namnest","namnest_back","arw","arw_back","nmmb","nmmb_back","hrrr","hrrr_back"]
model_labels = ["Nam Nest","Nam Nest -12h", "ARW", "ARW -12h", "NMMB", "NMMB -12h", "HRRR", "HRRR -6h"]
masks = ["namnest","namnest","arw","arw","nmmb","nmmb","hrrr","hrrr"]
aweasd_fils = glob.glob('/gpfs_backup/stormtrack/jtradfor/ensemble_data/rawdata/href/*aweasd*npy')

for fhour in range(0,31):
	print(fhour)
	plt.figure(figsize=(16,9))

	m = Basemap(projection='lcc',lat_0=5,lon_0=-100,llcrnrlon=-126,llcrnrlat=23,urcrnrlon=-63,urcrnrlat=50,resolution='l')
	shp_info = m.readshapefile('/gpfs_backup/stormtrack/jtradfor/ensemble_data/reference/st99_d00','states',drawbounds=False)

	ax = plt.gca()

	for nshape,seg in enumerate(m.states):
		poly = Polygon(seg,facecolor='white',edgecolor='white',zorder=1,linewidth=.5)
		poly2 = Polygon(seg,facecolor='none',edgecolor='black',zorder=3,linewidth=.5)
		ax.add_patch(poly)
		ax.add_patch(poly2)

	aweasds = []
	for k,model in enumerate(models):
		aweasd_fil = '/gpfs_backup/stormtrack/jtradfor/ensemble_data/rawdata/href/aweasd_%s.npy' % (model)
		aweasd_mask = '/gpfs_backup/stormtrack/jtradfor/ensemble_data/reference/%s_mask.npy' % (masks[k])
		temp_aweasd = np.load(aweasd_fil)[fhour]
		temp_mask = np.load(aweasd_mask)
		temp_aweasd[temp_mask] = np.nan
		aweasds.append(temp_aweasd/2.54)

	aweasds_copy = np.copy(aweasds) 
	aweasds_mean = np.mean(aweasds_copy,axis=0)
	aweasds_mean[aweasds_mean<0.25] = np.nan
	aweasds_mean[aweasds_mean>1000000] = np.nan
	aweasds_mean[aweasds_mean>24] = 24.0

	bounds = np.linspace(0,24,13)
	norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)

	im = m.imshow(aweasds_mean,zorder=2,norm=norm,aspect='equal',interpolation='none',vmin=0.0,vmax=24.0)
	cbar = plt.colorbar(im,fraction=0.023,ticks=[0,2,4,6,8,10,12,14,16,18,20,22,24])
	cbar.ax.yaxis.set_tick_params(color='w')
	cbar.ax.set_yticklabels([0,2,4,6,8,10,12,14,16,18,20,22,24],color='w')
	plt.box(False)
	meanfil = '/gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/uploads/outimages/href/%s%s00_aweasd_mean.png' % (datesub,str(fhour).zfill(2))
	plt.savefig(meanfil,facecolor='#101010',bbox_inches='tight',dpi=500)
	plt.close()



