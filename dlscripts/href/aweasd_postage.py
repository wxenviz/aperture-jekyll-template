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
cref_fils = glob.glob('/gpfs_backup/stormtrack/jtradfor/ensemble_data/rawdata/href/*aweasd*npy')

for fhour in range(0,31):

	fig, axes = plt.subplots(nrows=3,ncols=3,figsize=(16,9))
	plt.subplots_adjust(wspace=0.05,hspace=-.10)

	aweasd = []
	for k,model in enumerate(models):
		aweasd_fil = '/gpfs_backup/stormtrack/jtradfor/ensemble_data/rawdata/href/aweasd_%s.npy' % (model)
		aweasd_mask = '/gpfs_backup/stormtrack/jtradfor/ensemble_data/reference/%s_mask.npy' % (masks[k])
		temp_aweasd = np.load(aweasd_fil)[fhour]
		temp_mask = np.load(aweasd_mask)
		temp_aweasd[temp_mask] = np.nan
		aweasd.append(temp_aweasd)


	for i in range(0,8):

		ax = axes.flat[i]
		m = Basemap(projection='lcc',lat_0=5,lon_0=-100,llcrnrlon=-126,llcrnrlat=23,urcrnrlon=-63,urcrnrlat=50,resolution='l',ax=ax)
		shp_info = m.readshapefile('/gpfs_backup/stormtrack/jtradfor/ensemble_data/reference/st99_d00','states',drawbounds=False)

		for nshape,seg in enumerate(m.states):
			poly = Polygon(seg,facecolor='white',edgecolor='white',zorder=1,linewidth=.5)
			poly2 = Polygon(seg,facecolor='none',edgecolor='black',zorder=3,linewidth=.5)
			ax.add_patch(poly)
			ax.add_patch(poly2)
		aweasd[i] = aweasd[i]/2.54
		aweasd[i][aweasd[i]<0.25] = np.nan
		aweasd[i][aweasd[i]>1000000] = np.nan
		aweasd[i][aweasd[i]>24.] = 24.

		bounds = np.linspace(0,24,13)
		norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)

		im = m.imshow(aweasd[i],zorder=2,norm=norm,aspect='equal',interpolation='none',vmin=0.0,vmax=24.0)
		ax.text(.05,.05,model_labels[i],transform=ax.transAxes)

	axes[-1,-1].axis('off')
	cbar = fig.colorbar(im, ax=axes.ravel().tolist(),fraction=0.025,ticks=[0,2,4,6,8,10,12,14,16,18,20,22,24])
	cbar.ax.yaxis.set_tick_params(color='w')
	cbar.ax.set_yticklabels([0,2,4,6,8,10,12,14,16,18,20,22,24],color='w')
	postagefil = '/gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/uploads/outimages/href/%s%s00_aweasd_postage.png' % (datesub,str(fhour).zfill(2))
	plt.savefig(postagefil,facecolor='#101010',bbox_inches='tight',dpi=500)
	plt.close()



