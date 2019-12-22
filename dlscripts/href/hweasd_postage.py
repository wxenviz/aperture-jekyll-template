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
cref_fils = glob.glob('/gpfs_backup/stormtrack/jtradfor/ensemble_data/rawdata/href/*hweasd*npy')

for fhour in range(0,31):

	fig, axes = plt.subplots(nrows=3,ncols=3,figsize=(16,9))
	plt.subplots_adjust(wspace=0.05,hspace=-.10)

	hweasd = []
	for k,model in enumerate(models):
		hweasd_fil = '/gpfs_backup/stormtrack/jtradfor/ensemble_data/rawdata/href/hweasd_%s.npy' % (model)
		hweasd_mask = '/gpfs_backup/stormtrack/jtradfor/ensemble_data/reference/%s_mask.npy' % (masks[k])
		temp_hweasd = np.load(hweasd_fil)[fhour]
		temp_mask = np.load(hweasd_mask)
		temp_hweasd[temp_mask] = np.nan
		hweasd.append(temp_hweasd/2.54) 


	for i in range(0,8):

		ax = axes.flat[i]
		m = Basemap(projection='lcc',lat_0=5,lon_0=-100,llcrnrlon=-126,llcrnrlat=23,urcrnrlon=-63,urcrnrlat=50,resolution='l',ax=ax)
		shp_info = m.readshapefile('/gpfs_backup/stormtrack/jtradfor/ensemble_data/reference/st99_d00','states',drawbounds=False)

		for nshape,seg in enumerate(m.states):
			poly = Polygon(seg,facecolor='white',edgecolor='white',zorder=1,linewidth=.5)
			poly2 = Polygon(seg,facecolor='none',edgecolor='black',zorder=3,linewidth=.5)
			ax.add_patch(poly)
			ax.add_patch(poly2)

		hweasd[i][hweasd[i]<0.02] = np.nan
		hweasd[i][hweasd[i]>1000000] = np.nan
		hweasd[i][hweasd[i]>2.] = 2.0

		bounds = np.linspace(0,2,9)
		norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)

		im = m.imshow(hweasd[i],zorder=2,norm=norm,aspect='equal',interpolation='none',cmap='jet',vmin=0.0,vmax=2.0)
		ax.text(.05,.05,model_labels[i],transform=ax.transAxes)

	axes[-1,-1].axis('off')
	cbar = fig.colorbar(im, ax=axes.ravel().tolist(),fraction=0.025,ticks=[0,.2,.4,.6,.8,1.0,1.2,1.4,1.6,1.8,2.0])
	cbar.ax.yaxis.set_tick_params(color='w')
	cbar.ax.set_yticklabels([0,.2,.4,.6,.8,1.0,1.2,1.4,1.6,1.8,2.0],color='w')
	postagefil = '/gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/uploads/outimages/href/%s%s00_hweasd_postage.png' % (datesub,str(fhour).zfill(2))
	plt.savefig(postagefil,facecolor='#101010',bbox_inches='tight',dpi=500)
	plt.close()



