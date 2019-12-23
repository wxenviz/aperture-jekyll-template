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
apcp_fils = glob.glob('/gpfs_backup/stormtrack/jtradfor/ensemble_data/rawdata/href/*apcp*npy')

for fhour in range(0,31):
	plt.figure(figsize=(16,9))

	m = Basemap(projection='lcc',lat_0=5,lon_0=-100,llcrnrlon=-126,llcrnrlat=23,urcrnrlon=-63,urcrnrlat=50,resolution='l')
	shp_info = m.readshapefile('/gpfs_backup/stormtrack/jtradfor/ensemble_data/reference/st99_d00','states',drawbounds=False)

	ax = plt.gca()

	for nshape,seg in enumerate(m.states):
		poly = Polygon(seg,facecolor='white',edgecolor='white',zorder=1,linewidth=.5)
		poly2 = Polygon(seg,facecolor='none',edgecolor='black',zorder=3,linewidth=.5)
		ax.add_patch(poly)
		ax.add_patch(poly2)

	apcps = []
	for k,model in enumerate(models):
		apcp_fil = '/gpfs_backup/stormtrack/jtradfor/ensemble_data/rawdata/href/apcp_%s.npy' % (model)
		apcp_mask = '/gpfs_backup/stormtrack/jtradfor/ensemble_data/reference/%s_mask.npy' % (masks[k])
		temp_apcp = np.load(apcp_fil)[fhour]
		temp_mask = np.load(apcp_mask)
		temp_apcp[temp_mask] = np.nan
		temp_apcp[temp_apcp>1000000] = np.nan
		apcps.append(temp_apcp/25.4) 

	apcp_mean = np.mean(apcps,axis=0)
	apcps.append(apcp_mean)

	apcps_copy = np.copy(apcps)    
	
	apcps_sd = np.zeros_like(apcps_copy[0])

	for i in range(0,len(apcps_sd)):
		for j in range(0,len(apcps_sd[0])):
			c = 0
			for k in range(0,3):
				for l in range(0,3):
					if i%3==k and j%3==l:
						apcps_sd[i,j] = apcps_copy[c][i,j]
					c+=1
	apcps_sd[apcps_sd>1000000] = np.nan
	apcps_sd[apcps_sd>10] = 10.0
	apcps_sd[apcps_sd<=.1] = np.nan

	bounds = np.linspace(0,10,11)
	norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)

	im = m.imshow(apcps_sd,zorder=2,aspect='equal',norm=norm,cmap='jet',interpolation='none',vmin=0,vmax=10)
	cbar = plt.colorbar(im,fraction=0.023,ticks=[0,1,2,3,4,5,6,7,8,9,10])
	cbar.ax.yaxis.set_tick_params(color='w')
	cbar.ax.set_yticklabels([0,1,2,3,4,5,6,7,8,9,10],color='w')
	plt.box(False)
	sdfil = '/gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/uploads/outimages/href/%s%s00_apcp_ab.png' % (datesub,str(fhour).zfill(2))
	plt.savefig(sdfil,facecolor='#101010',bbox_inches='tight',dpi=800)
	plt.close()



