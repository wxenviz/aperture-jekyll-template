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

datesub = str(sys.argv[1]) 
levels = []
colors = []
with open('/gpfs_backup/stormtrack/jtradfor/ensemble_data/reference/N0Q_Color_Lookup.csv','r') as colorcsv:
	colorreader = csv.reader(colorcsv,delimiter=',')
	for line in colorreader:
		if float(line[1])>=0 and float(line[1])<=60:
			colorints = [int(i) for i in line[2:]]
			colors.append((colorints))
			levels.append(float(line[1]))  

colors = np.array(colors)/255.0
cmap1 = LinearSegmentedColormap.from_list("my_colormap",colors,N=len(levels),gamma=1.0) 
          
models = ["namnest","namnest_back","arw","arw_back","nmmb","nmmb_back","hrrr","hrrr_back"]
model_labels = ["Nam Nest","Nam Nest -12h", "ARW", "ARW -12h", "NMMB", "NMMB -12h", "HRRR", "HRRR -6h"]
masks = ["namnest","namnest","arw","arw","nmmb","nmmb","hrrr","hrrr"]
cref_fils = glob.glob('/gpfs_backup/stormtrack/jtradfor/ensemble_data/rawdata/href/*cref*npy')


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

	reflectivities = []
	for k,model in enumerate(models):
		cref_fil = '/gpfs_backup/stormtrack/jtradfor/ensemble_data/rawdata/href/cref_%s.npy' % (model)
		cref_mask = '/gpfs_backup/stormtrack/jtradfor/ensemble_data/reference/%s_mask.npy' % (masks[k])
		temp_cref = np.load(cref_fil)[fhour]
		temp_cref[temp_cref<0] = 0.0
		temp_mask = np.load(cref_mask)
		temp_cref[temp_mask] = np.nan
		temp_cref[temp_cref>1000000] = np.nan
		reflectivities.append(temp_cref) 

	rmean = np.mean(reflectivities,axis=0)
	reflectivities.append(rmean)
	reflectivities_copy = np.copy(reflectivities)	
	reflect_40 = []
	for c,reflect_member in enumerate(reflectivities_copy):
		reflect_member[reflect_member<20] = 0
		reflect_member[reflect_member>=20] = c+1
		reflect_40.append(reflect_member)

	reflect_pb_ab = np.zeros_like(reflect_40[0])

	for i in range(0,len(reflect_pb_ab)):
		for j in range(0,len(reflect_pb_ab[0])):
			c = 0
			for k in range(0,3):
				for l in range(0,3):
					if i%3==k and j%3==l:
						reflect_pb_ab[i,j] = reflect_40[c][i,j]
					c+=1

	reflect_pb_ab[reflect_pb_ab==0] = np.nan

	im = m.imshow(reflect_pb_ab,zorder=2,cmap='tab10',interpolation='none',vmin=1,vmax=10)
	plt.box(False)
	pbfil = '/gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/uploads/outimages/href/%s%s00_R_pbab.png' % (datesub,str(fhour).zfill(2))
	plt.savefig(pbfil,facecolor='#101010',bbox_inches='tight',dpi=500)
	plt.close()



