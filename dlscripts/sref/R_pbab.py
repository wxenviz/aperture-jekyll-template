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

reflectivities = np.load('/gpfs_backup/stormtrack/jtradfor/ensemble_data/rawdata/sref/%s_creflect.npy' % (forecasthoursub))
reflectivities_mask = np.load('/gpfs_backup/stormtrack/jtradfor/ensemble_data/reference/sref_arw_mask.npy')

reflect_40 = []
for c,reflect_member in enumerate(reflectivities):
	reflect_member[reflectivities_mask] = np.nan
	reflect_member[reflect_member<20] = np.nan
	reflect_member[reflect_member>1000000] = np.nan
	reflect_member[reflect_member>=20] = c+1
	reflect_40.append(reflect_member)


reflectivities_pb_sd = np.zeros_like(reflect_40[0])

for i in range(0,len(reflectivities_pb_sd)):
	for j in range(0,len(reflectivities_pb_sd[0])):
    		c = 0
    		for k in range(0,4):
			for l in range(0,4):
	    			if i%4==k and j%4==l:
					reflectivities_pb_sd[i,j] = reflect_40[c][i,j]
	    			c+=1
	
reflectivities_pb_sd[reflectivities_pb_sd==0]=np.nan

m.imshow(reflectivities_pb_sd,zorder=2,aspect='equal',interpolation='none',cmap='tab20',vmin=1,vmax=20)
plt.box(False)
pbsdfil = '/gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/uploads/outimages/sref/%s_R_pbab.png' % (forecasthoursub)
plt.savefig(pbsdfil,facecolor='#101010',bbox_inches='tight',dpi=800)
plt.close()


