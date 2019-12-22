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

os.chdir('/gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/')

def apcp(totalPrecips,outname):
    	plt.figure(figsize=(16,9))

   	m = Basemap(projection='lcc',lat_0=5,lon_0=-100,llcrnrlon=-126,llcrnrlat=23,urcrnrlon=-63,urcrnrlat=50,resolution='h')
    	shp_info = m.readshapefile('/gpfs_backup/stormtrack/jtradfor/ensemble_data/reference/st99_d00','states',drawbounds=False)
    
    	ax = plt.gca()
    
    	for nshape,seg in enumerate(m.states):
        	poly = Polygon(seg,facecolor='white',edgecolor='white',zorder=1,linewidth=1)
		poly2 = Polygon(seg,facecolor='none',edgecolor='black',zorder=3,linewidth=1)
        	ax.add_patch(poly)
		ax.add_patch(poly2)

	totalPrecip_mean = np.mean(totalPrecips,axis=0)
	totalPrecip_mean[totalPrecip_mean<0.25] = np.nan
	totalPrecip_mean[totalPrecip_mean>100] = 100
	totalPrecip_mean[-50:,:] = np.nan
	totalPrecip_mean[:50,:] = np.nan
	totalPrecip_mean[:,:50] = np.nan
	totalPrecip_mean[:,-50:] = np.nan
	m.imshow(totalPrecip_mean,zorder=2,interpolation='none',cmap='jet',vmin=0,vmax=100)
	plt.colorbar(fraction=0.023, pad=-0.02)
	plt.box(False)
	apcpfil = '/gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/uploads/outimages/hrrre/' + outname + '_apcp_mean.png'
	plt.savefig(apcpfil,facecolor='#101010',bbox_inches='tight',dpi=500)
	plt.close()

datapaths = glob.glob('/gpfs_backup/stormtrack/jtradfor/ensemble_data/rawdata/hrrre/*apcp*')

if len(datapaths)>0:

	latest = 0
	latestpath = datapaths[0]
	for datapath in datapaths:
		if int(os.path.basename(datapath)[9:11]) > latest:
			latest = int(os.path.basename(datapath)[9:11])
			latestpath = datapath
	fil = os.path.basename(latestpath)[:13]
	totalPrecips = np.load(latestpath)
	apcpfil = apcp(totalPrecips,fil)
