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

def hsnow(snowtotals,snowtotalsback,outname):
    	plt.figure(figsize=(16,9))

   	m = Basemap(projection='lcc',lat_0=5,lon_0=-100,llcrnrlon=-126,llcrnrlat=23,urcrnrlon=-63,urcrnrlat=50,resolution='h')
    	shp_info = m.readshapefile('/gpfs_backup/stormtrack/jtradfor/ensemble_data/reference/st99_d00','states',drawbounds=False)
    
    	ax = plt.gca()
    
    	for nshape,seg in enumerate(m.states):
        	poly = Polygon(seg,facecolor='white',edgecolor='white',zorder=1,linewidth=1)
		poly2 = Polygon(seg,facecolor='none',edgecolor='black',zorder=3,linewidth=1)
        	ax.add_patch(poly)
		ax.add_patch(poly2)

	hsnows = []
	for i in range(0,9):
		hsnows.append(snowtotals[i] - snowtotalsback[i])

	hsnow_mean = np.mean(hsnows,axis=0)
	hsnow_mean[hsnow_mean<0.0025] = np.nan
	hsnow_mean[hsnow_mean>0.05] = 0.05
	hsnow_mean[-50:,:] = np.nan
	hsnow_mean[:50,:] = np.nan
	hsnow_mean[:,:50] = np.nan
	hsnow_mean[:,-50:] = np.nan
	im = m.imshow(hsnow_mean,interpolation='none',cmap='Blues',vmin=0.,vmax=0.05)
	plt.colorbar(im,fraction=0.023, pad=-0.02)
	plt.box(False)
	hsnowfil = '/gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/uploads/outimages/hrrre/' + outname + '_hweasd_mean.png'
	plt.savefig(hsnowfil,facecolor='#101010',bbox_inches='tight',dpi=500)
	plt.close()


datapaths = glob.glob('/gpfs_backup/stormtrack/jtradfor/ensemble_data/rawdata/hrrre/*snowtotal*')

if len(datapaths)>0:
	latest = 0
	latestpath = datapaths[0]
	for datapath in datapaths:
		if int(os.path.basename(datapath)[9:11]) > latest:
			latest = int(os.path.basename(datapath)[9:11])
			latestpath = datapath
			back1hr = str(latest - 1).zfill(2)
			back1hrpath = os.path.dirname(latestpath) + '/' + os.path.basename(latestpath)[:9] + back1hr + os.path.basename(latestpath)[11:]
	fil = os.path.basename(latestpath)[:13]
	snowtotals = np.load(latestpath)
	filback = os.path.basename(back1hrpath)[:13]
	snowtotalsback = np.load(back1hrpath)
	snowfil = hsnow(snowtotals,snowtotalsback,fil)
