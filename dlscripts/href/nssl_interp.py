import ftplib
import glob
import subprocess as sp
import csv
import numpy as np
import netCDF4 as nc4
import pygrib as pg
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import scipy
import os
import sys
import re
import time
import subprocess as sp
import pickle

from mpl_toolkits.basemap import Basemap
from matplotlib.patches import Polygon
from matplotlib.colors import LinearSegmentedColormap
from scipy.spatial import Delaunay
from scipy.interpolate import LinearNDInterpolator
from shutil import copyfile
from datetime import datetime,timedelta

datesub = str(sys.argv[1])

Xrange = np.arange(-126,-63,0.025)
Yrange = np.arange(23,50,0.025)
[destmeshX,destmeshY] = np.meshgrid(Xrange,Yrange)
destpairs = np.zeros([destmeshX.shape[0]*destmeshX.shape[1],2])
destpairs[:,0] = destmeshX.flatten()
destpairs[:,1] = destmeshY.flatten()

m = Basemap(projection='lcc',lat_0=5,lon_0=-100,llcrnrlon=-126,llcrnrlat=23,urcrnrlon=-63,urcrnrlat=50,resolution='l')

X,Y = m(destmeshX,destmeshY)
nx = int((m.xmax-m.xmin)/3000.)+1; ny = int((m.ymax-m.ymin)/3000.)+1
crefs = []
brefs = []
apcps = []
aweasds = []
hpcps = []
hweasds = []

for fhour in range(0,31):
	fil = '/gpfs_backup/stormtrack/jtradfor/ensemble_data/rawdata/href/%s%s00_nssl.grib' % (datesub,str(fhour).zfill(2))
	backfil = fil[:69] + str(fhour-1).zfill(2) + fil[71:]
	grbs = pg.open(fil)
	lat,lon = grbs[1].latlons()

	if fhour==0:
		apcp = np.zeros_like(lat) 
		aweasd = np.zeros_like(lat)
		hpcp = np.zeros_like(lat)
		hweasd = np.zeros_like(lat)
		bref = np.zeros_like(lat) 
		cref = np.zeros_like(lat)
	else:
		grbsback = pg.open(backfil)
		for grb in grbs:
			if 'fcst time 0' in str(grb) and 'Total Precipitation' in str(grb):
				apcp = grb.values[:]
			elif 'Snow Fall water equivalent' in str(grb):
				aweasd = grb.values[:]
			elif 'level 1000' in str(grb):
				bref = grb.values[:]
			elif 'entireAtmosphere' in str(grb):
				cref = grb.values[:]
		apcp[apcp>1000000] = 0.0
		aweasd[aweasd>1000000] = 0.0
		for grbback in grbsback:
			if 'fcst time 0' in str(grbback) and 'Total Precipitation' in str(grbback):
				apcpback = grbback.values[:]
			elif 'Snow Fall water equivalent' in str(grbback):
				aweasdback = grbback.values[:]
		apcpback[apcpback>1000000] = 0.0
		aweasdback[aweasdback>1000000] - 0.0

		hpcp = apcp - apcpback
		hweasd = aweasd - aweasdback
				
	with open('/gpfs_backup/stormtrack/jtradfor/ensemble_data/rawdata/href/nssl.tri', 'rb') as nssl_tri:
		tri_nssl = pickle.load(nssl_tri)
	nssl_tri.close()

	for j,field in enumerate([bref,cref,apcp,hpcp,aweasd,hweasd]):
		field = field.flatten()
		interpolator = LinearNDInterpolator(tri_nssl,field)
		temp = interpolator(destmeshX,destmeshY)
		field = m.transform_scalar(temp,Xrange,Yrange,nx,ny,masked=True)
		field[field<=0] = np.nan

	crefs.append(cref)
	brefs.append(bref)
	apcps.append(apcp)
	aweasds.append(aweasd)
	hpcps.append(hpcp)
	hweasds.append(hweasd)

np.save('/gpfs_backup/stormtrack/jtradfor/ensemble_data/rawdata/href/cref_nssl.npy',crefs)
np.save('/gpfs_backup/stormtrack/jtradfor/ensemble_data/rawdata/href/bref_nssl.npy',brefs)
np.save('/gpfs_backup/stormtrack/jtradfor/ensemble_data/rawdata/href/apcp_nssl.npy',apcps)
np.save('/gpfs_backup/stormtrack/jtradfor/ensemble_data/rawdata/href/aweasd_nssl.npy',aweasds)
np.save('/gpfs_backup/stormtrack/jtradfor/ensemble_data/rawdata/href/hpcp_nssl.npy',hpcps)
np.save('/gpfs_backup/stormtrack/jtradfor/ensemble_data/rawdata/href/hweasd_nssl.npy',hweasds)
