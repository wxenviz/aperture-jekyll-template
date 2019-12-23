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
temps = []
for fhour in range(0,31):
	print(fhour)
	fil = '/gpfs_backup/stormtrack/jtradfor/ensemble_data/rawdata/href/%s%s00_namnest.grib2' % (datesub,str(fhour).zfill(2))
	backfil = fil[:69] + str(fhour-1).zfill(2) + fil[71:]
	grbs = pg.open(fil)
	lat,lon = grbs[1].latlons()

	
	if fhour==0:
		hpcp = np.zeros((1059,1799))
		apcp = np.zeros((1059,1799))
		hweasd = np.zeros((1059,1799))
		aweasd = np.zeros((1059,1799))
		bref = np.zeros((1059,1799))
		cref = np.zeros((1059,1799))
		temperature = np.zeros((1059,1799))
	else:
		grbsback = pg.open(backfil)
	
		if fhour%3==1:
			for grb in grbs:
				if 'Total Precipitation' in str(grb):
					hpcp = np.array(grb.values[:])
				elif 'Water equivalent' in str(grb) and 'instant' not in str(grb):
					hweasd = np.array(grb.values[:])
				elif 'level 1000 m' in str(grb):
					bref = grb.values[:]
				elif 'Maximum/Composite' in str(grb):
					cref = grb.values[:]
				elif 'temperature' in str(grb):
					temperature = grb.values[:]
		else:
			for grb in grbs:
				if 'Total Precipitation' in str(grb):
					temphpcp = np.array(grb.values[:])
				elif 'Water equivalent' in str(grb) and 'instant' not in str(grb):
					temphweasd = np.array(grb.values[:])
				elif 'level 1000 m' in str(grb):
					bref = grb.values[:]
				elif 'Maximum/Composite' in str(grb):
					cref = grb.values[:]
				elif 'temperature' in str(grb):
					temperature = grb.values[:]

			for grbback in grbsback:
				if 'Total Precipitation' in str(grbback):
					temphpcpback = np.array(grbback.values[:])
				elif 'Water equivalent' in str(grbback):
					temphweasdback = np.array(grbback.values[:])
			
			hpcp = temphpcp - temphpcpback
			hweasd = temphweasd - temphweasdback

		apcp = apcp + hpcp
		aweasd = aweasd + hweasd

	with open('/gpfs_backup/stormtrack/jtradfor/ensemble_data/reference/namnest.tri', 'rb') as namnest_tri:
		tri_namnest = pickle.load(namnest_tri)
	namnest_tri.close()

	cref_flatten = cref.flatten()
	interpolator = LinearNDInterpolator(tri_namnest,cref_flatten)
	temp = interpolator(destmeshX,destmeshY)
	cref_interp = m.transform_scalar(temp,Xrange,Yrange,nx,ny,masked=True)

	bref_flatten = bref.flatten()
	interpolator = LinearNDInterpolator(tri_namnest,bref_flatten)
	temp = interpolator(destmeshX,destmeshY)
	bref_interp = m.transform_scalar(temp,Xrange,Yrange,nx,ny,masked=True)

	apcp_flatten = apcp.flatten()
	interpolator = LinearNDInterpolator(tri_namnest,apcp_flatten)
	temp = interpolator(destmeshX,destmeshY)
	apcp_interp = m.transform_scalar(temp,Xrange,Yrange,nx,ny,masked=True)

	hpcp_flatten = hpcp.flatten()
	interpolator = LinearNDInterpolator(tri_namnest,hpcp_flatten)
	temp = interpolator(destmeshX,destmeshY)
	hpcp_interp = m.transform_scalar(temp,Xrange,Yrange,nx,ny,masked=True)

	aweasd_flatten = aweasd.flatten()
	interpolator = LinearNDInterpolator(tri_namnest,aweasd_flatten)
	temp = interpolator(destmeshX,destmeshY)
	aweasd_interp = m.transform_scalar(temp,Xrange,Yrange,nx,ny,masked=True)

	hweasd_flatten = hweasd.flatten()
	interpolator = LinearNDInterpolator(tri_namnest,hweasd_flatten)
	temp = interpolator(destmeshX,destmeshY)
	hweasd_interp = m.transform_scalar(temp,Xrange,Yrange,nx,ny,masked=True)

	temperature_flatten = temperature.flatten()
	interpolator = LinearNDInterpolator(tri_namnest,temperature_flatten)
	temp = interpolator(destmeshX,destmeshY)
	temperature_interp = m.transform_scalar(temp,Xrange,Yrange,nx,ny,masked=True)

	crefs.append(cref_interp.data)
	brefs.append(bref_interp.data)
	apcps.append(apcp_interp.data)
	aweasds.append(aweasd_interp.data)
	hpcps.append(hpcp_interp.data)
	hweasds.append(hweasd_interp.data)
	temps.append(temperature_interp.data)

np.save('/gpfs_backup/stormtrack/jtradfor/ensemble_data/rawdata/href/cref_namnest.npy',crefs)
np.save('/gpfs_backup/stormtrack/jtradfor/ensemble_data/rawdata/href/bref_namnest.npy',brefs)
np.save('/gpfs_backup/stormtrack/jtradfor/ensemble_data/rawdata/href/apcp_namnest.npy',apcps)
np.save('/gpfs_backup/stormtrack/jtradfor/ensemble_data/rawdata/href/aweasd_namnest.npy',aweasds)
np.save('/gpfs_backup/stormtrack/jtradfor/ensemble_data/rawdata/href/hpcp_namnest.npy',hpcps)
np.save('/gpfs_backup/stormtrack/jtradfor/ensemble_data/rawdata/href/hweasd_namnest.npy',hweasds)
np.save('/gpfs_backup/stormtrack/jtradfor/ensemble_data/rawdata/href/temp_namnest.npy',temps)
