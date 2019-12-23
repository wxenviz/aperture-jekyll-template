#!/usr/bin/env python
# coding: utf-8

# In[ ]:


"""
Created on Mon Apr 29 09:56:21 2019

@author: jacob
"""
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
import re

from mpl_toolkits.basemap import Basemap
from matplotlib.patches import Polygon
from matplotlib.colors import LinearSegmentedColormap
from scipy.spatial import Delaunay
from scipy.interpolate import LinearNDInterpolator
from shutil import copyfile

existsinall = np.load('/gpfs_backup/stormtrack/jtradfor/ensemble_data/rawdata/hrrre/existfils.npy')

if len(existsinall)>0:

	#Required mapping data for plotting
	levels = []
	colors = []
	with open('/gpfs_backup/stormtrack/jtradfor/ensemble_data/reference/N0Q_Color_Lookup.csv','r') as colorcsv:
	    colorreader = csv.reader(colorcsv,delimiter=',')
	    for line in colorreader:
		colorints = [int(i) for i in line[2:]]
		colors.append((colorints))
		levels.append(float(line[1]))  
	colors = np.array(colors)/255.0
	cmap1 = LinearSegmentedColormap.from_list("my_colormap",colors,N=len(levels),gamma=1.0)        

	curve = nc4.Dataset('/gpfs_backup/stormtrack/jtradfor/ensemble_data/reference/master.nc','r')
	lon_curve = np.flipud(curve.variables['gridlon_0'][:,:])
	lat_curve = np.flipud(curve.variables['gridlat_0'][:,:])
	orig = np.zeros([lon_curve.shape[0]*lon_curve.shape[1],2])
	orig[:,0] = lon_curve.flatten()
	orig[:,1] = lat_curve.flatten()
	tri = Delaunay(orig)

	Xrange = np.arange(-126,-63,0.05)
	Yrange = np.arange(23,50,0.05)
	[destmeshX,destmeshY] = np.meshgrid(Xrange,Yrange)
	destpairs = np.zeros([destmeshX.shape[0]*destmeshX.shape[1],2])
	destpairs[:,0] = destmeshX.flatten()
	destpairs[:,1] = destmeshY.flatten()


	m = Basemap(projection='lcc',lat_0=5,lon_0=-100,llcrnrlon=-126,llcrnrlat=23,urcrnrlon=-63,urcrnrlat=50,resolution='h')
	shp_info = m.readshapefile('/gpfs_backup/stormtrack/jtradfor/ensemble_data/reference/st99_d00','states',drawbounds=False)

	ax = plt.gca()

	for nshape,seg in enumerate(m.states):
	    poly = Polygon(seg,facecolor='none',edgecolor='black',zorder=1,linewidth=1)
	    ax.add_patch(poly)

	X,Y = m(destmeshX,destmeshY)
	nx = int((m.xmax-m.xmin)/3000.)+1; ny = int((m.ymax-m.ymin)/3000.)+1


	fil = existsinall[0]
	ensmembers = []
	for i in range(1,10):
		ensmembers.append('/gpfs_backup/stormtrack/jtradfor/ensemble_data/rawdata/hrrre/' + fil + '_' + str(i))
	    
	reflectivities = []
	basereflectivities = []
	basereflectivitiesdetect = []
	temps = []
	totalPrecips = []
	snowtotals = []
	for c,ensmember in enumerate(ensmembers):

		grbs = pg.open(ensmember)

		reflectivity = np.flipud(grbs[1].values[:])
		reflectivity[:,0:5] = -10.0
		reflectivity[0:5,:] = -10.0
		reflectivity[:,-5:] = -10.0
		reflectivity[-5:,:] = -10.0
		reflectivity = reflectivity.flatten()
		interpolator = LinearNDInterpolator(tri,reflectivity)
		reflectivityi = interpolator(destmeshX,destmeshY)
		reflectivityi[np.isnan(reflectivityi)] = -10.0
		reflect_interp = m.transform_scalar(reflectivityi,Xrange,Yrange,nx,ny,masked=True)
		reflectivities.append(reflect_interp)

		breflectivity = np.flipud(grbs[6].values[:])
		breflectivity[:,0:5] = -10.0
		breflectivity[0:5,:] = -10.0
		breflectivity[:,-5:] = -10.0
		breflectivity[-5:,:] = -10.0
		basereflectivitiesdetect.append(np.copy(breflectivity))
		breflectivity = breflectivity.flatten()
		interpolator = LinearNDInterpolator(tri,breflectivity)
		breflectivityi = interpolator(destmeshX,destmeshY)
		breflectivityi[np.isnan(breflectivityi)] = -10.0
		breflect_interp = m.transform_scalar(breflectivityi,Xrange,Yrange,nx,ny,masked=True)
		basereflectivities.append(breflect_interp)

		snowtotal = np.flipud(grbs[46].values[:])
		snowtotal[:,0:5] = 0.
		snowtotal[0:5,:] = 0.
		snowtotal[:,-5:] = 0.
		snowtotal[-5:,:] = 0.
		snowtotal = snowtotal.flatten()
		interpolator = LinearNDInterpolator(tri,snowtotal)
		snowtotali = interpolator(destmeshX,destmeshY)
		snowtotali[np.isnan(snowtotali)] = 0.0
		snowtotal_interp = m.transform_scalar(snowtotali,Xrange,Yrange,nx,ny,masked=True)
		snowtotals.append(snowtotal_interp)
		
		temp = np.flipud(grbs[70].values[:])
		temps.append(temp)

		fhour = os.path.basename(ensmember)[9:11]
		if fhour!='00':
			totalPrecip = np.flipud(grbs[154].values[:])
			totalPrecip[:,0:5] = 0.0
			totalPrecip[0:5,:] = 0.0
			totalPrecip[:,-5:] = 0.0
			totalPrecip[-5:,:] = 0.0
			totalPrecip = totalPrecip.flatten()
			interpolator = LinearNDInterpolator(tri,totalPrecip)
			totalPrecipi = interpolator(destmeshX,destmeshY)
			totalPrecipi[np.isnan(totalPrecipi)] = 0.0
			totalPrecip_interp = m.transform_scalar(totalPrecipi,Xrange,Yrange,nx,ny,masked=True)
			totalPrecips.append(totalPrecip_interp)
		else:
			totalPrecip_interp = np.zeros_like(reflect_interp)
			totalPrecips.append(totalPrecip_interp)
		

	for i in range(1,10):
		sp.call('rm /gpfs_backup/stormtrack/jtradfor/ensemble_data/rawdata/hrrre/' + fil + '_' + str(i),shell=True)
	sp.call('rm /gpfs_backup/stormtrack/jtradfor/ensemble_data/rawdata/hrrre/existfils.npy',shell=True)
    
	np.save('/gpfs_backup/stormtrack/jtradfor/ensemble_data/rawdata/hrrre/' + fil + '_creflect.npy',reflectivities)
	np.save('/gpfs_backup/stormtrack/jtradfor/ensemble_data/rawdata/hrrre/' + fil + '_breflect.npy',basereflectivities)
	np.save('/gpfs_backup/stormtrack/jtradfor/ensemble_data/rawdata/hrrre/' + fil + '_breflectdetect.npy',basereflectivitiesdetect)
	np.save('/gpfs_backup/stormtrack/jtradfor/ensemble_data/rawdata/hrrre/' + fil + '_temps.npy',temps)
	np.save('/gpfs_backup/stormtrack/jtradfor/ensemble_data/rawdata/hrrre/' + fil + '_apcp.npy',totalPrecips)
	np.save('/gpfs_backup/stormtrack/jtradfor/ensemble_data/rawdata/hrrre/' + fil + '_snowtotal.npy',snowtotals)
