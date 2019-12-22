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
import scipy
import os
import sys
import re
import time
import subprocess as sp

from mpl_toolkits.basemap import Basemap
from matplotlib.patches import Polygon
from matplotlib.colors import LinearSegmentedColormap
from scipy.spatial import Delaunay
from scipy.interpolate import LinearNDInterpolator
from shutil import copyfile
from datetime import datetime,timedelta

sp.call('/usr/bin/rm /gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/dlscripts/sref/*_err',shell=True)
sp.call('/usr/bin/rm /gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/dlscripts/sref/*_out',shell=True)
sp.call('/usr/bin/rm /gpfs_backup/stormtrack/jtradfor/ensemble_data/rawdata/sref/*',shell=True)

latlons = np.load('/gpfs_backup/stormtrack/jtradfor/ensemble_data/reference/sref_master.npy')
lats = latlons[0]
lons = latlons[1]

orig = np.zeros([lons.shape[0]*lons.shape[1],2])
orig[:,0] = lons.flatten()
orig[:,1] = lats.flatten()
tri = Delaunay(orig)

Xrange = np.arange(-126,-63,0.10)
Yrange = np.arange(23,50,0.10)
[destmeshX,destmeshY] = np.meshgrid(Xrange,Yrange)
destpairs = np.zeros([destmeshX.shape[0]*destmeshX.shape[1],2])
destpairs[:,0] = destmeshX.flatten()
destpairs[:,1] = destmeshY.flatten()

m = Basemap(projection='lcc',lat_0=5,lon_0=-100,llcrnrlon=-126,llcrnrlat=23,urcrnrlon=-63,urcrnrlat=50,resolution='h')
nx = int((m.xmax-m.xmin)/16000.)+1; ny = int((m.ymax-m.ymin)/16000.)+1


get_inv = '/gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/get_inv.pl'
get_grib = '/gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/get_grib.pl'
fields = "':(APCP|REFC|WEASD):'"

offset = 4 
timecheck = datetime.utcnow()
yrmnthday = (timecheck - timedelta(hours=offset)).strftime("%Y%m%d")
hr = ((timecheck - timedelta(hours=offset)).strftime("%H")).zfill(2)

newdate = (timecheck - timedelta(hours=offset)).strftime("%y%j%H") + "00"
print(newdate)
os.chdir('/gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/')
with open('/gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/sref.html') as f:
	s = f.read()
with open('/gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/sref.html', 'w') as f:
	s = re.sub(r'datei = \"\d{9}',r'datei = "' + newdate,s)
	f.write(s)
	f.close()
sp.call('git add -A',shell=True)
sp.call("git commit -m 'update html'",shell=True)
sp.call('git push origin master',shell=True)
print(newdate)
for fhour in range(0,90,3):
	outfil = newdate + str(fhour).zfill(2) + "00"
	outfilstart = newdate + "0000"
	outfilback = newdate + str(fhour - 3).zfill(2) + "00"
	print(fhour)
	fils = []	
	for core in ['arw','nmb']:
		for perturbation in ['n','p']:
			for i in range(1,5):
				griburl = "https://www.ftp.ncep.noaa.gov/data/nccf/com/sref/prod/sref.%s/%s/pgrb/sref_%s.t%sz.pgrb132.%s%s.f%s.grib2" % (yrmnthday,hr,core,hr,perturbation,i,str(fhour).zfill(2))
				idxurl = "https://www.ftp.ncep.noaa.gov/data/nccf/com/sref/prod/sref.%s/%s/pgrb/sref_%s.t%sz.pgrb132.%s%s.f%s.grib2.idx" % (yrmnthday,hr,core,hr,perturbation,i,str(fhour).zfill(2))
				outfilgrib = '/gpfs_backup/stormtrack/jtradfor/ensemble_data/rawdata/sref/%s_%s_%s%s.grib2' % (outfil, core, perturbation, i) 
				existfils = []
				k = 0
				while outfilgrib not in existfils:
					sp.call("%s %s | egrep %s | %s %s %s" % (get_inv, idxurl, fields, get_grib, griburl, outfilgrib),shell=True)
					existfils = glob.glob('/gpfs_backup/stormtrack/jtradfor/ensemble_data/rawdata/sref/*.grib2')

					if k>=1:
						time.sleep(30)
					if k>=5:
						break
					else:
						k+=1
					
				fils.append(outfilgrib)


	reflectivities = []
	totalprecips = []
	totalsnows = []

	for c,fil in enumerate(fils):

		grbs = pg.open(fil)
		grbs.seek(0)
		reflectivity = grbs[4].values[:]
		reflectivity = reflectivity.flatten()
		interpolator = LinearNDInterpolator(tri,reflectivity)
		reflectivityi = interpolator(destmeshX,destmeshY)
		reflect_interp = m.transform_scalar(reflectivityi,Xrange,Yrange,nx,ny,masked=False)
		reflectivities.append(reflect_interp)

		if(str(fhour).zfill(2)!='00'):
			
			if 'arw' in fil:
				totalprecip = grbs[2].values[:]
				totalprecip = totalprecip.flatten()
				interpolator = LinearNDInterpolator(tri,totalprecip)
				totalprecipi = interpolator(destmeshX,destmeshY)
				totalprecip_interp = m.transform_scalar(totalprecipi,Xrange,Yrange,nx,ny,masked=False)
				
				totalsnow = grbs[3].values[:]
				totalsnow = totalsnow.flatten()
				interpolator = LinearNDInterpolator(tri,totalsnow)
				totalsnowi = interpolator(destmeshX,destmeshY)
				totalsnow_interp = m.transform_scalar(totalsnowi,Xrange,Yrange,nx,ny,masked=False)
				
			else:
				precip3hr = grbs[2].values[:]
				precip3hr = precip3hr.flatten()
				interpolator = LinearNDInterpolator(tri,precip3hr)
				precip3hri = interpolator(destmeshX,destmeshY)
				precip3hr_interp = m.transform_scalar(precip3hri,Xrange,Yrange,nx,ny,masked=False)
				totalprecip_interp = (np.load('/gpfs_backup/stormtrack/jtradfor/ensemble_data/rawdata/sref/%s_apcp.npy' % (outfilback))[c])*25.4 + precip3hr_interp

				snow3hr = grbs[3].values[:]
				snow3hr = snow3hr.flatten()
				interpolator = LinearNDInterpolator(tri,snow3hr)
				snow3hri = interpolator(destmeshX,destmeshY)
				snow3hr_interp = m.transform_scalar(snow3hri,Xrange,Yrange,nx,ny,masked=False)
				snow3hr_interp = snow3hr_interp/2.54
				totalsnow_interp = (np.load('/gpfs_backup/stormtrack/jtradfor/ensemble_data/rawdata/sref/%s_asnow.npy' % (outfilback))[c])*2.54 + snow3hr_interp
		else:
			totalprecip_interp = np.zeros_like(reflect_interp)
			totalsnow_interp = np.zeros_like(reflect_interp)

		totalprecips.append((totalprecip_interp)/25.4)
		totalsnows.append((totalsnow_interp)/2.54)

	np.save('/gpfs_backup/stormtrack/jtradfor/ensemble_data/rawdata/sref/%s_creflect.npy' % (outfil),reflectivities)
	np.save('/gpfs_backup/stormtrack/jtradfor/ensemble_data/rawdata/sref/%s_apcp.npy' % (outfil),totalprecips)
	np.save('/gpfs_backup/stormtrack/jtradfor/ensemble_data/rawdata/sref/%s_asnow.npy' % (outfil),totalsnows)
	sp.call("sed -i -- 's,set datestr=.*,set datestr=\"%s\",g' /gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/dlscripts/sref/*_r.csh" % (outfil),shell=True)
	time.sleep(5)
	sp.call("bsub < /gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/dlscripts/sref/apcp_ab_r.csh",shell=True)
	sp.call("bsub < /gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/dlscripts/sref/apcp_mean_r.csh",shell=True)
	sp.call("bsub < /gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/dlscripts/sref/apcp_postage_r.csh",shell=True)
	sp.call("bsub < /gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/dlscripts/sref/aweasd_ab_r.csh",shell=True)
	sp.call("bsub < /gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/dlscripts/sref/aweasd_mean_r.csh",shell=True)	
	sp.call("bsub < /gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/dlscripts/sref/hpcp_ab_r.csh",shell=True)
	sp.call("bsub < /gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/dlscripts/sref/hpcp_postage_r.csh",shell=True)
	sp.call("bsub < /gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/dlscripts/sref/hpcp_mean_r.csh",shell=True)
	sp.call("bsub < /gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/dlscripts/sref/hweasd_ab_r.csh",shell=True)
	sp.call("bsub < /gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/dlscripts/sref/hweasd_mean_r.csh",shell=True)		
	sp.call("bsub < /gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/dlscripts/sref/R_pbab_r.csh",shell=True)
	sp.call("bsub < /gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/dlscripts/sref/R_pb_r.csh",shell=True)				
	sp.call("bsub < /gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/dlscripts/sref/R_mean_r.csh",shell=True)	
	sp.call("bsub < /gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/dlscripts/sref/R_ab_r.csh",shell=True)
	sp.call("bsub < /gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/dlscripts/sref/R_postage_r.csh",shell=True)
	
	time.sleep(180)
	sp.call('git add -A',shell=True)
	sp.call("git commit -m 'new images'",shell=True)
	sp.call('git push origin master',shell=True)










