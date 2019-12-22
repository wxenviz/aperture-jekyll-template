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
import pickle

from mpl_toolkits.basemap import Basemap
from matplotlib.patches import Polygon
from matplotlib.colors import LinearSegmentedColormap
from scipy.spatial import Delaunay
from scipy.interpolate import LinearNDInterpolator
from shutil import copyfile
from datetime import datetime,timedelta

sp.call('/usr/bin/rm /gpfs_backup/stormtrack/jtradfor/ensemble_data/rawdata/href/*.npy',shell=True)
sp.call('/usr/bin/rm /gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/dlscripts/href/*_err',shell=True)
sp.call('/usr/bin/rm /gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/dlscripts/href/*_out',shell=True)

os.chdir('/gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io')

sp.call('git ls-files --others --exclude-standard > /tmp/my_untracked_files',shell=True)
sp.call('mv .git/config /tmp/',shell=True)
sp.call('rm -rf .git',shell=True)
sp.call('git init',shell=True)
sp.call('git add .',shell=True)
sp.call('mv /tmp/config .git/',shell=True)
sp.call('cat /tmp/my_untracked_files | xargs -0 git rm --cached',shell=True)
sp.call("git commit -m 'Cleared History'",shell=True)
sp.call('git push -u --force origin master',shell=True)
time.sleep(60)
get_inv = '/gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/get_inv.pl'
get_grib = '/gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/get_grib.pl'
fields = "':(APCP|REFD|REFC|WEASD|TMIN):'"
hrrrfields = "':(APCP|REFD|REFC|WEASD|TMP:2 m above ground):'"

offset = 3 
timecheck = datetime.utcnow()
yrmnthday = (timecheck - timedelta(hours=offset)).strftime("%Y%m%d")
backyrmnthday = (timecheck - timedelta(hours=offset+12)).strftime("%Y%m%d")
back6yrmnthday = (timecheck - timedelta(hours=offset+6)).strftime("%Y%m%d")

hr = ((timecheck - timedelta(hours=offset)).strftime("%H")).zfill(2)
backhr = ((timecheck - timedelta(hours=offset+12)).strftime("%H")).zfill(2)
back6hr = ((timecheck - timedelta(hours=offset+6)).strftime("%H")).zfill(2)

newdate = (timecheck - timedelta(hours=offset)).strftime("%y%j%H") + "00"

with open('/gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/href.html') as f:
	s = f.read()
with open('/gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/href.html', 'w') as f:
	s = re.sub(r'datei = \"\d{9}',r'datei = "' + newdate,s)
	f.write(s)
	f.close()
sp.call('git add -A',shell=True)
sp.call("git commit -m 'update html'",shell=True)
sp.call('git push origin master',shell=True) 

print(newdate)
brefs = [[] for x in xrange(31)]
crefs = [[] for x in xrange(31)]
apcps = [[] for x in xrange(31)]
aweasds = [[] for x in xrange(31)]
hpcps = [[] for x in xrange(31)]
hweasds = [[] for x in xrange(31)]

for fhour in range(0,31,1):
	outfil = newdate + str(fhour).zfill(2) + "00"

	for core in ['arw','nmmb']:
		
		griburl = "https://www.ftp.ncep.noaa.gov/data/nccf/com/hiresw/prod/hiresw.%s/hiresw.t%sz.%s_2p5km.f%s.conus.grib2" % (yrmnthday,hr,core,str(fhour).zfill(2))
		idxurl = "https://www.ftp.ncep.noaa.gov/data/nccf/com/hiresw/prod/hiresw.%s/hiresw.t%sz.%s_2p5km.f%s.conus.grib2.idx" % (yrmnthday,hr,core,str(fhour).zfill(2))
		outfilgrib = '/gpfs_backup/stormtrack/jtradfor/ensemble_data/rawdata/href/%s_%s.grib2' % (outfil, core)
		print("%s %s | egrep %s | %s %s %s" % (get_inv, idxurl, fields, get_grib, griburl, outfilgrib))
		sp.call("%s %s | egrep %s | %s %s %s" % (get_inv, idxurl, fields, get_grib, griburl, outfilgrib),shell=True)

		griburlback = "https://www.ftp.ncep.noaa.gov/data/nccf/com/hiresw/prod/hiresw.%s/hiresw.t%sz.%s_2p5km.f%s.conus.grib2" % (backyrmnthday,backhr,core,str(fhour+12).zfill(2))
		idxurlback = "https://www.ftp.ncep.noaa.gov/data/nccf/com/hiresw/prod/hiresw.%s/hiresw.t%sz.%s_2p5km.f%s.conus.grib2.idx" % (backyrmnthday,backhr,core,str(fhour+12).zfill(2))
		outfilgribback = '/gpfs_backup/stormtrack/jtradfor/ensemble_data/rawdata/href/%s_%s_back.grib2' % (outfil, core) 
		sp.call("%s %s | egrep %s | %s %s %s" % (get_inv, idxurlback, fields, get_grib, griburlback, outfilgribback),shell=True)

	namnestgriburl = "https://www.ftp.ncep.noaa.gov/data/nccf/com/nam/prod/nam.%s/nam.t%sz.conusnest.hiresf%s.tm00.grib2" % (yrmnthday,hr,str(fhour).zfill(2))
	namnestidxurl = "https://www.ftp.ncep.noaa.gov/data/nccf/com/nam/prod/nam.%s/nam.t%sz.conusnest.hiresf%s.tm00.grib2.idx" % (yrmnthday,hr,str(fhour).zfill(2))
	namnestoutfilgrib = '/gpfs_backup/stormtrack/jtradfor/ensemble_data/rawdata/href/%s_namnest.grib2' % (outfil)
	sp.call("%s %s | egrep %s | %s %s %s" % (get_inv, namnestidxurl, fields, get_grib, namnestgriburl, namnestoutfilgrib),shell=True)

	namnestgriburlback = "https://www.ftp.ncep.noaa.gov/data/nccf/com/nam/prod/nam.%s/nam.t%sz.conusnest.hiresf%s.tm00.grib2" % (backyrmnthday,backhr,str(fhour+12).zfill(2))
	namnestidxurlback = "https://www.ftp.ncep.noaa.gov/data/nccf/com/nam/prod/nam.%s/nam.t%sz.conusnest.hiresf%s.tm00.grib2.idx" % (backyrmnthday,backhr,str(fhour+12).zfill(2))
	namnestoutfilgribback = '/gpfs_backup/stormtrack/jtradfor/ensemble_data/rawdata/href/%s_namnest_back.grib2' % (outfil)
	sp.call("%s %s | egrep %s | %s %s %s" % (get_inv, namnestidxurlback, fields, get_grib, namnestgriburlback, namnestoutfilgribback),shell=True)

	hrrrgriburl = "https://www.ftp.ncep.noaa.gov/data/nccf/com/hrrr/prod/hrrr.%s/conus/hrrr.t%sz.wrfsfcf%s.grib2" % (yrmnthday,hr,str(fhour).zfill(2))
	hrrridxurl = "https://www.ftp.ncep.noaa.gov/data/nccf/com/hrrr/prod/hrrr.%s/conus/hrrr.t%sz.wrfsfcf%s.grib2.idx" % (yrmnthday,hr,str(fhour).zfill(2))		
	hrrroutfilgrib = '/gpfs_backup/stormtrack/jtradfor/ensemble_data/rawdata/href/%s_hrrr.grib2' % (outfil)
	sp.call("%s %s | egrep %s | %s %s %s" % (get_inv, hrrridxurl, hrrrfields, get_grib, hrrrgriburl, hrrroutfilgrib),shell=True)

	hrrrgriburlback = "https://www.ftp.ncep.noaa.gov/data/nccf/com/hrrr/prod/hrrr.%s/conus/hrrr.t%sz.wrfsfcf%s.grib2" % (back6yrmnthday,back6hr,str(fhour+6).zfill(2))
	hrrridxurlback = "https://www.ftp.ncep.noaa.gov/data/nccf/com/hrrr/prod/hrrr.%s/conus/hrrr.t%sz.wrfsfcf%s.grib2.idx" % (back6yrmnthday,back6hr,str(fhour+6).zfill(2))		
	hrrroutfilgribback = '/gpfs_backup/stormtrack/jtradfor/ensemble_data/rawdata/href/%s_hrrr_back.grib2' % (outfil)
	sp.call("%s %s | egrep %s | %s %s %s" % (get_inv, hrrridxurlback, hrrrfields, get_grib, hrrrgriburlback, hrrroutfilgribback),shell=True)
	print hrrrgriburlback
	print hrrrgriburl



sp.call("sed -i -- 's,set datestr=.*,set datestr=\"%s\",g' /gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/dlscripts/href/*interp*.csh" % (newdate),shell=True)
sp.call('git pull',shell=True)
sp.call('git push origin master',shell=True)
time.sleep(5)

sp.call('bsub < /gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/dlscripts/href/namnest_interp.csh',shell=True)
sp.call('bsub < /gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/dlscripts/href/namnest_back_interp.csh',shell=True)
sp.call('bsub < /gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/dlscripts/href/arw_interp.csh',shell=True)
sp.call('bsub < /gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/dlscripts/href/arw_back_interp.csh',shell=True)
sp.call('bsub < /gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/dlscripts/href/nmmb_interp.csh',shell=True)
sp.call('bsub < /gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/dlscripts/href/nmmb_back_interp.csh',shell=True)
sp.call('bsub < /gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/dlscripts/href/hrrr_interp.csh',shell=True)
sp.call('bsub < /gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/dlscripts/href/hrrr_back_interp.csh',shell=True)
print datetime.now()
time.sleep(4200)
print datetime.now()
sp.call('/usr/bin/rm /gpfs_backup/stormtrack/jtradfor/ensemble_data/rawdata/href/*.grib2*',shell=True)

sp.call('/usr/bin/mv /gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/uploads/outimages/href/*.png /gpfs_backup/stormtrack/jtradfor/ensemble_data/old_href/' ,shell=True)

sp.call("sed -i -- 's,set datestr=.*,set datestr=\"%s\",g' /gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/dlscripts/href/*_r.csh" % (newdate),shell=True)
time.sleep(5)

sp.call('bsub < /gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/dlscripts/href/apcp_mean_r.csh',shell=True)
sp.call('bsub < /gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/dlscripts/href/apcp_postage_r.csh',shell=True)
sp.call('bsub < /gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/dlscripts/href/apcp_ab_r.csh',shell=True)

sp.call('bsub < /gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/dlscripts/href/hpcp_mean_r.csh',shell=True)
sp.call('bsub < /gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/dlscripts/href/hpcp_postage_r.csh',shell=True)
sp.call('bsub < /gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/dlscripts/href/hpcp_ab_r.csh',shell=True)

sp.call('bsub < /gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/dlscripts/href/aweasd_mean_r.csh',shell=True)
sp.call('bsub < /gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/dlscripts/href/aweasd_postage_r.csh',shell=True)
sp.call('bsub < /gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/dlscripts/href/aweasd_ab_r.csh',shell=True)

sp.call('bsub < /gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/dlscripts/href/hweasd_mean_r.csh',shell=True)
sp.call('bsub < /gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/dlscripts/href/hweasd_postage_r.csh',shell=True)
sp.call('bsub < /gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/dlscripts/href/hweasd_ab_r.csh',shell=True)

sp.call('bsub < /gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/dlscripts/href/R_mean_r.csh',shell=True)
sp.call('bsub < /gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/dlscripts/href/R_postage_r.csh',shell=True)
sp.call('bsub < /gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/dlscripts/href/R_ab_r.csh',shell=True)
sp.call('bsub < /gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/dlscripts/href/R_pb_r.csh',shell=True)
sp.call('bsub < /gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/dlscripts/href/R_pbab_r.csh',shell=True)

sp.call('bsub < /gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/dlscripts/href/BR_mean_r.csh',shell=True)
sp.call('bsub < /gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/dlscripts/href/BR_postage_r.csh',shell=True)
sp.call('bsub < /gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/dlscripts/href/BR_ab_r.csh',shell=True)
sp.call('bsub < /gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/dlscripts/href/BR_pb_r.csh',shell=True)
sp.call('bsub < /gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/dlscripts/href/BR_pbab_r.csh',shell=True)

sp.call('bsub < /gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/dlscripts/href/band_detection_mean_r.csh',shell=True)
sp.call('bsub < /gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/dlscripts/href/band_detection_postage_r.csh',shell=True)
sp.call('bsub < /gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/dlscripts/href/band_detection_probability_r.csh',shell=True)
sp.call('bsub < /gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/dlscripts/href/band_detection_probability2_r.csh',shell=True)

