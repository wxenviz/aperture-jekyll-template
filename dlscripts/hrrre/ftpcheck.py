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

os.chdir('/gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/')
#sys.stdout = open('/gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/dlscripts/hrrre/logfile','w')

sp.call('git pull',shell=True)

#Get latest downloaded initialization time, if it exists
try:
    dlfils = glob.glob('/gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/uploads/outimages/hrrre/*R_mean.png')
    latestfil = os.path.basename(dlfils[0])
    currentdate = datetime.datetime.strptime(latestfil[0:7],'%y%j%H')
except:
    currentdate = datetime.datetime.strptime('1500100','%y%j%H')

#Open FTP
ftp = ftplib.FTP('gsdftp.fsl.noaa.gov', 'anonymous', 'jtradfor@ncsu.edu')
ftp.cwd("/hrrre/conus/mem0001/wrftwo")
fils = []
ftp.retrlines('NLST', fils.append)

#Get latest FTP intialization time
ftplatest = datetime.datetime.strptime('1500100','%y%j%H')
ftplateststr = ""
for fil in fils:
    ftpcurrent = datetime.datetime.strptime(fil[0:7],'%y%j%H')
    if ftpcurrent>ftplatest:
        ftplatest = ftpcurrent
        ftplateststr = fil[0:7]
    newdate = fil[0:9]
#If it's a new initialization, delete old images
if ftplatest > currentdate:
    try:
        imagefils = glob.glob('/gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/uploads/outimages/hrrre/*.png')
        rawfils = glob.glob('/gpfs_backup/stormtrack/jtradfor/ensemble_data/rawdata/hrrre/*')
        with open('/gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/hrrre.html') as f:
            s = f.read()
        with open('/gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/hrrre.html', 'w') as f:
        	s = re.sub(r'datei = \"\d{9}',r'datei = "' + newdate,s)
        	f.write(s)
        	f.close()
        	sp.call('git add -A',shell=True)
        	sp.call("git commit -m 'update html'",shell=True)
        	sp.call('git push origin master',shell=True)

        if len(imagefils)>0:
           sp.call('mv /gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/uploads/outimages/hrrre/*.png /gpfs_backup/stormtrack/jtradfor/ensemble_data/old_outimages/',shell=True)
    	   sp.call('git add -A',shell=True)
           sp.call("git commit -m 'removed old images'",shell=True)
    	   sp.call('git push origin master',shell=True)

        if len(rawfils)>0:
    	   sp.call('rm /gpfs_backup/stormtrack/jtradfor/ensemble_data/rawdata/hrrre/*',shell=True)
    except: 
        pass
    dlfils = []
    
#Get all files for all members
allmembers = []
for i in range(1,10):
    member = "mem000" + str(i)
    ftp.cwd("/hrrre/conus/" + member + "/wrftwo")
    fils = []
    ftp.retrlines('NLST', fils.append)
    allmembers.append(fils)
existsinall = []
notexistsinall = []
#Go through member 1 files
for i in range(0,9):
    print(allmembers[i][-1])
for fil in allmembers[0]:
    #If file is from latest initialization and I haven't produced an image from it,
    #check if the file exists for all members. If it does, add to download list
    if '/gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/uploads/outimages/hrrre/' + fil + '_R_mean.png' not in dlfils and fil[:7]==ftplateststr:
	check = 0
        for i in range(1,9):
            if fil not in allmembers[i]:
                check = 1
        if check==0:
            existsinall.append(fil)
        else:
            notexistsinall.append(fil)

    else:
        continue
ftp.close()
existsinall = np.array(existsinall)
if len(existsinall)>0:
    np.save('/gpfs_backup/stormtrack/jtradfor/ensemble_data/rawdata/hrrre/existfils.npy',existsinall)

