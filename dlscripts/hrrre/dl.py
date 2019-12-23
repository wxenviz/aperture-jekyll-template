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
import time

from mpl_toolkits.basemap import Basemap
from matplotlib.patches import Polygon
from matplotlib.colors import LinearSegmentedColormap
from scipy.spatial import Delaunay
from scipy.interpolate import LinearNDInterpolator
from shutil import copyfile


existsinall = np.load('/gpfs_backup/stormtrack/jtradfor/ensemble_data/rawdata/hrrre/existfils.npy')
print(existsinall)
if len(existsinall)>0:
	ftp = ftplib.FTP('gsdftp.fsl.noaa.gov', 'anonymous', 'jtradfor@ncsu.edu')
	fil = existsinall[0]
	j=0
	for i in range(1,10):
	    while j<5:
	        try:
	    	    ftp.retrbinary("RETR /hrrre/conus/mem000" + str(i) + "/wrftwo/" + fil, open("/gpfs_backup/stormtrack/jtradfor/ensemble_data/rawdata/hrrre/" + fil + "_" + str(i),'wb').write)
	            break
		except:
		    print('trying again')
		    time.sleep(60)
		    ftp.retrbinary("RETR /hrrre/conus/mem000" + str(i) + "/wrftwo/" + fil, open("/gpfs_backup/stormtrack/jtradfor/ensemble_data/rawdata/hrrre/" + fil + "_" + str(i),'wb').write)
		    j+=1
	ftp.close()

