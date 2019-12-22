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

#sys.stdout = open('/gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/logfile','a+')

os.chdir('/gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/')

def hpcp_ab(totalPrecips,totalPrecipsback,outname):

    plt.figure(figsize=(16,9))

    m = Basemap(projection='lcc',lat_0=5,lon_0=-100,llcrnrlon=-126,llcrnrlat=23,urcrnrlon=-63,urcrnrlat=50,resolution='h')
    shp_info = m.readshapefile('/gpfs_backup/stormtrack/jtradfor/ensemble_data/reference/st99_d00','states',drawbounds=False)
    
    ax = plt.gca()
    
    for nshape,seg in enumerate(m.states):
        poly = Polygon(seg,facecolor='white',edgecolor='white',zorder=1,linewidth=1)
	poly2 = Polygon(seg,facecolor='none',edgecolor='black',zorder=3,linewidth=1)
        ax.add_patch(poly)
	ax.add_patch(poly2)

    totalPrecips_copy = np.copy(totalPrecips)
    totalPrecipsback_copy = np.copy(totalPrecipsback)    
    
    totalPrecips_ab = np.zeros_like(totalPrecips_copy[0])
    
    for i in range(0,len(totalPrecips_ab)):
        for j in range(0,len(totalPrecips_ab[0])):
            c = 0
            for k in range(0,3):
                for l in range(0,3):
                    if i%3==k and j%3==l:
                        totalPrecips_ab[i,j] = totalPrecips_copy[c][i,j] - totalPrecipsback_copy[c][i,j]
                    c+=1
                
    totalPrecips_ab[totalPrecips_ab<0.25] = np.nan
    totalPrecips_ab[totalPrecips_ab>25.0] = 25.0
    totalPrecips_ab[-50:,:] = np.nan
    totalPrecips_ab[:50,:] = np.nan
    totalPrecips_ab[:,:50] = np.nan
    totalPrecips_ab[:,-50:] = np.nan
    m.imshow(totalPrecips_ab,aspect='equal',interpolation='none',cmap='jet',vmin=0.0,vmax=25.0)
    plt.colorbar(fraction=0.023,pad=-0.02)
    plt.box(False)
    sdfil = '/gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/uploads/outimages/hrrre/' + outname + '_hpcp_ab.png'
    plt.savefig(sdfil,facecolor='#101010',bbox_inches='tight',dpi=800)
    plt.close()

datapaths = glob.glob('/gpfs_backup/stormtrack/jtradfor/ensemble_data/rawdata/hrrre/*apcp*')

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
	totalPrecips = np.load(latestpath)
	filback = os.path.basename(back1hrpath)[:13]
	totalPrecipsback = np.load(back1hrpath)
	hpcp_ab_fil = hpcp_ab(totalPrecips,totalPrecipsback,fil)

