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


def ens_BRmean(reflectivities,outname):
   
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
 
    plt.figure(figsize=(16,9))

    m = Basemap(projection='lcc',lat_0=5,lon_0=-100,llcrnrlon=-126,llcrnrlat=23,urcrnrlon=-63,urcrnrlat=50,resolution='h')
    shp_info = m.readshapefile('/gpfs_backup/stormtrack/jtradfor/ensemble_data/reference/st99_d00','states',drawbounds=False)
    
    ax = plt.gca()
    
    for nshape,seg in enumerate(m.states):
        poly = Polygon(seg,facecolor='white',edgecolor='white',zorder=1,linewidth=1)
	poly2 = Polygon(seg,facecolor='none',edgecolor='black',zorder=3,linewidth=1)
        ax.add_patch(poly)
	ax.add_patch(poly2)

    reflectivities_copy = np.copy(reflectivities) 
    reflect_mean = np.mean(reflectivities_copy,axis=0)
    reflect_mean[reflect_mean<-9.5] = np.nan
    reflect_mean[-50:,:] = np.nan
    reflect_mean[:50,:] = np.nan
    reflect_mean[:,:50] = np.nan
    reflect_mean[:,-50:] = np.nan
    m.imshow(reflect_mean,zorder=2,aspect='equal',interpolation='none',cmap=cmap1,vmin=-32.5,vmax=95.)
    cbar = plt.colorbar(fraction=0.023,pad=-0.02)
    cbar.ax.yaxis.set_tick_params(color='w')
    cbar.ax.set_yticklabels([-20,0,20,40,60,80],color='w')
    plt.box(False)
    meanfil = '/gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/uploads/outimages/hrrre/' + outname + '_BR_mean.png'
    plt.savefig(meanfil,facecolor='#101010',bbox_inches='tight',dpi=500)
    plt.close()


datapaths = glob.glob('/gpfs_backup/stormtrack/jtradfor/ensemble_data/rawdata/hrrre/*breflect.npy')
dateis = []

if len(datapaths)>0:

	latest = 0
	latestpath = datapaths[0]
	for datapath in datapaths:
	    print datapath
	    if int(os.path.basename(datapath)[9:11]) > latest:
		latest = int(os.path.basename(datapath)[9:11]) 
		latestpath = datapath

	fil = os.path.basename(latestpath)[0:13]
	reflectivities = np.load(latestpath)
	pbsdfil = ens_BRmean(reflectivities,fil)

