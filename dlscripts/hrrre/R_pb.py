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
sys.stdout = open('/gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/logfile','a+')
print 'starting rpb'

os.chdir('/gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/')


def pb(reflectivities,outname):
    
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
    reflect_pb = np.zeros_like(reflectivities_copy[0])
    for c,reflect_member in enumerate(reflectivities_copy):
        reflect_member[reflect_member<40] = 0
        reflect_member[reflect_member>=40] = c+1
        reflect_pb = np.max([reflect_pb,reflect_member],axis=0)
        
    reflect_pb[reflect_pb==0] = np.nan
    reflect_pb[-50:,:] = np.nan
    reflect_pb[:50,:] = np.nan
    reflect_pb[:,:50] = np.nan
    reflect_pb[:,-50:] = np.nan
        
    m.imshow(reflect_pb,zorder=2,cmap='tab10',interpolation='none',vmin=1,vmax=10)
    plt.box(False)
    pbfil = '/gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/uploads/outimages/hrrre/' + outname + '_R_pb.png'
    plt.savefig(pbfil,facecolor='#101010',bbox_inches='tight',dpi=500)
    plt.close()

datapaths = glob.glob('/gpfs_backup/stormtrack/jtradfor/ensemble_data/rawdata/hrrre/*creflect*')
dateis = []

if len(datapaths)>0:

	latest = 0
	latestpath = datapaths[0]
	for datapath in datapaths:
	    if int(os.path.basename(datapath)[9:11]) > latest:
		latest = int(os.path.basename(datapath)[9:11]) 
		latestpath = datapath

	fil = os.path.basename(latestpath)[0:13]
	reflectivities = np.load(latestpath)
	pbsdfil = pb(reflectivities,fil)
