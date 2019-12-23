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

#sys.stdout = open('/gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/logfile2','w')

os.chdir('/gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/')

def postage(apcp,outname):
   
    fig, axes = plt.subplots(nrows=3,ncols=3,figsize=(16,9))
    plt.subplots_adjust(wspace=0.05,hspace=-.10)
    for i in range(0,9):
	
	ax = axes.flat[i]
        m = Basemap(projection='lcc',lat_0=5,lon_0=-100,llcrnrlon=-126,llcrnrlat=23,urcrnrlon=-63,urcrnrlat=50,resolution='l',ax=ax)
        shp_info = m.readshapefile('/gpfs_backup/stormtrack/jtradfor/ensemble_data/reference/st99_d00','states',drawbounds=False)
    
        for nshape,seg in enumerate(m.states):
            poly = Polygon(seg,facecolor='white',edgecolor='white',zorder=1,linewidth=1)
	    poly2 = Polygon(seg,facecolor='none',edgecolor='black',zorder=3,linewidth=1)
            ax.add_patch(poly)
	    ax.add_patch(poly2)
	
        apcp[i][apcp[i]<0.25] = np.nan
	apcp[i][apcp[i]>100.] = 100.
        apcp[i][-50:,:] = np.nan
        apcp[i][:50,:] = np.nan
        apcp[i][:,:50] = np.nan
        apcp[i][:,-50:] = np.nan
        im = m.imshow(apcp[i],zorder=2,aspect='equal',interpolation='none',cmap='jet',vmin=0.0,vmax=100.0)

    fig.colorbar(im, ax=axes.ravel().tolist(),fraction=0.025)
    postagefil = '/gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/uploads/outimages/hrrre/' + outname + '_apcp_postage.png'
    plt.savefig(postagefil,facecolor='#101010',bbox_inches='tight',dpi=500)
    plt.close()

datapaths = glob.glob('/gpfs_backup/stormtrack/jtradfor/ensemble_data/rawdata/hrrre/*apcp*')

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
	pbsdfil = postage(reflectivities,fil)

f = open("/gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/dlscripts/hrrre/unlock/apcp_postage.unlock","w+")
