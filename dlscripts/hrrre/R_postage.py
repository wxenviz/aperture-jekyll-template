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

def postage(reflectivities,outname):
   
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
 
    fig, axes = plt.subplots(nrows=3,ncols=3,figsize=(16,9))
    plt.subplots_adjust(wspace=0.05,hspace=-0.10)
    for i in range(0,9):
	
	ax = axes.flat[i]
        m = Basemap(projection='lcc',lat_0=5,lon_0=-100,llcrnrlon=-126,llcrnrlat=23,urcrnrlon=-63,urcrnrlat=50,resolution='l',ax=ax)
        shp_info = m.readshapefile('/gpfs_backup/stormtrack/jtradfor/ensemble_data/reference/st99_d00','states',drawbounds=False)
    
        for nshape,seg in enumerate(m.states):
            poly = Polygon(seg,facecolor='white',edgecolor='white',zorder=1,linewidth=1)
	    poly2 = Polygon(seg,facecolor='none',edgecolor='black',zorder=3,linewidth=1)
            ax.add_patch(poly)
	    ax.add_patch(poly2)
	
        reflectivities[i][reflectivities[i]<-9.5] = np.nan
        reflectivities[i][-50:,:] = np.nan
        reflectivities[i][:50,:] = np.nan
        reflectivities[i][:,:50] = np.nan
        reflectivities[i][:,-50:] = np.nan
        im = m.imshow(reflectivities[i],zorder=2,aspect='equal',interpolation='none',cmap=cmap1,vmin=-32.5,vmax=95.)
    fig.colorbar(im, ax=axes.ravel().tolist(),fraction=0.025)
    postagefil = '/gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/uploads/outimages/hrrre/' + outname + '_R_postage.png'
    plt.savefig(postagefil,facecolor='#101010',bbox_inches='tight',dpi=500)
    plt.close()

datapaths = glob.glob('/gpfs_backup/stormtrack/jtradfor/ensemble_data/rawdata/hrrre/*creflect*')

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


f = open("/gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/dlscripts/hrrre/unlock/Rpostage.unlock","w+")
