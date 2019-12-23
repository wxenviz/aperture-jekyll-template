#!/usr/bin/env python
# coding: utf-8

# In[15]:


import subprocess as sp
import pygrib as pg
import numpy as np
import netCDF4 as nc4
import csv
import cv2
import os
import re

from geopy.distance import great_circle
from mpl_toolkits.basemap import Basemap, interp
from datetime import datetime, timedelta
from matplotlib.colors import LinearSegmentedColormap
from scipy import misc
from mpl_toolkits.basemap import Basemap
from matplotlib.patches import Polygon
import matplotlib.pyplot as plt
plt.switch_backend('agg')

sp.call("wget -P /gpfs_backup/stormtrack/jtradfor/ensemble_data/radardata/ https://mesonet.agron.iastate.edu/archive/data/2019/11/30/GIS/uscomp/n0q_201911301000.png",shell=True)

sp.call("wget -P /gpfs_backup/stormtrack/jtradfor/ensemble_data/radardata/ https://pando-rgw01.chpc.utah.edu/hrrr/sfc/20191130/hrrr.t10z.wrfsfcf00.grib2",shell=True)
grbs = pg.open('/gpfs_backup/stormtrack/jtradfor/ensemble_data/radardata/hrrr.t10z.wrfsfcf00.grib2')
temps = np.flipud(grbs[59].values[:])

temps[temps<=273.15] = 1
temps[temps>273.15]=0
    
    
curve = nc4.Dataset('/gpfs_backup/stormtrack/jtradfor/ensemble_data/reference/master.nc','r')
lon_curve = np.flipud(curve.variables['gridlon_0'][:,:])
lat_curve = np.flipud(curve.variables['gridlat_0'][:,:])

levels = []
colors = []
with open('/gpfs_backup/stormtrack/jtradfor/ensemble_data/reference/N0Q_Color_Lookup.csv','r') as colorcsv:
    colorreader = csv.reader(colorcsv,delimiter=',')
    for line in colorreader:
	if float(line[1])>=0 and float(line[1])<=60:
		colorints = [int(i) for i in line[2:]]
		colors.append((colorints))
		levels.append(float(line[1]))  
colors = np.array(colors)/255.0
cmap1 = LinearSegmentedColormap.from_list("my_colormap",colors,N=len(levels),gamma=1.0)

fildate = '201911301000'
fil = '/gpfs_backup/stormtrack/jtradfor/ensemble_data/radardata/n0q_%s.png' % (fildate)
lats = np.arange(23,50,0.005)
lons = np.arange(-126,-65,0.005)

outfilobj = datetime.strptime(fildate,'%Y%m%d%H%M')
outfil = outfilobj.strftime('%y%j%H%M')

reflectivity = np.flipud(np.array(misc.imread(fil,mode='P')).astype(float))
reflectivity = (reflectivity - 65) * .5
reflect_interp = interp(reflectivity,lons,lats,lon_curve,lat_curve,checkbounds=False,masked=True,order=1)
reflect_interp = reflect_interp.filled(fill_value=-32.5)
reflectivity = .5*np.round(reflect_interp/.5)

reflectivity_contour = np.copy(reflectivity)
kernel = np.ones((5,5),np.float32)/25.
reflectivity_contour = cv2.filter2D(reflectivity_contour,-1,kernel)

reflectivity_contour[reflectivity_contour>0]=255
reflectivity_contour[reflectivity_contour<0]=0
reflectivity_contour = np.array(reflectivity_contour).astype('uint8')

contours, _= cv2.findContours(reflectivity_contour,cv2.RETR_EXTERNAL,cv2.CHAIN_APPROX_SIMPLE)

overlay = np.zeros_like(lat_curve)
for contour in contours:

    if cv2.contourArea(contour)>=100 and cv2.contourArea(contour)<250000:

        mask = np.zeros_like(reflectivity).astype('uint8')
        cv2.drawContours(mask,[contour],0,255,-1)
        pixelpoints = cv2.findNonZero(mask)
        pixelpoints = pixelpoints.transpose()
        mask[temps==0] = 0
        if np.max(mask)==0:
            continue
        else:
            snowpixelpoints = cv2.findNonZero(mask)
            snowpixelpoints = snowpixelpoints.transpose()
            threshold = np.mean(reflectivity[snowpixelpoints[1],snowpixelpoints[0]]) + 1.25*np.std(reflectivity[snowpixelpoints[1],snowpixelpoints[0]])
            if threshold<15:
                threshold=15.0
            reflect_contour_2 = np.zeros_like(reflectivity)
            reflect_contour_2[pixelpoints[1],pixelpoints[0]] = reflectivity[pixelpoints[1],pixelpoints[0]]
            reflect_contour_2[reflect_contour_2>=threshold] = 255
            reflect_contour_2[reflect_contour_2<threshold] = 0
            reflect_contour_2 = np.array(reflect_contour_2).astype('uint8')

            interior_contours, _= cv2.findContours(reflect_contour_2,cv2.RETR_EXTERNAL,cv2.CHAIN_APPROX_SIMPLE)
            for interior_contour in interior_contours:
                if cv2.contourArea(interior_contour)>=100:
                    M = cv2.moments(interior_contour)
                    cx = int(M['m10']/M['m00'])
                    cy = int(M['m01']/M['m00'])
                    loncenter = lon_curve[cy,cx]
                    latcenter = lat_curve[cy,cx]

                    interior_mask = np.zeros_like(reflectivity).astype('uint8')
                    cv2.drawContours(interior_mask,[interior_contour],0,255,-1)
                    interiorpixelpoints = cv2.findNonZero(interior_mask)
                    interiorpixelpoints = interiorpixelpoints.transpose()

                    rect = cv2.minAreaRect(interior_contour)
                    angle = rect[2]
		    rectarea = rect[1][0]*rect[1][1]
                    box = cv2.boxPoints(rect)
                    box = np.int0(box)
                    avg1x = int(np.mean([box[0][0],box[1][0]]))
                    avg1y = int(np.mean([box[0][1],box[1][1]]))
                    avg2x = int(np.mean([box[1][0],box[2][0]]))
                    avg2y = int(np.mean([box[1][1],box[2][1]]))
                    avg3x = int(np.mean([box[2][0],box[3][0]]))
                    avg3y = int(np.mean([box[2][1],box[3][1]]))
                    avg4x = int(np.mean([box[3][0],box[0][0]]))
                    avg4y = int(np.mean([box[3][1],box[0][1]]))


                    width = np.min([great_circle((lat_curve[avg1y,avg1x],lon_curve[avg1y,avg1x]), (lat_curve[avg3y,avg3x],lon_curve[avg3y,avg3x])).kilometers,great_circle((lat_curve[avg2y,avg2x],lon_curve[avg2y,avg2x]), (lat_curve[avg4y,avg4x],lon_curve[avg4y,avg4x])).kilometers])
                    length = np.max([great_circle((lat_curve[avg1y,avg1x],lon_curve[avg1y,avg1x]), (lat_curve[avg3y,avg3x],lon_curve[avg3y,avg3x])).kilometers,great_circle((lat_curve[avg2y,avg2x],lon_curve[avg2y,avg2x]), (lat_curve[avg4y,avg4x],lon_curve[avg4y,avg4x])).kilometers])

                    if great_circle((lat_curve[avg1y,avg1x],lon_curve[avg1y,avg1x]), (lat_curve[avg3y,avg3x],lon_curve[avg3y,avg3x])).kilometers > great_circle((lat_curve[avg2y,avg2x],lon_curve[avg2y,avg2x]), (lat_curve[avg4y,avg4x],lon_curve[avg4y,avg4x])).kilometers:
                        endpt1lat = lat_curve[avg1y,avg1x]
                        endpt1lon = lon_curve[avg1y,avg1x]
                        endpt2lat = lat_curve[avg3y,avg3x]
                        endpt2lon = lon_curve[avg3y,avg3x]
                    else:
                        endpt1lat = lat_curve[avg2y,avg2x]
                        endpt1lon = lon_curve[avg2y,avg2x]
                        endpt2lat = lat_curve[avg4y,avg4x]
                        endpt2lon = lon_curve[avg4y,avg4x]
                        
                        
                    coldpts = np.sum(temps[interiorpixelpoints[1],interiorpixelpoints[0]])
                    allpts = len(interiorpixelpoints[0][0])
		    reflect_avg = np.mean(reflectivity[interiorpixelpoints[1],interiorpixelpoints[0]])
                    
		    if (length>=250 and length/width>=3.0 and coldpts/allpts>=0.50):
			print(threshold)
			solidity = (cv2.contourArea(interior_contour)/rectarea)
			realrectarea = length*width
			realarea = realrectarea*solidity
                        overlay[interiorpixelpoints[1],interiorpixelpoints[0]] = 1

      
fig1 = plt.figure(figsize=(16,9))
m = Basemap(projection='lcc',lat_0=5,lon_0=-100,llcrnrlon=-126,llcrnrlat=23,urcrnrlon=-65,urcrnrlat=50,resolution='h')
shp_info = m.readshapefile('/gpfs_backup/stormtrack/jtradfor/ensemble_data/reference/st99_d00','states',drawbounds=False)
X,Y = m(lon_curve,lat_curve)
ax = plt.gca()
for nshape,seg in enumerate(m.states):
   poly = Polygon(seg,facecolor='black',edgecolor='black',zorder=1,linewidth=1)
   poly2 = Polygon(seg,facecolor='none',edgecolor='#EFE2CB',zorder=3,linewidth=1)
   ax.add_patch(poly)
   ax.add_patch(poly2)

reflectivity[reflectivity<=0.0]=np.nan
im = m.contourf(X,Y,reflectivity,levels=levels,cmap=cmap1,zorder=2)
#m.contour(X,Y,overlay,colors='orange',linewidths=0.55,zorder=4)
#cbar = plt.colorbar(im,orientation="horizontal",pad=0.0,fraction=0.05,ticks=[0,10,20,30,40,50,60])
#cbar.ax.yaxis.set_tick_params(color='black')
#cbar.ax.set_yticklabels([0,10,20,30,40,50,60],color='black')
plt.box(False)

plt.savefig('/gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/uploads/outimages/radar/radar_%s_test.png' % (outfil),facecolor='#EFE2CB',bbox_inches='tight',dpi=500)
plt.close()

