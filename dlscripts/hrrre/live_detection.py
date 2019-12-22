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

os.chdir('/gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/')
sp.call('git pull',shell=True)


curtime = datetime.utcnow()

rndtime = curtime - timedelta(minutes=curtime.minute % 15,
                             seconds=curtime.second,
                             microseconds=curtime.microsecond)
year = str(rndtime.year)
month = str(rndtime.month).zfill(2)
day = str(rndtime.day).zfill(2)
hour = str(rndtime.hour).zfill(2)
minute = str(rndtime.minute).zfill(2)
fildate = year + month + day + hour + minute

hrrrtime = rndtime - timedelta(hours = 1)
hrrryear = str(hrrrtime.year)
hrrrmonth = str(hrrrtime.month).zfill(2)
hrrrday = str(hrrrtime.day).zfill(2)
hrrrhour = str(hrrrtime.hour).zfill(2)
hrrryrmnthday = hrrryear + hrrrmonth + hrrrday
hrrrfloor = datetime.strptime(hrrryear + hrrrmonth + hrrrday + hrrrhour,'%Y%m%d%H')
hrrrflooryear = str(hrrrfloor.year)
hrrrfloormonth = str(hrrrfloor.month).zfill(2)
hrrrfloorday = str(hrrrfloor.day).zfill(2)
hrrrfloorhour = str(hrrrfloor.hour).zfill(2)
hrrrflooryrmnthday = hrrrflooryear + hrrrfloormonth + hrrrfloorday

try:
    sp.call("wget -P /gpfs_backup/stormtrack/jtradfor/ensemble_data/radardata/ https://mesonet.agron.iastate.edu/archive/data/%s/%s/%s/GIS/uscomp/n0q_%s.png" % (year,month,day,fildate),shell=True)
except:
    rndtime = curtime - timedelta(minutes=15) - timedelta(minutes=curtime.minute % 15,
                             seconds=curtime.second,
                             microseconds=curtime.microsecond)
    year = str(rndtime.year)
    month = str(rndtime.month).zfill(2)
    day = str(rndtime.day).zfill(2)
    hour = str(rndtime.hour).zfill(2)
    minute = str(rndtime.minute).zfill(2)
    fildate = year + month + day + hour + minute
    
    hrrrtime = rndtime - timedelta(hours = 1)
    hrrryear = str(hrrrtime.year)
    hrrrmonth = str(hrrrtime.month).zfill(2)
    hrrrday = str(hrrrtime.day).zfill(2)
    hrrrhour = str(hrrrtime.hour).zfill(2)
    hrrryrmnthday = hrrryear + hrrrmonth + hrrrday
    hrrrfloor = datetime.strptime(hrrryear + hrrrmonth + hrrrday + hrrrhour,'%Y%m%d%H')
    hrrrflooryear = str(hrrrfloor.year)
    hrrrfloormonth = str(hrrrfloor.month).zfill(2)
    hrrrfloorday = str(hrrrfloor.day).zfill(2)
    hrrrfloorhour = str(hrrrfloor.hour).zfill(2)
    hrrrflooryrmnthday = hrrrflooryear + hrrrfloormonth + hrrrfloorday
    
    sp.call("wget -P /gpfs_backup/stormtrack/jtradfor/ensemble_data/radardata/ https://mesonet.agron.iastate.edu/archive/data/%s/%s/%s/GIS/uscomp/n0q_%s.png" % (year,month,day,fildate),shell=True)

subhr = int(((rndtime-hrrrfloor).seconds)/60.)
if (subhr==60):
    sp.call("wget -P /gpfs_backup/stormtrack/jtradfor/ensemble_data/radardata/ https://www.ftp.ncep.noaa.gov/data/nccf/com/hrrr/prod/hrrr.%s/conus/hrrr.t%sz.wrfsfcf01.grib2" % (hrrrflooryrmnthday,hrrrfloorhour),shell=True)
    grbs = pg.open('/gpfs_backup/stormtrack/jtradfor/ensemble_data/radardata/hrrr.t%sz.wrfsfcf01.grib2' % (hrrrfloorhour))
    temps = np.flipud(grbs[66].values[:])
elif (subhr==75):
    sp.call("wget -P /gpfs_backup/stormtrack/jtradfor/ensemble_data/radardata/ https://www.ftp.ncep.noaa.gov/data/nccf/com/hrrr/prod/hrrr.%s/conus/hrrr.t%sz.wrfsubhf02.grib2" % (hrrrflooryrmnthday,hrrrfloorhour),shell=True)
    grbs = pg.open('/gpfs_backup/stormtrack/jtradfor/ensemble_data/radardata/hrrr.t%sz.wrfsubhf02.grib2' % (hrrrfloorhour))
    temps = np.flipud(grbs[13].values[:])
elif (subhr==90):
    sp.call("wget -P /gpfs_backup/stormtrack/jtradfor/ensemble_data/radardata/ https://www.ftp.ncep.noaa.gov/data/nccf/com/hrrr/prod/hrrr.%s/conus/hrrr.t%sz.wrfsubhf02.grib2" % (hrrrflooryrmnthday,hrrrfloorhour),shell=True)
    grbs = pg.open('/gpfs_backup/stormtrack/jtradfor/ensemble_data/radardata/hrrr.t%sz.wrfsubhf02.grib2' % (hrrrfloorhour))
    temps = np.flipud(grbs[56].values[:])
elif (subhr==105):
    sp.call("wget -P /gpfs_backup/stormtrack/jtradfor/ensemble_data/radardata/ https://www.ftp.ncep.noaa.gov/data/nccf/com/hrrr/prod/hrrr.%s/conus/hrrr.t%sz.wrfsubhf02.grib2" % (hrrrflooryrmnthday,hrrrfloorhour),shell=True)
    grbs = pg.open('/gpfs_backup/stormtrack/jtradfor/ensemble_data/radardata/hrrr.t%sz.wrfsubhf02.grib2' % (hrrrfloorhour))
    temps = np.flipud(grbs[99].values[:])
    

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

bandstats = []
try:
	with open('/gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/dlscripts/hrrre/bandstats.csv','r') as statsfile:
		reader = csv.reader(statsfile,delimiter=',')
		for line in reader:
			bandstats.append(line)
	statsfile.close()
except:
	pass


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

    if cv2.contourArea(contour)>=100:

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
			solidity = (cv2.contourArea(interior_contour)/rectarea)
			print(solidity)
			realrectarea = length*width
			realarea = realrectarea*solidity
			print(realarea)
                        overlay[interiorpixelpoints[1],interiorpixelpoints[0]] = 1
			bandstats.append([str(outfil),str(latcenter),str(loncenter),str(endpt1lat),str(endpt1lon),str(endpt2lat),str(endpt2lon),str(length),str(width),str(realarea),str(threshold),str(reflect_avg)])

      
plt.figure(figsize=(16,9))
m = Basemap(projection='lcc',lat_0=5,lon_0=-100,llcrnrlon=-126,llcrnrlat=23,urcrnrlon=-65,urcrnrlat=50,resolution='h')
shp_info = m.readshapefile('/gpfs_backup/stormtrack/jtradfor/ensemble_data/reference/st99_d00','states',drawbounds=False)
X,Y = m(lon_curve,lat_curve)
ax = plt.gca()
for nshape,seg in enumerate(m.states):
   poly = Polygon(seg,facecolor='white',edgecolor='white',zorder=1,linewidth=1)
   poly2 = Polygon(seg,facecolor='none',edgecolor='black',zorder=3,linewidth=1)
   ax.add_patch(poly)
   ax.add_patch(poly2)

reflectivity[reflectivity<=0.0]=np.nan
im = m.contourf(X,Y,reflectivity,levels=levels,cmap=cmap1,zorder=2)
m.contour(X,Y,overlay,colors='orange',linewidths=0.50,zorder=4)
cbar = plt.colorbar(im,fraction=0.023,ticks=[0,10,20,30,40,50,60])
cbar.ax.yaxis.set_tick_params(color='w')
cbar.ax.set_yticklabels([0,10,20,30,40,50,60],color='w')
plt.box(False)


with open('/gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/dlscripts/hrrre/bandstats.csv','wb') as statsfile:
	csv.writer(statsfile).writerows(bandstats)
statsfile.close()

sp.call('/usr/bin/mv /gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/uploads/outimages/radar/*.png /gpfs_backup/stormtrack/jtradfor/ensemble_data/old_radar/',shell=True)
plt.savefig('/gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/uploads/outimages/radar/radar_%s.png' % (outfil),facecolor='#101010',bbox_inches='tight',dpi=500)
plt.close()
sp.call('rm /gpfs_backup/stormtrack/jtradfor/ensemble_data/radardata/*',shell=True)

with open('/gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/radar.html') as f:
	s = f.read()
with open('/gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/index.html', 'w') as f:
      	s = re.sub(r'datei = \"\d{9}',r'datei = "' + outfil,s)
        f.write(s)
        f.close()

sp.call('git add -A',shell=True)
sp.call("git commit -m 'radar update'",shell=True)
sp.call("git push origin master",shell=True)
