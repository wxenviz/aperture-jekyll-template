import subprocess as sp
import pygrib as pg
import numpy as np
import netCDF4 as nc4
import csv
import cv2
import os
import re
import glob

from geopy.distance import great_circle
from mpl_toolkits.basemap import Basemap, interp
from datetime import datetime, timedelta
from matplotlib.colors import LinearSegmentedColormap
from matplotlib import cm
from scipy import misc
from mpl_toolkits.basemap import Basemap
from matplotlib.patches import Polygon
import matplotlib.pyplot as plt
plt.switch_backend('agg')

os.chdir('/gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/')
sp.call('git pull',shell=True)

def detection(reflectivities,temps,outname):

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
		poly = Polygon(seg,facecolor='none',edgecolor='black',zorder=2,linewidth=1)
		ax.add_patch(poly)
	
	hrrr = nc4.Dataset('/gpfs_backup/stormtrack/jtradfor/ensemble_data/reference/master.nc','r')
	hrrr_lons = np.flipud(hrrr.variables['gridlon_0'][:,:]).astype('f4')
	hrrr_lats = np.flipud(hrrr.variables['gridlat_0'][:,:]).astype('f4')

	reflectivities_copy = np.copy(reflectivities)

	reflect_mean = np.mean(reflectivities_copy,axis=0)
	reflect_mean[reflect_mean<=-9.5] = np.nan
	reflect_mean[-50:,:] = np.nan
	reflect_mean[:50,:] = np.nan
	reflect_mean[:,:50] = np.nan
	reflect_mean[:,-50:] = np.nan

	X,Y = m(hrrr_lons,hrrr_lats)

	im = m.contourf(X,Y,reflect_mean,levels=levels,cmap=cmap1,zorder=1)

	
	for i,reflectivity_copy in enumerate(reflectivities_copy):

		temp = temps[i]
		temp[temp<=273.15] = 1
		temp[temp>273.15] = 0

		reflectivity_contour = np.copy(reflectivity_copy)
		reflectivity_contour[reflectivity_contour>0] = 255
		reflectivity_contour[reflectivity_contour<0] = 0
		reflectivity_contour = np.array(reflectivity_contour).astype('uint8')
		
		contours, _ = cv2.findContours(reflectivity_contour,cv2.RETR_EXTERNAL,cv2.CHAIN_APPROX_SIMPLE)

		overlay = np.zeros_like(hrrr_lats)
		for contour in contours:
			if cv2.contourArea(contour)>=100 and cv2.contourArea(contour)<250000:

				mask = np.zeros_like(reflectivity_copy)
				cv2.drawContours(mask,[contour],0,255,-1)
				pixelpoints = cv2.findNonZero(mask)
				pixelpoints = pixelpoints.transpose()
				mask[temp==0] = 0
				if np.max(mask)==0:
					continue
				else:
					snowpixelpoints = cv2.findNonZero(mask)
					snowpixelpoints = snowpixelpoints.transpose()
					threshold = np.mean(reflectivity_copy[snowpixelpoints[1],snowpixelpoints[0]]) + 1.25*np.std(reflectivity_copy[snowpixelpoints[1],snowpixelpoints[0]])
					if threshold < 0:
						threshold = 0
					reflect_contour_2 = np.zeros_like(reflectivity_copy)
					reflect_contour_2[pixelpoints[1],pixelpoints[0]] = reflectivity_copy[pixelpoints[1],pixelpoints[0]]
					reflect_contour_2[reflect_contour_2>=threshold] = 255
					reflect_contour_2[reflect_contour_2<threshold] = 0
					reflect_contour_2 = np.array(reflect_contour_2).astype('uint8')

					interior_contours,_ = cv2.findContours(reflect_contour_2,cv2.RETR_EXTERNAL,cv2.CHAIN_APPROX_SIMPLE)
					for interior_contour in interior_contours:
						if cv2.contourArea(interior_contour)>=100:
							M = cv2.moments(interior_contour)
							cx = int(M['m10']/M['m00'])
							cy = int(M['m01']/M['m00'])
							loncenter = hrrr_lons[cy,cx]
							latcenter = hrrr_lats[cy,cx]
	
							interior_mask = np.zeros_like(reflectivity_copy).astype('uint8')
							cv2.drawContours(interior_mask,[interior_contour],0,255,-1)
							interiorpixelpoints = cv2.findNonZero(interior_mask)
							interiorpixelpoints = interiorpixelpoints.transpose()

							rect = cv2.minAreaRect(interior_contour)
							angle = rect[2]
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

							
                    					width = np.min([great_circle((hrrr_lats[avg1y,avg1x],hrrr_lons[avg1y,avg1x]), (hrrr_lats[avg3y,avg3x],hrrr_lons[avg3y,avg3x])).kilometers,great_circle((hrrr_lats[avg2y,avg2x],hrrr_lons[avg2y,avg2x]), (hrrr_lats[avg4y,avg4x],hrrr_lons[avg4y,avg4x])).kilometers])
                    					length = np.max([great_circle((hrrr_lats[avg1y,avg1x],hrrr_lons[avg1y,avg1x]), (hrrr_lats[avg3y,avg3x],hrrr_lons[avg3y,avg3x])).kilometers,great_circle((hrrr_lats[avg2y,avg2x],hrrr_lons[avg2y,avg2x]), (hrrr_lats[avg4y,avg4x],hrrr_lons[avg4y,avg4x])).kilometers])

							if great_circle((hrrr_lats[avg1y,avg1x],hrrr_lons[avg1y,avg1x]), (hrrr_lats[avg3y,avg3x],hrrr_lons[avg3y,avg3x])).kilometers > great_circle((hrrr_lats[avg2y,avg2x],hrrr_lons[avg2y,avg2x]), (hrrr_lats[avg4y,avg4x],hrrr_lons[avg4y,avg4x])).kilometers:
                       						endpt1lat = hrrr_lats[avg1y,avg1x]
                        					endpt1lon = hrrr_lons[avg1y,avg1x]
                        					endpt2lat = hrrr_lats[avg3y,avg3x]
                        					endpt2lon = hrrr_lons[avg3y,avg3x]
                    					else:
                        					endpt1lat = hrrr_lats[avg2y,avg2x]
                        					endpt1lon = hrrr_lons[avg2y,avg2x]
                        					endpt2lat = hrrr_lats[avg4y,avg4x]
                        					endpt2lon = hrrr_lons[avg4y,avg4x]

							coldpts = np.sum(temp[interiorpixelpoints[1],interiorpixelpoints[0]])
                    					allpts = len(interiorpixelpoints[0][0])
                    
                    					if (length>=250 and length/width>=3.0 and coldpts/allpts>=0.5 and lonpts/allpts>=0.5 and radarpts/allpts>=0.5):
                        					overlay[interiorpixelpoints[1],interiorpixelpoints[0]] = 1

		m.contour(X,Y,overlay,colors=[(float(cm.tab10(i)[0]),float(cm.tab10(i)[1]),float(cm.tab10(i)[2]),float(cm.tab10(i)[3]))],zorder=3)

	cbar = plt.colorbar(im,fraction=0.023,pad=-0.02,ticks=[-20,0,20,40,60,80])
	plt.box(False)


	outfil = '/gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/uploads/outimages/hrrre/' + outname + '_snowbands.png'
	plt.savefig(outfil,bbox_inches='tight',dpi=500)
	plt.close()

datapaths = glob.glob('/gpfs_backup/stormtrack/jtradfor/ensemble_data/rawdata/hrrre/*breflectdetect*')
tempspath = glob.glob('/gpfs_backup/stormtrack/jtradfor/ensemble_data/rawdata/hrrre/*temps*')

if len(datapaths)>0:

	latest = 0
	latestpath = datapaths[0]
	latesttempspath = tempspath[0]
	for datapath in datapaths:
	    if int(os.path.basename(datapath)[9:11]) > latest:
		latest = int(os.path.basename(datapath)[9:11]) 
		latestpath = datapath
		latesttempspath = '/gpfs_backup/stormtrack/jtradfor/ensemble_data/rawdata/hrrre/' + os.path.basename(datapath)[:13] + '_temps.npy'
	fil = os.path.basename(latestpath)[0:13]
	reflectivities = np.load(latestpath)
	temps = np.load(latesttempspath)
	pbsdfil = detection(reflectivities,temps,fil)
