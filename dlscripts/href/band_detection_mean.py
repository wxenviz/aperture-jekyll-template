import subprocess as sp
import pygrib as pg
import numpy as np
import netCDF4 as nc4
import csv
import cv2
import os
import re
import glob
import sys

from geopy.distance import great_circle
from mpl_toolkits.basemap import Basemap, interp
from datetime import datetime, timedelta
from matplotlib.colors import LinearSegmentedColormap
from matplotlib import cm
from scipy import misc
from mpl_toolkits.basemap import Basemap
from matplotlib.patches import Polygon
from matplotlib import cm
from matplotlib.patches import Patch

import matplotlib.pyplot as plt
plt.switch_backend('agg')

dark2 = cm.get_cmap('Dark2')

datesub = str(sys.argv[1])
x = np.load('/gpfs_backup/stormtrack/jtradfor/ensemble_data/reference/basemapx.npy')
y = np.load('/gpfs_backup/stormtrack/jtradfor/ensemble_data/reference/basemapy.npy')

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

m = Basemap(projection='lcc',lat_0=5,lon_0=-100,llcrnrlon=-126,llcrnrlat=23,urcrnrlon=-63,urcrnrlat=50,resolution='l')
lons,lats = m(x,y,inverse=True)    

models = ["namnest","namnest_back","arw","arw_back","nmmb","nmmb_back","hrrr","hrrr_back"]
model_labels = ["Nam Nest","Nam Nest -12h", "ARW", "ARW -12h", "NMMB", "NMMB -12h", "HRRR", "HRRR -6h"]
masks = ["namnest","namnest","arw","arw","nmmb","nmmb","hrrr","hrrr"]

hrefstats = []
try:
	with open('/gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/dlscripts/href/hrefstats.csv','r') as statsfile:
		reader = csv.reader(statsfile,delimiter=',')
		for line in reader:
			hrefstats.append(line)
	statsfile.close()
except:
	pass

for fhour in range(0,31):
	print(fhour)
	start = datetime.now()
	reflectivities = []
	temps = []
	for k,model in enumerate(models):
		bref_fil = '/gpfs_backup/stormtrack/jtradfor/ensemble_data/rawdata/href/bref_%s.npy' % (model)
		temp_fil = '/gpfs_backup/stormtrack/jtradfor/ensemble_data/rawdata/href/temp_%s.npy' % (model)
		bref_mask = '/gpfs_backup/stormtrack/jtradfor/ensemble_data/reference/%s_mask.npy' % (masks[k])
		bref = np.load(bref_fil)[fhour]
		bref[bref<=0] = 0.0
		temp_mask = np.load(bref_mask)
		bref[temp_mask] = np.nan
		temp = np.load(temp_fil)[fhour]
		reflectivities.append(bref)
		temp[temp_mask] = np.nan
		temps.append(temp)

	reflectivities_copy = np.copy(reflectivities)
	reflect_mean = np.mean(reflectivities_copy,axis=0)

	plt.figure(figsize=(16,9))
	m = Basemap(projection='lcc',lat_0=5,lon_0=-100,llcrnrlon=-126,llcrnrlat=23,urcrnrlon=-63,urcrnrlat=50,resolution='l')	

	shp_info = m.readshapefile('/gpfs_backup/stormtrack/jtradfor/ensemble_data/reference/st99_d00','states',drawbounds=False)

	ax = plt.gca()

	for nshape,seg in enumerate(m.states):
		poly = Polygon(seg,facecolor='white',edgecolor='white',zorder=1,linewidth=.5)
		poly2 = Polygon(seg,facecolor='none',edgecolor='black',zorder=3,linewidth=.5)
		ax.add_patch(poly)
		ax.add_patch(poly2)

	reflect_mean[reflect_mean>1000000] = np.nan
	reflect_mean[reflect_mean>60.0] = 60.0
	reflect_mean[reflect_mean<=0] = np.nan


	im = m.contourf(x,y,reflect_mean,zorder=2,aspect='equal',interpolation='none',cmap=cmap1,levels=levels)

	overlay = []
	for i,reflectivity_copy in enumerate(reflectivities_copy):
		temp = temps[i]
		temp[temp<=273.15] = 1
		temp[temp>273.15] = 0


		reflectivity_contour = np.copy(reflectivity_copy)
		kernel = np.ones((5,5))/25.
		reflectivity_contour = cv2.filter2D(reflectivity_contour,-1,kernel)
		reflectivity_contour[reflectivity_contour>5] = 255
		reflectivity_contour[reflectivity_contour<5] = 0
		reflectivity_contour[np.isnan(reflectivity_contour)] = 0
		reflectivity_contour = np.array(reflectivity_contour).astype('uint8')
		
		contours, _ = cv2.findContours(reflectivity_contour,cv2.RETR_EXTERNAL,cv2.CHAIN_APPROX_SIMPLE)

		overlay = np.zeros_like(reflectivity_contour)
		for contour in contours:
			if cv2.contourArea(contour)>=100:

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
					if threshold < 15:
						threshold = 15
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
							loncenter = lons[cy,cx]
							latcenter = lats[cy,cx]
	
							interior_mask = np.zeros_like(reflectivity_copy).astype('uint8')
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

							try:	
								width = np.min([great_circle((lats[avg1y,avg1x],lons[avg1y,avg1x]), (lats[avg3y,avg3x],lons[avg3y,avg3x])).kilometers,great_circle((lats[avg2y,avg2x],lons[avg2y,avg2x]), (lats[avg4y,avg4x],lons[avg4y,avg4x])).kilometers])
								length = np.max([great_circle((lats[avg1y,avg1x],lons[avg1y,avg1x]), (lats[avg3y,avg3x],lons[avg3y,avg3x])).kilometers,great_circle((lats[avg2y,avg2x],lons[avg2y,avg2x]), (lats[avg4y,avg4x],lons[avg4y,avg4x])).kilometers])
								if great_circle((lats[avg1y,avg1x],lons[avg1y,avg1x]), (lats[avg3y,avg3x],lons[avg3y,avg3x])).kilometers > great_circle((lats[avg2y,avg2x],lons[avg2y,avg2x]), (lats[avg4y,avg4x],lons[avg4y,avg4x])).kilometers:
									endpt1lat = lats[avg1y,avg1x]
									endpt1lon = lons[avg1y,avg1x]
									endpt2lat = lats[avg3y,avg3x]
									endpt2lon = lons[avg3y,avg3x]
								else:
									endpt1lat = lats[avg2y,avg2x]
									endpt1lon = lons[avg2y,avg2x]
									endpt2lat = lats[avg4y,avg4x]
									endpt2lon = lons[avg4y,avg4x]

								reflect_avg = np.mean(reflectivity_copy[interiorpixelpoints[1],interiorpixelpoints[0]])
								coldpts = np.sum(temp[interiorpixelpoints[1],interiorpixelpoints[0]])
								allpts = len(interiorpixelpoints[0][0])
								
								if (length>=250 and length/width>=3.0 and coldpts/allpts>=0.5):
									solidity = (cv2.contourArea(interior_contour)/rectarea)
									print(solidity)
									realrectarea = length*width
									realarea = realrectarea*solidity
									print(rectarea,realarea,realrectarea)
									overlay[interiorpixelpoints[1],interiorpixelpoints[0]] = 1
									hrefstats.append([datesub + (str(fhour).zfill(2)) + '00',models[i],str(latcenter),str(loncenter),str(endpt1lat),str(endpt1lon),str(endpt2lat),str(endpt2lon),str(length),str(width),str(realarea),str(threshold),str(reflect_avg)])
							except:
								continue


		m.contour(x,y,overlay,colors=[dark2(i)],linewidths=0.75,zorder=3,label=models[i])

	cbar = plt.colorbar(im,fraction=0.025,ticks=[0,10,20,30,40,50,60])
	cbar.ax.yaxis.set_tick_params(color='w')
	cbar.ax.set_yticklabels([0,10,20,30,40,50,60],color='w')
	legend_elements = []
	for i in range(0,8):
		legend_element = Patch(facecolor=dark2(i),label=model_labels[i])
		legend_elements.append(legend_element)
	plt.legend(handles=legend_elements)
	plt.box(False)

	outfil = '/gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/uploads/outimages/href/%s%s00_snowbands.png' % (datesub,str(fhour).zfill(2))
	plt.savefig(outfil,facecolor='#101010',bbox_inches='tight',dpi=500)
	plt.close()

	print datetime.now() - start

with open('/gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/dlscripts/href/hrefstats.csv','wb') as statsfile:
	csv.writer(statsfile).writerows(hrefstats)

statsfile.close()

