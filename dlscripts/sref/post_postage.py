import glob
import subprocess as sp
import os
import sys

#sys.stdout = open('/gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/logfile','a+')

datapaths = glob.glob('/gpfs_backup/stormtrack/jtradfor/ensemble_data/rawdata/hrrre/*creflect*')
datapaths3 = glob.glob('/gpfs_backup/stormtrack/jtradfor/ensemble_data/rawdata/hrrre/*breflect*')
datapaths2 = glob.glob('/gpfs_backup/stormtrack/jtradfor/ensemble_data/rawdata/hrrre/*apcp*')

if len(datapaths)>0:
	
	earliest = 37
	earliestpath = datapaths[0]
	for datapath in datapaths:
	    if int(os.path.basename(datapath)[9:11]) < earliest:
		earliest = int(os.path.basename(datapath)[9:11])
		earliestpath = datapath
	sp.call('rm %s' % (earliestpath),shell=True)

if len(datapaths3)>0:
	
	earliest = 37
	earliestpath = datapaths3[0]
	for datapath in datapaths3:
	    if int(os.path.basename(datapath)[9:11]) < earliest:
		earliest = int(os.path.basename(datapath)[9:11])
		earliestpath = datapath
	sp.call('rm %s' % (earliestpath),shell=True)


if len(datapaths2)>0:
	
	earliest = 37
	earliestpath = datapaths2[0]
	print earliestpath
	for datapath in datapaths2:
	    if int(os.path.basename(datapath)[9:11]) < earliest:
		earliest = int(os.path.basename(datapath)[9:11])
		earliestpath = datapath
	sp.call('rm %s' % (earliestpath),shell=True)


os.chdir('/gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/')
sp.call('git add -A',shell=True)
sp.call("git commit -m 'new images'",shell=True)
sp.call('git push origin master',shell=True)

sp.call('rm /gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/dlscripts/unlock/*',shell=True)
