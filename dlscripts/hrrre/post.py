import glob
import subprocess as sp
import os
import sys

#sys.stdout = open('/gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/logfile','a+')

os.chdir('/gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/')
sp.call('git add -A',shell=True)
sp.call("git commit -m 'new images'",shell=True)
sp.call('git push origin master',shell=True)

