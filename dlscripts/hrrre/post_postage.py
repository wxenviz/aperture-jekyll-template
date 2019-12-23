import glob
import subprocess as sp
import os
import sys


os.chdir('/gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/')
sp.call('git add -A',shell=True)
sp.call("git commit -m 'new images'",shell=True)
sp.call('git push origin master',shell=True)

sp.call('rm /gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/dlscripts/hrrre/unlock/*',shell=True)
