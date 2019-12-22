#!/bin/csh
#BSUB -n 4
#BSUB -q mea716
#BSUB -R span[ptile=4]
#BSUB -W 12:15

source /gpfs_backup/stormtrack/jtradfor/aliases.csh  

set datestr="193552000"
python /gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/dlscripts/href/nssl_interp.py $datestr
