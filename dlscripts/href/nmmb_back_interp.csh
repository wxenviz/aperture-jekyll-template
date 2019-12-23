#!/bin/csh
#BSUB -o /gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/dlscripts/href/nmmb_back_interp_out
#BSUB -e /gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/dlscripts/href/nmmb_back_interp_err
#BSUB -n 4
#BSUB -q mea716
#BSUB -R span[ptile=4]
#BSUB -W 12:15

source /gpfs_backup/stormtrack/jtradfor/aliases.csh  

<<<<<<< HEAD
set datestr="193570000"
=======
set datestr="193561200"
>>>>>>> ed7c10900abed222533e501c7f20e345aa9f0518
python /gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/dlscripts/href/nmmb_back_interp.py $datestr
