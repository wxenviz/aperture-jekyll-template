#!/bin/csh
#BSUB -o /gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/dlscripts/href/aweasd_postage_out
#BSUB -e /gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/dlscripts/href/aweasd_postage_err
#BSUB -n 4
#BSUB -q mea716
#BSUB -R span[ptile=4]
#BSUB -W 12:15

source /usr/local/apps/mpich3/centos7/intelmpi2016.csh

unsetenv MPICH_NO_LOCAL

set datestr="193561200"
python /gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/dlscripts/href/aweasd_postage.py $datestr
