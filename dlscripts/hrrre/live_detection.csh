#!/bin/csh

source /usr/local/apps/mpich3/centos7/intelmpi2016.csh
source /gpfs_backup/stormtrack/jtradfor/aliases.csh

unsetenv MPICH_NO_LOCAL

python /gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/dlscripts/hrrre/live_detection.py
