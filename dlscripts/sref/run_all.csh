#!/bin/csh
source /gpfs_backup/stormtrack/jtradfor/aliases.csh

/usr/bin/rm /gpfs_backup/stormtrack/jtradfor/ensemble_data/rawdata/sref/*
wait
python /gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/dlscripts/sref/dl.py
