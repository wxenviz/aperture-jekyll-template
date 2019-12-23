#!/bin/csh
source /gpfs_backup/stormtrack/jtradfor/aliases.csh

if ( -e /gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/dlscripts/hrrre/unlock/ftpcheck.lock) then
	exit 1
endif

touch /gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/dlscripts/hrrre/unlock/ftpcheck.lock

python /gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/dlscripts/hrrre/ftpcheck.py

python /gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/dlscripts/hrrre/dl.py
wait
python /gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/dlscripts/hrrre/interpolate.py
wait
bsub < /gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/dlscripts/hrrre/BR_mean.csh
bsub < /gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/dlscripts/hrrre/R_mean.csh
bsub < /gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/dlscripts/hrrre/apcp_ab.csh
bsub < /gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/dlscripts/hrrre/apcp_mean.csh
bsub < /gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/dlscripts/hrrre/aweasd_ab.csh
bsub < /gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/dlscripts/hrrre/BR_ab.csh
bsub < /gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/dlscripts/hrrre/hpcp_ab.csh
bsub < /gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/dlscripts/hrrre/hpcp_mean.csh
bsub < /gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/dlscripts/hrrre/hweasd_ab.csh
bsub < /gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/dlscripts/hrrre/hweasd_mean.csh
bsub < /gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/dlscripts/hrrre/R_ab.csh
bsub < /gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/dlscripts/hrrre/R_pbab.csh
bsub < /gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/dlscripts/hrrre/R_pb.csh
bsub < /gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/dlscripts/hrrre/aweasd_mean.csh
bsub < /gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/dlscripts/hrrre/apcp_postage.csh
bsub < /gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/dlscripts/hrrre/BR_postage.csh
bsub < /gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/dlscripts/hrrre/hpcp_postage.csh
bsub < /gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/dlscripts/hrrre/R_postage.csh
wait
echo images_done

set flip = 0;
while ( $flip == 0)
	if ( -e /gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/dlscripts/hrrre/unlock/Rpostage.unlock && -e /gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/dlscripts/hrrre/unlock/apcp_postage.unlock && -e /gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/dlscripts/hrrre/unlock/BRpostage.unlock && -e /gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/dlscripts/hrrre/unlock/hpcp_postage.unlock) then
		set flip = 1
	else
		echo "sleeping"
		sleep 10
	endif
end
wait
echo "they all exist"

python /gpfs_backup/stormtrack/jtradfor/ensemble_data/wxenviz.github.io/dlscripts/hrrre/post_postage.py
echo fin
