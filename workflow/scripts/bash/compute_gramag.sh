base=/home/ROBARTS/rhaast/graham/scratch/B1corr
surfmorph_dir=$base/output/surfmorph

timepoints="orig corr"
tps=($timepoints)
tail -n +2 $base/participants.tsv | while read s ; do
    qmri_dir=$surfmorph_dir/qmri/${s}/anat
	for tp in $timepoints ; do
	    labels_dir=$surfmorph_dir/labels/${s}_${tp}/anat
		input_dir=$base/output/freesurfer_spm/${s}_${tp}.long.${s}
		if [ ! -f $input_dir/tmp/T1.nii.gz ] ; then
			mri_convert $input_dir/mri/T1.mgz $input_dir/tmp/T1.nii.gz -rl $labels_dir/${s}_${tp}_brain_T1w.nii.gz -nc
		fi
		if [ ! -f $input_dir/tmp/T1-gramag.nii.gz ] ; then
			wb_command -volume-gradient $input_dir/tmp/T1.nii.gz $input_dir/tmp/T1-gramag.nii.gz
		fi
	done
	if [ ! -f $qmri_dir/${s}_T1-gramag-diff.nii.gz ] ; then
	    mkdir -p $surfmorph_dir/qmri/${s}/anat/
	    corr_dir=$base/output/freesurfer_spm/${s}_${tps[1]}.long.${s}/tmp
	    orr_dir=$base/output/freesurfer_spm/${s}_${tps[0]}.long.${s}/tmp
		fslmaths $corr_dir/T1-gramag.nii.gz -sub $orr_dir/T1-gramag.nii.gz $qmri_dir/${s}_T1-gramag-diff.nii.gz
	fi
done
