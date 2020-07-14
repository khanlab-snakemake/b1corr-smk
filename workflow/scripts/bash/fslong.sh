#!usr/bin/bash
export FREESURFER_HOME=/project/6050199/rhaast/freesurfer
source $FREESURFER_HOME/SetUpFreeSurfer.sh

subj=$1

t1w_orig=$2
t1w_orig_base=`basename $t1w_orig .nii.gz`
t1w_corr=$3
t1w_corr_base=`basename $t1w_corr .nii.gz`

mask=$4
SUBJECTS_DIR=$5
EXPERT=$7

echo '----------------------- INPUT --------------------------'
echo 'Subject ID:          '$subj
echo 'MP2RAGE T1w orig:    '$t1w_orig
echo 'MP2RAGE T1w corr:    '$t1w_corr
echo 'Brainmask:           '$mask

echo '--------------------- FS (normal) ----------------------'

batch_file=$SUBJECTS_DIR/${subj}_fs_regular.txt
rm $batch_file

for tp in orig corr ; do
 subjid=${subj}_${tp}
 tpid_dir=$SUBJECTS_DIR/$subjid/mri/orig

 if [ ! -e $tpid_dir ] ; then
   mkdir -p $tpid_dir
 fi

 input=t1w_${tp}
 echo fscalc $mask mul ${!input} --o $tpid_dir/001.mgz
 echo recon-all -subjid $subjid -all -no-wsgcaatlas -parallel -openmp 4 -hires -expert ${EXPERT} >> $batch_file
done

#parallel -u -j 2 {} < $batch_file

echo '---------------------- FS (base) -----------------------'
recon-all -base $subj -tp ${subj}_orig -tp ${subj}_corr -all -parallel -openmp 8 -hires -expert ${EXPERT} 

echo '---------------------- FS (long) -----------------------'

batch_file=$SUBJECTS_DIR/${subj}_fs_long.txt
rm $batch_file

for tp in orig corr ; do
  subjid=${subj}_${tp}
  echo recon-all -long $subjid ${subj} -all -parallel -openmp 4 -hires -expert ${EXPERT} >> $batch_file
done

parallel -u -j 2 {} < $batch_file

if [[ ( -f $SUBJECTS_DIR/${subj}_corr.long.${subj}/scripts/recon-all.done && -f $SUBJECTS_DIR/${subj}_orig.long.${subj}/scripts/recon-all.done ) ]] ; then
    touch $6
fi

echo '------------------------ DONE --------------------------'
