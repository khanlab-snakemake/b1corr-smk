#!usr/bin/bash
export FREESURFER_HOME=/project/6050199/rhaast/freesurfer
source $FREESURFER_HOME/SetUpFreeSurfer.sh

subj=$1
uni_orig=$2
uni_corr=$3
t1_orig=$4
t1_corr=$5
b1=$6

out_dir=$7

export SUBJECTS_DIR=$8

modalities="uni t1 thickness"
timepoints="orig corr"
hemispheres="lh rh"

# Loop over different images, timepoints and hemispheres
for m in $modalities ; do
    for t in $timepoints ; do
        eval declare -n in_file="${m}_${t}"
        if [[ ( $m == "b1" ) && ( $t == "corr" ) ]] ; then
            continue 
        fi
        
        # Convert input images to .mgz and conform to FreeSurfer volume output
        mri_convert $in_file $SUBJECTS_DIR/${subj}_${t}.long.${subj}/tmp/${m}_${t}.mgz -rl $SUBJECTS_DIR/${subj}_${t}.long.${subj}/mri/orig.mgz -nc
        
        for h in $hemispheres ; do
            out_native_mgh="${out_dir}/${h}.${subj}_${m}_${t}.mgh"
            out_fsaverage_mgh="${out_dir}/${h}.${subj}_${m}_${t}.fsaverage.mgh"

            if [ $m == "thickness" ] ; then
                in_file=${SUBJECTS_DIR}/${subj}_${t}.long.${subj}/surf/${h}.thickness
                
                # Map directly to fsaverage in case of cortical thickness
                mri_surf2surf --srcsubject ${subj}_corr.long.${subj} --sval $in_file --trgsubject fsaverage --tval ${out_fsaverage_mgh} --hemi ${h} 
            else
            
                # Else map first onto subject surface mesh and then fsaverage
                mri_vol2surf --mov $in_file --regheader ${subj}_corr.long.${subj} --hemi $h --projfrac-avg 0.2 0.8 0.05 --o $out_native_mgh --cortex
                mri_surf2surf --srcsubject ${subj}_corr.long.${subj} --sval $out_native_mgh --trgsubject fsaverage --tval $out_fsaverage_mgh --hemi ${h}
            fi

            # Convert to gifti for connectome workbench compatibility
            out_fsaverage_gii="${out_dir}/${h}.${subj}_${m}_${t}.fsaverage.shape.gii"
            mri_convert $out_fsaverage_mgh $out_fsaverage_gii
        done 
    done
done

for h in $hemispheres ; do
    out_native_mgh="${out_dir}/${h}.${subj}_b1.mgh"
    out_fsaverage_mgh="${out_dir}/${h}.${subj}_b1.fsaverage.mgh"
    out_fsaverage_gii="${out_dir}/${h}.${subj}_b1.fsaverage.shape.gii"

    mri_vol2surf --mov $b1 --regheader ${subj}_corr.long.${subj} --hemi $h --projfrac-avg 0.2 0.8 0.05 --o $out_native_mgh --cortex
    mri_surf2surf --srcsubject ${subj}_corr.long.${subj} --sval $out_native_mgh --trgsubject fsaverage --tval $out_fsaverage_mgh --hemi ${h}
    mri_convert $out_fsaverage_mgh $out_fsaverage_gii
done
