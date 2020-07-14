#!usr/bin/bash

subj=$1

export out=`dirname $2`
export inv2=$3
export t1w=$4
export t1=$5

echo '-------------------BET brain extraction----------------------'
echo ""

fsl5.0-bet $inv2 $out/brain -v -R -m -f 0.20 -g 0

pushd $out

echo ""
echo '-----------Fine tune brain mask using CBS tools--------------'
echo ""

python3 - <<'EOF'
import os
import nighres
nighres.brain.mp2rage_dura_estimation(os.environ['inv2'],
                                      'brain_mask.nii.gz',
                                      background_distance=2.0,
                                      output_type='dura_prior',
                                      save_data=True,
                                      overwrite=True,
                                      output_dir='.',
                                      file_name='inv2')
EOF
popd

fsl5.0-fslmaths $out/inv2_dura-proba.nii.gz -thr 0.99 -uthr 1 -bin $out/inv2_dura-mask.nii.gz
fsl5.0-fslmaths $out/brain_mask.nii.gz -sub $out/inv2_dura-mask.nii.gz $out/brain_mask-filt.nii.gz

pushd $out
python3 - <<'EOF'
import os
import nighres
nighres.brain.mp2rage_dura_estimation(os.environ['inv2'],
                                      'brain_mask-filt.nii.gz',
                                      background_distance=2.0,
                                      output_type='bg_prior',
                                      save_data=True,
                                      overwrite=True,
                                      output_dir='.',
                                      file_name='inv2_bg')
EOF
popd

mean_inv2=`fsl5.0-fslstats $inv2 -m`
fsl5.0-fslmaths $inv2 -uthr $mean_inv2 -binv $out/inv2_mask.nii.gz

fsl5.0-fslmaths $out/inv2_dura-proba.nii.gz -thr 0.8 $out/inv2_dura-mask.nii.gz
fsl5.0-fslmaths $t1w -thr 350 -bin $out/t1w_csf-mask.nii.gz
fsl5.0-fslmaths $out/inv2_mask.nii.gz -sub $out/inv2_bg_dura-proba.nii.gz \
						-sub $out/inv2_dura-mask.nii.gz \
						-mul $out/t1w_csf-mask.nii.gz $out/inv2_mask_final.nii.gz

fsl5.0-fslmaths $out/inv2_mask_final.nii.gz -mul $out/brain_mask.nii.gz -thr 0 -mul $t1w $out/mp2rage_brain.nii.gz
fsl5.0-fslmaths $out/inv2_mask_final.nii.gz -mul $out/brain_mask.nii.gz -thr 0 $2

echo ""
echo '----------------------Saved and done-----------------------'
echo ""
