import numpy as np
import nibabel as nib

seg = nib.load(snakemake.input.seg)
seg_data = seg.get_fdata()

for r, (roi, labels) in enumerate(snakemake.params.structures.items()):
    binarized = np.zeros(seg_data.shape)
    binarized[np.isin(seg_data,labels)] = 1

    img = nib.Nifti1Image(binarized,affine=seg.affine,header=seg.header)
    nib.save(img,snakemake.output.seg[r])

print(snakemake.input.t1w)
t1w = nib.load(snakemake.input.t1w)
img = nib.Nifti1Image(t1w.get_fdata(), affine=t1w.affine, header=t1w.header) 
nib.save(img, snakemake.output.t1w)