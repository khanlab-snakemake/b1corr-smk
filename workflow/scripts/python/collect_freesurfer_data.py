import itertools
import numpy as np
import pandas as pd
import nibabel as nib
from os.path import dirname

# Subject info
subject = snakemake.wildcards.subject
if 'MSTRCHT' in subject:
    site='Maastricht'
elif 'NEW' in subject:
    site='Maastricht-2'
else:
    site='London'

# Some parameters for extracting the data
fs_dir = snakemake.params.subjects_dir
modalities = snakemake.params.modalities
timepoints = snakemake.params.timepoints
hemispheres = snakemake.params.hemispheres
rois = snakemake.params.structures

## Do surface work first
in_dir=dirname(snakemake.input.surface_maps[0])

columns = []
for combination in itertools.product(modalities,timepoints,hemispheres):

    # Define column names for output dataframe
    columns.append('{}-{}-{}'.format(combination[0],combination[1],combination[2]))
    
    # Load .shape.gii surface overlay
    in_file = '{}/{}.{}_{}_{}.fsaverage.shape.gii'.format(in_dir,combination[2],subject,combination[0],combination[1])
    gii = nib.load(in_file)
    if 'gii_data' in dir():
        gii_data = np.c_[gii_data, gii.darrays[0].data]
    else:
        gii_data = gii.darrays[0].data        

# Save to pandas pickle file
df_surface = pd.DataFrame(data=gii_data, columns=columns)
df_surface.to_pickle(snakemake.output.surface_data)

## Extract aseg volume data
data = []
for parc, timepoint in zip(snakemake.input.aseg,timepoints):
    in_aseg = nib.load(parc)
    aseg = in_aseg.get_fdata()

    in_t1 = nib.load('{0}/{1}_{2}.long.{1}/tmp/t1_{2}.mgz'.format(fs_dir,subject,timepoint))
    t1 = in_t1.get_fdata()  

    # Get volume (nr of voxels) and T1 per ROI
    for roi, labels in rois.items():
        for h, hh in enumerate(hemispheres):
            volume = (aseg == labels[h]).sum()
            t1_mean = np.nanmean(t1[aseg == labels[h]])
            data.append([subject,site,timepoint,roi,hh,volume,t1_mean])

    # Do the same for whole brain (non-zero aseg voxels)
    volume = (aseg != 0).sum()
    t1_mean = np.nanmean(t1[aseg > 0])
    data.append([subject,site,timepoint,'brain','both',volume,t1_mean])

# Convert to pandas dataframe and save to pandas pickle file
columns = ['subject','site','timepoint','roi','hemisphere','volume','t1']
df_aseg = pd.DataFrame(data=data, columns=columns)
df_aseg.to_pickle(snakemake.output.aseg_data)
