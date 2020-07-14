import itertools
import numpy as np
import pandas as pd
import nibabel as nib
from os.path import dirname

def bool_vec_dissimilarity(booldata1, booldata2, method):
    from scipy.spatial.distance import dice, jaccard
    methods = {'dice': dice, 'jaccard': jaccard}
    if not (np.any(booldata1) or np.any(booldata2)):
        return 0
    return 1 - methods[method](booldata1.flat, booldata2.flat)

def dice_jaccard(data1,data2):
    scale = 1.0

    jaccard = bool_vec_dissimilarity(data1 != 0, data2 != 0, method='jaccard')
    dice = 2.0 * jaccard / (jaccard + 1.0)
    volume1 = (scale * len(data1[data1 != 0]))
    volume2 = (scale * len(data2[data2 != 0]))        

    return [volume1, volume2, round(dice,5), round(jaccard,5)]

def find_border(data):
    from scipy.ndimage.morphology import binary_erosion
    eroded = binary_erosion(data)
    border = np.logical_and(data, np.logical_not(eroded))
    return border

def get_coordinates(data, affine):
    if len(data.shape) == 4:
        data = data[:, :, :, 0]
    indices = np.vstack(np.nonzero(data))
    indices = np.vstack((indices, np.ones(indices.shape[1])))
    coordinates = np.dot(affine, indices)
    return coordinates[:3, :]

def hausdorff(origdata1, origdata2, affine):
    from scipy.spatial.distance import cdist

    data1 = np.zeros(origdata1.shape)
    data1[origdata1 != 0] = 1
    
    data2 = np.zeros(origdata2.shape)
    data2[origdata2 != 0] = 1
    
    data1 = np.logical_not(
        np.logical_or(data1 == 0, np.isnan(data1)))
    data2 = np.logical_not(
        np.logical_or(data2 == 0, np.isnan(data2)))

    border1 = find_border(data1)
    border2 = find_border(data2)

    set1_coordinates = get_coordinates(border1, affine)
    set2_coordinates = get_coordinates(border2, affine)
    distances = cdist(set1_coordinates.T, set2_coordinates.T)
    mins = np.concatenate((np.amin(distances, axis=0),
                           np.amin(distances, axis=1)))

    return round(np.max(mins),5)

# Subject info
subject = snakemake.wildcards.subject
if 'MSTRCHT' in subject:
    site='Maastricht'
elif 'NEW' in subject:
    site='Maastricht-2'
else:
    site='London'

# Some parameters for extracting the data
timepoints = snakemake.params.timepoints
hemispheres = snakemake.params.hemispheres
rois = snakemake.params.structures

## Extract aseg volume data
data = []

# Get segmentation metrics per ROI
for roi, labels in rois.items():
    for h, hh in enumerate(hemispheres):
        in_orig = nib.load(snakemake.input[0])
        orig = in_orig.get_fdata()
        orig[orig != labels[h]] = 0

        in_corr = nib.load(snakemake.input[1])
        corr = in_corr.get_fdata()
        corr[corr != labels[h]] = 0

        overlap  = dice_jaccard(orig,corr)
        distance = hausdorff(orig,corr,in_orig.affine)

        data.append([subject,site,roi,hh,overlap[1]-overlap[0],overlap[2],overlap[3],distance])

# Convert to pandas dataframe and save to pandas pickle file
columns = ['subject','site','roi','hemisphere','difference','dice','jaccard','hausdorff']
df_aseg = pd.DataFrame(data=data, columns=columns)
df_aseg.to_pickle(snakemake.output[0])
