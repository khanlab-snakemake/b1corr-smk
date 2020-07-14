import numpy as np
import pandas as pd
from os import path
from pathlib import Path

## Combines the output from each individual subject into group- and modality-based
## datasets for analyses and visualization

# Do surface work first
subjects = []
for data in snakemake.input.surface_data:
    df = pd.read_pickle(data)
    
    base = path.basename(data)

    subject = base.split('_')[0]
    subjects.append(subject)

    group = 'Maastricht' if 'MSTRCHT' in subject else 'London'
  
    for column in df.columns:
        column_split = column.split('-')
        var_name = '{}_{}_{}_{}'.format(group,column_split[0],column_split[1],column_split[2])
        
        if var_name in dir():
            exec("{0} = np.c_[{0}, df['{1}'].to_numpy()]".format(var_name,column))
        else:
            exec("{0} = df['{1}'].to_numpy()".format(var_name,column))          

# Save to pandas pickle file
for output in snakemake.output.surface_data:
    base = Path(output).stem
    site = output.split('/')[-2]
    var_name = '{}_{}'.format(site,base)
    
    exec('df = pd.DataFrame(data={})'.format(var_name))
    df.to_pickle(output)


## Then extract aseg volume data. Simple
df = pd.read_pickle(snakemake.input.aseg_data[0])
for data in snakemake.input.aseg_data:
    tmp = pd.read_pickle(data)
    df = pd.concat([df, tmp], ignore_index=True)

# Save to pandas pickle file
df.to_pickle(snakemake.output.aseg_data)


## Finally, combine segmentation metrics. Simple
df = pd.read_pickle(snakemake.input.seg_data[0])
for data in snakemake.input.seg_data:
    tmp = pd.read_pickle(data)
    df = pd.concat([df, tmp], ignore_index=True)

# Save to pandas pickle file
df.to_pickle(snakemake.output.seg_data)
