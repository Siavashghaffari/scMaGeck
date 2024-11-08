import scanpy as sc
from scipy import io
import pandas as pd
import numpy as np


file_path = "/gstore/scratch/u/ghaffars/Dataset/sublib4/raw_qc.h5ad"
Name = 'sublib4'
Dataset_HOME = os.path.join('/gstore/scratch/u/ghaffars/scMaGeck', Name)
folderName=os.path.join(Dataset_HOME, f'matrix_files_'+Name)
meta_path = os.path.join(Dataset_HOME, f'metadata_'+Name+'.csv')


#Read the AnnData
bdata = sc.read(file_path)


with open(folderName + '/barcodes.tsv', 'w') as f:
    for item in bdata.obs_names:
        f.write(item + '\n')
        
with open(folderName + '/features.tsv', 'w') as f:
    for item in ['\t'.join([x,x,'Gene Expression']) for x in bdata.var_names]:
        f.write(item + '\n')
        
io.mmwrite(folderName +'/matrix', bdata.X.T)

bdata.obs.to_csv(meta_path)