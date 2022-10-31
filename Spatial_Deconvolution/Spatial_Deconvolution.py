##import packages
import numpy as np
import pandas as pd
import sys
import os
import time as tm
import pickle
from functools import partial
import scipy.stats as st
from scipy.stats import wasserstein_distance
import scipy.stats
import copy
from sklearn.model_selection import KFold
import pandas as pd
import multiprocessing
import matplotlib as mpl 
import matplotlib.pyplot as plt
import scanpy as sc
import warnings
import subprocess
import seaborn as sns
from sklearn.metrics import mean_squared_error
from scipy.spatial.distance import jensenshannon
from scipy.stats import pearsonr,ttest_ind,mannwhitneyu
import matplotlib
import time


### please in the SpatialBenchmarking dir.
os.chdir('~/Spatial_Deconvolution/')

### please add the your SpatialBenchmarking dir into the python path 
sys.path.append('./')

import Benchmarking.DeconvolutionSpot as DeconvolutionSpot




### for Cell2location, Stereoscope, Tangram, DestVI, you must have .h5ad files as input.
RNA_h5ad = '~/Spatial_Deconvolution/WAT_scRNA_MTX.h5ad'
Spatial_h5ad = '~/Spatial_Deconvolution/st_s' + sys.argv[1] + '_MTX.h5ad'  # data can be downloaded from https://data.mendeley.com/datasets/3bs5f8mvbs/1
celltype_key = 'celltype_new'  # key of the celltype in RNA_h5Seurat
output_path = '~/Spatial_Deconvolution/s'+ sys.argv[1]
if not os.path.exists(output_path):
    os.mkdir(output_path)

test = DeconvolutionSpot.Deconvolutions(RNA_h5ad = RNA_h5ad, Spatial_h5ad = Spatial_h5ad, celltype_key = celltype_key, output_path = output_path)
Methods = ['Cell2location', 'Stereoscope','Tangram', 'DestVI']
Result = test.Dencon(Methods)




### for RCTD and SPOTlight, you must have .h5seurat files as input.
RNA_h5Seurat = '~/Spatial_Deconvolution/WAT_scRNA_MTX.h5seurat'
Spatial_h5Seurat = '~/Spatial_Deconvolution/st_s' + sys.argv[1] + '_MTX.h5seurat'  # data can be downloaded from https://data.mendeley.com/datasets/3bs5f8mvbs/1 
celltype_key = 'celltype_new'  # key of the celltype in RNA_h5Seurat
output_path = '~/zjw/20220528WAT/cell_metabolism_spatial/S' + sys.argv[1]
python_path = '~/anaconda3/bin/python'
if not os.path.exists(output_path):
    os.mkdir(output_path)

test = DeconvolutionSpot.Deconvolutions(RNA_h5Seurat= RNA_h5Seurat, Spatial_h5Seurat = Spatial_h5Seurat, celltype_key = celltype_key, python_path = python_path, output_path = output_path)
Methods = ['RCTD', 'SPOTlight']
Result = test.Dencon(Methods)
