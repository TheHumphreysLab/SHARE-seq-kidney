import numpy as np
import pandas as pd
import scanpy as sc
import decoupler as dc
import seaborn as sns

sc.settings.verbosity = 3
sc.logging.print_header()
import matplotlib.pyplot as plt
from matplotlib import rcParams
import matplotlib.font_manager
rcParams['font.sans-serif']=['Arial']
sc.settings.set_figure_params(dpi=100, facecolor='white',fontsize=12)

raw_adata=sc.read('RNA_raw.h5ad') #this is the adata file for all RNA cells
adata = dc.get_pseudobulk(raw_adata, sample_col='patient_region_id', groups_col=None,  min_prop=0.05, min_smpls=2,use_raw=False)
sc.pp.filter_genes(adata, min_cells=adata.shape[0]) #the resulting file is n_obs × n_vars = 48 × 8836

sc.pp.normalize_total(adata, target_sum=1e6)
sc.pp.log1p(adata)

sc.tl.pca(adata, svd_solver='arpack')
sc.pl.pca_variance_ratio(adata, log=True)
sc.pp.neighbors(adata, n_neighbors=5,n_pcs=15,metric='cosine')
sc.tl.umap(adata, min_dist=0.5)
sc.pl.umap(adata,color='renal_region',palette=['#4C9150','#7A339E','#E0AB3D','#CC2114','black'],title='RNA gene expression pseudobulk\nclustering (kidney regions)')