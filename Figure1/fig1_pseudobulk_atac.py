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

raw_adata=sc.read('ATAC_raw.h5ad') #this is the adata file for all ATAC cells
adata = dc.get_pseudobulk(raw_adata, sample_col='patient_region_id', groups_col=None,  min_prop=0.05, min_smpls=2,use_raw=False)
sc.pp.filter_genes(adata, min_cells=adata.shape[0])

sc.pp.normalize_total(adata, target_sum=1e6)
sc.pp.log1p(adata)

sc.tl.pca(adata, svd_solver='arpack')
sc.pl.pca_variance_ratio(adata, log=True)
sc.pp.neighbors(adata, n_neighbors=8,n_pcs=20,metric='cosine')
sc.tl.umap(adata, min_dist=1)
sc.pl.umap(adata,color='renal_region_new',palette=['#4C9150','#7A339E','#E0AB3D','#CC2114','black'],title='ATAC gene activity pseudobulk\nclustering (kidney regions)')