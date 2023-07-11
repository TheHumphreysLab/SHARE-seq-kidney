import numpy as np
import pandas as pd
import scanpy as sc

sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
from matplotlib import rcParams
import matplotlib.font_manager
rcParams['font.sans-serif']=['Arial']
sc.settings.set_figure_params(dpi=100, facecolor='white',fontsize=12)
import matplotlib.pyplot as plt
import seaborn as sns

adata=sc.read('RNA_446267cells_raw.h5ad')

tal_use=[]
for i in adata.obs['celltype_2023']:
    if i in ['TAL1','TAL2','TAL3','tL-TAL']:
        tal_use.append(1)
    else:
        tal_use.append(0)
adata.obs['tal_use']=tal_use
adata2 = adata[adata.obs.tal_use == 1, :]

dis_use=[]
for i in adata2.obs['renal_region_new']:
    if i in ['RA','U']:
        dis_use.append(0)
    else:
        dis_use.append(1)
adata2.obs['dis_use']=dis_use
adata2 = adata2[adata2.obs.dis_use == 1, :]

sc.pp.filter_genes(adata2, min_cells=30)
adata2.var['mt'] = adata2.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata2, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
sc.pp.normalize_total(adata2)
sc.pp.log1p(adata2)
sc.pp.highly_variable_genes(adata2, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata2.raw = adata2
adata3 = adata2[:, adata2.var.highly_variable]
sc.pp.regress_out(adata3, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(adata3, max_value=10)
sc.tl.pca(adata3, svd_solver='arpack')

sc.external.pp.harmony_integrate(adata3, 'prep_date')
adata3.obsm['X_pca_first']=adata3.obsm['X_pca']
adata3.obsm['X_pca'] = adata3.obsm['X_pca_harmony']

sc.pp.neighbors(adata3, n_neighbors=30,n_pcs=30)
sc.tl.umap(adata3,min_dist=0.1)
sc.pl.umap(adata3, color=['renal_region_new'])