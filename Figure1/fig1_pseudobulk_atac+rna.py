import numpy as np
import pandas as pd
import scanpy as sc
sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
from matplotlib import rcParams
import matplotlib.font_manager
rcParams['font.sans-serif']=['Arial']
sc.settings.set_figure_params(dpi=100, facecolor='white',fontsize=12)

import anndata
from anndata import AnnData
import decoupler as dc

adata_atac=sc.read('pseudobulk_ATAC.h5ad')#pseudobulk data without normalization
adata_atac_index=[]
for i in adata_atac.obs.index:
    adata_atac_index.append(str(str('ATAC_')+str(i)))
adata_atac.obs.index=adata_atac_index

adata_rna=sc.read('pseudobulk_RNA.h5ad')#pseudobulk data without normalization
adata_rna_index=[]
for i in adata_rna.obs.index:
    adata_rna_index.append(str(str('RNA_')+str(i)))
adata_rna.obs.index=adata_rna_index

#combine ATAC and RNA seq data
adata_merge=anndata.concat([adata_rna,adata_atac],join="inner")

module=[]
for i in adata_merge.obs.index:
    if i[:3] == 'RNA':
        module.append('RNA')
    elif i[:4] == 'ATAC':
        module.append('ATAC')
adata_merge.obs['module']=module

sc.pp.filter_genes(adata_merge, min_cells=adata_merge.shape[0])
dc.plot_violins(adata_merge, log=True,figsize=(30, 5), dpi=200)#check normalization
sc.pp.normalize_total(adata_merge, target_sum=1e6)
sc.pp.log1p(adata_merge)
dc.plot_violins(adata_merge, title='Normalized counts per sample',figsize=(30,5), dpi=200)

sc.tl.pca(adata_merge, svd_solver='arpack')
sc.pl.pca_variance_ratio(adata_merge, log=True)

sc.pp.neighbors(adata_merge, n_neighbors=10,n_pcs=30,metric='cosine')
sc.tl.umap(adata_merge, min_dist=1.5)
sc.pl.umap(adata_merge,color=['renal_region','module'],palette=['#4C9150','#7A339E','#E0AB3D','#CC2114','black'])
sc.pl.umap(adata_merge,color='module',palette=['red','blue'],title='Pseudobulk clustering of combined\nRNA+ATAC data (modality)')