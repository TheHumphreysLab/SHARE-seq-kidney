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

#novaseqs2 file
adata1 = sc.read('share-nova-rna.h5ad')
new_index=[]
for i in adata1.obs.index:
    new_index.append(str(i+',B1'))
adata1.obs.index=new_index

#novaseqs4-1st flowcell
adata2=sc.read_10x_h5("RNA_all.hg19.gene.bc.matrices.h5")
new_index=[]
for i in adata2.obs.index:
    new_index.append(str(i+',B2'))
adata2.obs.index=new_index

#novaseqs4-2nd flowcell lane3
adata3=sc.read_10x_h5("RNAlane3.hg19.gene.bc.matrices.h5")
new_index=[]
for i in adata3.obs.index:
    new_index.append(str(i+',B3'))
adata3.obs.index=new_index

#novaseqs4-2nd flowcell lane4
adata4=sc.read_10x_h5("RNAlane4.hg19.gene.bc.matrices.h5")
new_index=[]
for i in adata4.obs.index:
    new_index.append(str(i+',B4'))
adata4.obs.index=new_index

#combine together
adata_merge=anndata.concat([adata1,adata2,adata3,adata4],join="outer")

new_index=[]
for i in adata_merge.obs.index:
    if i[-2:]=='B1':
        new_index.append(str(i.split('.')[0]+'.'+i.split('.')[1]+','+i.split('.')[2]+'.'+i.split('.')[3]+','+i.split('.')[4]+'.'+i.split('.')[5]+','+i.split('.')[6]+'.'+i.split('.')[7]))
    else:
        new_index.append(i)
adata_merge.obs.index=new_index

sc.pp.filter_cells(adata_merge, min_genes=200)

R1_id=[]
R2_id=[]
R3_id=[]
P1_id=[]
B_id=[]
for i in adata_merge.obs.index:
    R1_id.append(int(i.split(',')[0].split('.')[-1]))
    R2_id.append(int(i.split(',')[1].split('.')[-1]))
    R3_id.append(int(i.split(',')[2].split('.')[-1]))
    P1_id.append(int(i.split(',')[3].split('.')[-1]))
    B_id.append(int(i.split(',')[4][-1]))

adata_merge.obs['R1_id']=R1_id
adata_merge.obs['R2_id']=R2_id
adata_merge.obs['R3_id']=R3_id
adata_merge.obs['P1_id']=P1_id
adata_merge.obs['B_id']=B_id

barcode_error=[]
for i in range(len(adata_merge.obs)):
    if adata_merge.obs['R1_id'][i]>96 or adata_merge.obs['R2_id'][i]>96 or adata_merge.obs['R3_id'][i]>96:
        barcode_error.append(1)
    else:
        barcode_error.append(0)
        
adata_merge.obs['barcode_error']=barcode_error

adata = adata_merge[adata_merge.obs['barcode_error'] == 0, :]

adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
adata = adata[adata.obs.n_genes_by_counts < 5000, :]
sc.pp.filter_genes(adata, min_cells=30)

adata = adata[adata.obs.pct_counts_mt < 4, :]

#final QC
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_cells(adata, max_genes=5000)
sc.pp.filter_cells(adata, min_counts=300)
sc.pp.filter_cells(adata, max_counts=20000)

#doublet imputation
sc.external.pp.scrublet(adata,expected_doublet_rate=0.06,n_neighbors = 30)
adata=adata[adata.obs['doublet_score'] <= 0.2, :]
sc.pp.filter_genes(adata, min_cells=50)

sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata.raw = adata
adata2 = adata[:, adata.var.highly_variable]

sc.pp.regress_out(adata2, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(adata2, max_value=10)
sc.tl.pca(adata2, svd_solver='arpack')

#Run harmony
sc.external.pp.harmony_integrate(adata2, 'prep_date')
adata2.obsm['X_pca_first']=adata2.obsm['X_pca']
adata2.obsm['X_pca'] = adata2.obsm['X_pca_harmony']

sc.pp.neighbors(adata2, n_neighbors=30,n_pcs=30,metric='cosine')
sc.tl.umap(adata2,min_dist=0.1,maxiter=2000)#also cluster QC before this

#sc.tl.leiden(adata2,resolution=1,key_added='leiden1')

#figure production
sc.settings.set_figure_params(dpi=300, facecolor='white',fontsize=12)
sc.pl.umap(adata2, color=['celltype_2023'], alpha=0.1, size=2.5, title='446,267 Human Kidney Cells (RNA)')