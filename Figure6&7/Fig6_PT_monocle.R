setwd("/home/data/hli/SHARE_RNA_analysis2/PT_subclustering")
library(Seurat)
library(SeuratData)
library(SeuratDisk)

library(monocle3)
library(ggplot2)
library(dplyr)
novaseq <- LoadH5Seurat("RNA_PT_raw_21375X2110_metabolism.h5seurat")
pData <- data.frame(gene_short_name = colnames(novaseq@assays$RNA@data), 
                    celltype=novaseq@meta.data$celltype,
                    renal_region_new=novaseq@meta.data$renal_region_new,
                    patient_region_id_new=novaseq@meta.data$patient_region_id_new,
                    prep_date=novaseq@meta.data$prep_date,
                    row.names = colnames(novaseq@assays$RNA@data))

fData <- data.frame(gene_short_name = rownames(novaseq@assays$RNA@data),
                    row.names = rownames(novaseq@assays$RNA@data))

cds <- new_cell_data_set(as(novaseq@assays$RNA@data, "sparseMatrix"),
                         cell_metadata = pData,
                         gene_metadata = fData)

cds<-preprocess_cds(
  cds,
  method = "PCA",
  num_dim = 20)

#plot_pc_variance_explained(cds)

cds <- align_cds(cds, alignment_group = "prep_date")
cds <- reduce_dimension(cds,reduction_method='UMAP',core=8,preprocess_method='Aligned',
                        umap.min_dist = 0.001,umap.n_neighbors = 30)

plot_cells(cds, color_cells_by = "celltype",cell_size = 1,alpha=0.6)
#plot_cells(cds, color_cells_by = "prep_date")
cds <- cluster_cells(cds)
plot_cells(cds, color_cells_by = "partition")

cds <- learn_graph(cds,use_partition = FALSE, close_loop=FALSE)

p=plot_cells(cds,
           color_cells_by = "celltype",cell_size = 1, alpha=1,
           trajectory_graph_segment_size=2,
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,label_roots = FALSE)+
  theme(legend.position = "right",legend.text=element_text(size=12),legend.title=element_blank(),
        legend.key = element_rect(fill = "white"))+
  guides(color = guide_legend(override.aes = list(size = 7)))+
  scale_color_manual(values=c('#A5C940','#4F9900', '#013220'))+
  theme(axis.title = element_text(size = 15))

ggsave(filename = 'plots/21375X2110_monocle3_celltype.png', plot = p, width = 8, height = 7, units = "in", dpi = 300)

cds <- order_cells(cds) #choose PT
p=plot_cells(cds,
           color_cells_by = "pseudotime",cell_size = 1, alpha=0.6,
           trajectory_graph_segment_size=2,
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,label_roots = FALSE)+
  theme(axis.title = element_text(size = 15))

ggsave(filename = 'plots/21375X2110_monocle3_pseudotime.png', plot = p,
       width = 8, height = 7, units = "in", dpi = 300)

#save(cds, file='20230111_RNA_PT_raw_21375X2110_metabolism_monocle3.Rdata')


p=plot_genes_in_pseudotime(cds['FAAH2'],color_cells_by="celltype",min_expr=0.5,vertical_jitter=0.2)+
  theme(legend.position = "right",legend.text=element_text(size=12),legend.title=element_blank(),
        legend.key = element_rect(fill = "white"))+
  guides(color = guide_legend(override.aes = list(size = 7)))+
  scale_color_manual(values=c('#A5C940','#4F9900', '#013220'))+
  theme(axis.title = element_text(size = 15))+
  labs(y="Expression", x = "Pseudotime")
ggsave(filename = 'plots/21375X2110_monocle3_FAAH2.png', plot = p,
       width = 8, height = 7, units = "in", dpi = 300)