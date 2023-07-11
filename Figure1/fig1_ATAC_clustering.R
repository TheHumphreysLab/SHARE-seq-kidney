library(Signac)
library(Seurat)
library(GenomeInfoDb)
#BiocManager::install(c('BSgenome.Hsapiens.UCSC.hg19', 'EnsDb.Hsapiens.v75'))
library(EnsDb.Hsapiens.v75)
library(ggplot2)
library(patchwork)
set.seed(1234)

fragments<-CreateFragmentObject('../processed_data/ATAC.hg19.fragments.B123.sorted.tsv.gz')
features<-CallPeaks(fragments,format='BED',outdir = "../processed_data/",name='name',extsize=150,
                    macs2.path ='/home/users/hli/anaconda3/envs/py2/bin/macs2',cleanup=FALSE)

print("step3: FeatureMatrix")
print(Sys.time())
matrix<-FeatureMatrix(fragments,features,process_n = 100000)

chrom_assay <- CreateChromatinAssay(
  counts = matrix,
  genome = 'hg19',
  fragments = 'ATAC.hg19.fragments.B123.sorted.tsv.gz',
    min.cells = 1, ranges=features)

novaseq <- CreateSeuratObject(
  counts = chrom_assay,
  assay = 'peaks')

# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75)

# change to UCSC style since the data was mapped to hg19
seqlevelsStyle(annotations) <- 'UCSC'
# add the gene information to the object
Annotation(novaseq) <- annotations

# compute nucleosome signal score per cell
novaseq <- NucleosomeSignal(object = novaseq)

# compute TSS enrichment score per cell
novaseq <- TSSEnrichment(object = novaseq, fast = FALSE)

novaseq$high.tss <- ifelse(novaseq$TSS.enrichment > 1, 'High', 'Low')
TSSPlot(novaseq, group.by = 'high.tss') + NoLegend()

novaseq$nucleosome_group <- ifelse(novaseq$nucleosome_signal > 2.5, 'NS > 2.5', 'NS < 2.5')
FragmentHistogram(object = novaseq, group.by = 'nucleosome_group')

VlnPlot(object = novaseq,
  features = c('nCount_peaks','nFeature_peaks','TSS.enrichment', 'nucleosome_signal'),
  pt.size = 0.1,ncol =4)

total_fragments <- CountFragments('ATAC.hg19.fragments.B123.sorted.tsv.gz')
rownames(total_fragments)<-total_fragments$CB
head(total_fragments)

novaseq$fragments <- total_fragments[colnames(novaseq), "frequency_count"]
novaseq$reads_count <- total_fragments[colnames(novaseq), "reads_count"]
novaseq$mononucleosomal <- total_fragments[colnames(novaseq), "mononucleosomal"]
novaseq$nucleosome_free <- total_fragments[colnames(novaseq), "nucleosome_free"]

novaseq <- FRiP(
  object = novaseq,
  assay = 'peaks',
  total.fragments = 'fragments'
)
VlnPlot(object = novaseq,
  features = c('fragments','FRiP','reads_count'),
  pt.size = 0)

novaseq$blacklist_fraction <- FractionCountsInRegion(
  object = novaseq, 
  assay = 'peaks',
  regions = blacklist_hg19
)



print("First QC")
novaseq <- subset(x = novaseq,
  subset = nCount_peaks > 400 &
                  nCount_peaks < 50000 &
                  FRiP > 0.1 &
                  blacklist_fraction < 0.05 &
                  nucleosome_signal < 2.5 &
                  TSS.enrichment > 1)
length(novaseq@meta.data$orig.ident)

load("atac_seurat_object_raw_meta_table.Rdata")##here we load a metatable, since we want to add some RNA-seq metadata into the ATAC-seq data of the same cells.
metadataDF = meta_tbl_all[,c("sample_id", "batch", "atac_prep_date", "renal_region",
                             "n_genes_by_counts", "total_counts", "total_counts_mt", "pct_counts_mt",
                             "celltype4", "R1_id","rna_barcode",'doublet_score')]
rownames(metadataDF) = meta_tbl_all$atac_barcode

metadataDF2<-metadataDF[colnames(novaseq),]

novaseq <- AddMetaData(novaseq, metadata = metadataDF2)

##remove doublet > 0.2
print("Removing RNA doublets")
novaseq <- subset(x = novaseq, subset = doublet_score <= 0.2)
length(novaseq@meta.data$orig.ident)


##remove ids>96
print("Removing error cells with ids>96")
novaseq_rownames = data.frame(row_name = rownames(novaseq@meta.data))
test = lapply(novaseq_rownames$row_name, FUN = function(x){
  tmp = strsplit(x, split = ",")[[1]]
  tmp1 = (substr(tmp[1], 4,6))
  tmp2 = (substr(tmp[2], 4,6))
  tmp3 = (substr(tmp[3], 4,6))
  tmp4 = (substr(tmp[4], 4,5))
  c(tmp1, tmp2, tmp3, tmp4)
})

test2 = as.data.frame(do.call(rbind, test))
novaseq_rownames = cbind(novaseq_rownames, test2)
for (i in 2:5){
  novaseq_rownames[,i] = as.numeric(novaseq_rownames[,i])
}

novaseq_rownames$remove = FALSE
novaseq_rownames$remove = sapply(1:dim(novaseq_rownames)[1], FUN = function(x){
  if (novaseq_rownames[x, "V1"] > 96 | novaseq_rownames[x, "V2"] > 96 | 
      novaseq_rownames[x, "V3"] > 96 | novaseq_rownames[x, "V4"] > 96){
    return(TRUE)
  } else {
    return(FALSE)
  }
})

cell_barcode_keep = novaseq_rownames$row_name[!novaseq_rownames$remove]
novaseq = subset(novaseq, cells = cell_barcode_keep)
length(novaseq@meta.data$orig.ident)

print("RunTFIDF")
novaseq <- RunTFIDF(novaseq)
novaseq <- FindTopFeatures(novaseq, min.cutoff = 'q0')
novaseq <- RunSVD(novaseq)
#save(novaseq,file='SHARE_ATAC_analysis1/harmony_before.Rdata')

#Run Harmony
library(harmony)
novaseq <- RunHarmony(
  object = novaseq,
  group.by.vars = 'atac_prep_date',
  reduction = 'lsi',
  assay.use = 'peaks',
  project.dim = FALSE
)

#save(novaseq,file='SHARE_ATAC_analysis1/harmony.Rdata')

novaseq <- RunUMAP(object = novaseq, reduction = 'harmony', dims = 2:25, min.dist = 0.1, n.neighbors=30)

#novaseq <- FindClusters(object = novaseq, algorithm = 3)