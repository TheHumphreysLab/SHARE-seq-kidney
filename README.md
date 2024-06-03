## Transcriptomic, epigenomic and spatial metabolomic cell profiling redefines regional human kidney anatomy
<br>
This repository documents the scripts to generate data for our manuscript studying human kidney anatomical regions with SHARE-seq and imaging mass spectrometry (Cell Metabolism 2024). <br>
For citation (DOI: https://doi.org/10.1016/j.cmet.2024.02.015) (PMID: 38513647):
```
Li, H., Li, D., Ledru, N., Xuanyuan, Q., Wu, H., Asthana, A., Byers, L.N., Tullius, S.G., Orlando, G., Waikar, S.S. and Humphreys, B.D., 2024. Transcriptomic, epigenomic, and spatial metabolomic cell profiling redefines regional human kidney anatomy. Cell metabolism, 36(5), pp.1105-1125.
```
<br>
More information about **MALDIpy**, a package for single-cell analysis of MALDI-MS imaging mass spectrometry data, can be found in another GitHub repository: https://github.com/TheHumphreysLab/MALDIpy. <br>

Raw data (.fastq), processed data (count matrix .h5 files and fragment .bed files), sublibrary primer sequences and metadata of the SHARE-seq data have been deposited in NCBI’s Gene Expression Omnibus and are available through GEO Series accession number GSE234788. <br> <br>
A searchable database is available at our Kidney Interactive Transcriptomics (K.I.T.) website: http://humphreyslab.com/SingleCell/.<link> <br><br>
Pre-processing of raw fastq files was performed as previously described (Ma et al. 2020, Cell): https://github.com/masai1116/SHARE-seq-alignment and https://github.com/masai1116/SHARE-seq-alignmentV2/.<br>

### Descriptions

#### 1. Scripts for Figure 1<br>
Peusdobulk analysis for scRNA-seq<br>
Peusdobulk analysis for scATAC-seq<br>
Peusdobulk analysis for co-embedded scRNA-seq and scATAC-seq<br>
Quality control and single-cell clustering (scRNA-seq)<br>
.bed file preprocessing, peak calling and metadata inclusion (scATAC-seq)<br>
Quality control and single-cell clustering (scATAC-seq)<br>
Circular UMAP visualization with <a href="https://github.com/TheHumphreysLab/plot1cell">plot1cell</a><br>


#### 2. Scripts for Figure 2<br>
Single-cell and gene expression visualizations and correlation analysis (scRNA-seq)<br>
Single-cell and gene expression visualizations and correlation analysis (scATAC-seq)<br>
Calculating cell type proportions across samples or groups<br>
Calculating sample or group distributions for each cell type<br>
Calculating the consistency of cell cluster annotation between scRNA-seq and scATAC-seq<br>

#### 3. Scripts for Figure 3<br> 
Single-cell analysis of spatially resolved metabolomics data<br>
Single-cell metabolite feature visualization and tissue projection<br>

#### 4. Scripts for Figure 4<br>
Quality control and single-cell clustering for tL<br>
Quality control and single-cell clustering for TAL<br>
Quality control and single-cell clustering for distal tubular cells<br>
Gene expression visualization and calculating cell type proportions across samples or groups<br>

#### 5. Scripts for Figure 5<br>
Processing a total of 2111 metabolism associated genes only<br>
Quality control and single-cell clustering with metabolism genes<br>
Gene module scoring for metabolic profiles<br>
Data visualization<br>

#### 6. Scripts for Figure 6&7<br>
Acylcarnitine feature analysis in the spatially resolved metabolomics data<br>
Monocle3-based trajectory analysis for PT cells<br>
Marker visualization across predicted pseudotime<br>
Cinical data integration and regression analysis<br>

#### 7. Scripts for other multimodal analysis (multiome_code)<br>
Weighted Nearest Neighbour (WNN)<br>
chromVAR<br>
Cicero<br>
RENIN<br>


<br>
**************<br>


Find us on Twitter: 
<br/>
  <a href="https://twitter.com/HumphreysLab?ref_src=twsrc%5Etfw" class="twitter-follow-button" data-show-count="false"> **@HumphreysLab**</a>
<br/><br/>
