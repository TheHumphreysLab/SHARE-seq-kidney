library(Seurat)
library(SeuratDisk)
novaseq <- LoadH5Seurat("rna.h5seurat")#or load your Rdata
library(plot1cell)
library(plotrix)

novaseq<-SetIdent(novaseq,value = "celltype_2023")
###Check and see the meta data info on your Seurat object
colnames(novaseq@meta.data)  

###Prepare data for ploting

set.seed(1234)
cluster_colors<-c('#7D4729',
                           '#8a3e6a',
                           '#be658d',
                           '#F9CC72',
                           '#E2062B',
                           '#860111',
                           '#B4041E',
                           '#617A2E',
                           '#A57C00',
                           '#FF8933',
                           '#86DEBB',
                           '#6a3070',
                           '#4c2564',
                           '#0077BE',
                           '#00B5EB',
                           '#8DC71E',
                           '#69B41E',
                           '#013220',
                           '#E97E88',
                           '#F8D1CD',
                           '#E15566',
                           '#128394',
                           '#62CCCC',
                           '#046494',
                           '#C9F5E6',
                           '#936210',
                           '#5E2A0F', 
                           '#092092',
                           '#1C3BAC')

group_colors<-c('#4c9150', '#7a339e', '#e0ab3d', '#cc2114', '#000000')

plot_circlize2 <- function(
    data_plot,
    do.label = T,
    contour.levels = c(0.2,0.3),
    pt.size = 0.5,
    kde2d.n = 1000,
    contour.nlevels = 100,
    bg.color='#F9F2E4',
    col.use=NULL,
    label.cex = 0.5,
    repel=FALSE
) {
  data_plot %>%
    dplyr::group_by(Cluster) %>%
    summarise(x = median(x = x), y = median(x = y)) -> centers
  z <- MASS::kde2d(data_plot$x, data_plot$y, n=1000,lims = c(-1, 1, -1,1))
  celltypes<-names(table(data_plot$Cluster))
  cell_colors <- scales::hue_pal()(length(celltypes))
  if(!is.null(col.use)){
    cell_colors=col.use
    col_df<-data.frame(Cluster=celltypes, color2=col.use)
    cells_order<-rownames(data_plot)
    data_plot<-merge(data_plot, col_df, by="Cluster")
    rownames(data_plot)<-data_plot$cells
    data_plot<-data_plot[cells_order,]
    data_plot$Colors<-data_plot$color2
  }
  circos.clear()
  par(bg = bg.color)
  circos.par(cell.padding=c(0,0,0,0), track.margin=c(0.01,0),"track.height" = 0.01, gap.degree =c(rep(2, (length(celltypes)-1)),12),points.overflow.warning=FALSE)
  circos.initialize(sectors =  data_plot$Cluster, x = data_plot$x_polar2)
  circos.track(data_plot$Cluster, data_plot$x_polar2, y=data_plot$dim2, bg.border=NA,panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter,
                CELL_META$cell.ylim[2]+ mm_y(4.8),
                CELL_META$sector.index,
                cex=0.75, col = 'black', facing = "bending.inside", niceFacing = T)
    circos.axis(labels.cex = 0.38, col = 'black', labels.col =  'black')
  })
  for(i in 1:length(celltypes)){
    dd<-data_plot[data_plot$Cluster==celltypes[i],]
    circos.segments(x0 = min(dd$x_polar2), y0 = 0, x1 = max(dd$x_polar2), y1 = 0, col = cell_colors[i],  lwd=3, sector.index = celltypes[i])
  }
  contour(z, drawlabels=F, levels  =  c(0.2,0.4),col = '#d4b981', add=TRUE)
  text(x = 1, y=0.105, labels = "Cluster", cex = 0.65, col = 'black',srt=-90)
  points(data_plot$x,data_plot$y, pch = 19, col = alpha(data_plot$Colors,0.1), cex = pt.size);
  
  if(do.label){
    if(repel){
      textplot(x=centers$x, y=centers$y, words =  centers$Cluster,cex = label.cex, new = F,show.lines=F)
    } else {
      text(centers$x,centers$y, labels=centers$Cluster, cex = label.cex, col = 'black')
    }
  } 
}

add_track2 <- function(
    data_plot, 
    group, 
    track_num, 
    colors = NULL
){
  if(track_num<2){
    stop("The first track is the cluster track. Please change the track_num to a value greater than 1")
  }
  circos.track(data_plot$Cluster, data_plot$x_polar2, y=data_plot$dim2, bg.border=NA)
  celltypes<-names(table(data_plot$Cluster))
  group_names<-names(table(data_plot[,group]))
  if(is.null(colors)){
    col_group = scales::hue_pal()(length(group_names))
  } else {
    col_group = colors
  }
  for(i in 1:length(celltypes)) {
    data_plot_cl<-data_plot[data_plot$Cluster==celltypes[i],]
    dat_seg<-get_segment(data_plot_cl, group = group)
    dat_seg2<-c(dat_seg[-1]-1, nrow(data_plot_cl))
    scale_factor<-max(data_plot_cl$x_polar2)/nrow(data_plot_cl)
    dat_seg<-scale_factor*dat_seg
    dat_seg2<-scale_factor*dat_seg2
    circos.segments(x0 = dat_seg, y0 = 0, x1 = dat_seg2, y1 = 0, col = col_group, sector.index = celltypes[i], lwd=3)
  }
}


circ_data <- prepare_circlize_data(novaseq, scale = 0.7)
pdf(file = "rna_circular_use.pdf")
plot_circlize2(circ_data,do.label = F, pt.size = 0.01, col.use = cluster_colors ,bg.color = 'white', kde2d.n = 200, repel = T, label.cex = 0.6)
add_track2(circ_data, group = "renal_region_new", colors = group_colors, track_num = 2) ## can change it to one of the columns in the meta data of your seurat object
text(x = 0.955, y=0.1, labels = "Region", cex = 0.65, col = 'black',srt=-90)
arctext(x = "446,267 cells | RNA-seq", center = c(0, 0), radius = 0.9, middle = 7*pi/4 , clockwise = FALSE,
        cex = 1.2, stretch = 1)
dev.off()



#for ATAC:
load("SHARE_ATAC_analysis.Rdata")
cluster_colors <- c('#7D4729',
                             '#8a3e6a',
                             '#be658d',
                             '#F9CC72',
                             '#E2062B',
                             '#860111',
                             '#B4041E',
                             '#617A2E',
                             '#86DEBB',
                             '#6a3070',
                             '#4c2564',
                             '#00B5EB',
                             '#8DC71E',
                             '#013220',
                             '#F8D1CD',
                             '#128394',
                             '#046494',
                             '#C9F5E6',
                             '#936210',
                             '#5E2A0F',
                             '#092092')
                             
group_colors<-c('#4c9150', '#7a339e', '#e0ab3d', '#cc2114', '#000000')
plot_circlize2 <- function(
    data_plot,
    do.label = T,
    contour.levels = c(0.2,0.3),
    pt.size = 0.5,
    kde2d.n = 1000,
    contour.nlevels = 100,
    bg.color='#F9F2E4',
    col.use=NULL,
    label.cex = 0.5,
    repel=FALSE
) {
  data_plot %>%
    dplyr::group_by(Cluster) %>%
    summarise(x = median(x = x), y = median(x = y)) -> centers
  z <- MASS::kde2d(data_plot$x, data_plot$y, n=1000,lims = c(-1, 1, -1,1))
  celltypes<-names(table(data_plot$Cluster))
  cell_colors <- scales::hue_pal()(length(celltypes))
  if(!is.null(col.use)){
    cell_colors=col.use
    col_df<-data.frame(Cluster=celltypes, color2=col.use)
    cells_order<-rownames(data_plot)
    data_plot<-merge(data_plot, col_df, by="Cluster")
    rownames(data_plot)<-data_plot$cells
    data_plot<-data_plot[cells_order,]
    data_plot$Colors<-data_plot$color2
  }
  circos.clear()
  par(bg = bg.color)
  circos.par(cell.padding=c(0,0,0,0), track.margin=c(0.01,0),"track.height" = 0.01, gap.degree =c(rep(2, (length(celltypes)-1)),12),points.overflow.warning=FALSE)
  circos.initialize(sectors =  data_plot$Cluster, x = data_plot$x_polar2)
  circos.track(data_plot$Cluster, data_plot$x_polar2, y=data_plot$dim2, bg.border=NA,panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter,
                CELL_META$cell.ylim[2]+ mm_y(4.8),
                CELL_META$sector.index,
                cex=0.75, col = 'black', facing = "bending.inside", niceFacing = T)
    circos.axis(labels.cex = 0.4, col = 'black', labels.col =  'black')
  })
  for(i in 1:length(celltypes)){
    dd<-data_plot[data_plot$Cluster==celltypes[i],]
    circos.segments(x0 = min(dd$x_polar2), y0 = 0, x1 = max(dd$x_polar2), y1 = 0, col = cell_colors[i],  lwd=3, sector.index = celltypes[i])
  }
  contour(z, drawlabels=F, levels  =  c(0.2,0.1),col = '#d4b981', add=TRUE)
  text(x = 1, y=0.105, labels = "Cluster", cex = 0.65, col = 'black',srt=-90)
  points(data_plot$x,data_plot$y, pch = 19, col = alpha(data_plot$Colors,0.2), cex = pt.size);
  
  if(do.label){
    if(repel){
      textplot(x=centers$x, y=centers$y, words =  centers$Cluster,cex = label.cex, new = F,show.lines=F)
    } else {
      text(centers$x,centers$y, labels=centers$Cluster, cex = label.cex, col = 'black')
    }
  } 
}

add_track2 <- function(
    data_plot, 
    group, 
    track_num, 
    colors = NULL
){
  if(track_num<2){
    stop("The first track is the cluster track. Please change the track_num to a value greater than 1")
  }
  circos.track(data_plot$Cluster, data_plot$x_polar2, y=data_plot$dim2, bg.border=NA)
  celltypes<-names(table(data_plot$Cluster))
  group_names<-names(table(data_plot[,group]))
  if(is.null(colors)){
    col_group = scales::hue_pal()(length(group_names))
  } else {
    col_group = colors
  }
  for(i in 1:length(celltypes)) {
    data_plot_cl<-data_plot[data_plot$Cluster==celltypes[i],]
    dat_seg<-get_segment(data_plot_cl, group = group)
    dat_seg2<-c(dat_seg[-1]-1, nrow(data_plot_cl))
    scale_factor<-max(data_plot_cl$x_polar2)/nrow(data_plot_cl)
    dat_seg<-scale_factor*dat_seg
    dat_seg2<-scale_factor*dat_seg2
    circos.segments(x0 = dat_seg, y0 = 0, x1 = dat_seg2, y1 = 0, col = col_group, sector.index = celltypes[i], lwd=3)
  }
  #text(x = (1-0.03*(track_num-1)), y=0.1, labels = group, cex = 0.4, col = 'black',srt=-90)
}

novaseq<-SetIdent(novaseq,value = "celltype_2023")
###plot and save figures

circ_data <- prepare_circlize_data(novaseq, scale = 0.86)
#png(filename =  '20230604_test_circlize_plot.png', width = 6, height = 6,units = 'in', res = 300)
pdf(file = "atac_circular_use.pdf")
plot_circlize2(circ_data,do.label = F, pt.size = 0.01, col.use = cluster_colors ,bg.color = 'white', kde2d.n = 1000, repel = T, label.cex = 0.6)
add_track2(circ_data, group = "renal_region_new", colors = group_colors, track_num = 2) ## can change it to one of the columns in the meta data of your seurat object
text(x = 0.955, y=0.1, labels = "Region", cex = 0.65, col = 'black',srt=-90)
arctext(x = "401,875 cells | ATAC-seq", center = c(0, 0), radius = 0.9, middle = 7*pi/4 , clockwise = FALSE,
        cex = 1.2, stretch = 1)
dev.off()