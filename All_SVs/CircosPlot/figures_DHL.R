# WatchDHL
library(rstudioapi)

# Getting the path of your current open file
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))
print( getwd() )

cancer_genes = c("IGH","BCL2","MYC","BCL6","AICDA","IGK","IGL")

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 120)[1:n]
}


hg38.genes = read.delim("./hg38.genes.txt", stringsAsFactors = F)
hg38.genes = hg38.genes[!(hg38.genes$gene_type %in% c("lncRNA","processed_pseudogene","misc_RNA")),]
library(GenomicRanges)
hg38g = makeGRangesFromDataFrame(df = hg38.genes,keep.extra.columns = T)


sv_color_mapper <-function(x){
  y = switch(x,
             "BND" = gg_color_hue(6)[6],
             "DEL" = "blue",
             "INS" = gg_color_hue(6)[3],
             "INV" = gg_color_hue(6)[1],
             "DUP" = "orange")
  return(y)
}


annoteSVwithGene <-function(tmp_sv_bed, hg38g){
  names(tmp_sv_bed)[1:3] = c("chromosome","start","end")
  tmpgr = makeGRangesFromDataFrame(df = tmp_sv_bed[,1:3],keep.extra.columns = T)
  overlapGenes.df <- as.data.frame(distanceToNearest(tmpgr, hg38g, select = "all")) #returned are indexes
  overlapGenes.df = overlapGenes.df[overlapGenes.df$distance < 500,]
  overlapGenes.df = overlapGenes.df[!duplicated(overlapGenes.df$subjectHits),]
  
  res_sv_bed = hg38.genes[overlapGenes.df$subjectHits,]
  return(res_sv_bed)
}

######Figure1
datafile <- read.table("DHL.confident_SVs.txt",header = T)
other <- datafile[which(datafile$Type!="BND"),]
transloc <- datafile[which(datafile$Type=="BND"),]


sv2_start <- transloc[,c(1,2,3,8,10,11)]
sv2_end <- transloc[,c(5,6,7,8,10,11)]
sv2_col = unlist(lapply(transloc$Type,sv_color_mapper))

sv1_start <- other[,c(1,2,3,8,10,11)]
sv1_end <- other[,c(5,6,7,8,10,11)]
sv1_col = unlist(lapply(other$Type,sv_color_mapper))

genes <- unique(c(as.character(annoteSVwithGene(sv1_start,hg38g)$gene_name),as.character(annoteSVwithGene(sv2_start,hg38g)$gene_name),
                  as.character(annoteSVwithGene(sv1_end,hg38g)$gene_name),as.character(annoteSVwithGene(sv2_end,hg38g)$gene_name),
                  as.character(datafile$GeneB),as.character(datafile$GeneA)))
tmp_sv_bed <- hg38.genes[which(hg38.genes$gene_name %in% genes),]
sv_lab = tmp_sv_bed[,c(1:3,5)]
#now visualize with Circlize


annoteSVwithGene(sv1_start,hg38g)$gene_name

library(circlize)
cm = colorRamp2(c(-2,0,2),c("blue","grey95","red"),transparency = 0)

names(sv2_start)[1:3] = c("chromosome","start","end")
names(sv2_end)[1:3] = c("chromosome","start","end")
names(sv1_start)[1:3] = c("chromosome","start","end")
names(sv1_end)[1:3] = c("chromosome","start","end")

pdf(file = paste0("./SV_DHL.pdf"),width=10,height=8, onefile = T)

circos.clear()
chromosome.index = paste0("chr", c(1:22,"X","Y"))
cytoband = read.cytoband(species = "hg38")$df
circos.par(start.degree=0,gap.degree = 2)

circos.initializeWithIdeogram(cytoband, plotType = NULL, chromosome.index = chromosome.index)

if (nrow(tmp_sv_bed) > 0){
  if (nrow(sv_lab) > 0){
    circos.genomicLabels(sv_lab, labels.column = 4, side = "outside",cex = 1.2,font = par("font"),line_col = "#b3b3b3",
                         connection_height = convert_height(3, "mm"),labels_height = convert_height(0.4, "cm"),
                         col = ifelse(sv_lab$gene_name %in% cancer_genes,"red","black"))
  }}

circos.track(ylim = c(0, 0.5), panel.fun = function(x, y) {
  chr = substr(CELL_META$sector.index,4,5)
  xlim = CELL_META$xlim
  ylim = CELL_META$ylim
  circos.rect(xlim[1], -0.5, xlim[2], 0.7, col = "white")
  circos.text(mean(xlim), 0.2, chr, cex = 1.4, col = "black",
              facing = "inside", niceFacing = TRUE)
}, track.height = 0.06, bg.border = NA)

circos.genomicIdeogram(cytoband,track.height = 0.05)

circos.genomicLink(sv2_start, sv2_end, col = rand_color(1, transparency = 0.1),
                   border =   sv2_col, lwd = log2(transloc$support*2+1))
circos.genomicLink(sv1_start, sv1_end, col = rand_color(1, transparency = 0.1),
                   border =   sv1_col, lwd = log2(other$support*2+1),h=0.2)

library(stringr)


circos.clear()


dev.off() 
