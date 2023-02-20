# WatchDHL
library(rstudioapi)
library(stringr)
library(GenomicRanges)
library(circlize)

# Getting the path of your current open file
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))
print( getwd() )

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 150)[1:n]
}

cancer_genes = c("IgH","IgL","BCL2","MYC")

intra_genes = c("XXX")

fusion_genes = c("XXX")

hg38.genes = read.delim("hg38.genes.txt", stringsAsFactors = F)
#hg38.genes = hg38.genes[hg38.genes$gene_type %in% c("protein_coding","lincRNA","intergenic"),]

hg38g = makeGRangesFromDataFrame(df = hg38.genes,keep.extra.columns = T)

ncg = read.delim("IntOGen-DriverGenes.tsv", stringsAsFactors = F)

sv_color_mapper <-function(x){
  y = switch(x,
             "Both" = gg_color_hue(6)[3],
             "DHL" = gg_color_hue(6)[1])
  return(y)
}

sv = read.delim("BCL2_MYC.bedpe", stringsAsFactors = F)
sv = sv[which(sv$ID!="P3"),]
sv_raw = sv 
sv$theID = paste0(" the",sv$ID,":",sv$GENE_A,":",sv$GENE_B)
sv = sv[!duplicated(sv$theID),]


sv_start = sv[,c(1:3,7,8,10,11)]
sv_end = sv[,c(4:6,7,9,10,11)]
sv_col = unlist(lapply(sv$TYPE,sv_color_mapper))


annoteSVwithGene <-function(tmp_sv_bed, hg38g){
  names(tmp_sv_bed)[1:3] = c("chromosome","start","end")
  tmpgr = makeGRangesFromDataFrame(df = tmp_sv_bed[,1:3],keep.extra.columns = T)
  overlapGenes.df <- as.data.frame(distanceToNearest(tmpgr, hg38g, select = "all")) #returned are indexes
  overlapGenes.df = overlapGenes.df[overlapGenes.df$distance < 500,]
  overlapGenes.df = overlapGenes.df[!duplicated(overlapGenes.df$subjectHits),]
  
  res_sv_bed = hg38.genes[overlapGenes.df$subjectHits,]
  return(res_sv_bed)
}


#now visualize with Circlize
cm = colorRamp2(c(-2,0,2),c("blue","grey95","red"),transparency = 0)
names(sv_start)[1:3] = c("chromosome","start","end")
names(sv_end)[1:3] = c("chromosome","start","end")
sv_bed = rbind(sv_end[,c(1,2,3,7)],sv_start[,c(1,2,3,7)])

sv_annot = annoteSVwithGene(tmp_sv_bed = sv_bed, hg38g = hg38g)
sv_lab = sv_annot[,c(1:3,7)]

pdf(file = paste0("./chromosome_SV.pdf"),width=7,height=6, onefile = T)

circos.clear()
chromosome.index = paste0("chr", c(3,8,12,14,18,22))
cytoband = read.cytoband(species = "hg38")$df
circos.par(start.degree=2,gap.degree = 2)
circos.initializeWithIdeogram(cytoband, plotType = NULL, chromosome.index = chromosome.index)
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[2] + mm_y(2), CELL_META$sector.index, cex = 0.8, niceFacing = TRUE)
}, track.height = mm_h(1), cell.padding = c(0, 0, 0, 0), bg.border = NA)

highlight.chromosome(paste0("chr", c(3)), 
                     col = "#1b9e77", track.index = 1)
highlight.chromosome(paste0("chr", c(8)), 
                     col = "#7570b3", track.index = 1)
highlight.chromosome(paste0("chr", c(12)), 
                     col = "#e7298a", track.index = 1)
highlight.chromosome(paste0("chr", c(14)), 
                     col = "#fc8d62", track.index = 1)
highlight.chromosome(paste0("chr", c(18)), 
                     col = "#b3de69", track.index = 1)
highlight.chromosome(paste0("chr", c(22)), 
                     col = "#ffffb3", track.index = 1)

circos.genomicIdeogram(cytoband,track.height = 0.05)

#circos.genomicLink(sv_start, sv_end, col ="grey50", border =  sv_col, lwd = 1.6,h=0.6)
circos.genomicLink(sv_start, sv_end, col ="grey50", border =  sv_col, lwd = 1.6,h=0.8)
#circos.genomicLink(sv_start, sv_end, col ="grey50", border =  sv_col, lwd = 1.6,h=1)

dev.off() 
circos.clear()

 



#FIGURE2 

    
pdf(file = paste0("./CircosPlot_for_All_BCL2_MYC_SV_.pdf"),width=12,height=9, onefile = T)
    
rad51bRDS <-  read.table("SV.hg38.rds",header = T)
head(rad51bRDS)
circos.genomicInitialize(rad51bRDS, axis.labels.cex = 1,labels.cex = 1.5)
circos.track(ylim = c(0, 1),  bg.col = c("#FF000040", "#00FF0040", "#00FF0040", "#0000FF40","#bc80bd","#ff7f00","#80b1d3"),    bg.border = NA, track.height = 0.05)
sv2_start <- data.frame(sv_start$GENE_A,sv_start$start,sv_start$end,sv_start$theID,sv_start$TYPE)
sv2_end <- data.frame(sv_end$GENE_B,sv_end$start,sv_end$end,sv_end$theID,sv_end$TYPE)
colnames(sv2_start) = c("CHROM_A","START_A","END_A","ID","INFO")
colnames(sv2_end) = c("CHROM_B","START_B","END_B","ID","INFO")
circos.genomicLink(sv2_start, sv2_end, col = sv_col, border =  sv_col, lwd = 1.2,h=0.2)
circos.genomicLink(sv2_start, sv2_end, col = sv_col, border =  sv_col, lwd = 1.2,h=0.4)
circos.genomicLink(sv2_start, sv2_end, col = sv_col, border =  sv_col, lwd = 1.2,h=0.6)
circos.genomicLink(sv2_start, sv2_end, col = sv_col, border =  sv_col, lwd = 1.2,h=0.9)
n = max(tapply(rad51bRDS$transcript, rad51bRDS$gene, function(x) length(unique(x))))
circos.genomicTrack(rad51bRDS, ylim = c(0.5, n + 0.5), panel.fun = function(region, value, ...) {
    all_tx = unique(value$transcript)
    for(i in seq_along(all_tx)) {
           l = value$transcript == all_tx[i]
           current_tx_start = min(region[l, 1])
           current_tx_end = max(region[l, 2])
           circos.lines(c(current_tx_start, current_tx_end),  c(n - i + 1, n - i + 1), col = "grey40",lwd = 2)
           circos.genomicRect(region[l, , drop = FALSE], ytop = n - i + 1 + 0.4, 
                              ybottom = n - i + 1 - 0.4, col = "#254880", border = NA)
    }
}, bg.border = NA, track.height = 0.2)

circos.clear()

dev.off() 
    
    
    
    
  
  
  
  

