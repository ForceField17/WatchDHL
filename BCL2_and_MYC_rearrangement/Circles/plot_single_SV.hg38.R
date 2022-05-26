setwd("/Users/songdong/project/DLBCL/Basu/SV_hg38/")
cancer_genes = c("IDH1","TP53","CIC","FUBP1","EGFR","PTEN","CDK4","CDKN2A","PDGFRA","MET",
                "FGFR3","TACC3","MDM2","NF1","RB1","MTOR","MYC","MGMT","BCORL1","AKT1","AKT3",
                "SCAI","KIT","NRAS","SOS1","RRAS2","HRAS","KRAS","BCOR","CBL","SPTA1","BCL2","MYC")

intra_genes = c("XXX")

fusion_genes = c("XXX")

hg38.genes = read.delim("hg38.annotation_all.bed", stringsAsFactors = F)
#hg38.genes = hg38.genes[hg38.genes$gene_type %in% c("protein_coding","lincRNA","intergenic"),]
library(GenomicRanges)
hg38g = makeGRangesFromDataFrame(df = hg38.genes,keep.extra.columns = T)

ncg = read.delim("IntOGen-DriverGenes.tsv", stringsAsFactors = F)
#ncg = ncg[ncg$NCG6_oncogene>0 |ncg$NCG6_tsg>0,]

#pairs = read.delim("~/Documents/pangliomaevolution/SV/CGGA_tumor_pairs.txt", stringsAsFactors = F, header =T)

sv_filter <- function(bedpefile ){
  tmp = read.delim(bedpefile, stringsAsFactors = F)
  names(tmp)[1] = "CHROM_A"
  tmp = tmp[tmp$FILTER %in%  c("PASS","MGE10kb"),]
  #tmp = tmp[tmp$FILTER %in%  c("PASS","MGE10kb","MinSomaticScore"),]
  
  for (i in 1:nrow(tmp)){
    fmt = tmp$FORMAT[i] #format 
    nm = tmp[i,15] #normal
    tm = tmp[i,16] #tumor
    if (fmt == "PR:SR"){
      nml = unlist(strsplit(nm, ":")) #normal. before the ":" is PR, and after the ":" is SR
      nmpr = as.integer(unlist(strsplit(nml[1],","))) #normal, paired reads. before the "," is reference depth, while after "," is alternative depth
      nmsr = as.integer(unlist(strsplit(nml[2],","))) #normal, split reads. 
      
      tml = unlist(strsplit(tm, ":")) #tumor
      tmpr = as.integer(unlist(strsplit(tml[1],","))) #tumor, paired reads
      tmsr = as.integer(unlist(strsplit(tml[2],","))) #tumor, split reads
      
      nmprref = nmpr[1] #reference depth in NORMAL sample
      nmpralt = nmpr[2] #alternative depth in NORMAL sample
      tmprref = tmpr[1] #reference depth in TUMOR sample
      tmpralt = tmpr[2] #alternative depth in TUMOR sample
      
      nmsrref = nmsr[1] #reference depth in NORMAL sample
      nmsralt = nmsr[2] #alternative depth in NORMAL sample
      tmsrref = tmsr[1] #reference depth in TUMOR sample
      tmsralt = tmsr[2] #alternative depth in TUMOR sample
    } else if (fmt == "PR"){
      nml = nm #normal. before the ":" is PR, and after the ":" is SR
      nmpr = as.integer(unlist(strsplit(nml[1],","))) #normal, paired reads. before the "," is reference depth, while after "," is alternative depth
      nmsr = c(0,0) #normal, split reads, specified as 0
      
      tml = tm #tumor
      tmpr = as.integer(unlist(strsplit(tml[1],","))) #tumor, paired reads
      tmsr = c(0,0) #tumor, split reads, specified as 0
      
      nmprref = nmpr[1] #reference depth in NORMAL sample
      nmpralt = nmpr[2] #alternative depth in NORMAL sample
      tmprref = tmpr[1] #reference depth in TUMOR sample
      tmpralt = tmpr[2] #alternative depth in TUMOR sample
      
      nmsrref = nmsr[1] #reference depth in NORMAL sample
      nmsralt = nmsr[2] #alternative depth in NORMAL sample
      tmsrref = tmsr[1] #reference depth in TUMOR sample
      tmsralt = tmsr[2] #alternative depth in TUMOR sample
    }
    
    tmp$normal_paired_ref[i] = nmprref
    tmp$normal_paired_alt[i] = nmpralt
    tmp$tumor_paired_ref[i] = tmprref
    tmp$tumor_paired_alt[i] = tmpralt
    
    tmp$normal_split_ref[i] = nmsrref
    tmp$normal_split_alt[i] = nmsralt
    tmp$tumor_split_ref[i] = tmsrref
    tmp$tumor_split_alt[i] = tmsralt
  }
  
  tmp$NORMAL_SUPPORT = tmp$normal_paired_alt + tmp$normal_split_alt
  tmp$TUMOR_SUPPORT = tmp$tumor_paired_alt + tmp$tumor_split_alt
  
  
  tmp = tmp[tmp$TUMOR_SUPPORT >= 5 & tmp$tumor_split_alt >0 & tmp$NORMAL_SUPPORT == 0,]
  #tmp$ID = unlist(lapply(tmp$ID,FUN = function(x) paste(unlist(strsplit(x,":"))[1:7], collapse = ":")))
  #tmp = tmp[!duplicated(tmp$ID),]
  if(nrow(tmp)>0){
    if (substr(tmp$CHROM_A[1], start = 1, stop = 3) != "chr"){
      tmp$CHROM_A = paste0("chr",tmp$CHROM_A)
      tmp$CHROM_B = paste0("chr",tmp$CHROM_B)
    }}
  
  legal_chrs = paste0("chr", c(1:22,"X","Y"))
  tmp = tmp[tmp$CHROM_A %in% legal_chrs & tmp$CHROM_B %in% legal_chrs,]
  
  return(tmp)
}

sv_color_mapper <-function(x){
  y = switch(x,
             "DEL" = "#254880",
             "DUP" = "#FF073A",
             "INV" = "#189A52",
             "BND" = "#F1E952",
             "INS" = "orange")
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


  datafile = paste0("RAD51B.hg38.bedpe") 
  
#  pdf(file = paste0("RAD51B.pdf"),width=12,height=11, onefile = T)

  sv = sv_filter(datafile)
  sv_raw = sv #note the translocations are counted twice
  sv$ID = unlist(lapply(sv$ID,FUN = function(x) paste(unlist(strsplit(x,":"))[1:7], collapse = ":")))
  sv = sv[!duplicated(sv$ID),]

  sv1 = sv[!(sv$TYPE=="BND"),]
  sv2 = sv[(sv$TYPE=="BND"),]
  
  sv_start = sv1[,c(1:3,7,26,9)]
  sv_end = sv1[,c(4:6,7,26,10)]
  sv_col = unlist(lapply(sv1$TYPE,sv_color_mapper))
  
  sv2_start = sv2[,c(1:3,7,26,9)]
  sv2_end = sv2[,c(4:6,7,26,10)]
  sv2_col = unlist(lapply(sv2$TYPE,sv_color_mapper))
  
  #now visualize with Circlize
  
  sample_id = paste0("DLBCL")

  library(circlize)
   cm = colorRamp2(c(-2,0,2),c("blue","grey95","red"),transparency = 0)
    names(sv_start)[1:3] = c("chromosome","start","end")
    names(sv_end)[1:3] = c("chromosome","start","end")
    tmp_sv_bed = rbind( cbind(sv_start[!startsWith(sv_start$ID,prefix = "MantaBND"),],
                              sv_end[!startsWith(sv_end$ID,prefix = "MantaBND"),])[,c(1,2,9,4)])
    tmp_sv_bed = tmp_sv_bed[,1:3]
    #tmp_sv_bed = tmp_sv_bed[!startsWith(tmp_sv_bed$ID, "MantaINV"),1:3]
    if (nrow(tmp_sv_bed) > 0){
      sv_annot = annoteSVwithGene(tmp_sv_bed = tmp_sv_bed, hg38g = hg38g)
      sv_lab = sv_annot[(sv_annot$gene_name %in% ncg$symbol | sv_annot$gene_name %in% cancer_genes ),c(1:3,7)]
      sv_lab = sv_lab[!duplicated(sv_lab$gene_name),]
    }
    
    names(sv2_start)[1:3] = c("chromosome","start","end")
    names(sv2_end)[1:3] = c("chromosome","start","end")
    tmp_sv_bed2 = rbind(sv2_end[startsWith(sv2_end$ID,prefix = "MantaBND"),1:4],
                        sv2_start[startsWith(sv2_start$ID,prefix = "MantaBND"),1:4])
    if (nrow(tmp_sv_bed2) > 0){
      sv_annot2 = annoteSVwithGene(tmp_sv_bed = tmp_sv_bed2, hg38g = hg38g)
      #sv_lab2 = sv_annot2[,c(1:3,7)]
      sv_lab2 = sv_annot2[!(grepl("RP",sv_annot2$gene_name)),c(1:3,7)]
      sv_lab2 = sv_lab2[!duplicated(sv_lab2$gene_name),]
    }
    

    circos.clear()
    chromosome.index = paste0("chr", c(1,3,6,14))
    cytoband = read.cytoband(species = "hg38")$df
    circos.par(start.degree=2,gap.degree = 2)
    circos.initializeWithIdeogram(cytoband, plotType = NULL, 
                                  chromosome.index = chromosome.index)
    circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
      circos.text(CELL_META$xcenter, CELL_META$ylim[2] + mm_y(2), 
                  CELL_META$sector.index, cex = 0.8, niceFacing = TRUE)
    }, track.height = mm_h(1), cell.padding = c(0, 0, 0, 0), bg.border = NA)
    highlight.chromosome(paste0("chr", c(1)), 
                         col = "#1b9e77", track.index = 1)
    highlight.chromosome(paste0("chr", c(3)), 
                         col = "#7570b3", track.index = 1)
    highlight.chromosome(paste0("chr", c(6)), 
                         col = "#e7298a", track.index = 1)
    highlight.chromosome(paste0("chr", c(14)), 
                         col = "#fc8d62", track.index = 1)
    
    circos.genomicIdeogram(cytoband,track.height = 0.05)

    
    
    
    
    
    
    

    #circos.genomicDensity(covs[covs$seg.mean > 0.4,], col = c("#FF000080"), track.height = 0.1,window.size = 100000)
    #circos.genomicDensity(covs[covs$seg.mean < -0.4,], col = c("#0000FF80"), track.height = 0.1,window.size = 1000000)
    
    circos.genomicLink(sv_start, sv_end, col ="grey50", 
                       border =  sv_col, lwd = 1,h=0.07)
    
    my_col <- c("#66c2a5","#8da0cb","#8da0cb","#e78ac3")
    circos.genomicLink(sv2_start, sv2_end, col ="grey50", 
                       border =   my_col, lwd = 3)
    
    library(stringr)
    
    text(0, 0, id, cex = 1, col = "black")
    if (nrow(tmp_sv_bed2) > 0){
      if (nrow(sv_lab2) > 0){
        circos.genomicLabels(sv_lab2, labels.column = 4, side = "inside",cex = 0.6,font = par("font"),line_col = "#b3b3b3",
                             connection_height = convert_height(2, "mm"),labels_height = convert_height(0.2, "cm"),
                             col = ifelse(sv_lab2$gene_name %in% ncg$symbol,"red",ifelse(sv_lab2$gene_name %in% cancer_genes,"red",ifelse(str_split(sv_lab2$gene_name,";",simplify = T)[,1] %in% ncg$symbol,"blue",ifelse(str_split(sv_lab2$gene_name,";",simplify = T)[,2] %in% ncg$symbol,"blue",ifelse(grepl("IG",sv_lab2$gene_name),"blue","black"))))))
      }}
    
    circos.clear()
 #   legend("topright",fill = c("#254880","#FF073A","#189A52","#F1E952","orange"), legend = c("deletions","tandem duplications","inversions","translocations","insertions"), bty = "n",cex = 0.8)

 
#FIGURE2 
    df = data.frame(
      name  = c("TP53",  "TP63",    "TP73"),
      start = c(0, 189349205, 3569084),
      end   = c(7590856, 189615068, 3652765))
    circos.genomicInitialize(df)
    
    tp_family = readRDS(system.file(package = "circlize", "extdata", "tp_family_df.rds"))
    head(tp_family)
    
    rad51bRDS <-  read.table("RAD51B.hg38.rds",header = T)
    head(rad51bRDS)
    
    
    circos.genomicInitialize(rad51bRDS)
    circos.track(ylim = c(0, 1), 
                 bg.col = c("#FF000040", "#00FF0040", "#0000FF40","#254880","#FF073A","#189A52"), 
                 bg.border = NA, track.height = 0.05)
    
    n = max(tapply(rad51bRDS$transcript, rad51bRDS$gene, function(x) length(unique(x))))
    circos.genomicTrack(rad51bRDS, ylim = c(0.5, n + 0.5), 
                        panel.fun = function(region, value, ...) {
                          all_tx = unique(value$transcript)
                          for(i in seq_along(all_tx)) {
                            l = value$transcript == all_tx[i]
                            # for each transcript
                            current_tx_start = min(region[l, 1])
                            current_tx_end = max(region[l, 2])
                            circos.lines(c(current_tx_start, current_tx_end), 
                                         c(n - i + 1, n - i + 1), col = "#CCCCCC")
                            circos.genomicRect(region[l, , drop = FALSE], ytop = n - i + 1 + 0.4, 
                                               ybottom = n - i + 1 - 0.4, col = "orange", border = NA)
                          }
                        }, bg.border = NA, track.height = 0.4)
    
    my_col <- c("#66c2a5","#8da0cb","#8da0cb","#e78ac3")
    circos.genomicLink(sv2_start, sv2_end, col ="grey50", 
                       border =   my_col, lwd = 3)
    circos.clear()

    
    
    
  
  
  
  
  
 ## dev.off() 

