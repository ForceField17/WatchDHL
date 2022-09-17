# WatchDHL
library(rstudioapi)
library("ggplot2")
library("gridExtra")
library(grid)
# Getting the path of your current open file
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))
print( getwd() )

library(ggplot2)
library(gridExtra)
library(grid)
library(Gviz)
library(rtracklayer)
library(trackViewer)
library(VariantAnnotation)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

Mut <- read.table("Mut_DHL_PAX5.txt",header = F,sep = "\t")
#Mut <- Mut[which(Mut$V)]
Mut$Tag <-  paste0(Mut$V1,": ",Mut$V6)

SNP <- Mut$V3
sample.gr <- GRanges("chr9", IRanges(SNP, width=1, names=Mut$Tag))
sample.gr$color <- Mut$V5
sample.gr$border <- rep("gray30", length(SNP))


features <- GRanges("chr9", IRanges(c(37024550,37027250 ), 
                                    width=c(2700,0),
                                    names=c("block","5UTR")))

features$fill <- c("#e0e0e0", "#e0e0e0")
features$height <- c(0.015,0)

sample.gr.rot <- sample.gr
sample.gr.rot$label.parameter.rot <- 45
black <- gpar(col="black")
red <- gpar(col="red")
Mut$Color <- black
Mut$Color[which(Mut$V9=="Yes")] <- red
sample.gr.rot$label.parameter.gp <- Mut$Color

available.shapes <- c("circle", "square", "diamond", 
                      "triangle_point_up", "triangle_point_down")

sample.gr.rot$shape <- Mut$V4
pdf(file = paste0("./somatic_mutations_3kb_DHL_2022.pdf"),width=17,height=6, onefile = T)
sample.gr.rot$lwd <- 1.5
lolliplot(sample.gr.rot, features, legend=NA, txtSize=0.18,cex=1.4)

dev.off() 







