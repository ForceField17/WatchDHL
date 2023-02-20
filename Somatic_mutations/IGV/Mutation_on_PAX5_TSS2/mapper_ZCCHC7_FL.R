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



SNP <- c(37025150)
sample.gr <- GRanges("chr9", IRanges(SNP, width=1, names=c("P8: Del 1bp")))
sample.gr$color <- c("#fc8d62")
sample.gr$border <- rep("gray30", length(SNP))
sample.gr$label <- NULL#as.character(1:length(sample.gr))
#sample.gr$label.col <- ifelse(sample.gr$alpha>0.5, "white", "black")
sample.gr$label.size <- 2

#chr9:37,024,700-37,027,200


features <- GRanges("chr9", IRanges(c(37024550,37027250 ), 
                                    width=c(2700,0),
                                    names=c("block","5UTR")))

#chr9:37020399-37034510
#chr9:37020636-37020801
features$fill <- c("#e0e0e0", "#2166ac")
features$height <- c(0.015,0)

sample.gr.rot <- sample.gr
sample.gr.rot$label.parameter.rot <- 45
black <- gpar(col="black")
red <- gpar(col="red")
sample.gr.rot$label.parameter.gp <- list(black)

available.shapes <- c("circle", "square", "diamond", 
                      "triangle_point_up", "triangle_point_down")
#sample.gr.rot$shape <- sample(available.shapes, size = length(sample.gr), replace = TRUE)
#sample.gr.rot$fill
sample.gr.rot$shape <- c("triangle_point_up")

pdf(file = paste0("./ExtFig3d_somatic_mutations_3kb_FL.pdf"),width=17,height=6, onefile = T)
sample.gr.rot$lwd <- 1.5
lolliplot(sample.gr.rot, features, legend=NA, txtSize=0.18,cex=1.4)


dev.off() 




