# WatchDHL
library(rstudioapi)

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

available.shapes <- c("circle", "square", "diamond", "triangle_point_up", "triangle_point_down")

Anno <- read.table("gencode_annotation.txt",sep = "\t",header = T)
Anno1 <- Anno[which(Anno$chr == "chr18"),]
Anno2 <- Anno[which(Anno$chr == "chr8"),]
Anno3 <- Anno[which(Anno$chr == "chr14"),]

SVbp <- read.table("breakpoint_position.txt",sep = "\t",header = T)
SVbp$color <- gg_color_hue(6)[3]
SVbp$color[which(SVbp$Tumor == "DHL")] <- gg_color_hue(6)[1]
SVbp <- SVbp[which(SVbp$patient != "P3"),]


##### BCL2 locus
Loci1 <- SVbp[which(SVbp$Locus == "BCL2"),]

POS <- Loci1$position
NAME <- paste0(Loci1$Target," (",Loci1$patient,")")
COLOR <- Loci1$color

sample.gr <- GRanges("chr18", IRanges(POS, width=1, names=NAME))
sample.gr$color <- COLOR

sample.gr$border <- "gray30"
sample.gr$label <- NULL

START <- Anno1$start
END <- Anno1$end
MIN <- min(START,POS,END) - 5000
MAX <- max(START,POS,END) 
VECT <- c(MIN,START,MAX)
WIDTH <- c(0,END-START,5000)
LABEL <- c("",Anno1$type,"")

features <- GRanges("chr18", IRanges(VECT,  width=WIDTH,   names=LABEL ))

FILL <- rep(NA,length(LABEL))
FILL[which(LABEL != "")] <- "#2166ac"
features$fill <- FILL 
features$color <- FILL 

HEI <- rep(0,length(LABEL))
HEI[which(LABEL != "")] <- 0.03
features$height <- HEI


chr18 <- sample.gr
chr18$label.parameter.rot <- 90
chr18$lwd <- 1
chr18$shape <- "diamond"


##### MYC locus
Loci2 <- SVbp[which(SVbp$Locus == "MYC"),]

POS <- Loci2$position
NAME <- paste0(Loci2$Target," (",Loci2$patient,")")
COLOR <- Loci2$color

sample.gr <- GRanges("chr8", IRanges(POS, width=1, names=NAME))
sample.gr$color <- COLOR

sample.gr$border <- "gray30"
sample.gr$label <- NULL

START <- Anno2$start
END <- Anno2$end
MIN <- min(START,POS,END) - 5000
MAX <- max(START,POS,END) 
VECT <- c(MIN,START,MAX)
WIDTH <- c(0,END-START,5000)
LABEL <- c("",Anno2$type,"")

features2 <- GRanges("chr8", IRanges(VECT,  width=WIDTH,   names=LABEL ))

FILL <- rep(NA,length(LABEL))
FILL[which(LABEL != "")] <- "#2166ac"
features2$fill <- FILL 
features2$color <- FILL 

HEI <- rep(0,length(LABEL))
HEI[which(LABEL != "")] <- 0.03
features2$height <- HEI

chr8 <- sample.gr
chr8$label.parameter.rot <- 90
chr8$lwd <- 1
chr8$shape <- "diamond"

##### IgH locus
Loci3 <- SVbp[which(SVbp$Locus == "IgH"),]

POS <- Loci3$position
NAME <- paste0(Loci3$Target," (",Loci3$patient,")")
COLOR <- Loci3$color

sample.gr <- GRanges("chr14", IRanges(POS, width=1, names=NAME))
sample.gr$color <- COLOR

sample.gr$border <- "gray30"
sample.gr$label <- NULL

START <- Anno3$start
END <- Anno3$end
MIN <- min(START,POS,END) - 5000
MAX <- max(START,POS,END) 
VECT <- c(MIN,START,MAX)
WIDTH <- c(0,END-START,5000)
LABEL <- c("",Anno3$type,"")

features3 <- GRanges("chr14", IRanges(VECT,  width=WIDTH,   names=LABEL ))

FILL <- rep(NA,length(LABEL))
FILL[which(LABEL == "IgH (C)")] <- "#8dd3c7"
FILL[which(LABEL == "IgH (J)")] <- "#bc80bd"
FILL[which(LABEL == "IgH (D)")] <- "#fdb462"
FILL[which(LABEL == "IgH (V)")] <- "#80b1d3"
features3$fill <- FILL 
features3$color <- FILL 
HEI <- rep(0,length(LABEL))
HEI[which(LABEL != "")] <- 0.03
features3$height <- HEI

chr14 <- sample.gr
chr14$label.parameter.rot <- 90
chr14$lwd <- 1
chr14$shape <- "diamond"

##### LPP locus
Anno4 <- Anno[which(Anno$chr == "chr3"),]
Loci4 <- SVbp[which(SVbp$Locus == "LPP"),]

POS <- Loci4$position
NAME <- paste0(Loci4$Target," (",Loci4$patient,")")
COLOR <- Loci4$color

sample.gr <- GRanges("chr3", IRanges(POS, width=1, names=NAME))
sample.gr$color <- COLOR

sample.gr$border <- "gray30"
sample.gr$label <- NULL

START <- Anno4$start
END <- Anno4$end
MIN <- min(START,POS,END) - 5000
MAX <- max(START,POS,END) 
VECT <- c(MIN,START,MAX)
WIDTH <- c(0,END-START,5000)
LABEL <- c("",Anno4$type,"")

features4 <- GRanges("chr3", IRanges(VECT,  width=WIDTH,   names=LABEL ))

FILL <- rep(NA,length(LABEL))
FILL[which(LABEL == "LPP")] <- "#2166ac"
FILL[which(LABEL == "IgH (C)")] <- "#8dd3c7"
FILL[which(LABEL == "IgH (J)")] <- "#bc80bd"
FILL[which(LABEL == "IgH (D)")] <- "#fdb462"
features4$fill <- FILL 
features4$color <- FILL 
HEI <- rep(0,length(LABEL))
HEI[which(LABEL != "")] <- 0.03
features4$height <- HEI

chr3 <- sample.gr
chr3$label.parameter.rot <- 90
chr3$lwd <- 1
chr3$shape <- "diamond"



pdf(file = paste0("./Figure1B.pdf"),width=8,height=6, onefile = T)

lolliplot(chr18 , features,  legend=NULL, txtSize=1,cex=1,rescale = F,lwd=0.1)
lolliplot(chr8  , features2, legend=NULL, txtSize=1,cex=1,rescale = F,lwd=0.1)
lolliplot(chr14 , features3, legend=NULL, txtSize=1,cex=1,rescale = F,lwd=0.1)
lolliplot(chr3  , features4, legend=NULL, txtSize=1,cex=1,rescale = F,lwd=0.1)

dev.off() 











