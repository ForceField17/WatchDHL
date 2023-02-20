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

SVbp <- read.table("breakpoint_position.txt",sep = "\t",header = T)
SVbp$color <- gg_color_hue(6)[3]
SVbp$color[which(SVbp$Tumor == "DHL")] <- gg_color_hue(6)[1]
SVbp <- SVbp[which(SVbp$patient != "P3"),]


##### BCL2 locus
Loci1 <- SVbp[which(SVbp$Locus == "BCL2"),]
Anno1 <- Anno[which(Anno$locus == "BCL2"),]

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
LABEL <- c("",Anno1$locus,"")
PART <- c("",Anno1$part,"")

features <- GRanges("chr18", IRanges(VECT,  width=WIDTH,   names=LABEL ))

FILL <- rep(NA,length(LABEL))
FILL[which(LABEL != "")] <- "#2166ac"
features$fill <- FILL 
features$color <- FILL 

HEI <- rep(0,length(PART))
HEI[which(PART == "CDS")] <- 0.04
HEI[which(PART == "UTR")] <- 0.016
features$height <- HEI


chr18 <- sample.gr
chr18$label.parameter.rot <- 90
chr18$lwd <- 1
chr18$shape <- "diamond"


##### MYC locus
Loci2 <- SVbp[which(SVbp$Locus == "MYC"),]
Anno2 <- Anno[which(Anno$locus == "MYC" | Anno$locus == "PVT1"),]


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
LABEL <- c("",Anno2$locus,"")

features2 <- GRanges("chr8", IRanges(VECT,  width=WIDTH,   names=LABEL ))

FILL <- rep(NA,length(LABEL))
FILL[which(LABEL == "MYC")] <- "#2166ac"
FILL[which(LABEL == "PVT1")] <- "#de77ae"
features2$fill <- FILL 
features2$color <- FILL 

PART <- c("",Anno2$part,"")
HEI <- rep(0,length(PART))
HEI[which(PART == "CDS")] <- 0.04
HEI[which(PART == "UTR")] <- 0.016
HEI[which(PART == "transcript")] <- 0.01
features2$height <- HEI

chr8 <- sample.gr
chr8$label.parameter.rot <- 90
chr8$lwd <- 1
chr8$shape <- "diamond"

##### IgH locus
Loci3 <- SVbp[which(SVbp$Locus == "IgH"),]
Anno3 <- Anno[which(Anno$chr == "chr14"),]

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
LABEL <- c("",Anno3$locus,"")

features3 <- GRanges("chr14", IRanges(VECT,  width=WIDTH,   names=LABEL ))

FILL <- rep(NA,length(LABEL))
FILL[which(LABEL == "IgH (C)")] <- "#8dd3c7"
FILL[which(LABEL == "IgH (J)")] <- "#bc80bd"
FILL[which(LABEL == "IgH (D)")] <- "#fdb462"
FILL[which(LABEL == "IgH (V)")] <- "#d73027"
features3$fill <- FILL 
features3$color <- FILL 
HEI <- rep(0,length(LABEL))
HEI[which(LABEL != "")] <- 0.04
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
LABEL <- c("",Anno4$locus,"")

features4 <- GRanges("chr3", IRanges(VECT,  width=WIDTH,   names=LABEL ))

FILL <- rep(NA,length(LABEL))
FILL[which(LABEL == "LPP")] <- "#2166ac"
FILL[which(LABEL == "IgH (C)")] <- "#8dd3c7"
FILL[which(LABEL == "IgH (J)")] <- "#bc80bd"
FILL[which(LABEL == "IgH (D)")] <- "#fdb462"
features4$fill <- FILL 
features4$color <- FILL 
PART <- c("",Anno4$part,"")
HEI <- rep(0,length(PART))
HEI[which(PART == "CDS")] <- 0.04
HEI[which(PART == "UTR")] <- 0.016
features4$height <- HEI

chr3 <- sample.gr
chr3$label.parameter.rot <- 90
chr3$lwd <- 1
chr3$shape <- "diamond"

##### IgL locus
Loci5 <- SVbp[which(SVbp$Locus == "IgL"),]
Anno5 <- Anno[which(Anno$chr == "chr22"),]

POS <- Loci5$position
NAME <- paste0(Loci5$Target," (",Loci5$patient,")")
COLOR <- Loci5$color

sample.gr <- GRanges("chr22", IRanges(POS, width=1, names=NAME))
sample.gr$color <- COLOR

sample.gr$border <- "gray30"
sample.gr$label <- NULL

START <- Anno5$start
END <- Anno5$end
MIN <- min(START,POS,END) - 5000
MAX <- max(START,POS,END) 
VECT <- c(MIN,START,MAX)
WIDTH <- c(0,END-START,5000)
LABEL <- c("",Anno5$locus,"")

features5 <- GRanges("chr22", IRanges(VECT,  width=WIDTH,   names=LABEL ))

FILL <- rep(NA,length(LABEL))
FILL[which(LABEL == "IgL (C)")] <- "#8dd3c7"
FILL[which(LABEL == "IgL (J)")] <- "#bc80bd"
FILL[which(LABEL == "IgL (D)")] <- "#fdb462"
FILL[which(LABEL == "IgL (V)")] <- "#d73027"
features5$fill <- FILL 
features5$color <- FILL 
HEI <- rep(0,length(LABEL))
HEI[which(LABEL != "")] <- 0.04
features5$height <- HEI

chr22 <- sample.gr
chr22$label.parameter.rot <- 90
chr22$lwd <- 1
chr22$shape <- "diamond"


##### IRAG2 locus
Loci6 <- SVbp[which(SVbp$Locus == "IRAG2"),]
Anno6 <- Anno[which(Anno$locus == "IRAG2" ),]


POS <- Loci6$position
NAME <- paste0(Loci6$Target," (",Loci6$patient,")")
COLOR <- Loci6$color

sample.gr <- GRanges("chr12", IRanges(POS, width=1, names=NAME))
sample.gr$color <- COLOR

sample.gr$border <- "gray30"
sample.gr$label <- NULL

START <- Anno6$start
END <- Anno6$end
MIN <- min(START,POS,END) - 5000
MAX <- max(START,POS,END) 
VECT <- c(MIN,START,MAX)
WIDTH <- c(0,END-START,5000)
LABEL <- c("",Anno6$locus,"")

features6 <- GRanges("chr12", IRanges(VECT,  width=WIDTH,   names=LABEL ))

FILL <- rep(NA,length(LABEL))
FILL[which(LABEL == "IRAG2")] <- "#2166ac"
features6$fill <- FILL 
features6$color <- FILL 

PART <- c("",Anno6$part,"")
HEI <- rep(0,length(PART))
HEI[which(PART == "CDS")] <- 0.04
HEI[which(PART == "UTR")] <- 0.016
HEI[which(PART == "transcript")] <- 0.01
features6$height <- HEI

chr12 <- sample.gr
chr12$label.parameter.rot <- 90
chr12$lwd <- 1
chr12$shape <- "diamond"

##### CYRIB locus
Loci7 <- SVbp[which(SVbp$Locus == "CYRIB"),]
Anno7 <- Anno[which(Anno$locus == "CYRIB" ),]


POS <- Loci7$position
NAME <- paste0(Loci7$Target," (",Loci7$patient,")")
COLOR <- Loci7$color

sample.gr <- GRanges("chr8", IRanges(POS, width=1, names=NAME))
sample.gr$color <- COLOR

sample.gr$border <- "gray30"
sample.gr$label <- NULL

START <- Anno7$start
END <- Anno7$end
MIN <- min(START,POS,END) - 5000
MAX <- max(START,POS,END) 
VECT <- c(MIN,START,MAX)
WIDTH <- c(0,END-START,5000)
LABEL <- c("",Anno7$locus,"")

features7 <- GRanges("chr8", IRanges(VECT,  width=WIDTH,   names=LABEL ))

FILL <- rep(NA,length(LABEL))
FILL[which(LABEL == "CYRIB")] <- "#2166ac"
features7$fill <- FILL 
features7$color <- FILL 

PART <- c("",Anno7$part,"")
HEI <- rep(0,length(PART))
HEI[which(PART == "CDS")] <- 0.04
HEI[which(PART == "UTR")] <- 0.016
HEI[which(PART == "transcript")] <- 0.01
features7$height <- HEI

chr8.CYRIB <- sample.gr
chr8.CYRIB$label.parameter.rot <- 90
chr8.CYRIB$lwd <- 1
chr8.CYRIB$shape <- "diamond"

pdf(file = paste0("./All_MYC_BCL2_SV_related_breakpoints.pdf"),width=18,height=5, onefile = T)

lolliplot(chr18 ,      features,  legend=NULL, txtSize=1,cex=1,rescale = F,lwd=0.1)
lolliplot(chr8  ,      features2, legend=NULL, txtSize=1,cex=1,rescale = F,lwd=0.1)
lolliplot(chr14 ,      features3, legend=NULL, txtSize=1,cex=1,rescale = F,lwd=0.1)
lolliplot(chr3  ,      features4, legend=NULL, txtSize=1,cex=1,rescale = F,lwd=0.1)
lolliplot(chr22 ,      features5, legend=NULL, txtSize=1,cex=1,rescale = F,lwd=0.1)
lolliplot(chr12 ,      features6, legend=NULL, txtSize=1,cex=1,rescale = F,lwd=0.1)
lolliplot(chr8.CYRIB , features7, legend=NULL, txtSize=1,cex=1,rescale = F,lwd=0.1)

dev.off() 





