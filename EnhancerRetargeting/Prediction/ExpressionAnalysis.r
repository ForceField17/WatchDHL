# WatchDHL
library(rstudioapi)

# Getting the path of your current open file
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))
print( getwd() )
library(ggplot2)
library(DESeq2)
library(stringr)
library("pheatmap")
library("RColorBrewer")
library(fgsea)
library(ComplexHeatmap)
library(circlize)
library("SC3")
library("ggrepel")
library("preprocessCore")

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 150)[1:n]
}


Sample_features <- read.delim("./RNAseq/Info_11CellLines_discovery.txt",sep = "\t",header = T)
samples <- data.frame(Sample_features)
rownames(samples) <- Sample_features$case

counts1 <- read.csv("./RNAseq/AllCellLineRNA.txt", sep="", head=T, row.names = "CellLine",quote = "\"")
nrow(counts1)
counts <- na.omit(counts1)
counts <- counts[,which(colnames(counts) %in% rownames(samples))]
Name <- as.character(colnames(counts))
x.term <- as.character(trimws(Name, which = c("both", "left", "right"), whitespace = "[ \t\r\n]"))



####################
dds=DESeqDataSetFromMatrix(countData = counts,colData = samples,design = ~ Histology)
nrow(dds)
keep <- rowSums(counts(dds) >= 50) >= 3
dds <- dds[keep,]
nrow(dds)
vsd <- vst(dds, blind = FALSE)
head(assay(vsd), 3)
expMatrix <- assay(vsd)
expMatrix_new <- normalize.quantiles(expMatrix, copy = TRUE)
colnames(expMatrix_new) <- colnames(expMatrix)
rownames(expMatrix_new) <- rownames(expMatrix)
expMatrix <- expMatrix_new


####
genePairs <- read.table("final.TAD_gene_pairs.txt",sep = "\t",header = F)
colnames(genePairs) <- c("V1","V2","direction","distance")

data <- expMatrix
selected <- which(samples$tissue=="CellLine")

genePairs <- genePairs[which(genePairs$V1 %in% rownames(data) & genePairs$V2 %in% rownames(data)),]

pValue <- rep(1,nrow(genePairs))
pValue_LR <- rep(1,nrow(genePairs))
SpearmanR <- rep(1,nrow(genePairs))
MeanExpA <- rep(1,nrow(genePairs))
MeanExpB <- rep(1,nrow(genePairs))
for(i in 1:nrow(genePairs)){
  geneA <- data[which(rownames(data) == genePairs[i,1]),]
  geneB <- data[which(rownames(data) == genePairs[i,2]),]
  
    Sample <- samples
    Sample$Exp_geneA <- geneA
    Sample$Exp_geneB <- geneB
    Sample$Zscore_geneA <- (Sample$Exp_geneA-mean(Sample$Exp_geneA))/sd(Sample$Exp_geneA)
    Sample$Zscore_geneB <- (Sample$Exp_geneB-mean(Sample$Exp_geneB))/sd(Sample$Exp_geneB)
    
    a<-cor.test(Sample$Zscore_geneA,Sample$Zscore_geneB,method = "spearman")
    a$p.value
    a$estimate
    
    if( genePairs[i,1] =="PAX5" & genePairs[i,2] =="ZCCHC7"  ){
    new.plot<-ggplot() +theme_classic()
    new.plot <- new.plot+   geom_smooth(data = Sample, aes(x = Zscore_geneA, y = Zscore_geneB),method = "lm", se = F,formula= y ~ x,color="grey")
    new.plot <- new.plot+   geom_point(data = Sample, aes(x = Zscore_geneA, y = Zscore_geneB),color="#8dd3c7",alpha=0.9,size=3) +
      geom_text_repel(data =Sample ,aes(x=Zscore_geneA,y=Zscore_geneB,label=case),
                      stat = "identity",position = "identity",max.overlaps = getOption("ggrepel.max.overlaps", default = 30),force=1.5,
                      segment.size =0.2, size =2.4,color ='black')+
      geom_text(aes(x=-1.8,y=-2,label=paste0("pvalue = ",round(a$p.value,3))),color="blue",size=3)+
      geom_text(aes(x=-1.8,y=-2.23,label=paste0("RHO = ",round(a$estimate,3))),color="blue",size=3)
    new.plot<- new.plot +  theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(1,1,1,1),'lines'),legend.title=element_text(size=12,face='plain',color='black'),
                                 plot.title=element_text(size=14,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=12,face='bold'),
                                 legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position="right",legend.text=element_text(size=10,hjust=0,face='plain'),
                                 axis.text.x=element_text(size=12,angle=30,vjust=1,hjust=1,face='bold',color='black'),axis.text.y=element_text(size=12,face='bold',color='black'),
                                 axis.title.x=element_text(size=12,vjust=0,hjust=0.5,face='plain',color='black'),axis.title.y=element_text(size=12,face='plain',color='black'))
    new.plot<- new.plot+ylab(paste0("Expression of ",genePairs[i,2]," (Z-score)"))+xlab(paste0("Expression of ",genePairs[i,1]," (Z-score)"))
    
    new.plot<- new.plot+scale_y_continuous(expand=c(0,0),limits=c(-2.5,2.5),breaks = seq(-2,2,1)) 
    new.plot<- new.plot+scale_x_continuous(expand=c(0,0),limits=c(-2.5,2.5),breaks = seq(-2,2,1)) 
  
    new.plot<- new.plot+ggtitle(paste0(genePairs[i,1]," --- ",genePairs[i,2]))
    plot7<-cbind(ggplotGrob(new.plot),size="first")
    ggsave(file=paste0("./summary/ExtFig4a_linear_",genePairs[i,1],"_",genePairs[i,2],".pdf"), plot=plot7,bg = 'white', width = 10, height = 10, units = 'cm', dpi = 600)
    }
    
    pValue_LR[i] <- a$p.value
    SpearmanR[i] <- a$estimate
    MeanExpA[i] <-  mean(geneA)
    MeanExpB[i] <-  mean(geneB)
}

resultsAmp <-data.frame(genePairs$V1,genePairs$V2,MeanExpA,MeanExpB,SpearmanR,pValue_LR,genePairs$direction)
write.table(resultsAmp , file = "results_table.test.txt",sep = "\t",quote=F,row.names = F)


