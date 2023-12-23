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
library(rstatix)
library(ggpubr)
library(ggrepel)
library(preprocessCore)
library(clusterProfiler)
library(enrichplot)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 150)[1:n]
}

scale_rows = function(x){
  m = apply(x, 1, mean, na.rm = T)
  s = apply(x, 1, sd, na.rm = T)
  return((x - m) / s)
}

data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}


#raw data preprocessing
Sample_features <- read.delim("Data/Info_Nascent.txt",sep = "\t",header = T)
samples <- data.frame(Sample_features)
rownames(samples) <- Sample_features$case
################
counts1 <- read.csv("Data/Final_Nascent_protein_matrix.txt", sep="", head=T, row.names = "Gene",quote = "\"")
counts1 <- counts1[,which((colnames(counts1 ) %in% samples$case ))]

nrow(counts)
counts <- round(counts1/10000)
xxx <- counts[,c(4:6,10:12,16:18)]
keep <- rowSums(xxx == min(counts)) <= 4
counts <- counts[keep,]
nrow(counts)

dds=DESeqDataSetFromMatrix(countData = counts,colData = samples,design = ~ Group)

vsd <- vst(dds ,blind = TRUE)
head(assay(vsd), 3)
expMatrix_new <- assay(vsd)

all(rownames(samples) %in% colnames(expMatrix_new ))
all(colnames(expMatrix_new ) %in% rownames(samples))



#testing
data <- expMatrix_new
raw_matrix <- assay(dds)

Pvalue_OE_WT     <- rep(NA,nrow(data))
Pvalue_Repair_WT <- rep(NA,nrow(data))
FoldC_OE_WT      <- rep(NA,nrow(data))
FoldC_Repair_WT  <- rep(NA,nrow(data))
Rank_OE_WT      <- rep(NA,nrow(data))
Rank_Repair_WT  <- rep(NA,nrow(data))
Rank2_OE_WT      <- rep(NA,nrow(data))
Rank2_Repair_WT  <- rep(NA,nrow(data))
nVar1 <- 0
nVar2 <- 0
for(i in 1:nrow(data)){
  test <- data.frame(as.numeric(data[i,]), samples$Group,samples$Cell)
  colnames(test) <- c("value","Group","Cell")
  rownames(test) <- colnames(data)
  G1 <- test[which(test$Group=="nasWT"),1]
  G2 <- test[which(test$Group=="nasOE"),1]
  G3 <- test[which(test$Group=="nasRepair"),1]
  
  Var_OE_WT <- var.test(G2,G1)$p.value
  Var_Repair_WT <- var.test(G3,G1)$p.value
  
  if(Var_OE_WT > 0.05){
      Pvalue_OE_WT[i]  <- t.test(G2,G1,paired = F,var.equal = T)$p.value
      nVar1 <- nVar1 + 1
  }else{
      Pvalue_OE_WT[i]   <- t.test(G2,G1,paired = F,var.equal = F)$p.value
  }
  
  if(Var_Repair_WT > 0.05){
      Pvalue_Repair_WT[i] <- t.test(G3,G1,paired = F,var.equal = T)$p.value
      nVar2 <- nVar2 + 1
  }else{
      Pvalue_Repair_WT[i] <- t.test(G3,G1,paired = F,var.equal = F)$p.value
  }
  
  FoldC_OE_WT[i]      <-  mean(G2)-mean(G1)
  FoldC_Repair_WT[i]  <-  mean(G3)-mean(G1)

}

geneName <-rownames(data)

FDR_OE_WT <- p.adjust(Pvalue_OE_WT,method = "fdr")
FDR_Repair_WT <- p.adjust(Pvalue_Repair_WT,method = "fdr")

final_data <- data.frame(geneName,Pvalue_OE_WT,Pvalue_Repair_WT,FDR_OE_WT,FDR_Repair_WT,FoldC_OE_WT,FoldC_Repair_WT)



#########
TSG <- read.table("./lib/TS_list.txt",header=F)
ONCO <- read.table("./lib/OC_list.txt",header=F)
nascent <- read.table("figs/Statistical_two_way_ANOVA.txt",header = T,sep="\t")

plotTable <-  final_data[which(final_data$geneName %in% nascent$geneName[which(nascent$Rank_OPP>0 & nascent$FDR_OPP < 0.05)] ),]
plotTable_nonNas <- final_data[which(!(final_data$geneName %in% plotTable$geneName)),]

plotTable_change <- plotTable[( abs(plotTable$FoldC_OE_WT) > 0.5 | abs(plotTable$FoldC_Repair_WT) > 0.5) ,]

plotTable_Onco <- plotTable[(plotTable$geneName %in% ONCO$V1),]
plotTable_TSG <- plotTable[(plotTable$geneName %in% TSG$V1),]
plotTable_other <- plotTable[!(plotTable$geneName %in% ONCO$V1 | plotTable$geneName %in% TSG$V1),]


plotTable_sig1 <- plotTable_change[(plotTable_change$geneName %in% plotTable_Onco$geneName),]
plotTable_sig2 <- plotTable_change[(plotTable_change$geneName %in% plotTable_TSG$geneName),]
plotTable_sig <- rbind(plotTable_sig1,plotTable_sig2)

mySub2 <- ggplot()+theme_classic()
#mySub2 <- mySub2 + geom_point(data = plotTable_nonNas,   aes(x = FoldC_OE_WT, y = FoldC_Repair_WT), fill="black",alpha=0.4,size=0.5, shape = 22,stroke=0.2) 
mySub2 <- mySub2 + geom_point(data = plotTable_other ,   aes(x = FoldC_OE_WT, y = FoldC_Repair_WT), fill="grey",alpha=0.4,size=2, shape = 21,stroke=0.2) 
mySub2 <- mySub2 + geom_vline(xintercept = 0,linetype=2,color="grey40") + geom_hline(yintercept = 0,linetype=2,color="grey40") + geom_abline(slope = 1,linetype=2,color="grey40")
mySub2 <- mySub2 + geom_point(data = plotTable_Onco,      aes(x = FoldC_OE_WT, y = FoldC_Repair_WT), fill=gg_color_hue(3)[1],alpha=0.5,size=2, shape = 21,stroke=0.2) 
mySub2 <- mySub2 + geom_point(data = plotTable_TSG,    aes(x = FoldC_OE_WT, y = FoldC_Repair_WT), fill=gg_color_hue(3)[3],alpha=0.5,size=2, shape = 21,stroke=0.2) 

mySub2 <- mySub2 +geom_text_repel(data =plotTable,   aes(x = (FoldC_OE_WT), y = FoldC_Repair_WT,label = geneName), size =1,color ='black')

mySub2 <- mySub2 + xlab("Fold change (nascent OE / nascent WT)") + ylab("Fold change (nascent Mut / nascent WT)") +
  scale_y_continuous(expand=c(0,0),limits=c(-1.5,1.5),breaks=seq(-6,6,0.5)) + 
  scale_x_continuous(expand=c(0,0),limits=c(-1.5,1.5),breaks=seq(-6,6,0.5)) 
mySub2 <- mySub2 + theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,1,0.5,1),'lines'),
                         legend.key.width=unit(0.3,'cm'),legend.key.height=unit(1,'cm'),legend.position='right',legend.text=element_text(size=12),
                         axis.text.y=element_text(size=10,face='plain',color='black'),axis.text.x=element_text(size=10,face='plain',color='black'),
                         axis.title.x=element_text(size=12,face='plain',color='black'),axis.title.y=element_text(size=12,hjust=0.5,vjust=2,face='plain',color='black'))
figure_4<-rbind(ggplotGrob(mySub2),size="last")
ggsave(file="./figs/ExtFig10c_Highlighting_nascent_proteins.pdf", plot=figure_4,bg = 'white', width = 12, height = 11, units = 'cm', dpi = 600)



cor.test(~ FoldC_OE_WT + FoldC_Repair_WT,data=plotTable_TSG,method="spearman")
cor.test(~ FoldC_OE_WT + FoldC_Repair_WT,data=plotTable_Onco,method="spearman")
#######################
cor.test(~ FoldC_OE_WT + FoldC_Repair_WT,data=plotTable,method="pearson",exact =T, conf.level = 0.95)
cor.test(~ FoldC_OE_WT + FoldC_Repair_WT,data=plotTable,method="spearman",exact =T, conf.level = 0.95)
