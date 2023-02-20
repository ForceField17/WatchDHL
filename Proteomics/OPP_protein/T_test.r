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

#get CNV cases

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


#PCA
library(genefilter)
library(ggrepel)

expMatrix <- expMatrix_new
rv <- rowVars(expMatrix)
ntop <- 500
select <- order(rv, decreasing = TRUE)[seq_len(min(ntop,  length(rv)))]
pca <- prcomp(t(expMatrix[select, ]))
percentVar <- pca$sdev^2/sum(pca$sdev^2)
intgroup <- c("case",  "Cell",  "Group","Protein","Line")
if (!all(intgroup %in% colnames(samples))) {
  stop("the argument 'intgroup' should specify columns of colData(dds)")
}
intgroup.df <- as.data.frame(samples[, intgroup, drop = FALSE])
group <- if (length(intgroup) > 1){ factor(apply(intgroup.df, 1, paste, collapse = " : ")) }else{ samples[[intgroup]]}

d <- data.frame(pca$x,intgroup.df)

myPCA <- ggplot()+theme_classic()
myPCA <- myPCA + geom_point(data = d, aes(x =PC1, y = PC2,  fill = Line,shape = Protein),alpha=0.9,size = 3.6) + 
  xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) + scale_shape_manual(values=c(21,22))+
  scale_fill_manual(values = c(OE="#66c2a5",Repair="#fc8d62",WT="#6666F6"),guide = guide_legend(override.aes=list(shape=21)))+
  ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) + coord_fixed() 
myPCA <- myPCA + theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,1,0.5,1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                       text=element_text(size=12,face='plain',color='black'),legend.key.width=unit(1,'cm'),legend.key.height=unit(1,'cm'),legend.position='right',
                       legend.margin=margin(t=0.1,r=0.1,b=0,l=0.1,unit='cm'),legend.text=element_text(size=12,face='plain',color='black'),axis.text.y=element_text(size=12,face='plain',color='black'),
                       axis.text.x=element_text(size=12,face='plain',color='black'),axis.title.x=element_text(size=14,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))
myPCA <- myPCA + scale_y_continuous(expand=c(0,0),limits=c(-14,14),breaks = seq(-20,20,5)) +
                 scale_x_continuous(expand=c(0,0),limits=c(-14,14),breaks = seq(-20,20,5))
figure_1<-rbind(ggplotGrob(myPCA ),size="first")
ggsave(file="./figs/Fig5b_PCA.pdf", plot=figure_1,bg = 'white', width = 14, height = 10, units = 'cm', dpi = 600)
########


#testing
data <- expMatrix_new

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
  
  Rank_OE_WT[i] <- (mean(G2)-mean(G1)) / ( sd(G2)^2/length(G2) + sd(G1)^2/length(G1) )^0.5
  Rank_Repair_WT[i] <- (mean(G3)-mean(G1)) / ( sd(G3)^2/length(G3) + sd(G1)^2/length(G1) )^0.5

}

geneName <-rownames(expMatrix_new)
FDR_OE_WT <- p.adjust(Pvalue_OE_WT,method = "fdr")
FDR_Repair_WT <- p.adjust(Pvalue_Repair_WT,method = "fdr")

final_data <- data.frame(geneName,Pvalue_OE_WT,Pvalue_Repair_WT,FDR_OE_WT,FDR_Repair_WT,FoldC_OE_WT,FoldC_Repair_WT)
write.csv(final_data , file = "./figs/Statistical_T_test.csv")


####################
maker <- read.table("./lib/Drug_target.txt",header=T)
mydata <- t(expMatrix_new[which(rownames(expMatrix_new) %in% maker$Target),])
opp <- samples[which(samples$Protein=="nascent" & samples$Line != "OE"),]
mydata <- mydata[which(rownames(mydata) %in% opp$case),]

#heatmap
ann_colors = list(
  Group = c(nasOE="#7fbc41",nasRepair="#de77ae",nasWT="grey60")
)

mat <- expMatrix_new[which(rownames(expMatrix_new) %in%  maker$Target),which(colnames(expMatrix_new) %in% opp$case)]
anno_col <- as.data.frame(opp[, c("Group")])
colnames(anno_col) <- c("Group")

pdf(file = paste0("./figs/ExtFig15g_DrugTarget_heatmap.pdf"),width=4.2,height=5)
pheatmap(mat,cutree_rows=2,cluster_cols =F ,cluster_rows = T,#gaps_col=c(10),
         color = colorRampPalette(c('#29419E','#698DC9',"grey95",'#F67172','#DC2B18'))(50),
         annotation_col = anno_col,annotation_names_col=T,show_rownames = T,annotation_legend = T,
         scale = "row",annotation_colors = ann_colors,border_color = "grey20",
         clustering_distance_rows = "correlation",legend = TRUE, legend_breaks = NA,legend_labels = NA,
         clustering_distance_cols = "correlation", #"euclidean",#cellwidth = 15, #cellheight = 2,
         fontsize = 12,fontsize_row = 12, fontsize_col = 12, display_numbers = F,
         number_format = "%.0f", number_color = "grey6030", fontsize_number = 2)
dev.off()

Infos <-  final_data[which(final_data$geneName %in% maker$Target),]
write.table(Infos,file = "./figs/Pvalues_DrugTarget.txt",sep="\t",quote = F,row.names = F)
###################################################################


#################
TSG <- read.table("./lib/TS_list.txt",header=F)
ONCO <- read.table("./lib/OC_list.txt",header=F)

################# OE vs WT (OPP+)
plotTable <- final_data
plotTable_up <- plotTable[which( plotTable$FoldC_OE_WT  >=  0.4  & plotTable$FDR_OE_WT <  0.05 ) ,]
plotTable_down <- plotTable[which( plotTable$FoldC_OE_WT  <=  -0.4  & plotTable$FDR_OE_WT <  0.05 ) ,]
plotTable_other <- plotTable[which(!(plotTable$geneName %in% plotTable_up$geneName | plotTable$geneName %in% plotTable_down$geneName)) ,]

plotTable_Onco <- plotTable[(plotTable$geneName %in% ONCO$V1),]
plotTable_TSG <- plotTable[(plotTable$geneName %in% TSG$V1),]

plotTable_sig1 <- plotTable_up[(plotTable_up$geneName %in% plotTable_Onco$geneName ),]
plotTable_sig2 <- plotTable_down[(plotTable_down$geneName %in% plotTable_TSG$geneName),]
plotTable_sig <- rbind(plotTable_sig1,plotTable_sig2)

mySub2 <- ggplot()+theme_classic()
mySub2 <- mySub2 + geom_point(data = plotTable_other, aes(x =  (FoldC_OE_WT), y = -log10(FDR_OE_WT)),fill="grey",alpha=0.5,size=1.5, shape = 21,stroke=0.2) 
mySub2 <- mySub2 + geom_point(data = plotTable_up,    aes(x =  (FoldC_OE_WT), y = -log10(FDR_OE_WT)),fill=gg_color_hue(3)[1]  ,alpha=0.8,size=2, shape = 21,stroke=0.2) 
mySub2 <- mySub2 + geom_point(data = plotTable_down,    aes(x =  (FoldC_OE_WT), y = -log10(FDR_OE_WT)),fill=gg_color_hue(3)[3],alpha=0.8,size=2, shape = 21,stroke=0.2) 
#mySub2 <- mySub2 + geom_point(data = plotTable_sig1,    aes(x =   (FoldC_OE_WT), y = -log10(FDR_OE_WT)),fill=gg_color_hue(3)[1],alpha=0.9,size=2, shape = 21,stroke=0.2) 
#mySub2 <- mySub2 + geom_point(data = plotTable_sig2,    aes(x =   (FoldC_OE_WT), y = -log10(FDR_OE_WT)),fill=gg_color_hue(3)[3],alpha=0.9,size=2, shape = 21,stroke=0.2) 
mySub2 <- mySub2 + geom_vline(xintercept = 0.4,linetype=2,color="grey50") +
  geom_vline(xintercept = -0.4,linetype=2,color="grey50") + 
  geom_hline(yintercept = -log10(0.05),linetype=2,color="grey50") 
mySub2 <- mySub2 +geom_text(data =plotTable_up,   aes(x = (FoldC_OE_WT), y = -log10(FDR_OE_WT),label = geneName), size =0.3,color ='black')
mySub2 <- mySub2 +geom_text(data =plotTable_down,   aes(x = (FoldC_OE_WT), y = -log10(FDR_OE_WT),label = geneName), size =0.3,color ='black')

mySub2 <- mySub2 + xlab("Fold change log2( OE-OPP+ / WT-OPP+)") + ylab("-log10(adjusted p-value)") +
  scale_y_continuous(expand=c(0,0),limits=c(0,2.8),breaks=seq(-6,10,0.5)) + 
  scale_x_continuous(expand=c(0,0),limits=c(-1.77,1.3),breaks=seq(-6,6,0.4))
mySub2 <- mySub2 + theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,1,0.5,1),'lines'),plot.title=element_text(size=14,vjust=0.5,hjust=0.5,face='bold.italic'),
                         legend.key.width=unit(1.5,'cm'),legend.key.height=unit(0.5,'cm'),legend.position='top',
                         legend.margin=margin(t=0.1,r=0.1,b=0,l=0.1,unit='cm'),legend.text=element_text(size=12),axis.text.y=element_text(size=10,face='plain',color='black'),
                         axis.text.x=element_text(size=10,face='plain',color='black'),axis.title.x=element_text(size=12,face='plain',color='black'),axis.title.y=element_text(size=12,hjust=0.5,vjust=2,face='plain',color='black'))
figure_4<-rbind(ggplotGrob(mySub2),size="last")
ggsave(file="./figs/Fig5c_volcano_OE_vs_WT_OPP.pdf", plot=figure_4,bg = 'white', width = 12, height = 8, units = 'cm', dpi = 600)

a <- data.frame( c( nrow(plotTable_sig2), nrow(plotTable_TSG[which(!(plotTable_TSG$geneName %in% plotTable_down$geneName)),]) ),
                 c( nrow(plotTable_down[which(!(plotTable_down$geneName %in% plotTable_TSG$geneName)),]), nrow(plotTable[which(!(plotTable$geneName %in% plotTable_down$geneName | plotTable$geneName %in% plotTable_TSG$geneName)),])  ))
colnames(a) <- c("TSG","nonTSG")
rownames(a) <- c("Down","notDown")
a
fisher.test(a,alternative = "two.sided")

b <- data.frame( c( nrow(plotTable_sig1), nrow(plotTable_Onco[which(!(plotTable_Onco$geneName %in% plotTable_up$geneName)),]) ),
                 c( nrow(plotTable_up[which(!(plotTable_up$geneName %in% plotTable_Onco$geneName)),]), nrow(plotTable[which(!(plotTable$geneName %in% plotTable_up$geneName | plotTable$geneName %in% plotTable_Onco$geneName)),])  ))
colnames(b) <- c("ONCO","nonONCO")
rownames(b) <- c("Up","notUp")
b
fisher.test(b,alternative = "two.sided")



###################### PAX5-Mut vs WT (OPP+)
plotTable <- final_data
plotTable_up <- plotTable[which( plotTable$FoldC_Repair_WT  >=  0.6  & plotTable$FDR_Repair_WT <  0.05 ) ,]
plotTable_down <- plotTable[which( plotTable$FoldC_Repair_WT  <=  -0.6  & plotTable$FDR_Repair_WT <  0.05 ) ,]
plotTable_other <- plotTable[which(!(plotTable$geneName %in% plotTable_up$geneName | plotTable$geneName %in% plotTable_down$geneName)) ,]

plotTable_Onco <- plotTable[(plotTable$geneName %in% ONCO$V1),]
plotTable_TSG <- plotTable[(plotTable$geneName %in% TSG$V1),]

plotTable_sig1 <- plotTable_up[(plotTable_up$geneName %in% plotTable_Onco$geneName ),]
plotTable_sig2 <- plotTable_down[(plotTable_down$geneName %in% plotTable_TSG$geneName),]
plotTable_sig <- rbind(plotTable_sig1,plotTable_sig2)

mySub2 <- ggplot()+theme_classic()
mySub2 <- mySub2 + geom_point(data = plotTable_other, aes(x =  (FoldC_Repair_WT), y = -log10(FDR_Repair_WT)),fill="grey",alpha=0.5,size=1.5, shape = 21,stroke=0.2) 
mySub2 <- mySub2 + geom_point(data = plotTable_up,    aes(x =  (FoldC_Repair_WT), y = -log10(FDR_Repair_WT)),fill=gg_color_hue(3)[1]  ,alpha=0.8,size=2, shape = 21,stroke=0.2) 
mySub2 <- mySub2 + geom_point(data = plotTable_down,    aes(x =  (FoldC_Repair_WT), y = -log10(FDR_Repair_WT)),fill=gg_color_hue(3)[3],alpha=0.8,size=2, shape = 21,stroke=0.2) 
#mySub2 <- mySub2 + geom_point(data = plotTable_sig1,    aes(x =   (FoldC_Repair_WT), y = -log10(FDR_Repair_WT)),fill=gg_color_hue(3)[1],alpha=0.9,size=2, shape = 21,stroke=0.2) 
#mySub2 <- mySub2 + geom_point(data = plotTable_sig2,    aes(x =   (FoldC_Repair_WT), y = -log10(FDR_Repair_WT)),fill=gg_color_hue(3)[3],alpha=0.9,size=2, shape = 21,stroke=0.2) 
mySub2 <- mySub2 + geom_vline(xintercept = 0.6,linetype=2,color="grey50") +
                   geom_vline(xintercept = -0.6,linetype=2,color="grey50") + 
                   geom_hline(yintercept = -log10(0.05),linetype=2,color="grey50") 

mySub2 <- mySub2 +geom_text(data =plotTable_down,   aes(x = (FoldC_Repair_WT), y = -log10(FDR_Repair_WT),label = geneName), size =0.3,color ='black')
mySub2 <- mySub2 +geom_text(data =plotTable_up,   aes(x = (FoldC_Repair_WT), y = -log10(FDR_Repair_WT),label = geneName), size =0.3,color ='black')

mySub2 <- mySub2 + xlab("Fold change log2( Mut-OPP+ / WT-OPP+)") + ylab("-log10(adjusted p-value)") +
  scale_y_continuous(expand=c(0,0),limits=c(0,4),breaks=seq(-6,10,0.5)) + 
  scale_x_continuous(expand=c(0,0),limits=c(-2.8,1.85),breaks=seq(-6,6,0.6))
mySub2 <- mySub2 + theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,1,0.5,1),'lines'),plot.title=element_text(size=14,vjust=0.5,hjust=0.5,face='bold.italic'),
                         legend.key.width=unit(1.5,'cm'),legend.key.height=unit(0.5,'cm'),legend.position='top',
                         legend.margin=margin(t=0.1,r=0.1,b=0,l=0.1,unit='cm'),legend.text=element_text(size=12),axis.text.y=element_text(size=10,face='plain',color='black'),
                         axis.text.x=element_text(size=10,face='plain',color='black'),axis.title.x=element_text(size=12,face='plain',color='black'),axis.title.y=element_text(size=12,hjust=0.5,vjust=2,face='plain',color='black'))
figure_4<-rbind(ggplotGrob(mySub2),size="last")
ggsave(file="./figs/Fig5d_volcano_nascent_Repair_WT_up.pdf", plot=figure_4,bg = 'white', width = 12, height = 8, units = 'cm', dpi = 600)

a <- data.frame( c( nrow(plotTable_sig2), nrow(plotTable_TSG[which(!(plotTable_TSG$geneName %in% plotTable_down$geneName)),]) ),
                 c( nrow(plotTable_down[which(!(plotTable_down$geneName %in% plotTable_TSG$geneName)),]), nrow(plotTable[which(!(plotTable$geneName %in% plotTable_down$geneName | plotTable$geneName %in% plotTable_TSG$geneName)),])  ))
colnames(a) <- c("TSG","nonTSG")
rownames(a) <- c("Down","notDown")
a
fisher.test(a,alternative = "two.sided")

b <- data.frame( c( nrow(plotTable_sig1), nrow(plotTable_Onco[which(!(plotTable_Onco$geneName %in% plotTable_up$geneName)),]) ),
                 c( nrow(plotTable_up[which(!(plotTable_up$geneName %in% plotTable_Onco$geneName)),]), nrow(plotTable[which(!(plotTable$geneName %in% plotTable_up$geneName | plotTable$geneName %in% plotTable_Onco$geneName)),])  ))
colnames(b) <- c("ONCO","nonONCO")
rownames(b) <- c("Up","notUp")
b
fisher.test(b,alternative = "two.sided")

