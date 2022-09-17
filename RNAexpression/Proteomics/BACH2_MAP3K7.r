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

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 150)[1:n]
}


counts1 <- read.csv("./data/all_RNAexp_matrix_merged.2022April26.txt", sep="", head=T, row.names = "CellLine",quote = "\"")
nrow(counts1)


counts <- na.omit(counts1)

Name <- as.character(colnames(counts))
x.term <- as.character(trimws(Name, which = c("both", "left", "right"), whitespace = "[ \t\r\n]"))
Sample_features <- read.delim("./data/Info26.txt",sep = "\t",header = T)

samples <- data.frame(Sample_features)
rownames(samples) <- Sample_features$case

all(rownames(samples) %in% colnames(counts))


selected <- which(samples$Histology!="xxxxxxx")
#selected <- which(samples$tissue!="xxx" & samples$tissue!="CellLine")
selected <- which(samples$Site=="DLBCL" & samples$case!="RG126" & samples$tissue=="CellLine" & samples$Histology!="U")# & samples$Class=="Primary")
#selected <- which(samples$Histology=="GCB")
#selected <- which(samples$Histology!="xxx")
samples <- samples[selected,]
counts <- counts[,selected]
all(rownames(samples) %in% colnames(counts))


## create DESeq2 object
dds=DESeqDataSetFromMatrix(countData = counts,colData = samples,design = ~ BCL2)
nrow(dds)

#keep <- rowSums(counts(dds)) > 10
keep <- rowSums(counts(dds) >= 100) >= 3
dds <- dds[keep,]
nrow(dds)


#for n > 30
vsd <- vst(dds, blind = FALSE)
head(assay(vsd), 3)

#check
expMatrix <- assay(vsd)

library(preprocessCore)

expMatrix_new <- normalize.quantiles(expMatrix, copy = TRUE)
colnames(expMatrix_new) <- colnames(expMatrix)
rownames(expMatrix_new) <- rownames(expMatrix)

aaa <- colMeans(expMatrix)
bbb <- colMeans(expMatrix_new)
xxx_mean <- data.frame(aaa,bbb,samples$Type,samples$Gender)
distri_mean <- ggplot()+theme_classic()
distri_mean <- distri_mean + geom_point(data=xxx_mean,aes(x=aaa,y=bbb,color=samples.Type),alpha=0.7,size=3)+theme_bw()
distri_mean 

figure_1<-rbind(ggplotGrob(distri_mean),size="last")

ggsave(file="./temp/qn.pdf", plot=figure_1,bg = 'white', width = 24, height = 16, units = 'cm', dpi = 600)


####PCA
expMatrix <- expMatrix_new
rv <- rowVars(expMatrix)
ntop <- 500
select <- order(rv, decreasing = TRUE)[seq_len(min(ntop,  length(rv)))]
pca <- prcomp(t(expMatrix[select, ]))
percentVar <- pca$sdev^2/sum(pca$sdev^2)
intgroup <- c("Type","tissue", "Gender","Histology")
if (!all(intgroup %in% colnames(samples))) {
  stop("the argument 'intgroup' should specify columns of colData(dds)")
}
intgroup.df <- as.data.frame(samples[, intgroup, drop = FALSE])
group <- if (length(intgroup) > 1){ factor(apply(intgroup.df, 1, paste, collapse = " : ")) }else{ samples[[intgroup]]}

d <- data.frame(pca$x,intgroup.df)
highlight <- d[which(d$tissue!="Tumor"),]
myPCA <- ggplot()+theme_classic()
myPCA <- myPCA + geom_point(data = d, aes(x =PC1, y = PC2, shape = Type, fill = Histology),alpha=0.6,size = 3.6) + 
 # geom_text(data=highlight,aes(x=PC1,y=PC2,label = as.factor(rownames(highlight))),color="black",size=2.5 )+
  xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) + scale_shape_manual(values=c(21,22,23,24,25))+
  scale_fill_manual(values = gg_color_hue(6),guide = guide_legend(override.aes=list(shape=21)))+
  ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) + coord_fixed() 

myPCA <- myPCA + theme(panel.background=element_rect(fill='transparent',color='black'),plot.margin=unit(c(2,1,0.5,1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                       text=element_text(size=14,face='bold'),legend.key.width=unit(1,'cm'),legend.key.height=unit(1,'cm'),legend.position='right',
                       legend.margin=margin(t=0.1,r=0.1,b=0,l=0.1,unit='cm'),legend.text=element_text(size=14,face='bold.italic'),axis.text.y=element_text(size=14,face='bold',color='black'),
                       axis.text.x=element_text(size=14,face='bold',color='black'),axis.title.x=element_text(size=16,face='plain',color='black'),axis.title.y=element_text(size=16,hjust=0.5,vjust=2,face='plain',color='black'))
myPCA

myPCA <- ggplot()+theme_classic()
myPCA <- myPCA + geom_point(data = d, aes(x =PC1, y = PC2, fill = Histology, shape = Type),alpha=0.6,size = 3.6) + 
   geom_text(data=highlight,aes(x=PC1,y=PC2,label = as.factor(rownames(highlight))),color="black",size=2.5 )+
  xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) + scale_shape_manual(values=c(21,22,23,24,25))+
  scale_fill_manual(values = gg_color_hue(6),guide = guide_legend(override.aes=list(shape=21)))+
  ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) + coord_fixed() 
myPCA <- myPCA + theme(panel.background=element_rect(fill='transparent',color='black'),plot.margin=unit(c(2,1,0.5,1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                       text=element_text(size=14,face='bold'),legend.key.width=unit(1,'cm'),legend.key.height=unit(1,'cm'),legend.position='right',
                       legend.margin=margin(t=0.1,r=0.1,b=0,l=0.1,unit='cm'),legend.text=element_text(size=14,face='bold.italic'),axis.text.y=element_text(size=14,face='bold',color='black'),
                       axis.text.x=element_text(size=14,face='bold',color='black'),axis.title.x=element_text(size=16,face='plain',color='black'),axis.title.y=element_text(size=16,hjust=0.5,vjust=2,face='plain',color='black'))
myPCA
figure_1<-rbind(ggplotGrob(myPCA ),size="first")
ggsave(file="./temp/PCA.pdf", plot=figure_1,bg = 'white', width = 24, height = 22, units = 'cm', dpi = 600)
########




###########################heatmaps

gene <- c("ZCCHC7",	"AICDA")


pair_exp <- data.frame(t(expMatrix[which(rownames(expMatrix)==gene[1] | rownames(expMatrix)==gene[2]),]))

Sample <- samples
Sample$Exp_BACH2 <-  pair_exp[,colnames(pair_exp) == gene[1]] 
Sample$Exp_MAP3K7 <- pair_exp[,colnames(pair_exp) == gene[2]] 

#selected <- which(Sample$tissue!="xxx" & Sample$tissue=="CellLine" & Sample$case!="A3KAW" & Sample$case!="A4FUK" & Sample$case!="RG126" &  Sample$case!="U937")
#Sample <- Sample[selected,]

Sample$Zscore_MAP3K7 <- (Sample$Exp_MAP3K7-mean(Sample$Exp_MAP3K7))/sd(Sample$Exp_MAP3K7)
Sample$Zscore_BACH2 <- (Sample$Exp_BACH2-mean(Sample$Exp_BACH2))/sd(Sample$Exp_BACH2)

a<-cor.test(Sample$Zscore_MAP3K7,Sample$Zscore_BACH2,method = "spearman")
a$p.value
a$estimate

F1a.plot<-ggplot()+theme_classic()
F1a.plot<-F1a.plot+
  geom_point(data=Sample,aes(x=Zscore_BACH2,y=Zscore_MAP3K7,shape=ZCCHC7,color=Histology,fill=PAX5.noncoding),cex=1.8,alpha=.9)+
  #geom_text_repel(data =Sample ,aes(x=Zscore_BACH2,y=Zscore_MAP3K7,label=case),  stat = "identity",position = "identity",max.overlaps = getOption("ggrepel.max.overlaps", default = 30),force=1.5,     segment.size =0.2, size =2.4,color ='black')+
  geom_text(aes(x=2.45,y=2.8,label=paste0("pvalue = ",round(a$p.value,3))),color="blue",size=2.4)+
  geom_text(aes(x=2.45,y=2.65,label=paste0("SCC = ",round(a$estimate,3))),color="blue",size=2.4)
F1a.plot<-F1a.plot+theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,1,1,1),'lines'),legend.title=element_text(size=12,face='plain',color='black'),
                         plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=12,face='bold'),
                         legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position="right",legend.text=element_text(size=10,hjust=0,face='plain'),
                         axis.text.x=element_text(size=12,angle=30,vjust=1,hjust=1,face='bold',color='black'),axis.text.y=element_text(size=12,face='bold',color='black'),
                         axis.title.x=element_text(size=12,vjust=0,hjust=0.5,face='plain',color='black'),axis.title.y=element_text(size=12,face='plain',color='black'))
F1a.plot<-F1a.plot+ggtitle(NULL)+ylab(paste0("Expression of ",gene[2]))+xlab(paste0("Expression of ",gene[1]))
F1a.plot<-F1a.plot+scale_y_continuous(expand=c(0,0),limits=c(-3,3)) 
F1a.plot<-F1a.plot+scale_x_continuous(expand=c(0,0),limits=c(-3,3)) 
F1a.plot<-F1a.plot + scale_fill_manual(values=c(Mut="black",WT="transparent"),guide=guide_legend(override.aes=list(size=2,shape=21),nrow=2))
F1a.plot<-F1a.plot + scale_shape_manual(values=c(24,21))
F1a.plot<-F1a.plot + scale_color_manual(name="Histology",values=gg_color_hue(3),guide=guide_legend(override.aes=list(size=2,shape=21),nrow=3))
plotxxx<-cbind(ggplotGrob(F1a.plot),size="first")
ggsave(file=paste0("./figs_Morin/Linear_",gene[1],"_",gene[2],".pdf"), plot=plotxxx,bg = 'white', width = 17, height = 14, units = 'cm', dpi = 600)








###########################heatmaps
#testing
data <-expMatrix_new
xxx <- samples
samples <- xxx
#heatmap
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 80, c = 90)[1:n]
}

geneA <- data[which(rownames(data) == "ZCCHC7"),]
Sample <- samples
Sample$Exp_geneA <- geneA
Sample$Zscore_geneA <- (Sample$Exp_geneA-mean(Sample$Exp_geneA))/sd(Sample$Exp_geneA)
Sample$gp <- "ZCCHC7_High"
Sample$gp[which(Sample$Zscore_geneA > 0 & Sample$Zscore_geneA < 0)] <- "ZCCHC7_Medium"
Sample$gp[which(Sample$Zscore_geneA < 0 )] <- "ZCCHC7_Low"
samples <-  Sample[which(Sample$gp != "ZCCHC7_Medium"),]


ann_colors = list(
  Histology = c(ABC="#bebada",GCB="#8dd3c7",U="grey"),
  Markers = c(Upregulated="#FAFE01",DownRegulated="#1B00FF",Driver="black"),
  gp =   c(ZCCHC7_High=gg_color_hue(2)[1],ZCCHC7_Low=gg_color_hue(2)[2])
)

list1 <- read.table("./All_marker.txt",header = T)
list1 <- list1[which(list1$His != "others"),]
anno_row <- as.data.frame(list1$His)
colnames(anno_row) <- c("Markers")
rownames(anno_row) <- list1$Marker
MarkerList <- rownames(anno_row)

anno_col <- data.frame(samples[, c("gp","Histology")])

SubtypingMatrix <- expMatrix_new
topVarGenes <- which(rownames(SubtypingMatrix) %in% MarkerList)
mat <- SubtypingMatrix[topVarGenes,which(colnames(SubtypingMatrix) %in% rownames(samples))]

mat #<- mat - rowMeans(mat)
Mat <- mat[c(1,2,3,4,5,6,7,8,12,13,9,10,11),c(2,4,9,10,11,12,13,14,15,16,17,19,1,3,5,6,7,8,18,20)]
anno_col <- anno_col[c(2,4,9,10,11,12,13,14,15,16,17,19,1,3,5,6,7,8,18,20),]

pdf(file = paste0("./CellLine.pdf"),width=9,height=5, onefile = T)
pheatmap(Mat,cutree_rows=2,cluster_cols = F,cluster_rows = F,gaps_row=c(10),gaps_col=c(12),
         color = colorRampPalette(c("#1B00FF", "grey95","#FAFE01"))(30),
         annotation_col = anno_col,annotation_names_col=T,border_color=NA,
         annotation_names_row=F,show_rownames = T,annotation_legend = T,
         annotation_colors = ann_colors,scale = "row",annotation_row = anno_row,
         clustering_distance_cols = "correlation",
         clustering_distance_rows = "euclidean",#cellwidth = 10,# cellheight = 2,
         fontsize = 10,fontsize_row = 10, fontsize_col = 10, display_numbers = F,
         number_format = "%.0f", number_color = "grey60", fontsize_number = 8)
dev.off()
