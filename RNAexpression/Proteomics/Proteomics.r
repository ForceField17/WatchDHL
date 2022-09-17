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


counts1 <- read.csv("./data/FeatureCount_part3.txt", sep="", head=T, row.names = "CellLine",quote = "\"")
nrow(counts1)


counts <- na.omit(counts1)

Name <- as.character(colnames(counts))
x.term <- as.character(trimws(Name, which = c("both", "left", "right"), whitespace = "[ \t\r\n]"))
Sample_features <- read.delim("./data/Info_Morin.txt",sep = "\t",header = T)

samples <- data.frame(Sample_features)
rownames(samples) <- Sample_features$case

all(rownames(samples) %in% colnames(counts))


selected <- which(samples$Histology!="xxxxxxx")
#selected <- which(samples$tissue!="xxx" & samples$tissue!="CellLine")
#selected <- which(samples$Site=="DLBCL" & samples$case!="RG126" & samples$case!="RG003_1" & samples$case!="RG080_1" & samples$case!="RG078_1")# & samples$Class=="Primary")
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
expMatrix_new <- expMatrix
expMatrix_new <- normalize.quantiles(expMatrix, copy = TRUE)
colnames(expMatrix_new) <- colnames(expMatrix)
rownames(expMatrix_new) <- rownames(expMatrix)
#
#aaa <- colMeans(expMatrix)
#bbb <- colMeans(expMatrix_new)
#xxx_mean <- data.frame(aaa,bbb,samples$Type,samples$Gender)
#distri_mean <- ggplot()+theme_classic()
#distri_mean <- distri_mean + geom_point(data=xxx_mean,aes(x=aaa,y=bbb,color=samples.Type),alpha=0.7,size=3)+theme_bw()
#distri_mean 
#
#figure_1<-rbind(ggplotGrob(distri_mean),size="last")
#
#ggsave(file="./temp/qn.pdf", plot=figure_1,bg = 'white', width = 24, height = 16, units = 'cm', dpi = 600)
#

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
Sample$gp <- "ZCCHC7_2High"
Sample$gp[which(Sample$Zscore_geneA > 0 & Sample$Zscore_geneA < 0)] <- "ZCCHC7_Medium"
Sample$gp[which(Sample$Zscore_geneA < 0 )] <- "ZCCHC7_1Low"
samples <-  Sample[which(Sample$gp != "ZCCHC7_Medium"),]


ann_colors = list(
  Histology = c(ABC="#d01c8b",GCB="#4dac26",U="grey"),
  Markers = c(Upregulated="#FAFE01",DownRegulated="#1B00FF",Driver="white"),
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
Mat <- mat[c(1,3,5,7,9,10,11,12,13,14,8,4,2,6),samples$Order]
anno_col <- anno_col[samples$Order,]

pdf(file = paste0("./DLBCL.pdf"),width=16,height=4.5, onefile = T)
pheatmap(Mat,cutree_rows=2,cluster_cols = F,cluster_rows = F,gaps_row=c(10,13),gaps_col=c(46),
         color = colorRampPalette(c("#000099","#0000CC","#0000FF","#3333FF","#6666FF","#9999FF","#CCCCFF", "grey95","#FFCCCC","#FF9999","#FF6666","#FF3333", "#FF0000","#CC0000","#990000"))(50),
         annotation_col = anno_col,annotation_names_col=T,border_color=NA,
         annotation_names_row=F,show_rownames = T,annotation_legend = T,
         annotation_colors = ann_colors,scale = "row",annotation_row = anno_row,
         clustering_distance_cols = "correlation",
         clustering_distance_rows = "euclidean",#cellwidth = 10,# cellheight = 2,
         fontsize = 10,fontsize_row = 10, fontsize_col = 10, display_numbers = F,
         number_format = "%.0f", number_color = "grey60", fontsize_number = 8)
dev.off()






#################
UpGene <- c("PTTG1","TCL1A","RGS13","LAPTM5","THBS4")

geneA <- data[which(rownames(data) == "ZCCHC7"),]
Zscore_geneA <- (geneA-mean(geneA))/sd(geneA)
yyy <- data.frame(samples$case,rep("ZCCHC7",length(geneA)) ,Zscore_geneA ,samples$gp,rep("up",length(geneA)))
colnames(yyy) <- c("case","gene","Exp","Group","UpDown")

for(i in 1:length(UpGene)){
  geneA <- data[which(rownames(data) == UpGene[i]),]
  Zscore_geneA <- (geneA-mean(geneA))/sd(geneA)
  xxx <- data.frame(samples$case,rep(UpGene[i],length(geneA)) ,Zscore_geneA ,samples$gp,rep("up",length(geneA)))
  colnames(xxx) <- c("case","gene","Exp","Group","UpDown")
  yyy <- rbind(yyy,xxx)
}


xxx <- data.frame(c("aaa"),c("xxx"),c(-100) ,c("ZCCHC7_1Low"),c("up"))
colnames(xxx) <- c("case","gene","Exp","Group","UpDown")
yyy <- rbind(yyy,xxx)

xxx <- data.frame(c("aaa"),c("xxx"),c(-100) ,c("ZCCHC7_2High"),c("up"))
colnames(xxx) <- c("case","gene","Exp","Group","UpDown")
yyy <- rbind(yyy,xxx)

DownGene <- c("TRAF3","IRF5","BANK1","PKMYT1","CBFA2T3","ZMYM3","CTSD","MYSM1","PHLPP1","IKZF3")

for(i in 1:length(DownGene)){
  geneA <- data[which(rownames(data) == DownGene[i]),]
  Zscore_geneA <- (geneA-mean(geneA))/sd(geneA)
  xxx <- data.frame(samples$case,rep(DownGene[i],length(geneA)) ,Zscore_geneA ,samples$gp,rep("down",length(geneA)))
  colnames(xxx) <- c("case","gene","Exp","Group","UpDown")
  yyy <- rbind(yyy,xxx)
}

library(ggpubr)
#yyy <- yyy[which(yyy$gene!="ZCCHC7"),]
yyy$Order <- c(1:nrow(yyy))


GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, draw_group = function(self, data, ..., draw_quantiles = NULL){
  data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
  grp <- data[1,'group']
  newdata <- plyr::arrange(transform(data, x = if(grp%%2==1) xminv else xmaxv), if(grp%%2==1) y else -y)
  newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
  newdata[c(1,nrow(newdata)-1,nrow(newdata)), 'x'] <- round(newdata[1, 'x']) 
  if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
    stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <= 
                                              1))
    quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
    aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
    aesthetics$alpha <- rep(1, nrow(quantiles))
    both <- cbind(quantiles, aesthetics)
    quantile_grob <- GeomPath$draw_panel(both, ...)
    ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
  }
  else {
    ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
  }
})

geom_split_violin <- function (mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, position = position, show.legend = show.legend, inherit.aes = inherit.aes, params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}


NAN_plot <- ggplot(yyy,aes(x=reorder(gene,-Order),y=Exp,fill=Group))  + theme_classic()+
  geom_split_violin(trim = TRUE,width = 1.2,alpha=0.2) + 
  geom_boxplot(width = 0.25, notch = FALSE, notchwidth = .4, size=0.3,color="black",outlier.shape = NA, coef=0,alpha=0.6)
NAN_plot <-NAN_plot + theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,2,2,2),'lines'),
                            plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
                            legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='right',legend.text=element_text(size=14,hjust=0,face='bold'),
                            axis.text.x=element_text(size=12,angle=45,vjust=1,hjust=1,face='italic',color='black'),axis.text.y=element_text(size=12,face='plain',color='black'),
                            axis.title.x=element_text(size=14,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))
NAN_plot <- NAN_plot + scale_y_continuous(expand=c(0,0),limits=c(-3.6,3.6),breaks = seq(-10,10,1)) 
NAN_plot<- NAN_plot +ylab("RNA expression Z-score") +xlab(NULL)
NAN_plot <- NAN_plot + scale_fill_manual(name=NULL,values = c("#69E1E8",gg_color_hue(2)[1]))
NAN_plot <- NAN_plot +   stat_compare_means(label = "p.signif",method = "wilcox.test",aes(group = Group),label.y = c(3.2,3.2,3.2,3.2,3.2,3.2,3.2,3.2,3.2,3.2,3.2,3.2,3.2,3.2,3.2,3.2,3.2))


figure_2<-rbind(ggplotGrob(NAN_plot),size="last")

ggsave(file="All_gene.pdf", plot=figure_2,bg = 'white', width =32, height = 10, units = 'cm', dpi = 600)

