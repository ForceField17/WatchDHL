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


#selected <- which(samples$Histology!="xxxxxxx")
selected <- which(samples$Histology!="x" & samples$tissue!="Tumor" & samples$Type=="DHL" & samples$case!="RG047" )
#selected <- which(samples$tissue!="xxx" & samples$case!="RG126" & samples$tissue=="CellLine")
#selected <- which(samples$Histology=="GCB")
#selected <- which(samples$Histology!="xxx")
samples <- samples[selected,]
counts <- counts[,selected]
all(rownames(samples) %in% colnames(counts))


## create DESeq2 object
dds=DESeqDataSetFromMatrix(countData = counts,colData = samples,design = ~ ZCCHC7)
nrow(dds)

#keep <- rowSums(counts(dds)) > 10
keep <- rowSums(counts(dds) >= 50) >= 3
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






gene <- c("ZCCHC7",	"AICDA","PAX5")


pair_exp <- data.frame(t(expMatrix[which(rownames(expMatrix) %in% gene),]))

Sample <- samples
Sample$Exp_ZCCHC7 <-  pair_exp[,colnames(pair_exp) == gene[1]] 
Sample$Exp_AICDA <- pair_exp[,colnames(pair_exp) == gene[2]] 
Sample$Exp_PAX5 <- pair_exp[,colnames(pair_exp) == gene[3]] 

#selected <- which(Sample$Class!="xxx" & Sample$Class=="CellLine" & Sample$case!="A3KAW" & Sample$case!="A4FUK" & Sample$case!="RG126" &  Sample$case!="U937")
#Sample <- Sample[selected,]

Sample$Zscore_AICDA <- (Sample$Exp_AICDA-mean(Sample$Exp_AICDA))/sd(Sample$Exp_AICDA)
Sample$Zscore_ZCCHC7 <- (Sample$Exp_ZCCHC7-mean(Sample$Exp_ZCCHC7))/sd(Sample$Exp_ZCCHC7)
Sample$Zscore_PAX5 <- (Sample$Exp_PAX5-mean(Sample$Exp_PAX5))/sd(Sample$Exp_PAX5)

a<-cor.test(Sample$Zscore_AICDA,Sample$Zscore_ZCCHC7,method = "pearson")
a$p.value
a$estimate

F1a.plot<-ggplot()+theme_classic()
F1a.plot<-F1a.plot+
  geom_point(data=Sample,aes(x=Zscore_ZCCHC7,y=Zscore_AICDA,shape=Type,color=ZCCHC7),cex=1.8,alpha=.9)+
  #geom_text_repel(data =Sample ,aes(x=Zscore_ZCCHC7,y=Zscore_AICDA,label=case),  stat = "identity",position = "identity",max.overlaps = getOption("ggrepel.max.overlaps", default = 30),force=1.5,     segment.size =0.2, size =2.4,color ='black')+
  geom_text(aes(x=2.3,y=2.3,label=paste0("pvalue = ",round(a$p.value,3))),color="blue",size=2.4)+
  geom_text(aes(x=2.3,y=2.15,label=paste0("SCC = ",round(a$estimate,3))),color="blue",size=2.4)
F1a.plot<-F1a.plot+theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,1,1,1),'lines'),legend.title=element_text(size=12,face='plain',color='black'),
                         plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=12,face='bold'),
                         legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position="right",legend.text=element_text(size=10,hjust=0,face='plain'),
                         axis.text.x=element_text(size=12,angle=30,vjust=1,hjust=1,face='bold',color='black'),axis.text.y=element_text(size=12,face='bold',color='black'),
                         axis.title.x=element_text(size=12,vjust=0,hjust=0.5,face='plain',color='black'),axis.title.y=element_text(size=12,face='plain',color='black'))
F1a.plot<-F1a.plot+ggtitle(NULL)+ylab(paste0("Expression of ",gene[2]))+xlab(paste0("Expression of ",gene[1]))
F1a.plot<-F1a.plot+scale_y_continuous(expand=c(0,0),limits=c(-2,2.7)) 
F1a.plot<-F1a.plot+scale_x_continuous(expand=c(0,0),limits=c(-2,2.7)) 
F1a.plot<-F1a.plot + scale_fill_manual(values=c("black","red"),guide=guide_legend(override.aes=list(size=2,shape=21),nrow=2))
F1a.plot<-F1a.plot + scale_shape_manual(values=c(24,21,23))
F1a.plot<-F1a.plot + scale_color_manual(values=gg_color_hue(3),guide=guide_legend(override.aes=list(size=2,shape=21),nrow=3))
plotxxx<-cbind(ggplotGrob(F1a.plot),size="first")
ggsave(file=paste0("./Linear_",gene[1],"_",gene[2],".pdf"), plot=plotxxx,bg = 'white', width = 17, height = 14, units = 'cm', dpi = 600)




Sample$xxx <- "ZZZ"
Sample$xxx[which(Sample$case=="SUDHL6")] <- "other"

##################3333
NAN_plot <- ggplot() + theme_classic() 
NAN_plot <- NAN_plot + 
  geom_point(position=position_jitterdodge(jitter.width = 0.2),data=Sample,aes(x=reorder(ZCCHC7,-as.integer(as.factor(ZCCHC7))),y=Zscore_ZCCHC7,color=xxx),size=1.5,width=0.12,alpha=1)+
  geom_boxplot(data=Sample,aes(x=reorder(ZCCHC7,-as.integer(as.factor(ZCCHC7))),y=Zscore_ZCCHC7),width=0.4,alpha=0,size=0.4,outlier.size = 0.4,outlier.shape = NA)
NAN_plot <-NAN_plot + theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,2,2,2),'lines'),
                            plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
                            legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='right',legend.text=element_text(size=14,hjust=0,face='bold'),
                            axis.text.x=element_text(size=12,face='bold',color='black'),axis.text.y=element_text(size=12,face='bold',color='black'),
                            axis.title.x=element_text(size=14,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))
NAN_plot <- NAN_plot + scale_y_continuous(expand=c(0,0),limits=c(-1.8,3.2),breaks = seq(-10,10,1)) 
NAN_plot<- NAN_plot +ylab("Expression of ZCCHC7 (Z-score)") +xlab(NULL)
NAN_plot <- NAN_plot + scale_fill_manual(name=NULL,values = c("#bababa","transparent"))
NAN_plot <- NAN_plot + scale_color_manual(name=NULL,values = c("#e31a1c","#8dd3c7"))
the<-compare_means(p.adjust.method = "none",Zscore_ZCCHC7 ~ ZCCHC7,  data = Sample,paired = F )
#NAN_plot<- NAN_plot +stat_compare_means(comparisons = my_comparisons)
NAN_plot <- NAN_plot +stat_pvalue_manual(
  the, label = NULL, size = 3,linetype = 1,
  y.position = c( 2.9)
)

figure_2<-rbind(ggplotGrob(NAN_plot),size="last")

ggsave(file="ZCCHC7_XXX.pdf", plot=figure_2,bg = 'white', width =10, height = 10, units = 'cm', dpi = 600)


##################3333
NAN_plot <- ggplot() + theme_classic() 
NAN_plot <- NAN_plot + 
  geom_point(position=position_jitterdodge(jitter.width = 0.2),data=Sample,aes(x=reorder(ZCCHC7,-as.integer(as.factor(ZCCHC7))),y=Zscore_PAX5,color=xxx),size=1.5,width=0.12,alpha=1)+
  geom_boxplot(data=Sample,aes(x=reorder(ZCCHC7,-as.integer(as.factor(ZCCHC7))),y=Zscore_PAX5),width=0.4,alpha=0,size=0.4,outlier.size = 0.4,outlier.shape = NA)
NAN_plot <-NAN_plot + theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,2,2,2),'lines'),
                            plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
                            legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='right',legend.text=element_text(size=14,hjust=0,face='bold'),
                            axis.text.x=element_text(size=12,face='bold',color='black'),axis.text.y=element_text(size=12,face='bold',color='black'),
                            axis.title.x=element_text(size=14,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))
NAN_plot <- NAN_plot + scale_y_continuous(expand=c(0,0),limits=c(-1.8,3.2),breaks = seq(-10,10,1)) 
NAN_plot<- NAN_plot +ylab("Expression of PAX5 (Z-score)") +xlab(NULL)
NAN_plot <- NAN_plot + scale_fill_manual(name=NULL,values = c("#bababa","transparent"))
NAN_plot <- NAN_plot + scale_color_manual(name=NULL,values = c("#e31a1c","#8dd3c7"))
the<-compare_means(p.adjust.method = "none",Zscore_PAX5 ~ ZCCHC7,  data = Sample,paired = F )
#NAN_plot<- NAN_plot +stat_compare_means(comparisons = my_comparisons)
NAN_plot <- NAN_plot +stat_pvalue_manual(
  the, label = NULL, size = 3,linetype = 1,
  y.position = c( 2.9)
)

figure_2<-rbind(ggplotGrob(NAN_plot),size="last")

ggsave(file="PAX5_XXX.pdf", plot=figure_2,bg = 'white', width =10, height = 10, units = 'cm', dpi = 600)




