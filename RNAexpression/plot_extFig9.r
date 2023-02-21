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
library(circlize)
library(ggpubr)
library(ggbeeswarm)
library(preprocessCore)
library(ComplexHeatmap)
library(reshape2)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

#raw data preprocessing
Sample_features <- read.delim("./CellLine/Info_allCellLines.txt",sep = "\t",header = T)
samples <- data.frame(Sample_features)
rownames(samples) <- Sample_features$case
samples <- samples[which(samples$Batch %in% c("CCLE","GSE207388") & samples$Type != "Unknown" & samples$Disease=="DLBCL"),]

################
counts1 <- read.csv("./CellLine/RNA_allCellLines.txt", sep="\t", head=T, row.names = "Geneid")
nrow(counts1)
counts1 <- counts1[,which((colnames(counts1 ) %in% samples$case))]

dds=DESeqDataSetFromMatrix(countData = counts1,colData = samples,design = ~ Histology)
keep <- rowSums(counts(dds) >= 50) >= 5
dds <- dds[keep,]
nrow(dds)
vsd <- vst(dds ,blind = TRUE)
head(assay(vsd), 3)
expMatrix <- assay(vsd)
expMatrix_1 <- normalize.quantiles(expMatrix, copy = TRUE)
colnames(expMatrix_1) <- colnames(expMatrix)
rownames(expMatrix_1) <- rownames(expMatrix)

expMatrix_new <- rbind(expMatrix_1)

all(rownames(samples) %in% colnames(expMatrix_new ))
all(colnames(expMatrix_new ) %in% rownames(samples))

#testing
data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}
scale_rows = function(x){
  m = apply(x, 1, mean, na.rm = T)
  s = apply(x, 1, sd, na.rm = T)
  return((x - m) / s)
}

############3
CellLine <- read.delim("./CellLine/cells.txt",sep = "\t",header = T)

Sample <- samples[which(samples$Cell %in% CellLine$Cell),]
data <- expMatrix_new[,which(colnames(expMatrix_new) %in% Sample$case)]

gene <- c("PAX5","ZCCHC7","AICDA")
pair_exp <- data.frame(t(data[which(rownames(data) %in% gene),]))
Sample <- cbind(Sample,pair_exp)

CellLine$exp_ZCCHC7 <- 0
CellLine$exp_PAX5 <- 0
for(i in c(1:nrow(CellLine))){
    CellLine$exp_ZCCHC7[i] <- mean(Sample$ZCCHC7[which(Sample$Cell==CellLine$Cell[i])])
    CellLine$exp_PAX5[i] <- mean(Sample$PAX5[which(Sample$Cell==CellLine$Cell[i])])
}

##################################### All
Sample <- CellLine[which(!is.na(CellLine$exp_ZCCHC7)),]

Sample$ZCCHC7 <- (Sample$exp_ZCCHC7-mean(Sample$exp_ZCCHC7))/sd(Sample$exp_ZCCHC7)
Sample$PAX5 <- (Sample$exp_PAX5-mean(Sample$exp_PAX5))/sd(Sample$exp_PAX5)


Sample$Mutation <- "WT"
Sample$Mutation[which(  Sample$Pos6029>=5  |Sample$Pos6316>=5  |Sample$Pos6681>=5  |Sample$Pos6734>=5  | Sample$CNV == "Gain" | Sample$SV == "Fusion")] <- "Mut"

##################
NAN_plot <- ggplot(data=Sample,aes(x=reorder(Mutation,-as.integer(as.factor(Mutation))),y=ZCCHC7)) + theme_classic() 
NAN_plot <- NAN_plot + 
  geom_boxplot(data=Sample,aes(x=reorder(Mutation,-as.integer(as.factor(Mutation))),y=ZCCHC7),width = 0.4,size=0.5,fill="transparent")+
  geom_quasirandom(data=Sample,aes(x=reorder(Mutation,-as.integer(as.factor(Mutation))),y=ZCCHC7,color=Mutation,fill=Mutation),width = 0.25,size=2,alpha=0.7,stroke=0.8, varwidth = T)
NAN_plot <-NAN_plot + theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,2,2,2),'lines'),
                            plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
                            legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='right',legend.text=element_text(size=14,hjust=0,face='bold'),
                            axis.text.x=element_text(size=12,face='bold',color='black'),axis.text.y=element_text(size=12,face='bold',color='black'),
                            axis.title.x=element_text(size=14,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))
NAN_plot <- NAN_plot + scale_y_continuous(expand=c(0,0),limits=c(-2,3.45),breaks = seq(-10,20,1)) 
NAN_plot<- NAN_plot +ylab("Expression of ZCCHC7 (Z-score)") +xlab(NULL)
NAN_plot <- NAN_plot + scale_fill_manual(name=NULL,values = c(Mut="#d6604d",WT="#4393c3"))
NAN_plot <- NAN_plot + scale_color_manual(name=NULL,values =  c(Mut="#d6604d",WT="#4393c3"))

the<-compare_means(method = "wilcox.test",ZCCHC7 ~ Mutation,  data = Sample,paired = F )
NAN_plot <- NAN_plot +stat_pvalue_manual(
  the, label = NULL, size = 3,linetype = 1,
  y.position = c( 3.2)
)

figure_2<-rbind(ggplotGrob(NAN_plot),size="last")
ggsave(file="./ExtFig9/exp_ZCCHC7_Cell.pdf", plot=figure_2,bg = 'white', width =10, height = 9, units = 'cm', dpi = 600)


##################
NAN_plot <- ggplot(data=Sample,aes(x=reorder(Mutation,-as.integer(as.factor(Mutation))),y=PAX5)) + theme_classic() 
NAN_plot <- NAN_plot + 
  geom_boxplot(data=Sample,aes(x=reorder(Mutation,-as.integer(as.factor(Mutation))),y=PAX5),width = 0.4,size=0.5,fill="transparent",outlier.shape = NA)+
  geom_quasirandom(data=Sample,aes(x=reorder(Mutation,-as.integer(as.factor(Mutation))),y=PAX5,color=Mutation,fill=Mutation),width = 0.25,size=2,alpha=0.7,stroke=0.8, varwidth = T)
NAN_plot <-NAN_plot + theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,2,2,2),'lines'),
                            plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
                            legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='right',legend.text=element_text(size=14,hjust=0,face='bold'),
                            axis.text.x=element_text(size=12,face='bold',color='black'),axis.text.y=element_text(size=12,face='bold',color='black'),
                            axis.title.x=element_text(size=14,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))
NAN_plot <- NAN_plot + scale_y_continuous(expand=c(0,0),limits=c(-2.8,2.45),breaks = seq(-10,20,1)) 
NAN_plot<- NAN_plot +ylab("Expression of PAX5 (Z-score)") +xlab(NULL)
NAN_plot <- NAN_plot + scale_fill_manual(name=NULL,values = c(Mut="#d6604d",WT="#4393c3"))
NAN_plot <- NAN_plot + scale_color_manual(name=NULL,values =  c(Mut="#d6604d",WT="#4393c3"))
the<-compare_means(method = "wilcox.test",PAX5 ~ Mutation,  data = Sample,paired = F )
NAN_plot <- NAN_plot +stat_pvalue_manual(
  the, label = NULL, size = 3,linetype = 1,
  y.position = c( 2.2)
)

figure_2<-rbind(ggplotGrob(NAN_plot),size="last")
ggsave(file="./ExtFig9/exp_PAX5_Cell.pdf", plot=figure_2,bg = 'white', width =10, height = 9, units = 'cm', dpi = 600)




##########
Morethan2 <- c(37024996,37025363,37025364,37026028,37026029,37025178,37025179,37026315,37026316,37026678,37026680,37026681,37026733,37026734)
rank <- c(1:length(Morethan2))*-1
PosRank <- data.frame(Morethan2,rank)
colnames(PosRank ) <- c("Pos","Yrank")

IDrank <- Sample

Table <- read.table("./CellLine//matrix.txt",sep = "\t",header = T)
Table <- Table[which(Table$Depth!="."),]
Table <- Table[which(as.numeric(as.character(Table$Depth))>=5),]

######
case <- unique(IDrank$Cell)

gene.Matrix <- rep('N',length(case))
for( i in 2:length(Morethan2)){
  gene.Matrix <- cbind(gene.Matrix,rep('N',length(case)))
}

for(i in 1:length(case)){
  for(j in 1:length(Morethan2)){
    temp.P <- Table$sub[which(Table$Cell == as.character(case[i]) & Table$Pos == Morethan2[j] )]
    if(length(temp.P) == 1){
      gene.Matrix[i,j] <- as.character(temp.P)
    }
    if(length(temp.P) > 1){
      gene.Matrix[i,j] <- as.character(temp.P[1])
    }
    if("Ins" %in% temp.P){
      gene.Matrix[i,j] <- "Ins"
    }
  }
}

colnames(gene.Matrix) <- Morethan2
rownames(gene.Matrix) <- case
background <- melt(gene.Matrix)
colnames(background) <- c("Cell","Pos","sub")

new_table <- merge(background,PosRank )
new_table <- merge(new_table ,IDrank )

################3
IDrank <- IDrank[which(IDrank$Cell %in% case),]
IDrank$sample <- IDrank$Cell
######
F2a.pLot<-ggplot()+theme_classic()
F2a.pLot<-F2a.pLot+geom_tile(data = IDrank,aes(x=reorder(Cell,Xrank),y="Histology",fill=Histology),color='black',alpha=0.7,width=1,height=1,size=0.4,stat='identity')
F2a.pLot<-F2a.pLot+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(1,0,0,1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.68,'cm'),legend.position='non',axis.line.x = element_blank(),axis.line.y = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=16,face='bold.italic'), axis.text.y=element_text(size=12,face='plain',color='black'),axis.ticks.y = element_blank(),axis.ticks.x = element_blank(),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F2a.pLot<-F2a.pLot+scale_x_discrete(position = "bottom")+xlab(NULL)+ylab(NULL)+scale_fill_manual(name=NULL,values =c(GCB='#bf812d',ABC='#762a83'))
F2a.pLot

#######
F2b.pLot<-ggplot()+theme_classic()
F2b.pLot<-F2b.pLot+geom_tile(data = IDrank,aes(x=reorder(Cell,Xrank),y="ZCCHC7 mutation",fill=CNV),color='black',alpha=1,width=1,height=1,size=0.4,stat='identity')
IDrank_SV <- IDrank[which(IDrank$SV=="Fusion"),]
F2b.pLot<-F2b.pLot+geom_point(data = IDrank_SV,aes(x=reorder(Cell,Xrank),y="ZCCHC7 mutation"),shape=23,fill="black",color='black',alpha=1,size=3)
F2b.pLot<-F2b.pLot+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(0,0,0.2,1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.68,'cm'),legend.position='none',axis.line.x = element_blank(),axis.line.y = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=16,face='bold.italic'), axis.text.y=element_text(size=12,face='plain',color='black'),axis.ticks.y = element_blank(),axis.ticks.x = element_blank(),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F2b.pLot<-F2b.pLot+scale_x_discrete(position = "bottom")+xlab(NULL)+ylab(NULL)+scale_fill_manual(name=NULL,values =c(Gain=gg_color_hue(2)[1],WT="white",Unk="white"))
F2b.pLot


########################################################################################3333

tableA <- new_table[which(new_table$Pos >= 37026028 & new_table$Pos <= 37026029),]
F1b.pLot<-ggplot()+theme_classic()
F1b.pLot<-F1b.pLot+geom_tile(data = tableA,aes(x=reorder(Cell,Xrank),y=reorder(paste0("chr9:",Pos),Yrank),fill=sub),color='black',width=1,height=1,size=0.4,stat='identity')
F1b.pLot<-F1b.pLot+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(0,0,0,1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.68,'cm'),legend.position='none',axis.line.x = element_blank(),axis.line.y = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=16,face='bold.italic'), axis.text.y=element_text(size=12,face='plain',color='black'),axis.ticks.y = element_blank(),axis.ticks.x = element_blank(),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F1b.pLot<-F1b.pLot+scale_x_discrete(position = "bottom")+xlab(NULL)+ylab(NULL)
F1b.pLot<-F1b.pLot+scale_fill_manual(name=NULL,values =c(CT='#66c2a5',CG='#fc8d62' ,Del="#8da0cb",N="white"))
F1b.pLot
Part3 <- F1b.pLot


tableA <- new_table[which(new_table$Pos >= 37026315 & new_table$Pos <= 37026316),]
F1b.pLot<-ggplot()+theme_classic()
F1b.pLot<-F1b.pLot+geom_tile(data = tableA,aes(x=reorder(Cell,Xrank),y=reorder(paste0("chr9:",Pos),Yrank),fill=sub),color='black',width=1,height=1,size=0.4,stat='identity')
F1b.pLot<-F1b.pLot+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(0,0,0,1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.68,'cm'),legend.position='none',axis.line.x = element_blank(),axis.line.y = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=16,face='bold.italic'), axis.text.y=element_text(size=12,face='plain',color='black'),axis.ticks.y = element_blank(),axis.ticks.x = element_blank(),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F1b.pLot<-F1b.pLot+scale_x_discrete(position = "bottom")+xlab(NULL)+ylab(NULL)
F1b.pLot<-F1b.pLot+scale_fill_manual(name=NULL,values =c(CT='#66c2a5',CG='#fc8d62' ,Del="#8da0cb",N="white"))
F1b.pLot
Part32 <- F1b.pLot

tableA <- new_table[which(new_table$Pos >= 37026680 & new_table$Pos <= 37026681),]
F1b.pLot<-ggplot()+theme_classic()
F1b.pLot<-F1b.pLot+geom_tile(data = tableA,aes(x=reorder(Cell,Xrank),y=reorder(paste0("chr9:",Pos),Yrank),fill=sub),color='black',width=1,height=1,size=0.4,stat='identity')
F1b.pLot<-F1b.pLot+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(0,0,0,1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.68,'cm'),legend.position='none',axis.line.x = element_blank(),axis.line.y = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=16,face='bold.italic'), axis.text.y=element_text(size=12,face='plain',color='black'),axis.ticks.y = element_blank(),axis.ticks.x = element_blank(),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F1b.pLot<-F1b.pLot+scale_x_discrete(position = "bottom")+xlab(NULL)+ylab(NULL)
F1b.pLot<-F1b.pLot+scale_fill_manual(name=NULL,values =c(CT='#66c2a5',CG='#fc8d62',CA="#bc80bd",Del="#8da0cb",N="white"))
F1b.pLot
Part4 <- F1b.pLot

tableA <- new_table[which(new_table$Pos >= 37026733 & new_table$Pos <= 37026734),]
F1b.pLot<-ggplot()+theme_classic()
F1b.pLot<-F1b.pLot+geom_tile(data = tableA,aes(x=reorder(Cell,Xrank),y=reorder(paste0("chr9:",Pos),Yrank),fill=sub),color='black',width=1,height=1,size=0.4,stat='identity')
F1b.pLot<-F1b.pLot+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(0,0,0,1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.68,'cm'),legend.position='none',axis.line.x = element_blank(),axis.line.y = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=16,face='bold.italic'), axis.text.y=element_text(size=12,face='plain',color='black'),axis.ticks.y = element_blank(),axis.ticks.x = element_blank(),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F1b.pLot<-F1b.pLot+scale_x_discrete(position = "bottom")+xlab(NULL)+ylab(NULL)
F1b.pLot<-F1b.pLot+scale_fill_manual(name=NULL,values =c(CT='#66c2a5',CG='#fc8d62' ,Del="#8da0cb",N="white"))
F1b.pLot
Part5 <- F1b.pLot

tableA$ZCCHC7[which(tableA$ZCCHC7 < -2)] <- -2
tableA$PAX5[which(tableA$PAX5 < -2)] <- -2
tableA$ZCCHC7[which(tableA$ZCCHC7 >  2)] <-  2
tableA$PAX5[which(tableA$PAX5 >  2)] <-  2

F4c.plot<-ggplot()+theme_classic()
F4c.plot<-F4c.plot+geom_tile(data=tableA,aes(x=reorder(Cell,Xrank),y="ZCCHC7 expression",fill=ZCCHC7),color='black',width=1,height=1,size=0.4,stat='identity')
F4c.plot<-F4c.plot+scale_fill_gradientn(name=NULL,colours=c('#29419E','#698DC9',"grey95",'#F67172','#DC2B18'),limits=c(-2,2),breaks=seq(-3,3,1))
F4c.plot<-F4c.plot+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(0.1,0,0,1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.68,'cm'),legend.position='none',axis.line.x = element_blank(),axis.line.y = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=16,face='bold.italic'), axis.text.y=element_text(size=12,face='plain',color='black'),axis.ticks.y = element_blank(),axis.ticks.x = element_blank(),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F4c.plot<-F4c.plot+scale_x_discrete(position = "bottom")+xlab(NULL)+ylab(NULL)
F4c.plot
Part6 <- F4c.plot

F4c.plot<-ggplot()+theme_classic()
F4c.plot<-F4c.plot+geom_tile(data=tableA,aes(x=reorder(Cell,Xrank),y="PAX5 expression",fill=PAX5),color='black',width=1,height=1,size=0.4,stat='identity')
F4c.plot<-F4c.plot+scale_fill_gradientn(name=NULL,colours=colorRampPalette(c('#29419E','#698DC9',"grey95",'#F67172','#DC2B18'))(50),limits=c(-2,2),breaks=seq(-3,3,1))
F4c.plot<-F4c.plot+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(0,0,1,1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         legend.key.width=unit(1.2,'cm'),legend.key.height=unit(0.5,'cm'),legend.position='none',axis.line.x = element_blank(),axis.line.y = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=12,face='plain'), axis.text.y=element_text(size=12,face='plain',color='black'),axis.ticks.y = element_blank(),axis.ticks.x = element_blank(),
                         axis.text.x=element_text(size=12,angle=90,vjust=0.5,hjust=1,face='plain',color='black'),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F4c.plot<-F4c.plot+scale_x_discrete(position = "bottom")+xlab(NULL)+ylab(NULL)
F4c.plot
Part7 <- F4c.plot

Left <- rbind(ggplotGrob(F2a.pLot),ggplotGrob(F2b.pLot),ggplotGrob(Part3),ggplotGrob(Part32),ggplotGrob(Part4),ggplotGrob(Part5),ggplotGrob(Part6),ggplotGrob(Part7))


#raw data preprocessing
Sample_features <- read.delim("./PrimaryTumor/Info_allPrimaryTumors.txt",sep = "\t",header = T)
samples <- data.frame(Sample_features)
rownames(samples) <- Sample_features$case
samples <- samples[which((samples$Batch == "Morin")),]

################
counts1 <- read.csv("./PrimaryTumor/RNA_allPrimaryTumors.txt", sep="\t", head=T, row.names = "Geneid")
nrow(counts1)
counts1 <- counts1[,which((colnames(counts1 ) %in% samples$case))]

dds=DESeqDataSetFromMatrix(countData = counts1,colData = samples,design = ~ PAX5_TSS2)
keep <- rowSums(counts(dds) >= 50) >= 5
dds <- dds[keep,]
nrow(dds)
vsd <- vst(dds ,blind = TRUE)
head(assay(vsd), 3)
expMatrix <- assay(vsd)
expMatrix_1 <- normalize.quantiles(expMatrix, copy = TRUE)
colnames(expMatrix_1) <- colnames(expMatrix)
rownames(expMatrix_1) <- rownames(expMatrix)

expMatrix_new <- rbind(expMatrix_1)
all(rownames(samples) %in% colnames(expMatrix_new ))
all(colnames(expMatrix_new ) %in% rownames(samples))

##################################### 
CellLine <- read.delim("./PrimaryTumor/tumors.txt",sep = "\t",header = T)
CellLine <- CellLine[which( (CellLine$Data=="WGS" & 
                               CellLine$Pos6028_depth >= 10 & CellLine$Pos6029_depth >= 10 & 
                               CellLine$Pos6315_depth >= 10 & CellLine$Pos6316_depth >= 10 & 
                               CellLine$Pos6680_depth >= 10 & CellLine$Pos6681_depth >= 10 & 
                               CellLine$Pos6733_depth >= 10 & CellLine$Pos6734_depth >= 10 )  | 
                              CellLine$Hotspot>=5 | CellLine$Other3hotspot>=5 | CellLine$CNV == "Gain" | CellLine$SV == "Fusion"),]

Sample <- samples[which(samples$Cell %in% CellLine$Cell),]
data <- expMatrix_new[,which(colnames(expMatrix_new) %in% Sample$case)]

gene <- c("PAX5","ZCCHC7","AICDA")
pair_exp <- data.frame(t(data[which(rownames(data) %in% gene),]))
Sample <- cbind(Sample,pair_exp)

CellLine$exp_ZCCHC7 <- 0
CellLine$exp_PAX5 <- 0
for(i in c(1:nrow(CellLine))){
  CellLine$exp_ZCCHC7[i] <- mean(Sample$ZCCHC7[which(Sample$Cell==CellLine$Cell[i])])
  CellLine$exp_PAX5[i] <- mean(Sample$PAX5[which(Sample$Cell==CellLine$Cell[i])])
}

##################################### All
Sample <- CellLine[which(!is.na(CellLine$exp_ZCCHC7)),]

Sample$ZCCHC7 <- (Sample$exp_ZCCHC7-mean(Sample$exp_ZCCHC7))/sd(Sample$exp_ZCCHC7)
Sample$PAX5 <- (Sample$exp_PAX5-mean(Sample$exp_PAX5))/sd(Sample$exp_PAX5)

Sample$Mutation <- "WT"
Sample$Mutation[which( Sample$Hotspot>=5 | Sample$Other3hotspot>=5  | Sample$SV == "Fusion" | Sample$CNV == "Gain")] <- "Mut"

##################
NAN_plot <- ggplot(data=Sample,aes(x=reorder(Mutation,-as.integer(as.factor(Mutation))),y=ZCCHC7)) + theme_classic() 
NAN_plot <- NAN_plot + 
  geom_boxplot(data=Sample,aes(x=reorder(Mutation,-as.integer(as.factor(Mutation))),y=ZCCHC7),width = 0.4,size=0.5,fill="transparent")+
  geom_quasirandom(data=Sample,aes(x=reorder(Mutation,-as.integer(as.factor(Mutation))),y=ZCCHC7,color=Mutation,fill=Mutation),width = 0.25,size=2,alpha=0.7,stroke=0.8, varwidth = T)
NAN_plot <-NAN_plot + theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,2,2,2),'lines'),
                            plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
                            legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='right',legend.text=element_text(size=14,hjust=0,face='bold'),
                            axis.text.x=element_text(size=12,face='bold',color='black'),axis.text.y=element_text(size=12,face='bold',color='black'),
                            axis.title.x=element_text(size=14,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))
NAN_plot <- NAN_plot + scale_y_continuous(expand=c(0,0),limits=c(-3,2.4),breaks = seq(-10,20,1)) 
NAN_plot<- NAN_plot +ylab("Expression of ZCCHC7 (Z-score)") +xlab(NULL)
NAN_plot <- NAN_plot + scale_fill_manual(name=NULL,values = c(Mut="#d6604d",WT="#4393c3"))
NAN_plot <- NAN_plot + scale_color_manual(name=NULL,values =  c(Mut="#d6604d",WT="#4393c3"))
the<-compare_means(method = "wilcox.test",ZCCHC7 ~ Mutation,  data = Sample,paired = F )

NAN_plot <- NAN_plot +stat_pvalue_manual(
  the, label = NULL, size = 3,linetype = 1,
  y.position = c( 2.1)
)

figure_2<-rbind(ggplotGrob(NAN_plot),size="last")
ggsave(file="./ExtFig9/exp_ZCCHC7_Tumor.pdf", plot=figure_2,bg = 'white', width =10, height = 9, units = 'cm', dpi = 600)


NAN_plot <- ggplot(data=Sample,aes(x=reorder(Mutation,-as.integer(as.factor(Mutation))),y=PAX5)) + theme_classic() 
NAN_plot <- NAN_plot + 
  geom_boxplot(data=Sample,aes(x=reorder(Mutation,-as.integer(as.factor(Mutation))),y=PAX5),width = 0.4,size=0.5,fill="transparent",outlier.shape = NA)+
  geom_quasirandom(data=Sample,aes(x=reorder(Mutation,-as.integer(as.factor(Mutation))),y=PAX5,color=Mutation,fill=Mutation),width = 0.25,size=2,alpha=0.7,stroke=0.8, varwidth = T)
NAN_plot <-NAN_plot + theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,2,2,2),'lines'),
                            plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
                            legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='right',legend.text=element_text(size=14,hjust=0,face='bold'),
                            axis.text.x=element_text(size=12,face='bold',color='black'),axis.text.y=element_text(size=12,face='bold',color='black'),
                            axis.title.x=element_text(size=14,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))
NAN_plot <- NAN_plot + scale_y_continuous(expand=c(0,0),limits=c(-3,2.4),breaks = seq(-10,20,1)) 
NAN_plot<- NAN_plot +ylab("Expression of PAX5 (Z-score)") +xlab(NULL)
NAN_plot <- NAN_plot + scale_fill_manual(name=NULL,values = c(Mut="#d6604d",WT="#4393c3"))
NAN_plot <- NAN_plot + scale_color_manual(name=NULL,values =  c(Mut="#d6604d",WT="#4393c3"))
the<-compare_means(method = "wilcox.test",PAX5 ~ Mutation,  data = Sample,paired = F )
NAN_plot <- NAN_plot +stat_pvalue_manual(
  the, label = NULL, size = 3,linetype = 1,
  y.position = c( 2.1)
)

figure_2<-rbind(ggplotGrob(NAN_plot),size="last")
ggsave(file="./ExtFig9/exp_PAX5_Tumor.pdf", plot=figure_2,bg = 'white', width =10, height = 9, units = 'cm', dpi = 600)


##########################
Morethan2 <- c(37024996,37025363,37025364,37026028,37026029,37025178,37025179,37026315,37026316,37026678,37026680,37026681,37026733,37026734)
rank <- c(1:length(Morethan2))*-1
PosRank <- data.frame(Morethan2,rank)
colnames(PosRank ) <- c("Pos","Yrank")

IDrank <- Sample
xxx <- read.table("./PrimaryTumor/IDrank.txt",sep = "\t",header = T)

IDrank <- merge(IDrank,xxx,1,1)

Table <- read.table("./PrimaryTumor/matrix.txt",sep = "\t",header = T)
Table <- Table[which(Table$Depth!="."),]
Table <- Table[which(as.numeric(as.character(Table$Depth))>=5),]

######
case <- unique(IDrank$Cell)

gene.Matrix <- rep('N',length(case))
for( i in 2:length(Morethan2)){
  gene.Matrix <- cbind(gene.Matrix,rep('N',length(case)))
}

for(i in 1:length(case)){
  for(j in 1:length(Morethan2)){
    temp.P <- Table$sub[which(Table$Cell == as.character(case[i]) & Table$Pos == Morethan2[j] )]
    if(length(temp.P) == 1){
      gene.Matrix[i,j] <- as.character(temp.P)
    }
    if(length(temp.P) > 1){
      gene.Matrix[i,j] <- as.character(temp.P[1])
    }
    if("Ins" %in% temp.P){
      gene.Matrix[i,j] <- "Ins"
    }
  }
}

colnames(gene.Matrix) <- Morethan2
rownames(gene.Matrix) <- case
background <- melt(gene.Matrix)
colnames(background) <- c("Cell","Pos","sub")

new_table <- merge(background,PosRank )
new_table <- merge(new_table ,IDrank )


################3
IDrank <- IDrank[which(IDrank$Cell %in% case),]
IDrank$sample <- IDrank$Cell
######
F2a.pLot<-ggplot()+theme_classic()
F2a.pLot<-F2a.pLot+geom_tile(data = IDrank,aes(x=reorder(Cell,Xrank),y="Histology",fill=Histology),color='black',alpha=0.7,width=1,height=1,size=0.4,stat='identity')
F2a.pLot<-F2a.pLot+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(1,1,0,1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.68,'cm'),legend.position='non',axis.line.x = element_blank(),axis.line.y = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=16,face='bold.italic'), axis.text.y=element_blank(),axis.ticks.y = element_blank(),axis.ticks.x = element_blank(),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_blank())
F2a.pLot<-F2a.pLot+scale_x_discrete(position = "bottom")+xlab(NULL)+ylab(NULL)+scale_fill_manual(name=NULL,values =c(GCB='#bf812d',ABC='#762a83'))
F2a.pLot

#######
F2b.pLot<-ggplot()+theme_classic()
F2b.pLot<-F2b.pLot+geom_tile(data = IDrank,aes(x=reorder(Cell,Xrank),y="ZCCHC7 mutation",fill=CNV),color='black',alpha=1,width=1,height=1,size=0.4,stat='identity')
IDrank_SV <- IDrank[which(IDrank$SV=="Fusion"),]
F2b.pLot<-F2b.pLot+geom_point(data = IDrank_SV,aes(x=reorder(Cell,Xrank),y="ZCCHC7 mutation"),shape=23,fill="black",color='black',alpha=1,size=3)
F2b.pLot<-F2b.pLot+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(0,1,0.2,1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.68,'cm'),legend.position='none',axis.line.x = element_blank(),axis.line.y = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=16,face='bold.italic'), axis.text.y=element_blank(),axis.ticks.y = element_blank(),axis.ticks.x = element_blank(),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_blank())
F2b.pLot<-F2b.pLot+scale_x_discrete(position = "bottom")+xlab(NULL)+ylab(NULL)+scale_fill_manual(name=NULL,values =c(Gain=gg_color_hue(2)[1],WT="white",Unk="grey"))
F2b.pLot


########################################################################################3333

tableA <- new_table[which(new_table$Pos >= 37026028 & new_table$Pos <= 37026029),]
F1b.pLot<-ggplot()+theme_classic()
F1b.pLot<-F1b.pLot+geom_tile(data = tableA,aes(x=reorder(Cell,Xrank),y=reorder(paste0("chr9:",Pos),Yrank),fill=sub),color='black',width=1,height=1,size=0.4,stat='identity')
F1b.pLot<-F1b.pLot+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(0,1,0,1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.68,'cm'),legend.position='none',axis.line.x = element_blank(),axis.line.y = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=16,face='bold.italic'), axis.text.y=element_blank(),axis.ticks.y = element_blank(),axis.ticks.x = element_blank(),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_blank())
F1b.pLot<-F1b.pLot+scale_x_discrete(position = "bottom")+xlab(NULL)+ylab(NULL)
F1b.pLot<-F1b.pLot+scale_fill_manual(name=NULL,values =c(CT='#66c2a5',CG='#fc8d62' ,Del="#8da0cb",N="white"))
F1b.pLot
Part3 <- F1b.pLot

tableA <- new_table[which(new_table$Pos >= 37026315 & new_table$Pos <= 37026316),]
F1b.pLot<-ggplot()+theme_classic()
F1b.pLot<-F1b.pLot+geom_tile(data = tableA,aes(x=reorder(Cell,Xrank),y=reorder(paste0("chr9:",Pos),Yrank),fill=sub),color='black',width=1,height=1,size=0.4,stat='identity')
F1b.pLot<-F1b.pLot+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(0,1,0,1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.68,'cm'),legend.position='none',axis.line.x = element_blank(),axis.line.y = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=16,face='bold.italic'), axis.text.y=element_blank(),axis.ticks.y = element_blank(),axis.ticks.x = element_blank(),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_blank())
F1b.pLot<-F1b.pLot+scale_x_discrete(position = "bottom")+xlab(NULL)+ylab(NULL)
F1b.pLot<-F1b.pLot+scale_fill_manual(name=NULL,values =c(CT='#66c2a5',CG='#fc8d62' ,Del="#8da0cb",N="white"))
F1b.pLot
Part32 <- F1b.pLot

tableA <- new_table[which(new_table$Pos >= 37026680 & new_table$Pos <= 37026681),]
F1b.pLot<-ggplot()+theme_classic()
F1b.pLot<-F1b.pLot+geom_tile(data = tableA,aes(x=reorder(Cell,Xrank),y=reorder(paste0("chr9:",Pos),Yrank),fill=sub),color='black',width=1,height=1,size=0.4,stat='identity')
F1b.pLot<-F1b.pLot+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(0,1,0,1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.68,'cm'),legend.position='none',axis.line.x = element_blank(),axis.line.y = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=16,face='bold.italic'), axis.text.y=element_blank(),axis.ticks.y = element_blank(),axis.ticks.x = element_blank(),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_blank())
F1b.pLot<-F1b.pLot+scale_x_discrete(position = "bottom")+xlab(NULL)+ylab(NULL)
F1b.pLot<-F1b.pLot+scale_fill_manual(name=NULL,values =c(CT='#66c2a5',CG='#fc8d62',CA="#bc80bd",Del="#8da0cb",N="white"))
F1b.pLot
Part4 <- F1b.pLot

tableA <- new_table[which(new_table$Pos >= 37026733 & new_table$Pos <= 37026734),]
F1b.pLot<-ggplot()+theme_classic()
F1b.pLot<-F1b.pLot+geom_tile(data = tableA,aes(x=reorder(Cell,Xrank),y=reorder(paste0("chr9:",Pos),Yrank),fill=sub),color='black',width=1,height=1,size=0.4,stat='identity')
F1b.pLot<-F1b.pLot+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(0,1,0,1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.68,'cm'),legend.position='none',axis.line.x = element_blank(),axis.line.y = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=16,face='bold.italic'), axis.text.y=element_blank(),axis.ticks.y = element_blank(),axis.ticks.x = element_blank(),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_blank())
F1b.pLot<-F1b.pLot+scale_x_discrete(position = "bottom")+xlab(NULL)+ylab(NULL)
F1b.pLot<-F1b.pLot+scale_fill_manual(name=NULL,values =c(CT='#66c2a5',CG='#fc8d62' ,Del="#8da0cb",N="white"))
F1b.pLot
Part5 <- F1b.pLot

tableA$ZCCHC7[which(tableA$ZCCHC7 < -2)] <- -2
tableA$PAX5[which(tableA$PAX5 < -2)] <- -2
tableA$ZCCHC7[which(tableA$ZCCHC7 >  2)] <-  2
tableA$PAX5[which(tableA$PAX5 >  2)] <-  2


F4c.plot<-ggplot()+theme_classic()
F4c.plot<-F4c.plot+geom_tile(data=tableA,aes(x=reorder(Cell,Xrank),y="ZCCHC7 expression",fill=ZCCHC7),color='black',width=1,height=1,size=0.4,stat='identity')
F4c.plot<-F4c.plot+scale_fill_gradientn(name=NULL,colours=c('#29419E','#698DC9',"grey95",'#F67172','#DC2B18'),limits=c(-2,2),breaks=seq(-3,3,1))
F4c.plot<-F4c.plot+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(0.1,1,0,1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.68,'cm'),legend.position='none',axis.line.x = element_blank(),axis.line.y = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=16,face='bold.italic'), axis.text.y=element_blank(),axis.ticks.y = element_blank(),axis.ticks.x = element_blank(),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_blank())
F4c.plot<-F4c.plot+scale_x_discrete(position = "bottom")+xlab(NULL)+ylab(NULL)
F4c.plot
Part6 <- F4c.plot

F4c.plot<-ggplot()+theme_classic()
F4c.plot<-F4c.plot+geom_tile(data=tableA,aes(x=reorder(Cell,Xrank),y="PAX5 expression",fill=PAX5),color='black',width=1,height=1,size=0.4,stat='identity')
F4c.plot<-F4c.plot+scale_fill_gradientn(name=NULL,colours=colorRampPalette(c('#29419E','#698DC9',"grey95",'#F67172','#DC2B18'))(50),limits=c(-2,2),breaks=seq(-3,3,1))
F4c.plot<-F4c.plot+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(0,1,3,1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         legend.key.width=unit(1.2,'cm'),legend.key.height=unit(0.5,'cm'),legend.position='none',axis.line.x = element_blank(),axis.line.y = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=12,face='plain'), axis.text.y=element_blank(),axis.ticks.y = element_blank(),axis.ticks.x = element_blank(),
                         axis.text.x=element_text(size=12,angle=90,vjust=0.5,hjust=1,face='plain',color='black'),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_blank())
F4c.plot<-F4c.plot+scale_x_discrete(position = "bottom")+xlab(NULL)+ylab(NULL)
F4c.plot
Part7 <- F4c.plot

Right<- rbind(ggplotGrob(F2a.pLot),ggplotGrob(F2b.pLot),ggplotGrob(Part3),ggplotGrob(Part32),ggplotGrob(Part4),ggplotGrob(Part5),ggplotGrob(Part6),ggplotGrob(Part7))


figure<-cbind(Left,Right,size="last")
panels <- figure$layout$t[grep("panel", figure$layout$name)]

figure$widths[5]  <- unit(29/34,'null')
figure$widths[14] <- unit(24/34,'null')

figure$heights[panels][1] <- unit(1,'null')
figure$heights[panels][2] <- unit(1.12,'null')
figure$heights[panels][3] <- unit(2,'null')
figure$heights[panels][4] <- unit(2,'null')
figure$heights[panels][5] <- unit(2,'null')
figure$heights[panels][6] <- unit(2,'null')
figure$heights[panels][7] <- unit(1,'null')
figure$heights[panels][8] <- unit(1,'null')

figure$heights[panels][9] <- unit(1,'null')
figure$heights[panels][10] <- unit(1.12,'null')
figure$heights[panels][11] <- unit(2,'null')
figure$heights[panels][12] <- unit(2,'null')
figure$heights[panels][13] <- unit(2,'null')
figure$heights[panels][14] <- unit(2,'null')
figure$heights[panels][15] <- unit(1,'null')
figure$heights[panels][16] <- unit(1,'null')

#grid.draw(figure_1)
ggsave(file="./ExtFig9/Integated_landscape.pdf", plot=figure,bg = 'white', width = 28, height = 10, units = 'cm', dpi = 600)



