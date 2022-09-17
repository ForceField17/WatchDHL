# CELLOR
setwd("/Users/songdong/Dropbox/Dropbox/DLBCL/main_figures/landscape_noncoding/")
library(ggplot2)
library(gridExtra)
library(grid)


source("./lib/CELLO.R")

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}


interestedGene <- c("CD83","TCL1A","TMSB4X","VMP1","LTB","ZFP36L1","FAM186A","UHRF2","BIRC3","CFAP251","SOCS1","BMP7","ETS1","MIR142","BCR","BLK","SGK1","FOXO1","SEL1L3","MYC","FAM102A","FCHSD2","NCOA3","OSBPL10","PATL2","PIM1","S1PR2","LHPP","PALM2AKAP2","BTG2","MIR3681HG","MYO16","CXCR4","PAX5","ZCCHC7","IRF8","lnc-BCL6-7","RHOH","DTX1","ST6GAL1","CIITA","LPP","BACH2","DMD","BCL6","DIAPH2","BCL7A","IGK","FHIT","DISC1FP1","IMMP2L","IGL","BCL2","IGH")

#UST


Sample_features <- read.table("Genetic_profiles_new.txt",sep = "\t",header = T)
Sample_features[is.na(Sample_features)] <- "f_NA"
sample <- data.frame(Sample_features)
sample <- sample[which(sample$patient=='P8'),]
pat <- unique(sample$patient)
case <- unique(sample$case)
#raw data preprocessing

somatic <- read.delim('./Mut/final/Annotated.frequent.txt',header = T)
somatic <- somatic[which(somatic$CaseID %in% pat),]
FL.a <- somatic[which(somatic$LG_freq >= 5 & somatic$altdepth_LG >= 2),]
FL.a$ID <- "P8.FL.a"
FL.b <- somatic[which(as.numeric(as.character(somatic$LG2_freq)) >= 5 & as.numeric(as.character(somatic$altdepth_LG2)) >= 2),]
FL.b$ID <- "P8.FL.b"
DHL <- somatic[which(somatic$HG_freq >= 5 & somatic$altdepth_HG >= 2),]
DHL$ID <- "P8.DHL"

savi.table <- rbind(FL.a,FL.b,DHL)

savi.Sel <- savi.table[which(savi.table$ID %in% case),]

savi <- data.frame(savi.Sel$chr,savi.Sel$ID,savi.Sel$Gene_Name,savi.Sel$Effect_Impact,savi.Sel$Effect)
colnames(savi) <- c("chr","case","Gene_Name","Impact","Effect")

#analysis for selected genes
selGene <- interestedGene

gene.Matrix <- rep('N',length(case))
for( i in 2:length(selGene)){
  gene.Matrix <- cbind(gene.Matrix,rep('N',length(case)))
}

for(i in 1:length(case)){
  for(j in 1:length(selGene)){
    temp.P <- savi$Effect[which(savi$case == as.character(case[i]) & savi$Gene_Name == selGene[j] & savi$Impact!='HIGH' & savi$Impact!='MODERATE')]
    gene.Matrix[i,j] <- length(temp.P)
  }
}


colnames(gene.Matrix) <- selGene
rownames(gene.Matrix) <- case

mutation_gene_table <- gene.Matrix
#####


sample$order[which(sample$case=="P8.FL.a")] <- 1
sample$order[which(sample$case=="P8.FL.b")] <- 2
sample$order[which(sample$case=="P8.DHL")] <-  3
sample$order[which(sample$case=="P8.FL.a")] <- 1

plot_1b.data <- melt(mutation_gene_table)
colnames(plot_1b.data) <- c('ID','gene','value')


xxx<-as.integer(as.character(plot_1b.data$value))
mmm<-as.integer(as.character(plot_1b.data$value))

mmm[which(xxx == 0)] <-   "a0"
mmm[which(xxx == 1) ] <-  "b1"
mmm[which(xxx >= 2) ] <-  "c2"
mmm[which(xxx >= 5) ] <-  "c5"
mmm[which(xxx >= 10) ] <- "d10"
mmm[which(xxx >= 20) ] <- "e20"
mmm[which(xxx >= 30) ] <- "f30"
mmm[which(xxx >= 40) ] <- "g40"
plot_1b.data <- cbind(plot_1b.data,mmm)
new_table <- merge(plot_1b.data,sample,1,1)


#orderID<-c(1:nrow(plot_1b.data))
F1b.pLot<-ggplot()+theme_classic()
F1b.pLot<-F1b.pLot+geom_tile(data = new_table,aes(x=reorder(ID,order),y=gene,fill=mmm),color='black',width=1,height=1,size=0.5,stat='identity')

F1b.pLot<-F1b.pLot+scale_fill_manual(name="Mutation Number",values =c(a0='white',b1='#fee0d2',c2='#fcbba1',c5='#fb6a4a',d10='#ef3b2c',e20='#cb181d',f30='#a50f15',g40='#67000d'),labels=c(a0='0',b1='1',c2='>=2',c5='>=5',d10='>=10',e20='>=20',f30='>=30',g40='>=40'))
#scale_fill_gradientn(name="Number of mutation",colours=c('white','#fe9929','#de2d26','#67001f'),limits=c(0,1.01),breaks=c(0,0.5,1))
#scale_fill_manual(name="Mutation Number",values =c(a0='white',b1='#fee0d2',c2='#fcbba1',c5='#fb6a4a',d10='#ef3b2c',e20='#cb181d',f30='#a50f15',g40='#67000d'),labels=c(a0='0',b1='1',c2='>=2',c5='>=5',d10='>=10',e20='>=20',f30='>=30',g40='>=40'))
F1b.pLot<-F1b.pLot+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(0.1,0.1,0.1,0.1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.68,'cm'),legend.position='none',axis.ticks.x = element_blank(),axis.line.x = element_blank(),axis.line.y = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=16,face='bold.italic'),axis.ticks.y = element_blank(), axis.text.y=element_blank(),
                         axis.text.x=element_text(size=12,angle=45,vjust=1,hjust=1,face='bold',color='black'),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F1b.pLot<-F1b.pLot+scale_x_discrete(position = "bottom")+xlab(NULL)+ylab(NULL)
F1b.pLot
gene.scale1b = 1 + (ncol(mutation_gene_table) - 14)/14



yyy <- rep("Grade",length(sample$case))
data4<-data.frame(sample$case,yyy,sample$order,sample$Grade)
colnames(data4)<-c("xxx","yyy","zzz","his")

yyy <- rep("Gender",length(sample$case))
data5<-data.frame(sample$case,yyy,sample$order,sample$Gender)
colnames(data5)<-c("xxx","yyy","zzz","his")

yyy <- rep("BCL2.fusion",length(sample$case))
data7b<-data.frame(sample$case,yyy,sample$order,sample$BCL2)
colnames(data7b)<-c("xxx","yyy","zzz","his")

yyy <- rep("MYC.fusion",length(sample$case))
data9<-data.frame(sample$case,yyy,sample$order,sample$MYC)
colnames(data9)<-c("xxx","yyy","zzz","his")

F4b.pLot<-ggplot()+theme_classic()
F4b.pLot<-F4b.pLot+geom_tile(data = data4,aes(x=reorder(xxx,zzz),y=yyy,fill=his),color='black',width=1,height=1,size=0.5,stat='identity')
F4b.pLot<-F4b.pLot+scale_fill_manual(name=NULL,values=c(High=gg_color_hue(6)[1] ,Low=gg_color_hue(6)[3]),guide = guide_legend(override.aes=list(size=1),nrow=1),na.translate = F)
F4b.pLot<-F4b.pLot+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(1,0.1,0.1,0.1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.4,'cm'),legend.key.height=unit(0.4,'cm'),legend.position='none',axis.ticks.x = element_blank(),axis.line.x = element_blank(),axis.line.y = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=12,face='bold.italic'),axis.ticks.y = element_blank(), axis.text.y=element_blank(),legend.title=element_blank(),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F4b.pLot<-F4b.pLot+scale_x_discrete(position = "top")+xlab(NULL)+ylab(NULL)+guides(size=FALSE)
gene.scale4 = 1/14

F7b.pLot<-ggplot()+theme_classic()
F7b.pLot<-F7b.pLot+geom_tile(data = data7b,aes(x=reorder(xxx,zzz),y=yyy,fill=his),color='black',width=1,height=1,size=0.5,stat='identity')
F7b.pLot<-F7b.pLot+scale_fill_manual(name=NULL,values=c(Fusion='black',WT='white'),
                                     labels=c(Fusion='Fusion',WT='Wild type'),guide = guide_legend(override.aes=list(size=1),nrow=1),na.translate = F)
F7b.pLot<-F7b.pLot+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(0.1,0.1,0.1,0.1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.4,'cm'),legend.key.height=unit(0.4,'cm'),legend.position='none',axis.ticks.x = element_blank(),axis.line.x = element_blank(),axis.line.y = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=12,face='bold.italic'),axis.ticks.y = element_blank(), axis.text.y=element_blank(),legend.title=element_blank(),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F7b.pLot<-F7b.pLot+scale_x_discrete(position = "top")+xlab(NULL)+ylab(NULL)+guides(size=FALSE)
gene.scale7 = 1/14



F5b.pLot<-ggplot()+theme_classic()
F5b.pLot<-F5b.pLot+geom_tile(data = data5,aes(x=reorder(xxx,zzz),y=yyy,fill=his),color='black',width=1,height=1,size=0.5,stat='identity')
F5b.pLot<-F5b.pLot+scale_fill_manual(name=NULL,values=c(Male='#b3cde3',Female='#fccde5'),labels=c(Male='Male',Female='Female'),guide = guide_legend(override.aes=list(size=1),nrow=1),na.translate = F)
F5b.pLot<-F5b.pLot+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(0.1,0.1,0.1,0.1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.4,'cm'),legend.key.height=unit(0.4,'cm'),legend.position='none',axis.ticks.x = element_blank(),axis.line.x = element_blank(),axis.line.y = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=12,face='bold.italic'),axis.ticks.y = element_blank(), axis.text.y=element_blank(),legend.title=element_blank(),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F5b.pLot<-F5b.pLot+scale_x_discrete(position = "top")+xlab(NULL)+ylab(NULL)+guides(size=FALSE)
gene.scale5 = 1/14



F9b.pLot<-ggplot()+theme_classic()
F9b.pLot<-F9b.pLot+geom_tile(data = data9,aes(x=reorder(xxx,zzz),y=yyy,fill=his),color='black',width=1,height=1,size=0.5,stat='identity')
F9b.pLot<-F9b.pLot+scale_fill_manual(name=NULL,values=c(Fusion='black',WT='white'),
                                     labels=c(Fusion='Fusion',WT='Wild type'),guide = guide_legend(override.aes=list(size=1),nrow=1),na.translate = F)
F9b.pLot<-F9b.pLot+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(0.1,0.1,0.1,0.1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(1.1,'cm'),legend.key.height=unit(0.4,'cm'),legend.position='none',axis.ticks.x = element_blank(),axis.line.x = element_blank(),axis.line.y = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=12,face='bold'),axis.ticks.y = element_blank(), axis.text.y=element_blank(),legend.title=element_blank(),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F9b.pLot<-F9b.pLot+scale_x_discrete(position = "top")+xlab(NULL)+ylab(NULL)+guides(size=FALSE)
gene.scale11 = 1/14



###P1
sample <- data.frame(Sample_features)
sample <- sample[which(sample$patient=='P1'),]
pat <- unique(sample$patient)
case <- unique(sample$case)
#raw data preprocessing

somatic <- read.delim('./Mut/final/Annotated.frequent.txt',header = T)
somatic <- somatic[which(somatic$CaseID %in% pat),]
FL <- somatic[which(somatic$LG_freq >= 5 & somatic$altdepth_LG >= 2),]
FL$ID <- "P1.FL"
DHL <- somatic[which(somatic$HG_freq >= 5 & somatic$altdepth_HG >= 2),]
DHL$ID <- "P1.DHL"

savi.table <- rbind(FL,DHL)

savi.Sel <- savi.table[which(savi.table$ID %in% case),]

savi <- data.frame(savi.Sel$chr,savi.Sel$ID,savi.Sel$Gene_Name,savi.Sel$Effect_Impact,savi.Sel$Effect)
colnames(savi) <- c("chr","case","Gene_Name","Impact","Effect")

#analysis for selected genes
selGene <- interestedGene

gene.Matrix <- rep('N',length(case))
for( i in 2:length(selGene)){
  gene.Matrix <- cbind(gene.Matrix,rep('N',length(case)))
}

for(i in 1:length(case)){
  for(j in 1:length(selGene)){
    temp.P <- savi$Effect[which(savi$case == as.character(case[i]) & savi$Gene_Name == selGene[j] & savi$Impact!='HIGH' & savi$Impact!='MODERATE')]
    gene.Matrix[i,j] <- length(temp.P)
  }
}


colnames(gene.Matrix) <- selGene
rownames(gene.Matrix) <- case

mutation_gene_table <- gene.Matrix
#####


sample$order[which(sample$case=="P1.FL")] <- 1
sample$order[which(sample$case=="P1.DHL")] <-  2
sample$order[which(sample$case=="P1.FL")] <- 1

plot_1b.data <- melt(mutation_gene_table)
colnames(plot_1b.data) <- c('ID','gene','value')


xxx<-as.integer(as.character(plot_1b.data$value))
mmm<-as.integer(as.character(plot_1b.data$value))

mmm[which(xxx == 0)] <-   "a0"
mmm[which(xxx == 1) ] <-  "b1"
mmm[which(xxx >= 2) ] <-  "c2"
mmm[which(xxx >= 5) ] <-  "c5"
mmm[which(xxx >= 10) ] <- "d10"
mmm[which(xxx >= 20) ] <- "e20"
mmm[which(xxx >= 30) ] <- "f30"
mmm[which(xxx >= 40) ] <- "g40"
plot_1b.data <- cbind(plot_1b.data,mmm)
new_table <- merge(plot_1b.data,sample,1,1)

#orderID<-c(1:nrow(plot_1b.data))
F1b.ploT<-ggplot()+theme_classic()
F1b.ploT<-F1b.ploT+geom_tile(data = new_table,aes(x=reorder(ID,order),y=gene,fill=mmm),color='black',width=1,height=1,size=0.5,stat='identity')

F1b.ploT<-F1b.ploT+scale_fill_manual(name="Mutation Number",values =c(a0='white',b1='#fee0d2',c2='#fcbba1',c5='#fb6a4a',d10='#ef3b2c',e20='#cb181d',f30='#a50f15',g40='#67000d'),labels=c(a0='0',b1='1',c2='>=2',c5='>=5',d10='>=10',e20='>=20',f30='>=30',g40='>=40'))
F1b.ploT<-F1b.ploT+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(0.1,0.1,0.1,0.1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.68,'cm'),legend.position='none',axis.ticks.x = element_blank(),axis.line.x = element_blank(),axis.line.y = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=16,face='bold.italic'),axis.text.y=element_text(size=14,vjust=0.5,hjust=1,face='bold.italic',color='black'),axis.line = element_blank(),axis.ticks.y.left = element_line(),
                         axis.text.x=element_text(size=12,angle=45,vjust=1,hjust=1,face='bold',color='black'),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F1b.ploT<-F1b.ploT+scale_x_discrete(position = "bottom")+xlab(NULL)+ylab(NULL)
F1b.ploT
gene.scale1b = 1 + (ncol(mutation_gene_table) - 14)/14



yyy <- rep("Grade",length(sample$case))
data4<-data.frame(sample$case,yyy,sample$order,sample$Grade)
colnames(data4)<-c("xxx","yyy","zzz","his")

yyy <- rep("Gender",length(sample$case))
data5<-data.frame(sample$case,yyy,sample$order,sample$Gender)
colnames(data5)<-c("xxx","yyy","zzz","his")

yyy <- rep("Age",length(sample$case))
data6<-data.frame(sample$case,yyy,sample$order,sample$Age)
colnames(data6)<-c("xxx","yyy","zzz","his")

yyy <- rep("BCL2.fusion",length(sample$case))
data7b<-data.frame(sample$case,yyy,sample$order,sample$BCL2)
colnames(data7b)<-c("xxx","yyy","zzz","his")

yyy <- rep("MYC.fusion",length(sample$case))
data9<-data.frame(sample$case,yyy,sample$order,sample$MYC)
colnames(data9)<-c("xxx","yyy","zzz","his")

F4b.ploT<-ggplot()+theme_classic()
F4b.ploT<-F4b.ploT+geom_tile(data = data4,aes(x=reorder(xxx,zzz),y=yyy,fill=his),color='black',width=1,height=1,size=0.5,stat='identity')
F4b.ploT<-F4b.ploT+scale_fill_manual(name=NULL,values=c(High=gg_color_hue(6)[1] ,Low=gg_color_hue(6)[3]),guide = guide_legend(override.aes=list(size=1),nrow=1),na.translate = F)
F4b.ploT<-F4b.ploT+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(1,0.1,0.1,0.1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.4,'cm'),legend.key.height=unit(0.4,'cm'),legend.position='none',axis.ticks.x = element_blank(),axis.line.x = element_blank(),axis.line.y = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=12,face='bold.italic'),axis.text.y=element_text(size=14,vjust=0.5,hjust=1,face='bold.italic',color='black'),legend.title=element_blank(),axis.ticks.y = element_blank(),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F4b.ploT<-F4b.ploT+scale_x_discrete(position = "top")+xlab(NULL)+ylab(NULL)+guides(size=FALSE)
gene.scale4 = 1/14

F7b.ploT<-ggplot()+theme_classic()
F7b.ploT<-F7b.ploT+geom_tile(data = data7b,aes(x=reorder(xxx,zzz),y=yyy,fill=his),color='black',width=1,height=1,size=0.5,stat='identity')
F7b.ploT<-F7b.ploT+scale_fill_manual(name=NULL,values=c(Fusion='black',WT='white'),
                                     labels=c(Fusion='Fusion',WT='Wild type'),guide = guide_legend(override.aes=list(size=1),nrow=1),na.translate = F)
F7b.ploT<-F7b.ploT+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(0.1,0.1,0.1,0.1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.4,'cm'),legend.key.height=unit(0.4,'cm'),legend.position='none',axis.ticks.x = element_blank(),axis.line.x = element_blank(),axis.line.y = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=12,face='bold.italic'),axis.text.y=element_text(size=14,vjust=0.5,hjust=1,face='bold.italic',color='black'),legend.title=element_blank(),axis.ticks.y = element_blank(),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F7b.ploT<-F7b.ploT+scale_x_discrete(position = "top")+xlab(NULL)+ylab(NULL)+guides(size=FALSE)
gene.scale7 = 1/14



F5b.ploT<-ggplot()+theme_classic()
F5b.ploT<-F5b.ploT+geom_tile(data = data5,aes(x=reorder(xxx,zzz),y=yyy,fill=his),color='black',width=1,height=1,size=0.5,stat='identity')
F5b.ploT<-F5b.ploT+scale_fill_manual(name=NULL,values=c(Male='#b3cde3',Female='#fccde5'),labels=c(Male='Male',Female='Female'),guide = guide_legend(override.aes=list(size=1),nrow=1),na.translate = F)
F5b.ploT<-F5b.ploT+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(0.1,0.1,0.1,0.1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.4,'cm'),legend.key.height=unit(0.4,'cm'),legend.position='none',axis.ticks.x = element_blank(),axis.line.x = element_blank(),axis.line.y = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=12,face='bold.italic'),axis.text.y=element_text(size=14,vjust=0.5,hjust=1,face='bold.italic',color='black'),legend.title=element_blank(),axis.ticks.y = element_blank(),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F5b.ploT<-F5b.ploT+scale_x_discrete(position = "top")+xlab(NULL)+ylab(NULL)+guides(size=FALSE)
gene.scale5 = 1/14

F6b.ploT<-ggplot()+theme_classic()
F6b.ploT<-F6b.ploT+geom_tile(data = data6,aes(x=reorder(xxx,zzz),y=yyy,fill=as.numeric(as.character(his))),color='black',width=1,height=1,size=0.5,stat='identity')
F6b.ploT<-F6b.ploT+scale_fill_gradientn(name=NULL,colours=c('white','#f6e8c3','#dfc27d','#bf812d','#8c510a'),limits=c(50,80),breaks=seq(50,80,5))
F6b.ploT<-F6b.ploT+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(0.1,0.1,0.1,0.1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(1.1,'cm'),legend.key.height=unit(0.4,'cm'),legend.position='none',axis.ticks.x = element_blank(),axis.line.x = element_blank(),axis.line.y = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=12,face='bold'),axis.text.y=element_text(size=14,vjust=0.5,hjust=1,face='bold.italic',color='black'),legend.title=element_blank(),axis.ticks.y = element_blank(),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F6b.ploT<-F6b.ploT+scale_x_discrete(position = "top")+xlab(NULL)+ylab(NULL)+guides(size=FALSE)
gene.scale6 = 1/14




F9b.ploT<-ggplot()+theme_classic()
F9b.ploT<-F9b.ploT+geom_tile(data = data9,aes(x=reorder(xxx,zzz),y=yyy,fill=his),color='black',width=1,height=1,size=0.5,stat='identity')
F9b.ploT<-F9b.ploT+scale_fill_manual(name=NULL,values=c(Fusion='black',WT='white'),
                                     labels=c(Fusion='Fusion',WT='Wild type'),guide = guide_legend(override.aes=list(size=1),nrow=1),na.translate = F)
F9b.ploT<-F9b.ploT+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(0.1,0.1,0.1,0.1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(1.1,'cm'),legend.key.height=unit(0.4,'cm'),legend.position='none',axis.ticks.x = element_blank(),axis.line.x = element_blank(),axis.line.y = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=12,face='bold'),axis.text.y=element_text(size=14,vjust=0.5,hjust=1,face='bold.italic',color='black'),legend.title=element_blank(),axis.ticks.y = element_blank(),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F9b.ploT<-F9b.ploT+scale_x_discrete(position = "top")+xlab(NULL)+ylab(NULL)+guides(size=FALSE)
gene.scale11 = 1/14


##P2
sample <- data.frame(Sample_features)
sample <- sample[which(sample$patient=='P2'),]
pat <- unique(sample$patient)
case <- unique(sample$case)
#raw data preprocessing

somatic <- read.delim('./Mut/final/Annotated.frequent.txt',header = T)
somatic <- somatic[which(somatic$CaseID %in% pat),]
FL <- somatic[which(somatic$LG_freq >= 5 & somatic$altdepth_LG >= 2),]
FL$ID <- "P2.FL"
DHL <- somatic[which(somatic$HG_freq >= 5 & somatic$altdepth_HG >= 2),]
DHL$ID <- "P2.DHL"

savi.table <- rbind(FL,DHL)

savi.Sel <- savi.table[which(savi.table$ID %in% case),]

savi <- data.frame(savi.Sel$chr,savi.Sel$ID,savi.Sel$Gene_Name,savi.Sel$Effect_Impact,savi.Sel$Effect)
colnames(savi) <- c("chr","case","Gene_Name","Impact","Effect")

#analysis for selected genes
selGene <- interestedGene

gene.Matrix <- rep('N',length(case))
for( i in 2:length(selGene)){
  gene.Matrix <- cbind(gene.Matrix,rep('N',length(case)))
}

for(i in 1:length(case)){
  for(j in 1:length(selGene)){
    temp.P <- savi$Effect[which(savi$case == as.character(case[i]) & savi$Gene_Name == selGene[j] & savi$Impact!='HIGH' & savi$Impact!='MODERATE')]
    gene.Matrix[i,j] <- length(temp.P)
  }
}


colnames(gene.Matrix) <- selGene
rownames(gene.Matrix) <- case

mutation_gene_table <- gene.Matrix
#####


sample$order[which(sample$case=="P2.FL")] <- 1
sample$order[which(sample$case=="P2.DHL")] <-  2
sample$order[which(sample$case=="P2.FL")] <- 1

plot_1b.data <- melt(mutation_gene_table)
colnames(plot_1b.data) <- c('ID','gene','value')


xxx<-as.integer(as.character(plot_1b.data$value))
mmm<-as.integer(as.character(plot_1b.data$value))

mmm[which(xxx == 0)] <-   "a0"
mmm[which(xxx == 1) ] <-  "b1"
mmm[which(xxx >= 2) ] <-  "c2"
mmm[which(xxx >= 5) ] <-  "c5"
mmm[which(xxx >= 10) ] <- "d10"
mmm[which(xxx >= 20) ] <- "e20"
mmm[which(xxx >= 30) ] <- "f30"
mmm[which(xxx >= 40) ] <- "g40"
plot_1b.data <- cbind(plot_1b.data,mmm)
new_table <- merge(plot_1b.data,sample,1,1)

#orderID<-c(1:nrow(plot_1b.data))
F1b.plOT<-ggplot()+theme_classic()
F1b.plOT<-F1b.plOT+geom_tile(data = new_table,aes(x=reorder(ID,order),y=gene,fill=mmm),color='black',width=1,height=1,size=0.5,stat='identity')

F1b.plOT<-F1b.plOT+scale_fill_manual(name="Mutation Number",values =c(a0='white',b1='#fee0d2',c2='#fcbba1',c5='#fb6a4a',d10='#ef3b2c',e20='#cb181d',f30='#a50f15',g40='#67000d'),labels=c(a0='0',b1='1',c2='>=2',c5='>=5',d10='>=10',e20='>=20',f30='>=30',g40='>=40'))
F1b.plOT<-F1b.plOT+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(0.1,0.1,0.1,0.1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.68,'cm'),legend.position='none',axis.ticks.x = element_blank(),axis.line.x = element_blank(),axis.line.y = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=16,face='bold.italic'),axis.ticks.y = element_blank(), axis.text.y=element_blank(),
                         axis.text.x=element_text(size=12,angle=45,vjust=1,hjust=1,face='bold',color='black'),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F1b.plOT<-F1b.plOT+scale_x_discrete(position = "bottom")+xlab(NULL)+ylab(NULL)
F1b.plOT
gene.scale1b = 1 + (ncol(mutation_gene_table) - 14)/14



yyy <- rep("Grade",length(sample$case))
data4<-data.frame(sample$case,yyy,sample$order,sample$Grade)
colnames(data4)<-c("xxx","yyy","zzz","his")

yyy <- rep("Gender",length(sample$case))
data5<-data.frame(sample$case,yyy,sample$order,sample$Gender)
colnames(data5)<-c("xxx","yyy","zzz","his")

yyy <- rep("Age",length(sample$case))
data6<-data.frame(sample$case,yyy,sample$order,sample$Age)
colnames(data6)<-c("xxx","yyy","zzz","his")

yyy <- rep("BCL2.fusion",length(sample$case))
data7b<-data.frame(sample$case,yyy,sample$order,sample$BCL2)
colnames(data7b)<-c("xxx","yyy","zzz","his")

yyy <- rep("MYC.fusion",length(sample$case))
data9<-data.frame(sample$case,yyy,sample$order,sample$MYC)
colnames(data9)<-c("xxx","yyy","zzz","his")

F4b.plOT<-ggplot()+theme_classic()
F4b.plOT<-F4b.plOT+geom_tile(data = data4,aes(x=reorder(xxx,zzz),y=yyy,fill=his),color='black',width=1,height=1,size=0.5,stat='identity')
F4b.plOT<-F4b.plOT+scale_fill_manual(name=NULL,values=c(High=gg_color_hue(6)[1] ,Low=gg_color_hue(6)[3]),guide = guide_legend(override.aes=list(size=1),nrow=1),na.translate = F)
F4b.plOT<-F4b.plOT+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(1,0.1,0.1,0.1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.4,'cm'),legend.key.height=unit(0.4,'cm'),legend.position='none',axis.ticks.x = element_blank(),axis.line.x = element_blank(),axis.line.y = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=12,face='bold.italic'),axis.ticks.y = element_blank(), axis.text.y=element_blank(),legend.title=element_blank(),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F4b.plOT<-F4b.plOT+scale_x_discrete(position = "top")+xlab(NULL)+ylab(NULL)+guides(size=FALSE)
gene.scale4 = 1/14

F7b.plOT<-ggplot()+theme_classic()
F7b.plOT<-F7b.plOT+geom_tile(data = data7b,aes(x=reorder(xxx,zzz),y=yyy,fill=his),color='black',width=1,height=1,size=0.5,stat='identity')
F7b.plOT<-F7b.plOT+scale_fill_manual(name=NULL,values=c(Fusion='black',WT='white'),
                                     labels=c(Fusion='Fusion',WT='Wild type'),guide = guide_legend(override.aes=list(size=1),nrow=1),na.translate = F)
F7b.plOT<-F7b.plOT+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(0.1,0.1,0.1,0.1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.4,'cm'),legend.key.height=unit(0.4,'cm'),legend.position='none',axis.ticks.x = element_blank(),axis.line.x = element_blank(),axis.line.y = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=12,face='bold.italic'),axis.ticks.y = element_blank(), axis.text.y=element_blank(),legend.title=element_blank(),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F7b.plOT<-F7b.plOT+scale_x_discrete(position = "top")+xlab(NULL)+ylab(NULL)+guides(size=FALSE)
gene.scale7 = 1/14



F5b.plOT<-ggplot()+theme_classic()
F5b.plOT<-F5b.plOT+geom_tile(data = data5,aes(x=reorder(xxx,zzz),y=yyy,fill=his),color='black',width=1,height=1,size=0.5,stat='identity')
F5b.plOT<-F5b.plOT+scale_fill_manual(name=NULL,values=c(Male='#b3cde3',Female='#fccde5'),labels=c(Male='Male',Female='Female'),guide = guide_legend(override.aes=list(size=1),nrow=1),na.translate = F)
F5b.plOT<-F5b.plOT+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(0.1,0.1,0.1,0.1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.4,'cm'),legend.key.height=unit(0.4,'cm'),legend.position='none',axis.ticks.x = element_blank(),axis.line.x = element_blank(),axis.line.y = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=12,face='bold.italic'),axis.ticks.y = element_blank(), axis.text.y=element_blank(),legend.title=element_blank(),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F5b.plOT<-F5b.plOT+scale_x_discrete(position = "top")+xlab(NULL)+ylab(NULL)+guides(size=FALSE)
gene.scale5 = 1/14

F6b.plOT<-ggplot()+theme_classic()
F6b.plOT<-F6b.plOT+geom_tile(data = data6,aes(x=reorder(xxx,zzz),y=yyy,fill=as.numeric(as.character(his))),color='black',width=1,height=1,size=0.5,stat='identity')
F6b.plOT<-F6b.plOT+scale_fill_gradientn(name=NULL,colours=c('white','#f6e8c3','#dfc27d','#bf812d','#8c510a'),limits=c(50,80),breaks=seq(50,80,5))
F6b.plOT<-F6b.plOT+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(0.1,0.1,0.1,0.1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(1.1,'cm'),legend.key.height=unit(0.4,'cm'),legend.position='none',axis.ticks.x = element_blank(),axis.line.x = element_blank(),axis.line.y = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=12,face='bold'),axis.ticks.y = element_blank(), axis.text.y=element_blank(),legend.title=element_blank(),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F6b.plOT<-F6b.plOT+scale_x_discrete(position = "top")+xlab(NULL)+ylab(NULL)+guides(size=FALSE)
gene.scale6 = 1/14




F9b.plOT<-ggplot()+theme_classic()
F9b.plOT<-F9b.plOT+geom_tile(data = data9,aes(x=reorder(xxx,zzz),y=yyy,fill=his),color='black',width=1,height=1,size=0.5,stat='identity')
F9b.plOT<-F9b.plOT+scale_fill_manual(name=NULL,values=c(Fusion='black',WT='white'),
                                     labels=c(Fusion='Fusion',WT='Wild type'),guide = guide_legend(override.aes=list(size=1),nrow=1),na.translate = F)
F9b.plOT<-F9b.plOT+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(0.1,0.1,0.1,0.1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(1.1,'cm'),legend.key.height=unit(0.4,'cm'),legend.position='none',axis.ticks.x = element_blank(),axis.line.x = element_blank(),axis.line.y = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=12,face='bold'),axis.ticks.y = element_blank(), axis.text.y=element_blank(),legend.title=element_blank(),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F9b.plOT<-F9b.plOT+scale_x_discrete(position = "top")+xlab(NULL)+ylab(NULL)+guides(size=FALSE)
gene.scale11 = 1/14

###P3
sample <- data.frame(Sample_features)
sample <- sample[which(sample$patient=='P3'),]
pat <- unique(sample$patient)
case <- unique(sample$case)
#raw data preprocessing

somatic <- read.delim('./Mut/final/Annotated.frequent.txt',header = T)
somatic <- somatic[which(somatic$CaseID %in% pat),]
FL <- somatic[which(somatic$LG_freq >= 5 & somatic$altdepth_LG >= 2),]
FL$ID <- "P3.FL"
DHL <- somatic[which(somatic$HG_freq >= 5 & somatic$altdepth_HG >= 2),]
DHL$ID <- "P3.DHL"

savi.table <- rbind(FL,DHL)

savi.Sel <- savi.table[which(savi.table$ID %in% case),]

savi <- data.frame(savi.Sel$chr,savi.Sel$ID,savi.Sel$Gene_Name,savi.Sel$Effect_Impact,savi.Sel$Effect)
colnames(savi) <- c("chr","case","Gene_Name","Impact","Effect")

#analysis for selected genes
selGene <- interestedGene

gene.Matrix <- rep('N',length(case))
for( i in 2:length(selGene)){
  gene.Matrix <- cbind(gene.Matrix,rep('N',length(case)))
}

for(i in 1:length(case)){
  for(j in 1:length(selGene)){
    temp.P <- savi$Effect[which(savi$case == as.character(case[i]) & savi$Gene_Name == selGene[j] & savi$Impact!='HIGH' & savi$Impact!='MODERATE')]
    gene.Matrix[i,j] <- length(temp.P)
  }
}


colnames(gene.Matrix) <- selGene
rownames(gene.Matrix) <- case

mutation_gene_table <- gene.Matrix
#####


sample$order[which(sample$case=="P3.FL")] <- 1
sample$order[which(sample$case=="P3.DHL")] <-  2
sample$order[which(sample$case=="P3.FL")] <- 1

plot_1b.data <- melt(mutation_gene_table)
colnames(plot_1b.data) <- c('ID','gene','value')


xxx<-as.integer(as.character(plot_1b.data$value))
mmm<-as.integer(as.character(plot_1b.data$value))

mmm[which(xxx == 0)] <-   "a0"
mmm[which(xxx == 1) ] <-  "b1"
mmm[which(xxx >= 2) ] <-  "c2"
mmm[which(xxx >= 5) ] <-  "c5"
mmm[which(xxx >= 10) ] <- "d10"
mmm[which(xxx >= 20) ] <- "e20"
mmm[which(xxx >= 30) ] <- "f30"
mmm[which(xxx >= 40) ] <- "g40"
plot_1b.data <- cbind(plot_1b.data,mmm)
new_table <- merge(plot_1b.data,sample,1,1)

#orderID<-c(1:nrow(plot_1b.data))
F1b.pLOT<-ggplot()+theme_classic()
F1b.pLOT<-F1b.pLOT+geom_tile(data = new_table,aes(x=reorder(ID,order),y=gene,fill=mmm),color='black',width=1,height=1,size=0.5,stat='identity')

F1b.pLOT<-F1b.pLOT+scale_fill_manual(name="Mutation Number",values =c(a0='white',b1='#fee0d2',c2='#fcbba1',c5='#fb6a4a',d10='#ef3b2c',e20='#cb181d',f30='#a50f15',g40='#67000d'),labels=c(a0='0',b1='1',c2='>=2',c5='>=5',d10='>=10',e20='>=20',f30='>=30',g40='>=40'))
F1b.pLOT<-F1b.pLOT+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(0.1,0.1,0.1,0.1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.68,'cm'),legend.position='none',axis.ticks.x = element_blank(),axis.line.x = element_blank(),axis.line.y = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=16,face='bold.italic'),axis.ticks.y = element_blank(), axis.text.y=element_blank(),
                         axis.text.x=element_text(size=12,angle=45,vjust=1,hjust=1,face='bold',color='black'),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F1b.pLOT<-F1b.pLOT+scale_x_discrete(position = "bottom")+xlab(NULL)+ylab(NULL)
F1b.pLOT
gene.scale1b = 1 + (ncol(mutation_gene_table) - 14)/14



yyy <- rep("Grade",length(sample$case))
data4<-data.frame(sample$case,yyy,sample$order,sample$Grade)
colnames(data4)<-c("xxx","yyy","zzz","his")

yyy <- rep("Gender",length(sample$case))
data5<-data.frame(sample$case,yyy,sample$order,sample$Gender)
colnames(data5)<-c("xxx","yyy","zzz","his")

yyy <- rep("Age",length(sample$case))
data6<-data.frame(sample$case,yyy,sample$order,sample$Age)
colnames(data6)<-c("xxx","yyy","zzz","his")

yyy <- rep("BCL2.fusion",length(sample$case))
data7b<-data.frame(sample$case,yyy,sample$order,sample$BCL2)
colnames(data7b)<-c("xxx","yyy","zzz","his")

yyy <- rep("MYC.fusion",length(sample$case))
data9<-data.frame(sample$case,yyy,sample$order,sample$MYC)
colnames(data9)<-c("xxx","yyy","zzz","his")

F4b.pLOT<-ggplot()+theme_classic()
F4b.pLOT<-F4b.pLOT+geom_tile(data = data4,aes(x=reorder(xxx,zzz),y=yyy,fill=his),color='black',width=1,height=1,size=0.5,stat='identity')
F4b.pLOT<-F4b.pLOT+scale_fill_manual(name=NULL,values=c(High=gg_color_hue(6)[1] ,Low=gg_color_hue(6)[3]),guide = guide_legend(override.aes=list(size=1),nrow=1),na.translate = F)
F4b.pLOT<-F4b.pLOT+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(1,0.1,0.1,0.1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.4,'cm'),legend.key.height=unit(0.4,'cm'),legend.position='none',axis.ticks.x = element_blank(),axis.line.x = element_blank(),axis.line.y = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=12,face='bold.italic'),axis.ticks.y = element_blank(), axis.text.y=element_blank(),legend.title=element_blank(),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F4b.pLOT<-F4b.pLOT+scale_x_discrete(position = "top")+xlab(NULL)+ylab(NULL)+guides(size=FALSE)
gene.scale4 = 1/14

F7b.pLOT<-ggplot()+theme_classic()
F7b.pLOT<-F7b.pLOT+geom_tile(data = data7b,aes(x=reorder(xxx,zzz),y=yyy,fill=his),color='black',width=1,height=1,size=0.5,stat='identity')
F7b.pLOT<-F7b.pLOT+scale_fill_manual(name=NULL,values=c(Fusion='black',WT='white'),
                                     labels=c(Fusion='Fusion',WT='Wild type'),guide = guide_legend(override.aes=list(size=1),nrow=1),na.translate = F)
F7b.pLOT<-F7b.pLOT+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(0.1,0.1,0.1,0.1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.4,'cm'),legend.key.height=unit(0.4,'cm'),legend.position='none',axis.ticks.x = element_blank(),axis.line.x = element_blank(),axis.line.y = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=12,face='bold.italic'),axis.ticks.y = element_blank(), axis.text.y=element_blank(),legend.title=element_blank(),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F7b.pLOT<-F7b.pLOT+scale_x_discrete(position = "top")+xlab(NULL)+ylab(NULL)+guides(size=FALSE)
gene.scale7 = 1/14



F5b.pLOT<-ggplot()+theme_classic()
F5b.pLOT<-F5b.pLOT+geom_tile(data = data5,aes(x=reorder(xxx,zzz),y=yyy,fill=his),color='black',width=1,height=1,size=0.5,stat='identity')
F5b.pLOT<-F5b.pLOT+scale_fill_manual(name=NULL,values=c(Male='#b3cde3',Female='#fccde5'),labels=c(Male='Male',Female='Female'),guide = guide_legend(override.aes=list(size=1),nrow=1),na.translate = F)
F5b.pLOT<-F5b.pLOT+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(0.1,0.1,0.1,0.1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.4,'cm'),legend.key.height=unit(0.4,'cm'),legend.position='none',axis.ticks.x = element_blank(),axis.line.x = element_blank(),axis.line.y = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=12,face='bold.italic'),axis.ticks.y = element_blank(), axis.text.y=element_blank(),legend.title=element_blank(),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F5b.pLOT<-F5b.pLOT+scale_x_discrete(position = "top")+xlab(NULL)+ylab(NULL)+guides(size=FALSE)
gene.scale5 = 1/14

F6b.pLOT<-ggplot()+theme_classic()
F6b.pLOT<-F6b.pLOT+geom_tile(data = data6,aes(x=reorder(xxx,zzz),y=yyy,fill=as.numeric(as.character(his))),color='black',width=1,height=1,size=0.5,stat='identity')
F6b.pLOT<-F6b.pLOT+scale_fill_gradientn(name=NULL,colours=c('white','#f6e8c3','#dfc27d','#bf812d','#8c510a'),limits=c(50,80),breaks=seq(50,80,5))
F6b.pLOT<-F6b.pLOT+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(0.1,0.1,0.1,0.1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(1.1,'cm'),legend.key.height=unit(0.4,'cm'),legend.position='none',axis.ticks.x = element_blank(),axis.line.x = element_blank(),axis.line.y = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=12,face='bold'),axis.ticks.y = element_blank(), axis.text.y=element_blank(),legend.title=element_blank(),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F6b.pLOT<-F6b.pLOT+scale_x_discrete(position = "top")+xlab(NULL)+ylab(NULL)+guides(size=FALSE)
gene.scale6 = 1/14




F9b.pLOT<-ggplot()+theme_classic()
F9b.pLOT<-F9b.pLOT+geom_tile(data = data9,aes(x=reorder(xxx,zzz),y=yyy,fill=his),color='black',width=1,height=1,size=0.5,stat='identity')
F9b.pLOT<-F9b.pLOT+scale_fill_manual(name=NULL,values=c(Fusion='black',WT='white'),
                                     labels=c(Fusion='Fusion',WT='Wild type'),guide = guide_legend(override.aes=list(size=1),nrow=1),na.translate = F)
F9b.pLOT<-F9b.pLOT+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(0.1,0.1,0.1,0.1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(1.1,'cm'),legend.key.height=unit(0.4,'cm'),legend.position='none',axis.ticks.x = element_blank(),axis.line.x = element_blank(),axis.line.y = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=12,face='bold'),axis.ticks.y = element_blank(), axis.text.y=element_blank(),legend.title=element_blank(),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F9b.pLOT<-F9b.pLOT+scale_x_discrete(position = "top")+xlab(NULL)+ylab(NULL)+guides(size=FALSE)
gene.scale11 = 1/14


###P4
sample <- data.frame(Sample_features)
sample <- sample[which(sample$patient=='P4'),]
pat <- unique(sample$patient)
case <- unique(sample$case)
#raw data preprocessing

somatic <- read.delim('./Mut/final/Annotated.frequent.txt',header = T)
somatic <- somatic[which(somatic$CaseID %in% pat),]
FL <- somatic[which(somatic$LG_freq >= 5 & somatic$altdepth_LG >= 2),]
FL$ID <- "P4.FL"
DHL <- somatic[which(somatic$HG_freq >= 5 & somatic$altdepth_HG >= 2),]
DHL$ID <- "P4.DHL"

savi.table <- rbind(FL,DHL)

savi.Sel <- savi.table[which(savi.table$ID %in% case),]

savi <- data.frame(savi.Sel$chr,savi.Sel$ID,savi.Sel$Gene_Name,savi.Sel$Effect_Impact,savi.Sel$Effect)
colnames(savi) <- c("chr","case","Gene_Name","Impact","Effect")

#analysis for selected genes
selGene <- interestedGene

gene.Matrix <- rep('N',length(case))
for( i in 2:length(selGene)){
  gene.Matrix <- cbind(gene.Matrix,rep('N',length(case)))
}

for(i in 1:length(case)){
  for(j in 1:length(selGene)){
    temp.P <- savi$Effect[which(savi$case == as.character(case[i]) & savi$Gene_Name == selGene[j] & savi$Impact!='HIGH' & savi$Impact!='MODERATE')]
    gene.Matrix[i,j] <- length(temp.P)
  }
}


colnames(gene.Matrix) <- selGene
rownames(gene.Matrix) <- case

mutation_gene_table <- gene.Matrix
#####


sample$order[which(sample$case=="P4.FL")] <- 1
sample$order[which(sample$case=="P4.DHL")] <-  2
sample$order[which(sample$case=="P4.FL")] <- 1

plot_1b.data <- melt(mutation_gene_table)
colnames(plot_1b.data) <- c('ID','gene','value')


xxx<-as.integer(as.character(plot_1b.data$value))
mmm<-as.integer(as.character(plot_1b.data$value))

mmm[which(xxx == 0)] <-   "a0"
mmm[which(xxx == 1) ] <-  "b1"
mmm[which(xxx >= 2) ] <-  "c2"
mmm[which(xxx >= 5) ] <-  "c5"
mmm[which(xxx >= 10) ] <- "d10"
mmm[which(xxx >= 20) ] <- "e20"
mmm[which(xxx >= 30) ] <- "f30"
mmm[which(xxx >= 40) ] <- "g40"
plot_1b.data <- cbind(plot_1b.data,mmm)
new_table <- merge(plot_1b.data,sample,1,1)

#orderID<-c(1:nrow(plot_1b.data))
F1b.PLOT<-ggplot()+theme_classic()
F1b.PLOT<-F1b.PLOT+geom_tile(data = new_table,aes(x=reorder(ID,order),y=gene,fill=mmm),color='black',width=1,height=1,size=0.5,stat='identity')

F1b.PLOT<-F1b.PLOT+scale_fill_manual(name="Mutation Number",values =c(a0='white',b1='#fee0d2',c2='#fcbba1',c5='#fb6a4a',d10='#ef3b2c',e20='#cb181d',f30='#a50f15',g40='#67000d'),labels=c(a0='0',b1='1',c2='>=2',c5='>=5',d10='>=10',e20='>=20',f30='>=30',g40='>=40'))
F1b.PLOT<-F1b.PLOT+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(0.1,0.1,0.1,0.1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.68,'cm'),legend.position='none',axis.ticks.x = element_blank(),axis.line.x = element_blank(),axis.line.y = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=16,face='bold.italic'),axis.ticks.y = element_blank(), axis.text.y=element_blank(),
                         axis.text.x=element_text(size=12,angle=45,vjust=1,hjust=1,face='bold',color='black'),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F1b.PLOT<-F1b.PLOT+scale_x_discrete(position = "bottom")+xlab(NULL)+ylab(NULL)
F1b.PLOT
gene.scale1b = 1 + (ncol(mutation_gene_table) - 14)/14



yyy <- rep("Grade",length(sample$case))
data4<-data.frame(sample$case,yyy,sample$order,sample$Grade)
colnames(data4)<-c("xxx","yyy","zzz","his")

yyy <- rep("Gender",length(sample$case))
data5<-data.frame(sample$case,yyy,sample$order,sample$Gender)
colnames(data5)<-c("xxx","yyy","zzz","his")

yyy <- rep("Age",length(sample$case))
data6<-data.frame(sample$case,yyy,sample$order,sample$Age)
colnames(data6)<-c("xxx","yyy","zzz","his")

yyy <- rep("BCL2.fusion",length(sample$case))
data7b<-data.frame(sample$case,yyy,sample$order,sample$BCL2)
colnames(data7b)<-c("xxx","yyy","zzz","his")

yyy <- rep("MYC.fusion",length(sample$case))
data9<-data.frame(sample$case,yyy,sample$order,sample$MYC)
colnames(data9)<-c("xxx","yyy","zzz","his")

F4b.PLOT<-ggplot()+theme_classic()
F4b.PLOT<-F4b.PLOT+geom_tile(data = data4,aes(x=reorder(xxx,zzz),y=yyy,fill=his),color='black',width=1,height=1,size=0.5,stat='identity')
F4b.PLOT<-F4b.PLOT+scale_fill_manual(name=NULL,values=c(High=gg_color_hue(6)[1] ,Low=gg_color_hue(6)[3]),guide = guide_legend(override.aes=list(size=1),nrow=1),na.translate = F)
F4b.PLOT<-F4b.PLOT+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(1,0.1,0.1,0.1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.4,'cm'),legend.key.height=unit(0.4,'cm'),legend.position='none',axis.ticks.x = element_blank(),axis.line.x = element_blank(),axis.line.y = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=12,face='bold.italic'),axis.ticks.y = element_blank(), axis.text.y=element_blank(),legend.title=element_blank(),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F4b.PLOT<-F4b.PLOT+scale_x_discrete(position = "top")+xlab(NULL)+ylab(NULL)+guides(size=FALSE)
gene.scale4 = 1/14

F7b.PLOT<-ggplot()+theme_classic()
F7b.PLOT<-F7b.PLOT+geom_tile(data = data7b,aes(x=reorder(xxx,zzz),y=yyy,fill=his),color='black',width=1,height=1,size=0.5,stat='identity')
F7b.PLOT<-F7b.PLOT+scale_fill_manual(name=NULL,values=c(Fusion='black',WT='white'),
                                     labels=c(Fusion='Fusion',WT='Wild type'),guide = guide_legend(override.aes=list(size=1),nrow=1),na.translate = F)
F7b.PLOT<-F7b.PLOT+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(0.1,0.1,0.1,0.1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.4,'cm'),legend.key.height=unit(0.4,'cm'),legend.position='none',axis.ticks.x = element_blank(),axis.line.x = element_blank(),axis.line.y = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=12,face='bold.italic'),axis.ticks.y = element_blank(), axis.text.y=element_blank(),legend.title=element_blank(),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F7b.PLOT<-F7b.PLOT+scale_x_discrete(position = "top")+xlab(NULL)+ylab(NULL)+guides(size=FALSE)
gene.scale7 = 1/14



F5b.PLOT<-ggplot()+theme_classic()
F5b.PLOT<-F5b.PLOT+geom_tile(data = data5,aes(x=reorder(xxx,zzz),y=yyy,fill=his),color='black',width=1,height=1,size=0.5,stat='identity')
F5b.PLOT<-F5b.PLOT+scale_fill_manual(name=NULL,values=c(Male='#b3cde3',Female='#fccde5'),labels=c(Male='Male',Female='Female'),guide = guide_legend(override.aes=list(size=1),nrow=1),na.translate = F)
F5b.PLOT<-F5b.PLOT+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(0.1,0.1,0.1,0.1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.4,'cm'),legend.key.height=unit(0.4,'cm'),legend.position='none',axis.ticks.x = element_blank(),axis.line.x = element_blank(),axis.line.y = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=12,face='bold.italic'),axis.ticks.y = element_blank(), axis.text.y=element_blank(),legend.title=element_blank(),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F5b.PLOT<-F5b.PLOT+scale_x_discrete(position = "top")+xlab(NULL)+ylab(NULL)+guides(size=FALSE)
gene.scale5 = 1/14

F6b.PLOT<-ggplot()+theme_classic()
F6b.PLOT<-F6b.PLOT+geom_tile(data = data6,aes(x=reorder(xxx,zzz),y=yyy,fill=as.numeric(as.character(his))),color='black',width=1,height=1,size=0.5,stat='identity')
F6b.PLOT<-F6b.PLOT+scale_fill_gradientn(name=NULL,colours=c('white','#f6e8c3','#dfc27d','#bf812d','#8c510a'),limits=c(50,80),breaks=seq(50,80,5))
F6b.PLOT<-F6b.PLOT+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(0.1,0.1,0.1,0.1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(1.1,'cm'),legend.key.height=unit(0.4,'cm'),legend.position='none',axis.ticks.x = element_blank(),axis.line.x = element_blank(),axis.line.y = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=12,face='bold'),axis.ticks.y = element_blank(), axis.text.y=element_blank(),legend.title=element_blank(),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F6b.PLOT<-F6b.PLOT+scale_x_discrete(position = "top")+xlab(NULL)+ylab(NULL)+guides(size=FALSE)
gene.scale6 = 1/14




F9b.PLOT<-ggplot()+theme_classic()
F9b.PLOT<-F9b.PLOT+geom_tile(data = data9,aes(x=reorder(xxx,zzz),y=yyy,fill=his),color='black',width=1,height=1,size=0.5,stat='identity')
F9b.PLOT<-F9b.PLOT+scale_fill_manual(name=NULL,values=c(Fusion='black',WT='white'),
                                     labels=c(Fusion='Fusion',WT='Wild type'),guide = guide_legend(override.aes=list(size=1),nrow=1),na.translate = F)
F9b.PLOT<-F9b.PLOT+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(0.1,0.1,0.1,0.1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(1.1,'cm'),legend.key.height=unit(0.4,'cm'),legend.position='none',axis.ticks.x = element_blank(),axis.line.x = element_blank(),axis.line.y = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=12,face='bold'),axis.ticks.y = element_blank(), axis.text.y=element_blank(),legend.title=element_blank(),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F9b.PLOT<-F9b.PLOT+scale_x_discrete(position = "top")+xlab(NULL)+ylab(NULL)+guides(size=FALSE)
gene.scale11 = 1/14

##P5
sample <- data.frame(Sample_features)
sample <- sample[which(sample$patient=='P5'),]
pat <- unique(sample$patient)
case <- unique(sample$case)
#raw data preprocessing

somatic <- read.delim('./Mut/final/Annotated.frequent.txt',header = T)
somatic <- somatic[which(somatic$CaseID %in% pat),]
FL <- somatic[which(somatic$LG_freq >= 5 & somatic$altdepth_LG >= 2),]
FL$ID <- "P5.FL"
DHL <- somatic[which(somatic$HG_freq >= 5 & somatic$altdepth_HG >= 2),]
DHL$ID <- "P5.DHL"

savi.table <- rbind(FL,DHL)

savi.Sel <- savi.table[which(savi.table$ID %in% case),]

savi <- data.frame(savi.Sel$chr,savi.Sel$ID,savi.Sel$Gene_Name,savi.Sel$Effect_Impact,savi.Sel$Effect)
colnames(savi) <- c("chr","case","Gene_Name","Impact","Effect")

#analysis for selected genes
selGene <- interestedGene

gene.Matrix <- rep('N',length(case))
for( i in 2:length(selGene)){
  gene.Matrix <- cbind(gene.Matrix,rep('N',length(case)))
}

for(i in 1:length(case)){
  for(j in 1:length(selGene)){
    temp.P <- savi$Effect[which(savi$case == as.character(case[i]) & savi$Gene_Name == selGene[j] & savi$Impact!='HIGH' & savi$Impact!='MODERATE')]
    gene.Matrix[i,j] <- length(temp.P)
  }
}


colnames(gene.Matrix) <- selGene
rownames(gene.Matrix) <- case

mutation_gene_table <- gene.Matrix
#####


sample$order[which(sample$case=="P5.FL")] <- 1
sample$order[which(sample$case=="P5.DHL")] <-  2
sample$order[which(sample$case=="P5.FL")] <- 1

plot_1b.data <- melt(mutation_gene_table)
colnames(plot_1b.data) <- c('ID','gene','value')


xxx<-as.integer(as.character(plot_1b.data$value))
mmm<-as.integer(as.character(plot_1b.data$value))

mmm[which(xxx == 0)] <-   "a0"
mmm[which(xxx == 1) ] <-  "b1"
mmm[which(xxx >= 2) ] <-  "c2"
mmm[which(xxx >= 5) ] <-  "c5"
mmm[which(xxx >= 10) ] <- "d10"
mmm[which(xxx >= 20) ] <- "e20"
mmm[which(xxx >= 30) ] <- "f30"
mmm[which(xxx >= 40) ] <- "g40"
plot_1b.data <- cbind(plot_1b.data,mmm)
new_table <- merge(plot_1b.data,sample,1,1)

#orderID<-c(1:nrow(plot_1b.data))
F1b.PLOt<-ggplot()+theme_classic()
F1b.PLOt<-F1b.PLOt+geom_tile(data = new_table,aes(x=reorder(ID,order),y=gene,fill=mmm),color='black',width=1,height=1,size=0.5,stat='identity')

F1b.PLOt<-F1b.PLOt+scale_fill_manual(name="Mutation Number",values =c(a0='white',b1='#fee0d2',c2='#fcbba1',c5='#fb6a4a',d10='#ef3b2c',e20='#cb181d',f30='#a50f15',g40='#67000d'),labels=c(a0='0',b1='1',c2='>=2',c5='>=5',d10='>=10',e20='>=20',f30='>=30',g40='>=40'))
F1b.PLOt<-F1b.PLOt+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(0.1,0.1,0.1,0.1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.68,'cm'),legend.position='none',axis.ticks.x = element_blank(),axis.line.x = element_blank(),axis.line.y = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=16,face='bold.italic'),axis.ticks.y = element_blank(), axis.text.y=element_blank(),
                         axis.text.x=element_text(size=12,angle=45,vjust=1,hjust=1,face='bold',color='black'),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F1b.PLOt<-F1b.PLOt+scale_x_discrete(position = "bottom")+xlab(NULL)+ylab(NULL)
F1b.PLOt
gene.scale1b = 1 + (ncol(mutation_gene_table) - 14)/14



yyy <- rep("Grade",length(sample$case))
data4<-data.frame(sample$case,yyy,sample$order,sample$Grade)
colnames(data4)<-c("xxx","yyy","zzz","his")

yyy <- rep("Gender",length(sample$case))
data5<-data.frame(sample$case,yyy,sample$order,sample$Gender)
colnames(data5)<-c("xxx","yyy","zzz","his")

yyy <- rep("Age",length(sample$case))
data6<-data.frame(sample$case,yyy,sample$order,sample$Age)
colnames(data6)<-c("xxx","yyy","zzz","his")

yyy <- rep("BCL2.fusion",length(sample$case))
data7b<-data.frame(sample$case,yyy,sample$order,sample$BCL2)
colnames(data7b)<-c("xxx","yyy","zzz","his")

yyy <- rep("MYC.fusion",length(sample$case))
data9<-data.frame(sample$case,yyy,sample$order,sample$MYC)
colnames(data9)<-c("xxx","yyy","zzz","his")

F4b.PLOt<-ggplot()+theme_classic()
F4b.PLOt<-F4b.PLOt+geom_tile(data = data4,aes(x=reorder(xxx,zzz),y=yyy,fill=his),color='black',width=1,height=1,size=0.5,stat='identity')
F4b.PLOt<-F4b.PLOt+scale_fill_manual(name=NULL,values=c(High=gg_color_hue(6)[1] ,Low=gg_color_hue(6)[3]),guide = guide_legend(override.aes=list(size=1),nrow=1),na.translate = F)
F4b.PLOt<-F4b.PLOt+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(1,0.1,0.1,0.1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.4,'cm'),legend.key.height=unit(0.4,'cm'),legend.position='none',axis.ticks.x = element_blank(),axis.line.x = element_blank(),axis.line.y = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=12,face='bold.italic'),axis.ticks.y = element_blank(), axis.text.y=element_blank(),legend.title=element_blank(),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F4b.PLOt<-F4b.PLOt+scale_x_discrete(position = "top")+xlab(NULL)+ylab(NULL)+guides(size=FALSE)
gene.scale4 = 1/14

F7b.PLOt<-ggplot()+theme_classic()
F7b.PLOt<-F7b.PLOt+geom_tile(data = data7b,aes(x=reorder(xxx,zzz),y=yyy,fill=his),color='black',width=1,height=1,size=0.5,stat='identity')
F7b.PLOt<-F7b.PLOt+scale_fill_manual(name=NULL,values=c(Fusion='black',WT='white'),
                                     labels=c(Fusion='Fusion',WT='Wild type'),guide = guide_legend(override.aes=list(size=1),nrow=1),na.translate = F)
F7b.PLOt<-F7b.PLOt+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(0.1,0.1,0.1,0.1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.4,'cm'),legend.key.height=unit(0.4,'cm'),legend.position='none',axis.ticks.x = element_blank(),axis.line.x = element_blank(),axis.line.y = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=12,face='bold.italic'),axis.ticks.y = element_blank(), axis.text.y=element_blank(),legend.title=element_blank(),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F7b.PLOt<-F7b.PLOt+scale_x_discrete(position = "top")+xlab(NULL)+ylab(NULL)+guides(size=FALSE)
gene.scale7 = 1/14



F5b.PLOt<-ggplot()+theme_classic()
F5b.PLOt<-F5b.PLOt+geom_tile(data = data5,aes(x=reorder(xxx,zzz),y=yyy,fill=his),color='black',width=1,height=1,size=0.5,stat='identity')
F5b.PLOt<-F5b.PLOt+scale_fill_manual(name=NULL,values=c(Male='#b3cde3',Female='#fccde5'),labels=c(Male='Male',Female='Female'),guide = guide_legend(override.aes=list(size=1),nrow=1),na.translate = F)
F5b.PLOt<-F5b.PLOt+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(0.1,0.1,0.1,0.1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.4,'cm'),legend.key.height=unit(0.4,'cm'),legend.position='none',axis.ticks.x = element_blank(),axis.line.x = element_blank(),axis.line.y = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=12,face='bold.italic'),axis.ticks.y = element_blank(), axis.text.y=element_blank(),legend.title=element_blank(),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F5b.PLOt<-F5b.PLOt+scale_x_discrete(position = "top")+xlab(NULL)+ylab(NULL)+guides(size=FALSE)
gene.scale5 = 1/14

F6b.PLOt<-ggplot()+theme_classic()
F6b.PLOt<-F6b.PLOt+geom_tile(data = data6,aes(x=reorder(xxx,zzz),y=yyy,fill=as.numeric(as.character(his))),color='black',width=1,height=1,size=0.5,stat='identity')
F6b.PLOt<-F6b.PLOt+scale_fill_gradientn(name=NULL,colours=c('white','#f6e8c3','#dfc27d','#bf812d','#8c510a'),limits=c(50,80),breaks=seq(50,80,5))
F6b.PLOt<-F6b.PLOt+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(0.1,0.1,0.1,0.1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(1.1,'cm'),legend.key.height=unit(0.4,'cm'),legend.position='none',axis.ticks.x = element_blank(),axis.line.x = element_blank(),axis.line.y = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=12,face='bold'),axis.ticks.y = element_blank(), axis.text.y=element_blank(),legend.title=element_blank(),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F6b.PLOt<-F6b.PLOt+scale_x_discrete(position = "top")+xlab(NULL)+ylab(NULL)+guides(size=FALSE)
gene.scale6 = 1/14




F9b.PLOt<-ggplot()+theme_classic()
F9b.PLOt<-F9b.PLOt+geom_tile(data = data9,aes(x=reorder(xxx,zzz),y=yyy,fill=his),color='black',width=1,height=1,size=0.5,stat='identity')
F9b.PLOt<-F9b.PLOt+scale_fill_manual(name=NULL,values=c(Fusion='black',WT='white'),
                                     labels=c(Fusion='Fusion',WT='Wild type'),guide = guide_legend(override.aes=list(size=1),nrow=1),na.translate = F)
F9b.PLOt<-F9b.PLOt+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(0.1,0.1,0.1,0.1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(1.1,'cm'),legend.key.height=unit(0.4,'cm'),legend.position='none',axis.ticks.x = element_blank(),axis.line.x = element_blank(),axis.line.y = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=12,face='bold'),axis.ticks.y = element_blank(), axis.text.y=element_blank(),legend.title=element_blank(),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F9b.PLOt<-F9b.PLOt+scale_x_discrete(position = "top")+xlab(NULL)+ylab(NULL)+guides(size=FALSE)
gene.scale11 = 1/14


##P6
sample <- data.frame(Sample_features)
sample <- sample[which(sample$patient=='P6'),]
pat <- unique(sample$patient)
case <- unique(sample$case)
#raw data preprocessing

somatic <- read.delim('./Mut/final/Annotated.frequent.txt',header = T)
somatic <- somatic[which(somatic$CaseID %in% pat),]
FL <- somatic[which(somatic$LG_freq >= 5 & somatic$altdepth_LG >= 2),]
FL$ID <- "P6.FL"
DHL <- somatic[which(somatic$HG_freq >= 5 & somatic$altdepth_HG >= 2),]
DHL$ID <- "P6.DHL"

savi.table <- rbind(FL,DHL)

savi.Sel <- savi.table[which(savi.table$ID %in% case),]

savi <- data.frame(savi.Sel$chr,savi.Sel$ID,savi.Sel$Gene_Name,savi.Sel$Effect_Impact,savi.Sel$Effect)
colnames(savi) <- c("chr","case","Gene_Name","Impact","Effect")

#analysis for selected genes
selGene <- interestedGene

gene.Matrix <- rep('N',length(case))
for( i in 2:length(selGene)){
  gene.Matrix <- cbind(gene.Matrix,rep('N',length(case)))
}

for(i in 1:length(case)){
  for(j in 1:length(selGene)){
    temp.P <- savi$Effect[which(savi$case == as.character(case[i]) & savi$Gene_Name == selGene[j] & savi$Impact!='HIGH' & savi$Impact!='MODERATE')]
    gene.Matrix[i,j] <- length(temp.P)
  }
}


colnames(gene.Matrix) <- selGene
rownames(gene.Matrix) <- case

mutation_gene_table <- gene.Matrix
#####


sample$order[which(sample$case=="P6.FL")] <- 1
sample$order[which(sample$case=="P6.DHL")] <-  2
sample$order[which(sample$case=="P6.FL")] <- 1

plot_1b.data <- melt(mutation_gene_table)
colnames(plot_1b.data) <- c('ID','gene','value')


xxx<-as.integer(as.character(plot_1b.data$value))
mmm<-as.integer(as.character(plot_1b.data$value))

mmm[which(xxx == 0)] <-   "a0"
mmm[which(xxx == 1) ] <-  "b1"
mmm[which(xxx >= 2) ] <-  "c2"
mmm[which(xxx >= 5) ] <-  "c5"
mmm[which(xxx >= 10) ] <- "d10"
mmm[which(xxx >= 20) ] <- "e20"
mmm[which(xxx >= 30) ] <- "f30"
mmm[which(xxx >= 40) ] <- "g40"
plot_1b.data <- cbind(plot_1b.data,mmm)
new_table <- merge(plot_1b.data,sample,1,1)

#orderID<-c(1:nrow(plot_1b.data))
F1b.PLoT<-ggplot()+theme_classic()
F1b.PLoT<-F1b.PLoT+geom_tile(data = new_table,aes(x=reorder(ID,order),y=gene,fill=mmm),color='black',width=1,height=1,size=0.5,stat='identity')

F1b.PLoT<-F1b.PLoT+scale_fill_manual(name="Mutation Number",values =c(a0='white',b1='#fee0d2',c2='#fcbba1',c5='#fb6a4a',d10='#ef3b2c',e20='#cb181d',f30='#a50f15',g40='#67000d'),labels=c(a0='0',b1='1',c2='>=2',c5='>=5',d10='>=10',e20='>=20',f30='>=30',g40='>=40'))
F1b.PLoT<-F1b.PLoT+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(0.1,0.1,0.1,0.1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.68,'cm'),legend.position='none',axis.ticks.x = element_blank(),axis.line.x = element_blank(),axis.line.y = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=16,face='bold.italic'),axis.ticks.y = element_blank(), axis.text.y=element_blank(),
                         axis.text.x=element_text(size=12,angle=45,vjust=1,hjust=1,face='bold',color='black'),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F1b.PLoT<-F1b.PLoT+scale_x_discrete(position = "bottom")+xlab(NULL)+ylab(NULL)
F1b.PLoT
gene.scale1b = 1 + (ncol(mutation_gene_table) - 14)/14



yyy <- rep("Grade",length(sample$case))
data4<-data.frame(sample$case,yyy,sample$order,sample$Grade)
colnames(data4)<-c("xxx","yyy","zzz","his")

yyy <- rep("Gender",length(sample$case))
data5<-data.frame(sample$case,yyy,sample$order,sample$Gender)
colnames(data5)<-c("xxx","yyy","zzz","his")

yyy <- rep("Age",length(sample$case))
data6<-data.frame(sample$case,yyy,sample$order,sample$Age)
colnames(data6)<-c("xxx","yyy","zzz","his")

yyy <- rep("BCL2.fusion",length(sample$case))
data7b<-data.frame(sample$case,yyy,sample$order,sample$BCL2)
colnames(data7b)<-c("xxx","yyy","zzz","his")

yyy <- rep("MYC.fusion",length(sample$case))
data9<-data.frame(sample$case,yyy,sample$order,sample$MYC)
colnames(data9)<-c("xxx","yyy","zzz","his")

F4b.PLoT<-ggplot()+theme_classic()
F4b.PLoT<-F4b.PLoT+geom_tile(data = data4,aes(x=reorder(xxx,zzz),y=yyy,fill=his),color='black',width=1,height=1,size=0.5,stat='identity')
F4b.PLoT<-F4b.PLoT+scale_fill_manual(name=NULL,values=c(High=gg_color_hue(6)[1] ,Low=gg_color_hue(6)[3]),guide = guide_legend(override.aes=list(size=1),nrow=1),na.translate = F)
F4b.PLoT<-F4b.PLoT+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(1,0.1,0.1,0.1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.4,'cm'),legend.key.height=unit(0.4,'cm'),legend.position='none',axis.ticks.x = element_blank(),axis.line.x = element_blank(),axis.line.y = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=12,face='bold.italic'),axis.ticks.y = element_blank(), axis.text.y=element_blank(),legend.title=element_blank(),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F4b.PLoT<-F4b.PLoT+scale_x_discrete(position = "top")+xlab(NULL)+ylab(NULL)+guides(size=FALSE)
gene.scale4 = 1/14

F7b.PLoT<-ggplot()+theme_classic()
F7b.PLoT<-F7b.PLoT+geom_tile(data = data7b,aes(x=reorder(xxx,zzz),y=yyy,fill=his),color='black',width=1,height=1,size=0.5,stat='identity')
F7b.PLoT<-F7b.PLoT+scale_fill_manual(name=NULL,values=c(Fusion='black',WT='white'),
                                     labels=c(Fusion='Fusion',WT='Wild type'),guide = guide_legend(override.aes=list(size=1),nrow=1),na.translate = F)
F7b.PLoT<-F7b.PLoT+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(0.1,0.1,0.1,0.1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.4,'cm'),legend.key.height=unit(0.4,'cm'),legend.position='none',axis.ticks.x = element_blank(),axis.line.x = element_blank(),axis.line.y = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=12,face='bold.italic'),axis.ticks.y = element_blank(), axis.text.y=element_blank(),legend.title=element_blank(),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F7b.PLoT<-F7b.PLoT+scale_x_discrete(position = "top")+xlab(NULL)+ylab(NULL)+guides(size=FALSE)
gene.scale7 = 1/14



F5b.PLoT<-ggplot()+theme_classic()
F5b.PLoT<-F5b.PLoT+geom_tile(data = data5,aes(x=reorder(xxx,zzz),y=yyy,fill=his),color='black',width=1,height=1,size=0.5,stat='identity')
F5b.PLoT<-F5b.PLoT+scale_fill_manual(name=NULL,values=c(Male='#b3cde3',Female='#fccde5'),labels=c(Male='Male',Female='Female'),guide = guide_legend(override.aes=list(size=1),nrow=1),na.translate = F)
F5b.PLoT<-F5b.PLoT+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(0.1,0.1,0.1,0.1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.4,'cm'),legend.key.height=unit(0.4,'cm'),legend.position='none',axis.ticks.x = element_blank(),axis.line.x = element_blank(),axis.line.y = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=12,face='bold.italic'),axis.ticks.y = element_blank(), axis.text.y=element_blank(),legend.title=element_blank(),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F5b.PLoT<-F5b.PLoT+scale_x_discrete(position = "top")+xlab(NULL)+ylab(NULL)+guides(size=FALSE)
gene.scale5 = 1/14

F6b.PLoT<-ggplot()+theme_classic()
F6b.PLoT<-F6b.PLoT+geom_tile(data = data6,aes(x=reorder(xxx,zzz),y=yyy,fill=as.numeric(as.character(his))),color='black',width=1,height=1,size=0.5,stat='identity')
F6b.PLoT<-F6b.PLoT+scale_fill_gradientn(name=NULL,colours=c('white','#f6e8c3','#dfc27d','#bf812d','#8c510a'),limits=c(50,80),breaks=seq(50,80,5))
F6b.PLoT<-F6b.PLoT+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(0.1,0.1,0.1,0.1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(1.1,'cm'),legend.key.height=unit(0.4,'cm'),legend.position='none',axis.ticks.x = element_blank(),axis.line.x = element_blank(),axis.line.y = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=12,face='bold'),axis.ticks.y = element_blank(), axis.text.y=element_blank(),legend.title=element_blank(),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F6b.PLoT<-F6b.PLoT+scale_x_discrete(position = "top")+xlab(NULL)+ylab(NULL)+guides(size=FALSE)
gene.scale6 = 1/14




F9b.PLoT<-ggplot()+theme_classic()
F9b.PLoT<-F9b.PLoT+geom_tile(data = data9,aes(x=reorder(xxx,zzz),y=yyy,fill=his),color='black',width=1,height=1,size=0.5,stat='identity')
F9b.PLoT<-F9b.PLoT+scale_fill_manual(name=NULL,values=c(Fusion='black',WT='white'),
                                     labels=c(Fusion='Fusion',WT='Wild type'),guide = guide_legend(override.aes=list(size=1),nrow=1),na.translate = F)
F9b.PLoT<-F9b.PLoT+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(0.1,0.1,0.1,0.1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(1.1,'cm'),legend.key.height=unit(0.4,'cm'),legend.position='none',axis.ticks.x = element_blank(),axis.line.x = element_blank(),axis.line.y = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=12,face='bold'),axis.ticks.y = element_blank(), axis.text.y=element_blank(),legend.title=element_blank(),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F9b.PLoT<-F9b.PLoT+scale_x_discrete(position = "top")+xlab(NULL)+ylab(NULL)+guides(size=FALSE)
gene.scale11 = 1/14


###P7
sample <- data.frame(Sample_features)
sample <- sample[which(sample$patient=='P7'),]
pat <- unique(sample$patient)
case <- unique(sample$case)
#raw data preprocessing

somatic <- read.delim('./Mut/final/Annotated.frequent.txt',header = T)
somatic <- somatic[which(somatic$CaseID %in% pat),]
FL <- somatic[which(somatic$LG_freq >= 5 & somatic$altdepth_LG >= 2),]
FL$ID <- "P7.FL"
DHL <- somatic[which(somatic$HG_freq >= 5 & somatic$altdepth_HG >= 2),]
DHL$ID <- "P7.DHL"

savi.table <- rbind(FL,DHL)

savi.Sel <- savi.table[which(savi.table$ID %in% case),]

savi <- data.frame(savi.Sel$chr,savi.Sel$ID,savi.Sel$Gene_Name,savi.Sel$Effect_Impact,savi.Sel$Effect)
colnames(savi) <- c("chr","case","Gene_Name","Impact","Effect")

#analysis for selected genes
selGene <- interestedGene

gene.Matrix <- rep('N',length(case))
for( i in 2:length(selGene)){
  gene.Matrix <- cbind(gene.Matrix,rep('N',length(case)))
}

for(i in 1:length(case)){
  for(j in 1:length(selGene)){
    temp.P <- savi$Effect[which(savi$case == as.character(case[i]) & savi$Gene_Name == selGene[j] & savi$Impact!='HIGH' & savi$Impact!='MODERATE')]
    gene.Matrix[i,j] <- length(temp.P)
  }
}


colnames(gene.Matrix) <- selGene
rownames(gene.Matrix) <- case

mutation_gene_table <- gene.Matrix
#####


sample$order[which(sample$case=="P7.FL")] <- 1
sample$order[which(sample$case=="P7.DHL")] <-  2
sample$order[which(sample$case=="P7.FL")] <- 1

plot_1b.data <- melt(mutation_gene_table)
colnames(plot_1b.data) <- c('ID','gene','value')


xxx<-as.integer(as.character(plot_1b.data$value))
mmm<-as.integer(as.character(plot_1b.data$value))

mmm[which(xxx == 0)] <-   "a0"
mmm[which(xxx == 1) ] <-  "b1"
mmm[which(xxx >= 2) ] <-  "c2"
mmm[which(xxx >= 5) ] <-  "c5"
mmm[which(xxx >= 10) ] <- "d10"
mmm[which(xxx >= 20) ] <- "e20"
mmm[which(xxx >= 30) ] <- "f30"
mmm[which(xxx >= 40) ] <- "g40"
plot_1b.data <- cbind(plot_1b.data,mmm)
new_table <- merge(plot_1b.data,sample,1,1)

#orderID<-c(1:nrow(plot_1b.data))
F1b.PlOT<-ggplot()+theme_classic()
F1b.PlOT<-F1b.PlOT+geom_tile(data = new_table,aes(x=reorder(ID,order),y=gene,fill=mmm),color='black',width=1,height=1,size=0.5,stat='identity')

F1b.PlOT<-F1b.PlOT+scale_fill_manual(name="Mutation Number",values =c(a0='white',b1='#fee0d2',c2='#fcbba1',c5='#fb6a4a',d10='#ef3b2c',e20='#cb181d',f30='#a50f15',g40='#67000d'),labels=c(a0='0',b1='1',c2='>=2',c5='>=5',d10='>=10',e20='>=20',f30='>=30',g40='>=40'))
F1b.PlOT<-F1b.PlOT+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(0.1,0.1,0.1,0.1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.68,'cm'),legend.position='none',axis.ticks.x = element_blank(),axis.line.x = element_blank(),axis.line.y = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=16,face='bold.italic'),axis.ticks.y = element_blank(), axis.text.y=element_blank(),
                         axis.text.x=element_text(size=12,angle=45,vjust=1,hjust=1,face='bold',color='black'),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F1b.PlOT<-F1b.PlOT+scale_x_discrete(position = "bottom")+xlab(NULL)+ylab(NULL)
F1b.PlOT
gene.scale1b = 1 + (ncol(mutation_gene_table) - 14)/14



yyy <- rep("Grade",length(sample$case))
data4<-data.frame(sample$case,yyy,sample$order,sample$Grade)
colnames(data4)<-c("xxx","yyy","zzz","his")

yyy <- rep("Gender",length(sample$case))
data5<-data.frame(sample$case,yyy,sample$order,sample$Gender)
colnames(data5)<-c("xxx","yyy","zzz","his")

yyy <- rep("Age",length(sample$case))
data6<-data.frame(sample$case,yyy,sample$order,sample$Age)
colnames(data6)<-c("xxx","yyy","zzz","his")

yyy <- rep("BCL2.fusion",length(sample$case))
data7b<-data.frame(sample$case,yyy,sample$order,sample$BCL2)
colnames(data7b)<-c("xxx","yyy","zzz","his")

yyy <- rep("MYC.fusion",length(sample$case))
data9<-data.frame(sample$case,yyy,sample$order,sample$MYC)
colnames(data9)<-c("xxx","yyy","zzz","his")

F4b.PlOT<-ggplot()+theme_classic()
F4b.PlOT<-F4b.PlOT+geom_tile(data = data4,aes(x=reorder(xxx,zzz),y=yyy,fill=his),color='black',width=1,height=1,size=0.5,stat='identity')
F4b.PlOT<-F4b.PlOT+scale_fill_manual(name=NULL,values=c(High=gg_color_hue(6)[1] ,Low=gg_color_hue(6)[3]),guide = guide_legend(override.aes=list(size=1),nrow=1),na.translate = F)
F4b.PlOT<-F4b.PlOT+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(1,0.1,0.1,0.1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.4,'cm'),legend.key.height=unit(0.4,'cm'),legend.position='none',axis.ticks.x = element_blank(),axis.line.x = element_blank(),axis.line.y = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=12,face='bold.italic'),axis.ticks.y = element_blank(), axis.text.y=element_blank(),legend.title=element_blank(),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F4b.PlOT<-F4b.PlOT+scale_x_discrete(position = "top")+xlab(NULL)+ylab(NULL)+guides(size=FALSE)
gene.scale4 = 1/14

F7b.PlOT<-ggplot()+theme_classic()
F7b.PlOT<-F7b.PlOT+geom_tile(data = data7b,aes(x=reorder(xxx,zzz),y=yyy,fill=his),color='black',width=1,height=1,size=0.5,stat='identity')
F7b.PlOT<-F7b.PlOT+scale_fill_manual(name=NULL,values=c(Fusion='black',WT='white'),
                                     labels=c(Fusion='Fusion',WT='Wild type'),guide = guide_legend(override.aes=list(size=1),nrow=1),na.translate = F)
F7b.PlOT<-F7b.PlOT+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(0.1,0.1,0.1,0.1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.4,'cm'),legend.key.height=unit(0.4,'cm'),legend.position='none',axis.ticks.x = element_blank(),axis.line.x = element_blank(),axis.line.y = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=12,face='bold.italic'),axis.ticks.y = element_blank(), axis.text.y=element_blank(),legend.title=element_blank(),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F7b.PlOT<-F7b.PlOT+scale_x_discrete(position = "top")+xlab(NULL)+ylab(NULL)+guides(size=FALSE)
gene.scale7 = 1/14



F5b.PlOT<-ggplot()+theme_classic()
F5b.PlOT<-F5b.PlOT+geom_tile(data = data5,aes(x=reorder(xxx,zzz),y=yyy,fill=his),color='black',width=1,height=1,size=0.5,stat='identity')
F5b.PlOT<-F5b.PlOT+scale_fill_manual(name=NULL,values=c(Male='#b3cde3',Female='#fccde5'),labels=c(Male='Male',Female='Female'),guide = guide_legend(override.aes=list(size=1),nrow=1),na.translate = F)
F5b.PlOT<-F5b.PlOT+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(0.1,0.1,0.1,0.1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.4,'cm'),legend.key.height=unit(0.4,'cm'),legend.position='none',axis.ticks.x = element_blank(),axis.line.x = element_blank(),axis.line.y = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=12,face='bold.italic'),axis.ticks.y = element_blank(), axis.text.y=element_blank(),legend.title=element_blank(),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F5b.PlOT<-F5b.PlOT+scale_x_discrete(position = "top")+xlab(NULL)+ylab(NULL)+guides(size=FALSE)
gene.scale5 = 1/14

F6b.PlOT<-ggplot()+theme_classic()
F6b.PlOT<-F6b.PlOT+geom_tile(data = data6,aes(x=reorder(xxx,zzz),y=yyy,fill=as.numeric(as.character(his))),color='black',width=1,height=1,size=0.5,stat='identity')
F6b.PlOT<-F6b.PlOT+scale_fill_gradientn(name=NULL,colours=c('white','#f6e8c3','#dfc27d','#bf812d','#8c510a'),limits=c(50,80),breaks=seq(50,80,5))
F6b.PlOT<-F6b.PlOT+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(0.1,0.1,0.1,0.1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(1.1,'cm'),legend.key.height=unit(0.4,'cm'),legend.position='none',axis.ticks.x = element_blank(),axis.line.x = element_blank(),axis.line.y = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=12,face='bold'),axis.ticks.y = element_blank(), axis.text.y=element_blank(),legend.title=element_blank(),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F6b.PlOT<-F6b.PlOT+scale_x_discrete(position = "top")+xlab(NULL)+ylab(NULL)+guides(size=FALSE)
gene.scale6 = 1/14




F9b.PlOT<-ggplot()+theme_classic()
F9b.PlOT<-F9b.PlOT+geom_tile(data = data9,aes(x=reorder(xxx,zzz),y=yyy,fill=his),color='black',width=1,height=1,size=0.5,stat='identity')
F9b.PlOT<-F9b.PlOT+scale_fill_manual(name=NULL,values=c(Fusion='black',WT='white'),
                                     labels=c(Fusion='Fusion',WT='Wild type'),guide = guide_legend(override.aes=list(size=1),nrow=1),na.translate = F)
F9b.PlOT<-F9b.PlOT+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(0.1,0.1,0.1,0.1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(1.1,'cm'),legend.key.height=unit(0.4,'cm'),legend.position='none',axis.ticks.x = element_blank(),axis.line.x = element_blank(),axis.line.y = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=12,face='bold'),axis.ticks.y = element_blank(), axis.text.y=element_blank(),legend.title=element_blank(),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F9b.PlOT<-F9b.PlOT+scale_x_discrete(position = "top")+xlab(NULL)+ylab(NULL)+guides(size=FALSE)
gene.scale11 = 1/14


#P9
sample <- data.frame(Sample_features)
sample <- sample[which(sample$patient=='P9'),]
pat <- unique(sample$patient)
case <- unique(sample$case)
#raw data preprocessing

somatic <- read.delim('./Mut/final/Annotated.frequent.txt',header = T)
somatic <- somatic[which(somatic$CaseID %in% pat),]
FL <- somatic[which(somatic$LG_freq >= 5 & somatic$altdepth_LG >= 2),]
FL$ID <- "P9.FL"
DHL <- somatic[which(somatic$HG_freq >= 5 & somatic$altdepth_HG >= 2),]
DHL$ID <- "P9.DHL"

savi.table <- rbind(FL,DHL)

savi.Sel <- savi.table[which(savi.table$ID %in% case),]

savi <- data.frame(savi.Sel$chr,savi.Sel$ID,savi.Sel$Gene_Name,savi.Sel$Effect_Impact,savi.Sel$Effect)
colnames(savi) <- c("chr","case","Gene_Name","Impact","Effect")

#analysis for selected genes
selGene <- interestedGene

gene.Matrix <- rep('N',length(case))
for( i in 2:length(selGene)){
  gene.Matrix <- cbind(gene.Matrix,rep('N',length(case)))
}

for(i in 1:length(case)){
  for(j in 1:length(selGene)){
    temp.P <- savi$Effect[which(savi$case == as.character(case[i]) & savi$Gene_Name == selGene[j] & savi$Impact!='HIGH' & savi$Impact!='MODERATE')]
    gene.Matrix[i,j] <- length(temp.P)
  }
}


colnames(gene.Matrix) <- selGene
rownames(gene.Matrix) <- case

mutation_gene_table <- gene.Matrix
#####


sample$order[which(sample$case=="P9.FL")] <- 1
sample$order[which(sample$case=="P9.DHL")] <-  2
sample$order[which(sample$case=="P9.FL")] <- 1

plot_1b.data <- melt(mutation_gene_table)
colnames(plot_1b.data) <- c('ID','gene','value')


xxx<-as.integer(as.character(plot_1b.data$value))
mmm<-as.integer(as.character(plot_1b.data$value))

mmm[which(xxx == 0)] <-   "a0"
mmm[which(xxx == 1) ] <-  "b1"
mmm[which(xxx >= 2) ] <-  "c2"
mmm[which(xxx >= 5) ] <-  "c5"
mmm[which(xxx >= 10) ] <- "d10"
mmm[which(xxx >= 20) ] <- "e20"
mmm[which(xxx >= 30) ] <- "f30"
mmm[which(xxx >= 40) ] <- "g40"
plot_1b.data <- cbind(plot_1b.data,mmm)
new_table <- merge(plot_1b.data,sample,1,1)

#orderID<-c(1:nrow(plot_1b.data))
F1b.PloT<-ggplot()+theme_classic()
F1b.PloT<-F1b.PloT+geom_tile(data = new_table,aes(x=reorder(ID,order),y=gene,fill=mmm),color='black',width=1,height=1,size=0.5,stat='identity')

F1b.PloT<-F1b.PloT+scale_fill_manual(name="Mutation Number",values =c(a0='white',b1='#fee0d2',c2='#fcbba1',c5='#fb6a4a',d10='#ef3b2c',e20='#cb181d',f30='#a50f15',g40='#67000d'),labels=c(a0='0',b1='1',c2='>=2',c5='>=5',d10='>=10',e20='>=20',f30='>=30',g40='>=40'))
F1b.PloT<-F1b.PloT+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(0.1,0.1,0.1,0.1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.68,'cm'),legend.position='none',axis.ticks.x = element_blank(),axis.line.x = element_blank(),axis.line.y = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=16,face='bold.italic'),axis.ticks.y = element_blank(), axis.text.y=element_blank(),
                         axis.text.x=element_text(size=12,angle=45,vjust=1,hjust=1,face='bold',color='black'),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F1b.PloT<-F1b.PloT+scale_x_discrete(position = "bottom")+xlab(NULL)+ylab(NULL)
F1b.PloT
gene.scale1b = 1 + (ncol(mutation_gene_table) - 14)/14



yyy <- rep("Grade",length(sample$case))
data4<-data.frame(sample$case,yyy,sample$order,sample$Grade)
colnames(data4)<-c("xxx","yyy","zzz","his")

yyy <- rep("Gender",length(sample$case))
data5<-data.frame(sample$case,yyy,sample$order,sample$Gender)
colnames(data5)<-c("xxx","yyy","zzz","his")

yyy <- rep("Age",length(sample$case))
data6<-data.frame(sample$case,yyy,sample$order,sample$Age)
colnames(data6)<-c("xxx","yyy","zzz","his")

yyy <- rep("BCL2.fusion",length(sample$case))
data7b<-data.frame(sample$case,yyy,sample$order,sample$BCL2)
colnames(data7b)<-c("xxx","yyy","zzz","his")

yyy <- rep("MYC.fusion",length(sample$case))
data9<-data.frame(sample$case,yyy,sample$order,sample$MYC)
colnames(data9)<-c("xxx","yyy","zzz","his")

F4b.PloT<-ggplot()+theme_classic()
F4b.PloT<-F4b.PloT+geom_tile(data = data4,aes(x=reorder(xxx,zzz),y=yyy,fill=his),color='black',width=1,height=1,size=0.5,stat='identity')
F4b.PloT<-F4b.PloT+scale_fill_manual(name=NULL,values=c(High=gg_color_hue(6)[1] ,Low=gg_color_hue(6)[3]),guide = guide_legend(override.aes=list(size=1),nrow=1),na.translate = F)
F4b.PloT<-F4b.PloT+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(1,0.1,0.1,0.1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.4,'cm'),legend.key.height=unit(0.4,'cm'),legend.position='none',axis.ticks.x = element_blank(),axis.line.x = element_blank(),axis.line.y = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=12,face='bold.italic'),axis.ticks.y = element_blank(), axis.text.y=element_blank(),legend.title=element_blank(),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F4b.PloT<-F4b.PloT+scale_x_discrete(position = "top")+xlab(NULL)+ylab(NULL)+guides(size=FALSE)
gene.scale4 = 1/14

F7b.PloT<-ggplot()+theme_classic()
F7b.PloT<-F7b.PloT+geom_tile(data = data7b,aes(x=reorder(xxx,zzz),y=yyy,fill=his),color='black',width=1,height=1,size=0.5,stat='identity')
F7b.PloT<-F7b.PloT+scale_fill_manual(name=NULL,values=c(Fusion='black',WT='white'),
                                     labels=c(Fusion='Fusion',WT='Wild type'),guide = guide_legend(override.aes=list(size=1),nrow=1),na.translate = F)
F7b.PloT<-F7b.PloT+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(0.1,0.1,0.1,0.1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.4,'cm'),legend.key.height=unit(0.4,'cm'),legend.position='none',axis.ticks.x = element_blank(),axis.line.x = element_blank(),axis.line.y = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=12,face='bold.italic'),axis.ticks.y = element_blank(), axis.text.y=element_blank(),legend.title=element_blank(),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F7b.PloT<-F7b.PloT+scale_x_discrete(position = "top")+xlab(NULL)+ylab(NULL)+guides(size=FALSE)
gene.scale7 = 1/14



F5b.PloT<-ggplot()+theme_classic()
F5b.PloT<-F5b.PloT+geom_tile(data = data5,aes(x=reorder(xxx,zzz),y=yyy,fill=his),color='black',width=1,height=1,size=0.5,stat='identity')
F5b.PloT<-F5b.PloT+scale_fill_manual(name=NULL,values=c(Male='#b3cde3',Female='#fccde5'),labels=c(Male='Male',Female='Female'),guide = guide_legend(override.aes=list(size=1),nrow=1),na.translate = F)
F5b.PloT<-F5b.PloT+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(0.1,0.1,0.1,0.1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.4,'cm'),legend.key.height=unit(0.4,'cm'),legend.position='none',axis.ticks.x = element_blank(),axis.line.x = element_blank(),axis.line.y = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=12,face='bold.italic'),axis.ticks.y = element_blank(), axis.text.y=element_blank(),legend.title=element_blank(),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F5b.PloT<-F5b.PloT+scale_x_discrete(position = "top")+xlab(NULL)+ylab(NULL)+guides(size=FALSE)
gene.scale5 = 1/14

F6b.PloT<-ggplot()+theme_classic()
F6b.PloT<-F6b.PloT+geom_tile(data = data6,aes(x=reorder(xxx,zzz),y=yyy,fill=as.numeric(as.character(his))),color='black',width=1,height=1,size=0.5,stat='identity')
F6b.PloT<-F6b.PloT+scale_fill_gradientn(name=NULL,colours=c('white','#f6e8c3','#dfc27d','#bf812d','#8c510a'),limits=c(50,80),breaks=seq(50,80,5))
F6b.PloT<-F6b.PloT+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(0.1,0.1,0.1,0.1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(1.1,'cm'),legend.key.height=unit(0.4,'cm'),legend.position='none',axis.ticks.x = element_blank(),axis.line.x = element_blank(),axis.line.y = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=12,face='bold'),axis.ticks.y = element_blank(), axis.text.y=element_blank(),legend.title=element_blank(),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F6b.PloT<-F6b.PloT+scale_x_discrete(position = "top")+xlab(NULL)+ylab(NULL)+guides(size=FALSE)
gene.scale6 = 1/14




F9b.PloT<-ggplot()+theme_classic()
F9b.PloT<-F9b.PloT+geom_tile(data = data9,aes(x=reorder(xxx,zzz),y=yyy,fill=his),color='black',width=1,height=1,size=0.5,stat='identity')
F9b.PloT<-F9b.PloT+scale_fill_manual(name=NULL,values=c(Fusion='black',WT='white'),
                                     labels=c(Fusion='Fusion',WT='Wild type'),guide = guide_legend(override.aes=list(size=1),nrow=1),na.translate = F)
F9b.PloT<-F9b.PloT+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(0.1,0.1,0.1,0.1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(1.1,'cm'),legend.key.height=unit(0.4,'cm'),legend.position='none',axis.ticks.x = element_blank(),axis.line.x = element_blank(),axis.line.y = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=12,face='bold'),axis.ticks.y = element_blank(), axis.text.y=element_blank(),legend.title=element_blank(),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F9b.PloT<-F9b.PloT+scale_x_discrete(position = "top")+xlab(NULL)+ylab(NULL)+guides(size=FALSE)
gene.scale11 = 1/14




figure_1<-cbind(
          rbind(ggplotGrob(F5b.ploT),ggplotGrob(F4b.ploT),ggplotGrob(F1b.ploT)),
          rbind(ggplotGrob(F5b.plOT),ggplotGrob(F4b.plOT),ggplotGrob(F1b.plOT)),
          rbind(ggplotGrob(F5b.pLOT),ggplotGrob(F4b.pLOT),ggplotGrob(F1b.pLOT)),
          rbind(ggplotGrob(F5b.PLOT),ggplotGrob(F4b.PLOT),ggplotGrob(F1b.PLOT)),
          rbind(ggplotGrob(F5b.PLOt),ggplotGrob(F4b.PLOt),ggplotGrob(F1b.PLOt)),
          rbind(ggplotGrob(F5b.PLoT),ggplotGrob(F4b.PLoT),ggplotGrob(F1b.PLoT)),
          rbind(ggplotGrob(F5b.PlOT),ggplotGrob(F4b.PlOT),ggplotGrob(F1b.PlOT)),
          rbind(ggplotGrob(F5b.pLot),ggplotGrob(F4b.pLot),ggplotGrob(F1b.pLot)),
          rbind(ggplotGrob(F5b.PloT),ggplotGrob(F4b.PloT),ggplotGrob(F1b.PloT)),
          size="last")

panels <- figure_1$layout$t[grep("panel", figure_1$layout$name)]

figure_1$widths[5]  <- unit(2/4,'null')
figure_1$widths[14] <- unit(2/4,'null')
figure_1$widths[23] <- unit(2/4,'null')
figure_1$widths[32] <- unit(2/4,'null')
figure_1$widths[41] <- unit(2/4,'null')
figure_1$widths[50] <- unit(2/4,'null')
figure_1$widths[59] <- unit(2/4,'null')
figure_1$widths[68] <- unit(3/4,'null')
figure_1$widths[77] <- unit(2/4,'null')

figure_1$heights[panels][1] <- unit(gene.scale4,'null')
figure_1$heights[panels][2] <- unit(gene.scale4,'null')
figure_1$heights[panels][3] <- unit(gene.scale1b,'null')

figure_1$heights[panels][4] <- unit(gene.scale4,'null')
figure_1$heights[panels][5] <- unit(gene.scale4,'null')
figure_1$heights[panels][6] <- unit(gene.scale1b,'null')

figure_1$heights[panels][7] <- unit(gene.scale4,'null')
figure_1$heights[panels][8] <- unit(gene.scale4,'null')
figure_1$heights[panels][9] <- unit(gene.scale1b,'null')

figure_1$heights[panels][10] <- unit(gene.scale4,'null')
figure_1$heights[panels][11] <- unit(gene.scale4,'null')
figure_1$heights[panels][12] <- unit(gene.scale1b,'null')

figure_1$heights[panels][13] <- unit(gene.scale4,'null')
figure_1$heights[panels][14] <- unit(gene.scale4,'null')
figure_1$heights[panels][15] <- unit(gene.scale1b,'null')

figure_1$heights[panels][16] <- unit(gene.scale4,'null')
figure_1$heights[panels][17] <- unit(gene.scale4,'null')
figure_1$heights[panels][18] <- unit(gene.scale1b,'null')

figure_1$heights[panels][19] <- unit(gene.scale4,'null')
figure_1$heights[panels][20] <- unit(gene.scale4,'null')
figure_1$heights[panels][21] <- unit(gene.scale1b,'null')

figure_1$heights[panels][22] <- unit(gene.scale4,'null')
figure_1$heights[panels][23] <- unit(gene.scale4,'null')
figure_1$heights[panels][24] <- unit(gene.scale1b,'null')

figure_1$heights[panels][25] <- unit(gene.scale4,'null')
figure_1$heights[panels][26] <- unit(gene.scale4,'null')
figure_1$heights[panels][27] <- unit(gene.scale1b,'null')







#grid.draw(figure_1)
ggsave(file="landscape_noncoding_v4.pdf", plot=figure_1,bg = 'white', width = 20, height = 38, units = 'cm', dpi = 600)



