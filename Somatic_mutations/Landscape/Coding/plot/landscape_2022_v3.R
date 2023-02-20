# WatchDHL
library(rstudioapi)

# Getting the path of your current open file
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))
print( getwd() )

library(ggplot2)
library(gridExtra)
library(grid)


#source("./lib/CELLO.R")

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}


interestedGene <- c("ZCCHC7","CDKN2A/B","MDM2","IKZF3","EP300","SF3B1","EZH2","MYD88","EBF1","SOCS1",
                    "RFTN1","POU2AF1","LTB","EEF1A1","HIST1H1B","HIST1H2AC","HIST1H1C","BCL7A","SGK1","SLC30A4","DTX1",
                    "MED12","BCR","CCND3","FOXO1","FAM102A","IRF8","BMP7","PIM1","PCDHA11","TP53","TNFRSF14","CREBBP","KMT2D","MYC","BCL2")


CNV_gene <- c("CDKN2A/B","MDM2","IKZF3","ZCCHC7")



Sample_features <- read.table("Genetic_profiles_new.txt",sep = "\t",header = T)
Sample_features[is.na(Sample_features)] <- "f_NA"
sample <- data.frame(Sample_features)
sample <- sample[which(sample$patient=='P8'),]
pat <- unique(sample$patient)
case <- unique(sample$case)
#raw data preprocessing

somatic <- read.delim('./Mut/Annotated_Coding.txt',header = T)
somatic <- somatic[which(somatic$CaseID %in% pat),]
FL.a <- somatic[which(somatic$T1_freq >= 5 & somatic$HighQread_T1 >= 2),]
FL.a$ID <- "P8.FL.a"
FL.b <- somatic[which(as.numeric(as.character(somatic$T2_freq)) >= 10 & as.numeric(as.character(somatic$HighQread_T2)) >= 2),]
FL.b$ID <- "P8.FL.b"
DHL <- somatic[which(somatic$T3_freq >= 10 & somatic$HighQread_T3 >= 2),]
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
    temp.P <- savi$Effect[which(savi$case == as.character(case[i]) & savi$Gene_Name == selGene[j] & savi$Impact!='LOW' & savi$Impact!='MODIFIER')]
    if( any(grepl(pattern = 'stop',x = temp.P,fixed=T) | grepl(pattern = 'nonsense',x = temp.P,fixed=T) ) ){
      gene.Matrix[i,j] <- 'StopGained'
    }
    else if( any(grepl(pattern = 'frameshift_variant',x = temp.P,fixed=T)) ){
      gene.Matrix[i,j] <- 'Frameshift'
    }
    else if( any( grepl(pattern = 'splice',x = temp.P,fixed=T) | grepl(pattern = 'splice_acceptor',x = temp.P,fixed=T) ) ){
      gene.Matrix[i,j] <- 'SpliceSite'
    }
    else if( any( grepl(pattern = 'inframe',x = temp.P,fixed=T) ) ){
      gene.Matrix[i,j] <- 'InframeIndel'
    }
    else if( any( grepl(pattern = 'missense_variant', x = temp.P,fixed=T)) ){
      gene.Matrix[i,j] <- 'Missense'
    }
    else if( any( grepl(pattern = 'start', x = temp.P,fixed=T)) ){
      gene.Matrix[i,j] <- 'StartLost'
    }
    else{}
  }
}


colnames(gene.Matrix) <- selGene
rownames(gene.Matrix) <- case

mutation_gene_table <- gene.Matrix
#####


#sample$order[which(sample$case=="P8.FL.a")] <- 1
sample$order[which(sample$case=="P8.FL.b")] <- 2
sample$order[which(sample$case=="P8.DHL")] <-  3
sample$order[which(sample$case=="P8.FL.a")] <- 1

plot_1b.data <- melt(mutation_gene_table)
colnames(plot_1b.data) <- c('ID','gene','value')

FocalCNV <- read.delim('./Mut/CNV_matrix.txt',header = T)
FocalCNV <- FocalCNV[which(as.character(FocalCNV$Gene) %in% interestedGene),]
xxx <- FocalCNV[,which(colnames(FocalCNV) %in% case)]

rownames(xxx) <- FocalCNV$Gene
NEW_CNV <- melt(t(xxx))
colnames(NEW_CNV) <- c('ID','gene','value')

theAmp <- NEW_CNV[which(NEW_CNV$value >= 0.5 ),]
theDel <- NEW_CNV[which(NEW_CNV$value <= -0.5 ),]


new_table <- merge(plot_1b.data,sample,1,1)

SV <- read.delim('./Mut/SV.txt',header = T)
theSV <- SV[which(SV$ID %in% case),]

LOH <- read.delim('./Mut/LOH_matrix.txt',header = T)
LOH <- LOH[which(as.character(LOH$Gene) %in% interestedGene),]
xxx <- LOH[,which(colnames(LOH) %in% case)]
rownames(xxx) <- LOH$Gene
NEW_LOH <- melt(t(xxx))
colnames(NEW_LOH) <- c('ID','gene','value')
theLOH <- NEW_LOH[which(NEW_LOH$value == 1 ),]

#orderID<-c(1:nrow(plot_1b.data))
F1b.pLot<-ggplot()+theme_classic()
F1b.pLot<-F1b.pLot+geom_tile(data = new_table,aes(x=reorder(ID,order),y=gene,fill=value),color='black',width=1,height=1,size=0.5,stat='identity')
F1b.pLot<-F1b.pLot+geom_point(data = theSV,aes(x=ID,y=gene),shape=23,fill='black',size=3)
F1b.pLot<-F1b.pLot+geom_point(data = theAmp,aes(x=ID,y=gene),shape=24,fill='red',size=3)
F1b.pLot<-F1b.pLot+geom_point(data = theDel,aes(x=ID,y=gene),shape=25,fill='blue',size=3)
F1b.pLot<-F1b.pLot+geom_text(data = theLOH,aes(x=ID,y=gene),label="L",color='white',size=4)
F1b.pLot<-F1b.pLot+scale_fill_manual(name=NULL,values=c(N='white',Missense='#66c2a5',InframeIndel='#a6d854',Frameshift='#fc8d62',SpliceSite=gg_color_hue(6)[5],StopGained='#df65b0',StartLost='#762a83'),labels=c(N='None',Missense='Missense',Inframe='In frame indel',Frameshift='Frame shift', SpliceSite='Splice site',StopGained='Stop gained',StartLost='Start lost'),guide = guide_legend(override.aes=list(size=1),nrow=7),na.translate = F)
F1b.pLot<-F1b.pLot+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(0.1,0.1,0.1,0.1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.68,'cm'),legend.position='none',axis.ticks.x = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=16,face='bold.italic'),axis.ticks.y = element_blank(), axis.text.y=element_blank(),axis.line = element_blank(),
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

F4b.pLot<-ggplot()+theme_classic()
F4b.pLot<-F4b.pLot+geom_tile(data = data4,aes(x=reorder(xxx,zzz),y=yyy,fill=his),color='black',width=1,height=1,size=0.5,stat='identity')
F4b.pLot<-F4b.pLot+scale_fill_manual(name=NULL,values=c(High=gg_color_hue(6)[1] ,Low=gg_color_hue(6)[3]),guide = guide_legend(override.aes=list(size=1),nrow=1),na.translate = F)
F4b.pLot<-F4b.pLot+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(1,0.1,0,0.1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.4,'cm'),legend.key.height=unit(0.4,'cm'),legend.position='none',axis.ticks.x = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=12,face='bold.italic'),axis.ticks.y = element_blank(), axis.text.y=element_blank(),legend.title=element_blank(),axis.line = element_blank(),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F4b.pLot<-F4b.pLot+scale_x_discrete(position = "top")+xlab(NULL)+ylab(NULL)+guides(size=FALSE)
gene.scale4 = 1/14


F5b.pLot<-ggplot()+theme_classic()
F5b.pLot<-F5b.pLot+geom_tile(data = data5,aes(x=reorder(xxx,zzz),y=yyy,fill=his),color='black',width=1,height=1,size=0.5,stat='identity')
F5b.pLot<-F5b.pLot+scale_fill_manual(name=NULL,values=c(Male='#b3cde3',Female='#fccde5'),labels=c(Male='Male',Female='Female'),guide = guide_legend(override.aes=list(size=1),nrow=1),na.translate = F)
F5b.pLot<-F5b.pLot+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(0,0.1,0.1,0.1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.4,'cm'),legend.key.height=unit(0.4,'cm'),legend.position='none',axis.ticks.x = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=12,face='bold.italic'),axis.ticks.y = element_blank(), axis.text.y=element_blank(),legend.title=element_blank(),axis.line = element_blank(),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F5b.pLot<-F5b.pLot+scale_x_discrete(position = "top")+xlab(NULL)+ylab(NULL)+guides(size=FALSE)
gene.scale5 = 1/14





###P1
sample <- data.frame(Sample_features)
sample <- sample[which(sample$patient=='P1'),]
pat <- unique(sample$patient)
case <- unique(sample$case)
#raw data preprocessing

somatic <- read.delim('./Mut/Annotated_Coding.txt',header = T)
somatic <- somatic[which(somatic$CaseID %in% pat),]
FL <- somatic[which(somatic$T1_freq >= 10 & somatic$HighQread_T1 >= 2),]
FL$ID <- "P1.FL"
DHL <- somatic[which(somatic$T3_freq >= 10 & somatic$HighQread_T3 >= 2),]
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
    temp.P <- savi$Effect[which(savi$case == as.character(case[i]) & savi$Gene_Name == selGene[j] & savi$Impact!='LOW' & savi$Impact!='MODIFIER')]
    if( any(grepl(pattern = 'stop',x = temp.P,fixed=T) | grepl(pattern = 'nonsense',x = temp.P,fixed=T) ) ){
      gene.Matrix[i,j] <- 'StopGained'
    }
    else if( any(grepl(pattern = 'frameshift_variant',x = temp.P,fixed=T)) ){
      gene.Matrix[i,j] <- 'Frameshift'
    }
    else if( any( grepl(pattern = 'splice',x = temp.P,fixed=T) | grepl(pattern = 'splice_acceptor',x = temp.P,fixed=T) ) ){
      gene.Matrix[i,j] <- 'SpliceSite'
    }
    else if( any( grepl(pattern = 'inframe',x = temp.P,fixed=T) ) ){
      gene.Matrix[i,j] <- 'InframeIndel'
    }
    else if( any( grepl(pattern = 'missense_variant', x = temp.P,fixed=T)) ){
      gene.Matrix[i,j] <- 'Missense'
    }
    else if( any( grepl(pattern = 'start', x = temp.P,fixed=T)) ){
      gene.Matrix[i,j] <- 'StartLost'
    }
    else{}
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

FocalCNV <- read.delim('./Mut/CNV_matrix.txt',header = T)
FocalCNV <- FocalCNV[which(as.character(FocalCNV$Gene) %in% interestedGene),]
xxx <- FocalCNV[,which(colnames(FocalCNV) %in% case)]
rownames(xxx) <- FocalCNV$Gene
NEW_CNV <- melt(t(xxx))
colnames(NEW_CNV) <- c('ID','gene','value')

theAmp <- NEW_CNV[which(NEW_CNV$value >= 0.5 ),]
theDel <- NEW_CNV[which(NEW_CNV$value <= -0.5 ),]


new_table <- merge(plot_1b.data,sample,1,1)

SV <- read.delim('./Mut/SV.txt',header = T)
theSV <- SV[which(SV$ID %in% case),]

LOH <- read.delim('./Mut/LOH_matrix.txt',header = T)
LOH <- LOH[which(as.character(LOH$Gene) %in% interestedGene),]
xxx <- LOH[,which(colnames(LOH) %in% case)]
rownames(xxx) <- LOH$Gene
NEW_LOH <- melt(t(xxx))
colnames(NEW_LOH) <- c('ID','gene','value')
theLOH <- NEW_LOH[which(NEW_LOH$value == 1 ),]

#orderID<-c(1:nrow(plot_1b.data))
F1b.ploT<-ggplot()+theme_classic()
F1b.ploT<-F1b.ploT+geom_tile(data = new_table,aes(x=reorder(ID,order),y=gene,fill=value),color='black',width=1,height=1,size=0.5,stat='identity')
F1b.ploT<-F1b.ploT+geom_point(data = theSV,aes(x=ID,y=gene),shape=23,fill='black',size=3)
F1b.ploT<-F1b.ploT+geom_point(data = theAmp,aes(x=ID,y=gene),shape=24,fill='red',size=3)
F1b.ploT<-F1b.ploT+geom_point(data = theDel,aes(x=ID,y=gene),shape=25,fill='blue',size=3)
F1b.ploT<-F1b.ploT+geom_text(data = theLOH,aes(x=ID,y=gene),label="L",color='white',size=4)
F1b.ploT<-F1b.ploT+scale_fill_manual(name=NULL,values=c(N='white',Missense='#66c2a5',InframeIndel='#a6d854',Frameshift='#fc8d62',SpliceSite=gg_color_hue(6)[5],StopGained='#df65b0',StartLost='#762a83'),labels=c(N='None',Missense='Missense',Inframe='In frame indel',Frameshift='Frame shift', SpliceSite='Splice site',StopGained='Stop gained',StartLost='Start lost'),guide = guide_legend(override.aes=list(size=1),nrow=7),na.translate = F)
F1b.ploT<-F1b.ploT+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(0.1,0.1,0.1,1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.68,'cm'),legend.position='none',axis.ticks.x = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=16,face='bold.italic'),axis.text.y=element_text(size=14,vjust=0.5,hjust=1,face='italic',color='black'),axis.line = element_blank(),#axis.ticks.y.left = element_line(),
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



F4b.ploT<-ggplot()+theme_classic()
F4b.ploT<-F4b.ploT+geom_tile(data = data4,aes(x=reorder(xxx,zzz),y=yyy,fill=his),color='black',width=1,height=1,size=0.5,stat='identity')
F4b.ploT<-F4b.ploT+scale_fill_manual(name=NULL,values=c(High=gg_color_hue(6)[1] ,Low=gg_color_hue(6)[3]),guide = guide_legend(override.aes=list(size=1),nrow=1),na.translate = F)
F4b.ploT<-F4b.ploT+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(1,0.1,0,0.1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.4,'cm'),legend.key.height=unit(0.4,'cm'),legend.position='none',axis.ticks.x = element_blank(),axis.line = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=12,face='bold.italic'),axis.text.y=element_text(size=14,vjust=0.5,hjust=1,face='plain',color='black'),legend.title=element_blank(),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F4b.ploT<-F4b.ploT+scale_x_discrete(position = "top")+xlab(NULL)+ylab(NULL)+guides(size=FALSE)
gene.scale4 = 1/14



F5b.ploT<-ggplot()+theme_classic()
F5b.ploT<-F5b.ploT+geom_tile(data = data5,aes(x=reorder(xxx,zzz),y=yyy,fill=his),color='black',width=1,height=1,size=0.5,stat='identity')
F5b.ploT<-F5b.ploT+scale_fill_manual(name=NULL,values=c(Male='#b3cde3',Female='#fccde5'),labels=c(Male='Male',Female='Female'),guide = guide_legend(override.aes=list(size=1),nrow=1),na.translate = F)
F5b.ploT<-F5b.ploT+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(0,0.1,0.1,0.1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.4,'cm'),legend.key.height=unit(0.4,'cm'),legend.position='none',axis.ticks.x = element_blank(),axis.line = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=12,face='bold.italic'),axis.text.y=element_text(size=14,vjust=0.5,hjust=1,face='plain',color='black'),legend.title=element_blank(),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F5b.ploT<-F5b.ploT+scale_x_discrete(position = "top")+xlab(NULL)+ylab(NULL)+guides(size=FALSE)
gene.scale5 = 1/14




##P2
sample <- data.frame(Sample_features)
sample <- sample[which(sample$patient=='P2'),]
pat <- unique(sample$patient)
case <- unique(sample$case)
#raw data preprocessing

somatic <- read.delim('./Mut/Annotated_Coding.txt',header = T)
somatic <- somatic[which(somatic$CaseID %in% pat),]
FL <- somatic[which(somatic$T1_freq >= 5 & somatic$HighQread_T1 >= 2),]
FL$ID <- "P2.FL"
DHL <- somatic[which(somatic$T3_freq >= 10 & somatic$HighQread_T3 >= 2),]
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
    temp.P <- savi$Effect[which(savi$case == as.character(case[i]) & savi$Gene_Name == selGene[j] & savi$Impact!='LOW' & savi$Impact!='MODIFIER')]
    if( any(grepl(pattern = 'stop',x = temp.P,fixed=T) | grepl(pattern = 'nonsense',x = temp.P,fixed=T) ) ){
      gene.Matrix[i,j] <- 'StopGained'
    }
    else if( any(grepl(pattern = 'frameshift_variant',x = temp.P,fixed=T)) ){
      gene.Matrix[i,j] <- 'Frameshift'
    }
    else if( any( grepl(pattern = 'splice',x = temp.P,fixed=T) | grepl(pattern = 'splice_acceptor',x = temp.P,fixed=T) ) ){
      gene.Matrix[i,j] <- 'SpliceSite'
    }
    else if( any( grepl(pattern = 'inframe',x = temp.P,fixed=T) ) ){
      gene.Matrix[i,j] <- 'InframeIndel'
    }
    else if( any( grepl(pattern = 'missense_variant', x = temp.P,fixed=T)) ){
      gene.Matrix[i,j] <- 'Missense'
    }
    else if( any( grepl(pattern = 'start', x = temp.P,fixed=T)) ){
      gene.Matrix[i,j] <- 'StartLost'
    }
    else{}
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

FocalCNV <- read.delim('./Mut/CNV_matrix.txt',header = T)
FocalCNV <- FocalCNV[which(as.character(FocalCNV$Gene) %in% interestedGene),]
xxx <- FocalCNV[,which(colnames(FocalCNV) %in% case)]
rownames(xxx) <- FocalCNV$Gene
NEW_CNV <- melt(t(xxx))
colnames(NEW_CNV) <- c('ID','gene','value')

theAmp <- NEW_CNV[which(NEW_CNV$value >= 0.5 ),]
theDel <- NEW_CNV[which(NEW_CNV$value <= -0.5 ),]


new_table <- merge(plot_1b.data,sample,1,1)

SV <- read.delim('./Mut/SV.txt',header = T)
theSV <- SV[which(SV$ID %in% case),]

LOH <- read.delim('./Mut/LOH_matrix.txt',header = T)
LOH <- LOH[which(as.character(LOH$Gene) %in% interestedGene),]
xxx <- LOH[,which(colnames(LOH) %in% case)]
rownames(xxx) <- LOH$Gene
NEW_LOH <- melt(t(xxx))
colnames(NEW_LOH) <- c('ID','gene','value')
theLOH <- NEW_LOH[which(NEW_LOH$value == 1 ),]

#orderID<-c(1:nrow(plot_1b.data))
F1b.plOT<-ggplot()+theme_classic()
F1b.plOT<-F1b.plOT+geom_tile(data = new_table,aes(x=reorder(ID,order),y=gene,fill=value),color='black',width=1,height=1,size=0.5,stat='identity')
F1b.plOT<-F1b.plOT+geom_point(data = theSV,aes(x=ID,y=gene),shape=23,fill='black',size=3)
F1b.plOT<-F1b.plOT+geom_point(data = theAmp,aes(x=ID,y=gene),shape=24,fill='red',size=3)
F1b.plOT<-F1b.plOT+geom_point(data = theDel,aes(x=ID,y=gene),shape=25,fill='blue',size=3)
F1b.plOT<-F1b.plOT+geom_text(data = theLOH,aes(x=ID,y=gene),label="L",color='white',size=4)
F1b.plOT<-F1b.plOT+scale_fill_manual(name=NULL,values=c(N='white',Missense='#66c2a5',InframeIndel='#a6d854',Frameshift='#fc8d62',SpliceSite=gg_color_hue(6)[5],StopGained='#df65b0',StartLost='#762a83'),labels=c(N='None',Missense='Missense',Inframe='In frame indel',Frameshift='Frame shift', SpliceSite='Splice site',StopGained='Stop gained',StartLost='Start lost'),guide = guide_legend(override.aes=list(size=1),nrow=7),na.translate = F)
F1b.plOT<-F1b.plOT+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(0.1,0.1,0.1,0.1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.68,'cm'),legend.position='none',axis.ticks.x = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=16,face='bold.italic'),axis.ticks.y = element_blank(), axis.text.y=element_blank(),axis.line = element_blank(),
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



F4b.plOT<-ggplot()+theme_classic()
F4b.plOT<-F4b.plOT+geom_tile(data = data4,aes(x=reorder(xxx,zzz),y=yyy,fill=his),color='black',width=1,height=1,size=0.5,stat='identity')
F4b.plOT<-F4b.plOT+scale_fill_manual(name=NULL,values=c(High=gg_color_hue(6)[1] ,Low=gg_color_hue(6)[3]),guide = guide_legend(override.aes=list(size=1),nrow=1),na.translate = F)
F4b.plOT<-F4b.plOT+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(1,0.1,0,0.1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.4,'cm'),legend.key.height=unit(0.4,'cm'),legend.position='none',axis.ticks.x = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=12,face='bold.italic'),axis.ticks.y = element_blank(), axis.text.y=element_blank(),legend.title=element_blank(),axis.line = element_blank(),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F4b.plOT<-F4b.plOT+scale_x_discrete(position = "top")+xlab(NULL)+ylab(NULL)+guides(size=FALSE)
gene.scale4 = 1/14




F5b.plOT<-ggplot()+theme_classic()
F5b.plOT<-F5b.plOT+geom_tile(data = data5,aes(x=reorder(xxx,zzz),y=yyy,fill=his),color='black',width=1,height=1,size=0.5,stat='identity')
F5b.plOT<-F5b.plOT+scale_fill_manual(name=NULL,values=c(Male='#b3cde3',Female='#fccde5'),labels=c(Male='Male',Female='Female'),guide = guide_legend(override.aes=list(size=1),nrow=1),na.translate = F)
F5b.plOT<-F5b.plOT+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(0,0.1,0.1,0.1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.4,'cm'),legend.key.height=unit(0.4,'cm'),legend.position='none',axis.ticks.x = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=12,face='bold.italic'),axis.ticks.y = element_blank(), axis.text.y=element_blank(),legend.title=element_blank(),axis.line = element_blank(),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F5b.plOT<-F5b.plOT+scale_x_discrete(position = "top")+xlab(NULL)+ylab(NULL)+guides(size=FALSE)
gene.scale5 = 1/14




###P4
sample <- data.frame(Sample_features)
sample <- sample[which(sample$patient=='P4'),]
pat <- unique(sample$patient)
case <- unique(sample$case)
#raw data preprocessing

somatic <- read.delim('./Mut/Annotated_Coding.txt',header = T)
somatic <- somatic[which(somatic$CaseID %in% pat),]
FL <- somatic[which(somatic$T1_freq >= 10 & somatic$HighQread_T1 >= 2),]
FL$ID <- "P4.FL"
DHL <- somatic[which(somatic$T3_freq >= 10 & somatic$HighQread_T3 >= 2),]
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
    temp.P <- savi$Effect[which(savi$case == as.character(case[i]) & savi$Gene_Name == selGene[j] & savi$Impact!='LOW' & savi$Impact!='MODIFIER')]
    if( any(grepl(pattern = 'stop',x = temp.P,fixed=T) | grepl(pattern = 'nonsense',x = temp.P,fixed=T) ) ){
      gene.Matrix[i,j] <- 'StopGained'
    }
    else if( any(grepl(pattern = 'frameshift_variant',x = temp.P,fixed=T)) ){
      gene.Matrix[i,j] <- 'Frameshift'
    }
    else if( any( grepl(pattern = 'splice',x = temp.P,fixed=T) | grepl(pattern = 'splice_acceptor',x = temp.P,fixed=T) ) ){
      gene.Matrix[i,j] <- 'SpliceSite'
    }
    else if( any( grepl(pattern = 'inframe',x = temp.P,fixed=T) ) ){
      gene.Matrix[i,j] <- 'InframeIndel'
    }
    else if( any( grepl(pattern = 'missense_variant', x = temp.P,fixed=T)) ){
      gene.Matrix[i,j] <- 'Missense'
    }
    else if( any( grepl(pattern = 'start', x = temp.P,fixed=T)) ){
      gene.Matrix[i,j] <- 'StartLost'
    }
    else{}
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

FocalCNV <- read.delim('./Mut/CNV_matrix.txt',header = T)
FocalCNV <- FocalCNV[which(as.character(FocalCNV$Gene) %in% interestedGene),]
xxx <- FocalCNV[,which(colnames(FocalCNV) %in% case)]
rownames(xxx) <- FocalCNV$Gene
NEW_CNV <- melt(t(xxx))
colnames(NEW_CNV) <- c('ID','gene','value')

theAmp <- NEW_CNV[which(NEW_CNV$value >= 0.5 ),]
theDel <- NEW_CNV[which(NEW_CNV$value <= -0.5 ),]


new_table <- merge(plot_1b.data,sample,1,1)

SV <- read.delim('./Mut/SV.txt',header = T)
theSV <- SV[which(SV$ID %in% case),]

LOH <- read.delim('./Mut/LOH_matrix.txt',header = T)
LOH <- LOH[which(as.character(LOH$Gene) %in% interestedGene),]
xxx <- LOH[,which(colnames(LOH) %in% case)]
rownames(xxx) <- LOH$Gene
NEW_LOH <- melt(t(xxx))
colnames(NEW_LOH) <- c('ID','gene','value')
theLOH <- NEW_LOH[which(NEW_LOH$value == 1 ),]

#orderID<-c(1:nrow(plot_1b.data))
F1b.PLOT<-ggplot()+theme_classic()
F1b.PLOT<-F1b.PLOT+geom_tile(data = new_table,aes(x=reorder(ID,order),y=gene,fill=value),color='black',width=1,height=1,size=0.5,stat='identity')
F1b.PLOT<-F1b.PLOT+geom_point(data = theSV,aes(x=ID,y=gene),shape=23,fill='black',size=3)
F1b.PLOT<-F1b.PLOT+geom_point(data = theAmp,aes(x=ID,y=gene),shape=24,fill='red',size=3)
F1b.PLOT<-F1b.PLOT+geom_point(data = theDel,aes(x=ID,y=gene),shape=25,fill='blue',size=3)+geom_text(data = theLOH,aes(x=ID,y=gene),label="L",color='white',size=4)
F1b.PLOT<-F1b.PLOT+scale_fill_manual(name=NULL,values=c(N='white',Missense='#66c2a5',InframeIndel='#a6d854',Frameshift='#fc8d62',SpliceSite=gg_color_hue(6)[5],StopGained='#df65b0',StartLost='#762a83'),labels=c(N='None',Missense='Missense',Inframe='In frame indel',Frameshift='Frame shift', SpliceSite='Splice site',StopGained='Stop gained',StartLost='Start lost'),guide = guide_legend(override.aes=list(size=1),nrow=7),na.translate = F)
F1b.PLOT<-F1b.PLOT+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(0.1,0.1,0.1,0.1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.68,'cm'),legend.position='none',axis.ticks.x = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=16,face='bold.italic'),axis.ticks.y = element_blank(), axis.text.y=element_blank(),axis.line = element_blank(),
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



F4b.PLOT<-ggplot()+theme_classic()
F4b.PLOT<-F4b.PLOT+geom_tile(data = data4,aes(x=reorder(xxx,zzz),y=yyy,fill=his),color='black',width=1,height=1,size=0.5,stat='identity')
F4b.PLOT<-F4b.PLOT+scale_fill_manual(name=NULL,values=c(High=gg_color_hue(6)[1] ,Low=gg_color_hue(6)[3]),guide = guide_legend(override.aes=list(size=1),nrow=1),na.translate = F)
F4b.PLOT<-F4b.PLOT+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(1,0.1,0,0.1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.4,'cm'),legend.key.height=unit(0.4,'cm'),legend.position='none',axis.ticks.x = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=12,face='bold.italic'),axis.ticks.y = element_blank(), axis.text.y=element_blank(),legend.title=element_blank(),axis.line = element_blank(),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F4b.PLOT<-F4b.PLOT+scale_x_discrete(position = "top")+xlab(NULL)+ylab(NULL)+guides(size=FALSE)
gene.scale4 = 1/14




F5b.PLOT<-ggplot()+theme_classic()
F5b.PLOT<-F5b.PLOT+geom_tile(data = data5,aes(x=reorder(xxx,zzz),y=yyy,fill=his),color='black',width=1,height=1,size=0.5,stat='identity')
F5b.PLOT<-F5b.PLOT+scale_fill_manual(name=NULL,values=c(Male='#b3cde3',Female='#fccde5'),labels=c(Male='Male',Female='Female'),guide = guide_legend(override.aes=list(size=1),nrow=1),na.translate = F)
F5b.PLOT<-F5b.PLOT+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(0,0.1,0.1,0.1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.4,'cm'),legend.key.height=unit(0.4,'cm'),legend.position='none',axis.ticks.x = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=12,face='bold.italic'),axis.ticks.y = element_blank(), axis.text.y=element_blank(),legend.title=element_blank(),axis.line = element_blank(),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F5b.PLOT<-F5b.PLOT+scale_x_discrete(position = "top")+xlab(NULL)+ylab(NULL)+guides(size=FALSE)
gene.scale5 = 1/14

##P5
sample <- data.frame(Sample_features)
sample <- sample[which(sample$patient=='P5'),]
pat <- unique(sample$patient)
case <- unique(sample$case)
#raw data preprocessing

somatic <- read.delim('./Mut/Annotated_Coding.txt',header = T)
somatic <- somatic[which(somatic$CaseID %in% pat),]
FL <- somatic[which(somatic$T1_freq >= 10 & somatic$HighQread_T1 >= 2),]
FL$ID <- "P5.FL"
DHL <- somatic[which(somatic$T3_freq >= 10 & somatic$HighQread_T3 >= 2),]
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
    temp.P <- savi$Effect[which(savi$case == as.character(case[i]) & savi$Gene_Name == selGene[j] & savi$Impact!='LOW' & savi$Impact!='MODIFIER')]
    if( any(grepl(pattern = 'stop',x = temp.P,fixed=T) | grepl(pattern = 'nonsense',x = temp.P,fixed=T) ) ){
      gene.Matrix[i,j] <- 'StopGained'
    }
    else if( any(grepl(pattern = 'frameshift_variant',x = temp.P,fixed=T)) ){
      gene.Matrix[i,j] <- 'Frameshift'
    }
    else if( any( grepl(pattern = 'splice',x = temp.P,fixed=T) | grepl(pattern = 'splice_acceptor',x = temp.P,fixed=T) ) ){
      gene.Matrix[i,j] <- 'SpliceSite'
    }
    else if( any( grepl(pattern = 'inframe',x = temp.P,fixed=T) ) ){
      gene.Matrix[i,j] <- 'InframeIndel'
    }
    else if( any( grepl(pattern = 'missense_variant', x = temp.P,fixed=T)) ){
      gene.Matrix[i,j] <- 'Missense'
    }
    else if( any( grepl(pattern = 'start', x = temp.P,fixed=T)) ){
      gene.Matrix[i,j] <- 'StartLost'
    }
    else{}
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

FocalCNV <- read.delim('./Mut/CNV_matrix.txt',header = T)
FocalCNV <- FocalCNV[which(as.character(FocalCNV$Gene) %in% interestedGene),]
xxx <- FocalCNV[,which(colnames(FocalCNV) %in% case)]
rownames(xxx) <- FocalCNV$Gene
NEW_CNV <- melt(t(xxx))
colnames(NEW_CNV) <- c('ID','gene','value')

theAmp <- NEW_CNV[which(NEW_CNV$value >= 0.5 ),]
theDel <- NEW_CNV[which(NEW_CNV$value <= -0.5 ),]


new_table <- merge(plot_1b.data,sample,1,1)

SV <- read.delim('./Mut/SV.txt',header = T)
theSV <- SV[which(SV$ID %in% case),]

LOH <- read.delim('./Mut/LOH_matrix.txt',header = T)
LOH <- LOH[which(as.character(LOH$Gene) %in% interestedGene),]
xxx <- LOH[,which(colnames(LOH) %in% case)]
rownames(xxx) <- LOH$Gene
NEW_LOH <- melt(t(xxx))
colnames(NEW_LOH) <- c('ID','gene','value')
theLOH <- NEW_LOH[which(NEW_LOH$value == 1 ),]

#orderID<-c(1:nrow(plot_1b.data))
F1b.PLOt<-ggplot()+theme_classic()
F1b.PLOt<-F1b.PLOt+geom_tile(data = new_table,aes(x=reorder(ID,order),y=gene,fill=value),color='black',width=1,height=1,size=0.5,stat='identity')
F1b.PLOt<-F1b.PLOt+geom_point(data = theSV,aes(x=ID,y=gene),shape=23,fill='black',size=3)
F1b.PLOt<-F1b.PLOt+geom_point(data = theAmp,aes(x=ID,y=gene),shape=24,fill='red',size=3)
F1b.PLOt<-F1b.PLOt+geom_point(data = theDel,aes(x=ID,y=gene),shape=25,fill='blue',size=3)+geom_text(data = theLOH,aes(x=ID,y=gene),label="L",color='white',size=4)
F1b.PLOt<-F1b.PLOt+scale_fill_manual(name=NULL,values=c(N='white',Missense='#66c2a5',InframeIndel='#a6d854',Frameshift='#fc8d62',SpliceSite=gg_color_hue(6)[5],StopGained='#df65b0',StartLost='#762a83'),labels=c(N='None',Missense='Missense',Inframe='In frame indel',Frameshift='Frame shift', SpliceSite='Splice site',StopGained='Stop gained',StartLost='Start lost'),guide = guide_legend(override.aes=list(size=1),nrow=7),na.translate = F)
F1b.PLOt<-F1b.PLOt+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(0.1,0.1,0.1,0.1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.68,'cm'),legend.position='none',axis.ticks.x = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=16,face='bold.italic'),axis.ticks.y = element_blank(), axis.text.y=element_blank(),axis.line = element_blank(),
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


F4b.PLOt<-ggplot()+theme_classic()
F4b.PLOt<-F4b.PLOt+geom_tile(data = data4,aes(x=reorder(xxx,zzz),y=yyy,fill=his),color='black',width=1,height=1,size=0.5,stat='identity')
F4b.PLOt<-F4b.PLOt+scale_fill_manual(name=NULL,values=c(High=gg_color_hue(6)[1] ,Low=gg_color_hue(6)[3]),guide = guide_legend(override.aes=list(size=1),nrow=1),na.translate = F)
F4b.PLOt<-F4b.PLOt+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(1,0.1,0,0.1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.4,'cm'),legend.key.height=unit(0.4,'cm'),legend.position='none',axis.ticks.x = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=12,face='bold.italic'),axis.ticks.y = element_blank(), axis.text.y=element_blank(),legend.title=element_blank(),axis.line = element_blank(),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F4b.PLOt<-F4b.PLOt+scale_x_discrete(position = "top")+xlab(NULL)+ylab(NULL)+guides(size=FALSE)
gene.scale4 = 1/14


F5b.PLOt<-ggplot()+theme_classic()
F5b.PLOt<-F5b.PLOt+geom_tile(data = data5,aes(x=reorder(xxx,zzz),y=yyy,fill=his),color='black',width=1,height=1,size=0.5,stat='identity')
F5b.PLOt<-F5b.PLOt+scale_fill_manual(name=NULL,values=c(Male='#b3cde3',Female='#fccde5'),labels=c(Male='Male',Female='Female'),guide = guide_legend(override.aes=list(size=1),nrow=1),na.translate = F)
F5b.PLOt<-F5b.PLOt+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(0,0.1,0.1,0.1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.4,'cm'),legend.key.height=unit(0.4,'cm'),legend.position='none',axis.ticks.x = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=12,face='bold.italic'),axis.ticks.y = element_blank(), axis.text.y=element_blank(),legend.title=element_blank(),axis.line = element_blank(),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F5b.PLOt<-F5b.PLOt+scale_x_discrete(position = "top")+xlab(NULL)+ylab(NULL)+guides(size=FALSE)
gene.scale5 = 1/14


##P6
sample <- data.frame(Sample_features)
sample <- sample[which(sample$patient=='P6'),]
pat <- unique(sample$patient)
case <- unique(sample$case)
#raw data preprocessing

somatic <- read.delim('./Mut/Annotated_Coding.txt',header = T)
somatic <- somatic[which(somatic$CaseID %in% pat),]
FL <- somatic[which(somatic$T1_freq >= 10 & somatic$HighQread_T1 >= 2),]
FL$ID <- "P6.FL"
DHL <- somatic[which(somatic$T3_freq >= 10 & somatic$HighQread_T3 >= 2),]
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
    temp.P <- savi$Effect[which(savi$case == as.character(case[i]) & savi$Gene_Name == selGene[j] & savi$Impact!='LOW' & savi$Impact!='MODIFIER')]
    if( any(grepl(pattern = 'stop',x = temp.P,fixed=T) | grepl(pattern = 'nonsense',x = temp.P,fixed=T) ) ){
      gene.Matrix[i,j] <- 'StopGained'
    }
    else if( any(grepl(pattern = 'frameshift_variant',x = temp.P,fixed=T)) ){
      gene.Matrix[i,j] <- 'Frameshift'
    }
    else if( any( grepl(pattern = 'splice',x = temp.P,fixed=T) | grepl(pattern = 'splice_acceptor',x = temp.P,fixed=T) ) ){
      gene.Matrix[i,j] <- 'SpliceSite'
    }
    else if( any( grepl(pattern = 'inframe',x = temp.P,fixed=T) ) ){
      gene.Matrix[i,j] <- 'InframeIndel'
    }
    else if( any( grepl(pattern = 'missense_variant', x = temp.P,fixed=T)) ){
      gene.Matrix[i,j] <- 'Missense'
    }
    else if( any( grepl(pattern = 'start', x = temp.P,fixed=T)) ){
      gene.Matrix[i,j] <- 'StartLost'
    }
    else{}
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

FocalCNV <- read.delim('./Mut/CNV_matrix.txt',header = T)
FocalCNV <- FocalCNV[which(as.character(FocalCNV$Gene) %in% interestedGene),]
xxx <- FocalCNV[,which(colnames(FocalCNV) %in% case)]
rownames(xxx) <- FocalCNV$Gene
NEW_CNV <- melt(t(xxx))
colnames(NEW_CNV) <- c('ID','gene','value')

theAmp <- NEW_CNV[which(NEW_CNV$value >= 0.5 ),]
theDel <- NEW_CNV[which(NEW_CNV$value <= -0.5 ),]


new_table <- merge(plot_1b.data,sample,1,1)
SV <- read.delim('./Mut/SV.txt',header = T)
theSV <- SV[which(SV$ID %in% case),]

LOH <- read.delim('./Mut/LOH_matrix.txt',header = T)
LOH <- LOH[which(as.character(LOH$Gene) %in% interestedGene),]
xxx <- LOH[,which(colnames(LOH) %in% case)]
rownames(xxx) <- LOH$Gene
NEW_LOH <- melt(t(xxx))
colnames(NEW_LOH) <- c('ID','gene','value')
theLOH <- NEW_LOH[which(NEW_LOH$value == 1 ),]

#orderID<-c(1:nrow(plot_1b.data))
F1b.PLoT<-ggplot()+theme_classic()
F1b.PLoT<-F1b.PLoT+geom_tile(data = new_table,aes(x=reorder(ID,order),y=gene,fill=value),color='black',width=1,height=1,size=0.5,stat='identity')
F1b.PLoT<-F1b.PLoT+geom_point(data = theSV,aes(x=ID,y=gene),shape=23,fill='black',size=3)
F1b.PLoT<-F1b.PLoT+geom_point(data = theAmp,aes(x=ID,y=gene),shape=24,fill='red',size=3)
F1b.PLoT<-F1b.PLoT+geom_point(data = theDel,aes(x=ID,y=gene),shape=25,fill='blue',size=3)+geom_text(data = theLOH,aes(x=ID,y=gene),label="L",color='white',size=4)
F1b.PLoT<-F1b.PLoT+scale_fill_manual(name=NULL,values=c(N='white',Missense='#66c2a5',InframeIndel='#a6d854',Frameshift='#fc8d62',SpliceSite=gg_color_hue(6)[5],StopGained='#df65b0',StartLost='#762a83'),labels=c(N='None',Missense='Missense',Inframe='In frame indel',Frameshift='Frame shift', SpliceSite='Splice site',StopGained='Stop gained',StartLost='Start lost'),guide = guide_legend(override.aes=list(size=1),nrow=7),na.translate = F)
F1b.PLoT<-F1b.PLoT+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(0.1,0.1,0.1,0.1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.68,'cm'),legend.position='none',axis.ticks.x = element_blank(),axis.line = element_blank(),
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


F4b.PLoT<-ggplot()+theme_classic()
F4b.PLoT<-F4b.PLoT+geom_tile(data = data4,aes(x=reorder(xxx,zzz),y=yyy,fill=his),color='black',width=1,height=1,size=0.5,stat='identity')
F4b.PLoT<-F4b.PLoT+scale_fill_manual(name=NULL,values=c(High=gg_color_hue(6)[1] ,Low=gg_color_hue(6)[3]),guide = guide_legend(override.aes=list(size=1),nrow=1),na.translate = F)
F4b.PLoT<-F4b.PLoT+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(1,0.1,0,0.1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.4,'cm'),legend.key.height=unit(0.4,'cm'),legend.position='none',axis.ticks.x = element_blank(),axis.line = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=12,face='bold.italic'),axis.ticks.y = element_blank(), axis.text.y=element_blank(),legend.title=element_blank(),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F4b.PLoT<-F4b.PLoT+scale_x_discrete(position = "top")+xlab(NULL)+ylab(NULL)+guides(size=FALSE)
gene.scale4 = 1/14




F5b.PLoT<-ggplot()+theme_classic()
F5b.PLoT<-F5b.PLoT+geom_tile(data = data5,aes(x=reorder(xxx,zzz),y=yyy,fill=his),color='black',width=1,height=1,size=0.5,stat='identity')
F5b.PLoT<-F5b.PLoT+scale_fill_manual(name=NULL,values=c(Male='#b3cde3',Female='#fccde5'),labels=c(Male='Male',Female='Female'),guide = guide_legend(override.aes=list(size=1),nrow=1),na.translate = F)
F5b.PLoT<-F5b.PLoT+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(0,0.1,0.1,0.1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.4,'cm'),legend.key.height=unit(0.4,'cm'),legend.position='none',axis.ticks.x = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=12,face='bold.italic'),axis.ticks.y = element_blank(), axis.text.y=element_blank(),legend.title=element_blank(),axis.line = element_blank(),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F5b.PLoT<-F5b.PLoT+scale_x_discrete(position = "top")+xlab(NULL)+ylab(NULL)+guides(size=FALSE)
gene.scale5 = 1/14


###P7
sample <- data.frame(Sample_features)
sample <- sample[which(sample$patient=='P7'),]
pat <- unique(sample$patient)
case <- unique(sample$case)
#raw data preprocessing

somatic <- read.delim('./Mut/Annotated_Coding.txt',header = T)
somatic <- somatic[which(somatic$CaseID %in% pat),]
FL <- somatic[which(somatic$T1_freq >= 10 & somatic$HighQread_T1 >= 2),]
FL$ID <- "P7.FL"
DHL <- somatic[which(somatic$T3_freq >= 10 & somatic$HighQread_T3 >= 2),]
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
    temp.P <- savi$Effect[which(savi$case == as.character(case[i]) & savi$Gene_Name == selGene[j] & savi$Impact!='LOW' & savi$Impact!='MODIFIER')]
    if( any(grepl(pattern = 'stop',x = temp.P,fixed=T) | grepl(pattern = 'nonsense',x = temp.P,fixed=T) ) ){
      gene.Matrix[i,j] <- 'StopGained'
    }
    else if( any(grepl(pattern = 'frameshift_variant',x = temp.P,fixed=T)) ){
      gene.Matrix[i,j] <- 'Frameshift'
    }
    else if( any( grepl(pattern = 'splice',x = temp.P,fixed=T) | grepl(pattern = 'splice_acceptor',x = temp.P,fixed=T) ) ){
      gene.Matrix[i,j] <- 'SpliceSite'
    }
    else if( any( grepl(pattern = 'inframe',x = temp.P,fixed=T) ) ){
      gene.Matrix[i,j] <- 'InframeIndel'
    }
    else if( any( grepl(pattern = 'missense_variant', x = temp.P,fixed=T)) ){
      gene.Matrix[i,j] <- 'Missense'
    }
    else if( any( grepl(pattern = 'start', x = temp.P,fixed=T)) ){
      gene.Matrix[i,j] <- 'StartLost'
    }
    else{}
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

FocalCNV <- read.delim('./Mut/CNV_matrix.txt',header = T)
FocalCNV <- FocalCNV[which(as.character(FocalCNV$Gene) %in% interestedGene),]
xxx <- FocalCNV[,which(colnames(FocalCNV) %in% case)]
rownames(xxx) <- FocalCNV$Gene
NEW_CNV <- melt(t(xxx))
colnames(NEW_CNV) <- c('ID','gene','value')

theAmp <- NEW_CNV[which(NEW_CNV$value >= 0.5 ),]
theDel <- NEW_CNV[which(NEW_CNV$value <= -0.5 ),]


new_table <- merge(plot_1b.data,sample,1,1)
SV <- read.delim('./Mut/SV.txt',header = T)
theSV <- SV[which(SV$ID %in% case),]

LOH <- read.delim('./Mut/LOH_matrix.txt',header = T)
LOH <- LOH[which(as.character(LOH$Gene) %in% interestedGene),]
xxx <- LOH[,which(colnames(LOH) %in% case)]
rownames(xxx) <- LOH$Gene
NEW_LOH <- melt(t(xxx))
colnames(NEW_LOH) <- c('ID','gene','value')
theLOH <- NEW_LOH[which(NEW_LOH$value == 1 ),]

#orderID<-c(1:nrow(plot_1b.data))
F1b.PlOT<-ggplot()+theme_classic()
F1b.PlOT<-F1b.PlOT+geom_tile(data = new_table,aes(x=reorder(ID,order),y=gene,fill=value),color='black',width=1,height=1,size=0.5,stat='identity')
F1b.PlOT<-F1b.PlOT+geom_point(data = theSV,aes(x=ID,y=gene),shape=23,fill='black',size=3)
F1b.PlOT<-F1b.PlOT+geom_point(data = theAmp,aes(x=ID,y=gene),shape=24,fill='red',size=3)
F1b.PlOT<-F1b.PlOT+geom_point(data = theDel,aes(x=ID,y=gene),shape=25,fill='blue',size=3)+geom_text(data = theLOH,aes(x=ID,y=gene),label="L",color='white',size=4)
F1b.PlOT<-F1b.PlOT+scale_fill_manual(name=NULL,values=c(N='white',Missense='#66c2a5',InframeIndel='#a6d854',Frameshift='#fc8d62',SpliceSite=gg_color_hue(6)[5],StopGained='#df65b0',StartLost='#762a83'),labels=c(N='None',Missense='Missense',Inframe='In frame indel',Frameshift='Frame shift', SpliceSite='Splice site',StopGained='Stop gained',StartLost='Start lost'),guide = guide_legend(override.aes=list(size=1),nrow=7),na.translate = F)
F1b.PlOT<-F1b.PlOT+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(0.1,0.1,0.1,0.1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.68,'cm'),legend.position='none',axis.ticks.x = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=16,face='bold.italic'),axis.ticks.y = element_blank(), axis.text.y=element_blank(),axis.line = element_blank(),
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



F4b.PlOT<-ggplot()+theme_classic()
F4b.PlOT<-F4b.PlOT+geom_tile(data = data4,aes(x=reorder(xxx,zzz),y=yyy,fill=his),color='black',width=1,height=1,size=0.5,stat='identity')
F4b.PlOT<-F4b.PlOT+scale_fill_manual(name=NULL,values=c(High=gg_color_hue(6)[1] ,Low=gg_color_hue(6)[3]),guide = guide_legend(override.aes=list(size=1),nrow=1),na.translate = F)
F4b.PlOT<-F4b.PlOT+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(1,0.1,0,0.1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.4,'cm'),legend.key.height=unit(0.4,'cm'),legend.position='none',axis.ticks.x = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=12,face='bold.italic'),axis.ticks.y = element_blank(), axis.text.y=element_blank(),legend.title=element_blank(),axis.line = element_blank(),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F4b.PlOT<-F4b.PlOT+scale_x_discrete(position = "top")+xlab(NULL)+ylab(NULL)+guides(size=FALSE)
gene.scale4 = 1/14


F5b.PlOT<-ggplot()+theme_classic()
F5b.PlOT<-F5b.PlOT+geom_tile(data = data5,aes(x=reorder(xxx,zzz),y=yyy,fill=his),color='black',width=1,height=1,size=0.5,stat='identity')
F5b.PlOT<-F5b.PlOT+scale_fill_manual(name=NULL,values=c(Male='#b3cde3',Female='#fccde5'),labels=c(Male='Male',Female='Female'),guide = guide_legend(override.aes=list(size=1),nrow=1),na.translate = F)
F5b.PlOT<-F5b.PlOT+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(0,0.1,0.1,0.1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.4,'cm'),legend.key.height=unit(0.4,'cm'),legend.position='none',axis.ticks.x = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=12,face='bold.italic'),axis.ticks.y = element_blank(), axis.text.y=element_blank(),legend.title=element_blank(),axis.line = element_blank(),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F5b.PlOT<-F5b.PlOT+scale_x_discrete(position = "top")+xlab(NULL)+ylab(NULL)+guides(size=FALSE)
gene.scale5 = 1/14



#P9
sample <- data.frame(Sample_features)
sample <- sample[which(sample$patient=='P9'),]
pat <- unique(sample$patient)
case <- unique(sample$case)
#raw data preprocessing

somatic <- read.delim('./Mut/Annotated_Coding.txt',header = T)
somatic <- somatic[which(somatic$CaseID %in% pat),]
FL <- somatic[which(somatic$T1_freq >= 10 & somatic$HighQread_T1 >= 2),]
FL$ID <- "P9.FL"
DHL <- somatic[which(somatic$T3_freq >= 10 & somatic$HighQread_T3 >= 2),]
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
    temp.P <- savi$Effect[which(savi$case == as.character(case[i]) & savi$Gene_Name == selGene[j] & savi$Impact!='LOW' & savi$Impact!='MODIFIER')]
    if( any(grepl(pattern = 'stop',x = temp.P,fixed=T) | grepl(pattern = 'nonsense',x = temp.P,fixed=T) ) ){
      gene.Matrix[i,j] <- 'StopGained'
    }
    else if( any(grepl(pattern = 'frameshift_variant',x = temp.P,fixed=T)) ){
      gene.Matrix[i,j] <- 'Frameshift'
    }
    else if( any( grepl(pattern = 'splice',x = temp.P,fixed=T) | grepl(pattern = 'splice_acceptor',x = temp.P,fixed=T) ) ){
      gene.Matrix[i,j] <- 'SpliceSite'
    }
    else if( any( grepl(pattern = 'inframe',x = temp.P,fixed=T) ) ){
      gene.Matrix[i,j] <- 'InframeIndel'
    }
    else if( any( grepl(pattern = 'missense_variant', x = temp.P,fixed=T)) ){
      gene.Matrix[i,j] <- 'Missense'
    }
    else if( any( grepl(pattern = 'start', x = temp.P,fixed=T)) ){
      gene.Matrix[i,j] <- 'StartLost'
    }
    else{}
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

FocalCNV <- read.delim('./Mut/CNV_matrix.txt',header = T)
FocalCNV <- FocalCNV[which(as.character(FocalCNV$Gene) %in% interestedGene),]
xxx <- FocalCNV[,which(colnames(FocalCNV) %in% case)]
rownames(xxx) <- FocalCNV$Gene
NEW_CNV <- melt(t(xxx))
colnames(NEW_CNV) <- c('ID','gene','value')

theAmp <- NEW_CNV[which(NEW_CNV$value >= 0.5 ),]
theDel <- NEW_CNV[which(NEW_CNV$value <= -0.5 ),]


new_table <- merge(plot_1b.data,sample,1,1)

SV <- read.delim('./Mut/SV.txt',header = T)
theSV <- SV[which(SV$ID %in% case),]

LOH <- read.delim('./Mut/LOH_matrix.txt',header = T)
LOH <- LOH[which(as.character(LOH$Gene) %in% interestedGene),]
xxx <- LOH[,which(colnames(LOH) %in% case)]
rownames(xxx) <- LOH$Gene
NEW_LOH <- melt(t(xxx))
colnames(NEW_LOH) <- c('ID','gene','value')
theLOH <- NEW_LOH[which(NEW_LOH$value == 1 ),]

#orderID<-c(1:nrow(plot_1b.data))
F1b.PloT<-ggplot()+theme_classic()
F1b.PloT<-F1b.PloT+geom_tile(data = new_table,aes(x=reorder(ID,order),y=gene,fill=value),color='black',width=1,height=1,size=0.5,stat='identity')
F1b.PloT<-F1b.PloT+geom_point(data = theSV,aes(x=ID,y=gene),shape=23,fill='black',size=3)
F1b.PloT<-F1b.PloT+geom_point(data = theAmp,aes(x=ID,y=gene),shape=24,fill='red',size=3)
F1b.PloT<-F1b.PloT+geom_point(data = theDel,aes(x=ID,y=gene),shape=25,fill='blue',size=3)+geom_text(data = theLOH,aes(x=ID,y=gene),label="L",color='white',size=4)
F1b.PloT<-F1b.PloT+scale_fill_manual(name=NULL,values=c(N='white',Missense='#66c2a5',InframeIndel='#a6d854',Frameshift='#fc8d62',SpliceSite=gg_color_hue(6)[5],StopGained='#df65b0',StartLost='#762a83'),labels=c(N='None',Missense='Missense',Inframe='In frame indel',Frameshift='Frame shift', SpliceSite='Splice site',StopGained='Stop gained',StartLost='Start lost'),guide = guide_legend(override.aes=list(size=1),nrow=7),na.translate = F)
F1b.PloT<-F1b.PloT+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(0.1,1,1,0.1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.68,'cm'),legend.position='none',axis.ticks.x = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=16,face='bold.italic'),axis.ticks.y = element_blank(), axis.text.y=element_blank(),axis.line = element_blank(),
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


F4b.PloT<-ggplot()+theme_classic()
F4b.PloT<-F4b.PloT+geom_tile(data = data4,aes(x=reorder(xxx,zzz),y=yyy,fill=his),color='black',width=1,height=1,size=0.5,stat='identity')
F4b.PloT<-F4b.PloT+scale_fill_manual(name=NULL,values=c(High=gg_color_hue(6)[1] ,Low=gg_color_hue(6)[3]),guide = guide_legend(override.aes=list(size=1),nrow=1),na.translate = F)
F4b.PloT<-F4b.PloT+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(1,1,0,0.1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.4,'cm'),legend.key.height=unit(0.4,'cm'),legend.position='none',axis.ticks.x = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=12,face='bold.italic'),axis.ticks.y = element_blank(), axis.text.y=element_blank(),legend.title=element_blank(),axis.line = element_blank(),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F4b.PloT<-F4b.PloT+scale_x_discrete(position = "top")+xlab(NULL)+ylab(NULL)+guides(size=FALSE)
gene.scale4 = 1/14


F5b.PloT<-ggplot()+theme_classic()
F5b.PloT<-F5b.PloT+geom_tile(data = data5,aes(x=reorder(xxx,zzz),y=yyy,fill=his),color='black',width=1,height=1,size=0.5,stat='identity')
F5b.PloT<-F5b.PloT+scale_fill_manual(name=NULL,values=c(Male='#b3cde3',Female='#fccde5'),labels=c(Male='Male',Female='Female'),guide = guide_legend(override.aes=list(size=1),nrow=1),na.translate = F)
F5b.PloT<-F5b.PloT+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(0,1,0.1,0.1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.4,'cm'),legend.key.height=unit(0.4,'cm'),legend.position='none',axis.ticks.x = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=12,face='bold.italic'),axis.ticks.y = element_blank(), axis.text.y=element_blank(),legend.title=element_blank(),axis.line = element_blank(),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F5b.PloT<-F5b.PloT+scale_x_discrete(position = "top")+xlab(NULL)+ylab(NULL)+guides(size=FALSE)
gene.scale5 = 1/14



figure_1<-cbind(
          rbind(ggplotGrob(F4b.ploT),ggplotGrob(F5b.ploT),ggplotGrob(F1b.ploT)),
          rbind(ggplotGrob(F4b.plOT),ggplotGrob(F5b.plOT),ggplotGrob(F1b.plOT)),
          rbind(ggplotGrob(F4b.PLOT),ggplotGrob(F5b.PLOT),ggplotGrob(F1b.PLOT)),
          rbind(ggplotGrob(F4b.PLOt),ggplotGrob(F5b.PLOt),ggplotGrob(F1b.PLOt)),
          rbind(ggplotGrob(F4b.PLoT),ggplotGrob(F5b.PLoT),ggplotGrob(F1b.PLoT)),
          rbind(ggplotGrob(F4b.PlOT),ggplotGrob(F5b.PlOT),ggplotGrob(F1b.PlOT)),
          rbind(ggplotGrob(F4b.pLot),ggplotGrob(F5b.pLot),ggplotGrob(F1b.pLot)),
          rbind(ggplotGrob(F4b.PloT),ggplotGrob(F5b.PloT),ggplotGrob(F1b.PloT)),
          size="last")

panels <- figure_1$layout$t[grep("panel", figure_1$layout$name)]

figure_1$widths[5]  <- unit(2/4,'null')
figure_1$widths[14] <- unit(2/4,'null')
figure_1$widths[23] <- unit(2/4,'null')
figure_1$widths[32] <- unit(2/4,'null')
figure_1$widths[41] <- unit(2/4,'null')
figure_1$widths[50] <- unit(2/4,'null')
figure_1$widths[59] <- unit(3/4,'null')
figure_1$widths[68] <- unit(2/4,'null')


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




#grid.draw(figure_1)
ggsave(file="landscape_coding.pdf", plot=figure_1,bg = 'white', width = 16, height = 24, units = 'cm', dpi = 600)



