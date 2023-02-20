# WatchDHL
library(rstudioapi)

# Getting the path of your current open file
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))
print( getwd() )

library(ggplot2)
library(gridExtra)
library(grid)
library(reshape2) 


gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}


interestedGene <- c("IGH","BCL2","IGL","IMMP2L","FHIT","BCL6","BCL7A","DIAPH2","CIITA","BACH2",
                    "ROBO2","IGK","PAX5","RHOH","DMD","GRID2","DTX1","LPP","BCL6-SE","IRF8",
                    "ST6GAL1","CXCR4","BTG2","ZCCHC7","MYO16","S1PR2","PALM2AKAP2","PIM1","NCOA3",
                    "FOXO1","SGK1","SEL1L3","OSBPL10","MYC","BLK","WNK2","MSH6-SE","MIR3681HG",
                    "KLHL1","ETS1","CFAP251","SUGCT","MIR142","TCL1A","LHPP","UHRF2","RFTN1",
                    "MYO1E","IL4R","BIRC3","VMP1","RUBCNL","MCPH1","POU2AF1","LTB","IKZF3","CD83","MYBL1","BTG1","BMP7")


Sample_features <- read.table("Genetic_profiles_new.txt",sep = "\t",header = T)
Sample_features[is.na(Sample_features)] <- "f_NA"
sample <- data.frame(Sample_features)
sample <- sample[which(sample$patient=='P8'),]
pat <- unique(sample$patient)
case <- unique(sample$case)
#raw data preprocessing

somatic <- read.delim('./Mut/annotated.frequent.txt',header = T)
somatic <- somatic[which(somatic$CaseID %in% pat),]
FL.a <- somatic[which(somatic$T1_freq >= 10 & somatic$altdepth_T1 >= 5),]
FL.a$ID <- "P8.FL.a"
FL.b <- somatic[which(as.numeric(as.character(somatic$T2_freq)) >= 20 & as.numeric(as.character(somatic$altdepth_T2)) >= 5),]
FL.b$ID <- "P8.FL.b"
DHL <- somatic[which(somatic$T3_freq >= 20 & somatic$altdepth_T3 >= 5),]
DHL$ID <- "P8.DHL"

savi.table <- rbind(FL.a,FL.b,DHL)

savi.Sel <- savi.table[which(savi.table$ID %in% case),]

savi <- data.frame(savi.Sel$chr,savi.Sel$ID,savi.Sel$GeneID,savi.Sel$Effect_Impact,savi.Sel$Effect)
colnames(savi) <- c("chr","case","GeneID","Impact","Effect")

#analysis for selected genes
selGene <- interestedGene

gene.Matrix <- rep('N',length(case))
for( i in 2:length(selGene)){
  gene.Matrix <- cbind(gene.Matrix,rep('N',length(case)))
}

for(i in 1:length(case)){
  for(j in 1:length(selGene)){
    temp.P <- savi$Effect[which(savi$case == as.character(case[i]) & savi$GeneID == selGene[j] & savi$Impact!='HIGH' & savi$Impact!='MODERATE')]
    gene.Matrix[i,j] <- length(temp.P)
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


xxx<-as.integer(as.character(plot_1b.data$value))
mmm<-as.integer(as.character(plot_1b.data$value))

mmm[which(xxx == 0)] <-   "a0"
mmm[which(xxx == 1) ] <-  "b1"
mmm[which(xxx >= 2) ] <-  "c2"
mmm[which(xxx >= 4) ] <-  "c5"
mmm[which(xxx >= 8) ] <- "d10"
mmm[which(xxx >= 16) ] <- "e20"
mmm[which(xxx >= 32) ] <- "f30"
mmm[which(xxx >= 64) ] <- "g40"
plot_1b.data <- cbind(plot_1b.data,mmm)
new_table <- merge(plot_1b.data,sample,1,1)


#orderID<-c(1:nrow(plot_1b.data))
F1b.pLot<-ggplot()+theme_classic()
F1b.pLot<-F1b.pLot+geom_tile(data = new_table,aes(x=reorder(ID,-order),y=gene,fill=mmm),color='black',width=1,height=1,size=0.5,stat='identity')

F1b.pLot<-F1b.pLot+scale_fill_manual(name="Mutation Number",values =c(a0='white',b1='#fee0d2',c2='#fcbba1',c5='#fb6a4a',d10='#ef3b2c',e20='#cb181d',f30='#a50f15',g40='#67000d'),labels=c(a0='0',b1='1',c2='>=2',c5='>=5',d10='>=10',e20='>=20',f30='>=30',g40='>=40'))
F1b.pLot<-F1b.pLot+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(0,1,0,1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.68,'cm'),legend.position='none',axis.ticks.x = element_blank(),axis.line.x = element_blank(),axis.line.y = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=16,face='bold.italic'),axis.ticks.y = element_blank(), axis.text.y=element_text(size=12,vjust=0.5,hjust=1,face='plain',color='black'),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F1b.pLot<-F1b.pLot+scale_x_discrete(position = "bottom")+xlab(NULL)+ylab(NULL)+coord_flip()
F1b.pLot
gene.scale1b = 1 + (ncol(mutation_gene_table) - 14)/14



yyy <- rep("Grade",length(sample$case))
data4<-data.frame(sample$case,yyy,sample$order,sample$Grade)
colnames(data4)<-c("xxx","yyy","zzz","his")



F4b.pLot<-ggplot()+theme_classic()
F4b.pLot<-F4b.pLot+geom_tile(data = data4,aes(x=reorder(xxx,zzz),y=yyy,fill=his),color='black',width=1,height=1,size=0.5,stat='identity')
F4b.pLot<-F4b.pLot+scale_fill_manual(name=NULL,values=c(High=gg_color_hue(6)[1] ,Low=gg_color_hue(6)[3]),guide = guide_legend(override.aes=list(size=1),nrow=1),na.translate = F)
F4b.pLot<-F4b.pLot+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(1,0.1,0.1,0.1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.4,'cm'),legend.key.height=unit(0.4,'cm'),legend.position='none',axis.ticks.x = element_blank(),axis.line.x = element_blank(),axis.line.y = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=12,face='bold.italic'),axis.ticks.y = element_blank(), axis.text.y=element_blank(),legend.title=element_blank(),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F4b.pLot<-F4b.pLot+scale_x_discrete(position = "top")+xlab(NULL)+ylab(NULL)+guides(size=FALSE)
gene.scale4 = 1/14


###P1
sample <- data.frame(Sample_features)
sample <- sample[which(sample$patient=='P1'),]
pat <- unique(sample$patient)
case <- unique(sample$case)
#raw data preprocessing

somatic <- read.delim('./Mut/annotated.frequent.txt',header = T)
somatic <- somatic[which(somatic$CaseID %in% pat),]
FL <- somatic[which(somatic$T1_freq >= 20 & somatic$altdepth_T1 >= 5),]
FL$ID <- "P1.FL"
DHL <- somatic[which(somatic$T3_freq >= 20 & somatic$altdepth_T3 >= 5),]
DHL$ID <- "P1.DHL"

savi.table <- rbind(FL,DHL)

savi.Sel <- savi.table[which(savi.table$ID %in% case),]

savi <- data.frame(savi.Sel$chr,savi.Sel$ID,savi.Sel$GeneID,savi.Sel$Effect_Impact,savi.Sel$Effect)
colnames(savi) <- c("chr","case","GeneID","Impact","Effect")

#analysis for selected genes
selGene <- interestedGene

gene.Matrix <- rep('N',length(case))
for( i in 2:length(selGene)){
  gene.Matrix <- cbind(gene.Matrix,rep('N',length(case)))
}

for(i in 1:length(case)){
  for(j in 1:length(selGene)){
    temp.P <- savi$Effect[which(savi$case == as.character(case[i]) & savi$GeneID == selGene[j] & savi$Impact!='HIGH' & savi$Impact!='MODERATE')]
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
mmm[which(xxx >= 4) ] <-  "c5"
mmm[which(xxx >= 8) ] <- "d10"
mmm[which(xxx >= 16) ] <- "e20"
mmm[which(xxx >= 32) ] <- "f30"
mmm[which(xxx >= 64) ] <- "g40"
plot_1b.data <- cbind(plot_1b.data,mmm)
new_table <- merge(plot_1b.data,sample,1,1)

#orderID<-c(1:nrow(plot_1b.data))
F1b.ploT<-ggplot()+theme_classic()
F1b.ploT<-F1b.ploT+geom_tile(data = new_table,aes(x=reorder(ID,-order),y=gene,fill=mmm),color='black',width=1,height=1,size=0.5,stat='identity')

F1b.ploT<-F1b.ploT+scale_fill_manual(name="Mutation Number",values =c(a0='white',b1='#fee0d2',c2='#fcbba1',c5='#fb6a4a',d10='#ef3b2c',e20='#cb181d',f30='#a50f15',g40='#67000d'),labels=c(a0='0',b1='1',c2='>=2',c5='>=5',d10='>=10',e20='>=20',f30='>=30',g40='>=40'))
F1b.ploT<-F1b.ploT+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(0,1,0,1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.68,'cm'),legend.position='none',axis.ticks.x = element_blank(),axis.line.x = element_blank(),axis.line.y = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=16,face='bold.italic'),axis.text.y=element_text(size=12,vjust=0.5,hjust=1,face='plain',color='black'),axis.line = element_blank(),axis.ticks.y = element_blank(),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F1b.ploT<-F1b.ploT+scale_x_discrete(position = "bottom")+xlab(NULL)+ylab(NULL)+coord_flip()
F1b.ploT
gene.scale1b = 1 + (ncol(mutation_gene_table) - 14)/14



yyy <- rep("Grade",length(sample$case))
data4<-data.frame(sample$case,yyy,sample$order,sample$Grade)
colnames(data4)<-c("xxx","yyy","zzz","his")



F4b.ploT<-ggplot()+theme_classic()
F4b.ploT<-F4b.ploT+geom_tile(data = data4,aes(x=reorder(xxx,zzz),y=yyy,fill=his),color='black',width=1,height=1,size=0.5,stat='identity')
F4b.ploT<-F4b.ploT+scale_fill_manual(name=NULL,values=c(High=gg_color_hue(6)[1] ,Low=gg_color_hue(6)[3]),guide = guide_legend(override.aes=list(size=1),nrow=1),na.translate = F)
F4b.ploT<-F4b.ploT+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(1,0.1,0.1,0.1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.4,'cm'),legend.key.height=unit(0.4,'cm'),legend.position='none',axis.ticks.x = element_blank(),axis.line.x = element_blank(),axis.line.y = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=12,face='bold.italic'),axis.text.y=element_text(size=14,vjust=0.5,hjust=1,face='bold.italic',color='black'),legend.title=element_blank(),axis.ticks.y = element_blank(),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F4b.ploT<-F4b.ploT+scale_x_discrete(position = "top")+xlab(NULL)+ylab(NULL)+guides(size=FALSE)
gene.scale4 = 1/14



##P2
sample <- data.frame(Sample_features)
sample <- sample[which(sample$patient=='P2'),]
pat <- unique(sample$patient)
case <- unique(sample$case)
#raw data preprocessing

somatic <- read.delim('./Mut/annotated.frequent.txt',header = T)
somatic <- somatic[which(somatic$CaseID %in% pat),]
FL <- somatic[which(somatic$T1_freq >= 10 & somatic$altdepth_T1 >= 5),]
FL$ID <- "P2.FL"
DHL <- somatic[which(somatic$T3_freq >= 20 & somatic$altdepth_T3 >= 5),]
DHL$ID <- "P2.DHL"

savi.table <- rbind(FL,DHL)

savi.Sel <- savi.table[which(savi.table$ID %in% case),]

savi <- data.frame(savi.Sel$chr,savi.Sel$ID,savi.Sel$GeneID,savi.Sel$Effect_Impact,savi.Sel$Effect)
colnames(savi) <- c("chr","case","GeneID","Impact","Effect")

#analysis for selected genes
selGene <- interestedGene

gene.Matrix <- rep('N',length(case))
for( i in 2:length(selGene)){
  gene.Matrix <- cbind(gene.Matrix,rep('N',length(case)))
}

for(i in 1:length(case)){
  for(j in 1:length(selGene)){
    temp.P <- savi$Effect[which(savi$case == as.character(case[i]) & savi$GeneID == selGene[j] & savi$Impact!='HIGH' & savi$Impact!='MODERATE')]
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
mmm[which(xxx >= 4) ] <-  "c5"
mmm[which(xxx >= 8) ] <- "d10"
mmm[which(xxx >= 16) ] <- "e20"
mmm[which(xxx >= 32) ] <- "f30"
mmm[which(xxx >= 64) ] <- "g40"
plot_1b.data <- cbind(plot_1b.data,mmm)
new_table <- merge(plot_1b.data,sample,1,1)

#orderID<-c(1:nrow(plot_1b.data))
F1b.plOT<-ggplot()+theme_classic()
F1b.plOT<-F1b.plOT+geom_tile(data = new_table,aes(x=reorder(ID,-order),y=gene,fill=mmm),color='black',width=1,height=1,size=0.5,stat='identity')

F1b.plOT<-F1b.plOT+scale_fill_manual(name="Mutation Number",values =c(a0='white',b1='#fee0d2',c2='#fcbba1',c5='#fb6a4a',d10='#ef3b2c',e20='#cb181d',f30='#a50f15',g40='#67000d'),labels=c(a0='0',b1='1',c2='>=2',c5='>=5',d10='>=10',e20='>=20',f30='>=30',g40='>=40'))
F1b.plOT<-F1b.plOT+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(0,1,0,1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.68,'cm'),legend.position='none',axis.ticks.x = element_blank(),axis.line.x = element_blank(),axis.line.y = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=16,face='bold.italic'),axis.ticks.y = element_blank(), axis.text.y=element_text(size=12,vjust=0.5,hjust=1,face='plain',color='black'),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F1b.plOT<-F1b.plOT+scale_x_discrete(position = "bottom")+xlab(NULL)+ylab(NULL)+coord_flip()
F1b.plOT
gene.scale1b = 1 + (ncol(mutation_gene_table) - 14)/14



yyy <- rep("Grade",length(sample$case))
data4<-data.frame(sample$case,yyy,sample$order,sample$Grade)
colnames(data4)<-c("xxx","yyy","zzz","his")



F4b.plOT<-ggplot()+theme_classic()
F4b.plOT<-F4b.plOT+geom_tile(data = data4,aes(x=reorder(xxx,zzz),y=yyy,fill=his),color='black',width=1,height=1,size=0.5,stat='identity')
F4b.plOT<-F4b.plOT+scale_fill_manual(name=NULL,values=c(High=gg_color_hue(6)[1] ,Low=gg_color_hue(6)[3]),guide = guide_legend(override.aes=list(size=1),nrow=1),na.translate = F)
F4b.plOT<-F4b.plOT+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(1,0.1,0.1,0.1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.4,'cm'),legend.key.height=unit(0.4,'cm'),legend.position='none',axis.ticks.x = element_blank(),axis.line.x = element_blank(),axis.line.y = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=12,face='bold.italic'),axis.ticks.y = element_blank(), axis.text.y=element_blank(),legend.title=element_blank(),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F4b.plOT<-F4b.plOT+scale_x_discrete(position = "top")+xlab(NULL)+ylab(NULL)+guides(size=FALSE)
gene.scale4 = 1/14



##P6
sample <- data.frame(Sample_features)
sample <- sample[which(sample$patient=='P6'),]
pat <- unique(sample$patient)
case <- unique(sample$case)
#raw data preprocessing

somatic <- read.delim('./Mut/annotated.frequent.txt',header = T)
somatic <- somatic[which(somatic$CaseID %in% pat),]
FL <- somatic[which(somatic$T1_freq >= 20 & somatic$altdepth_T1 >= 5),]
FL$ID <- "P6.FL"
DHL <- somatic[which(somatic$T3_freq >= 20 & somatic$altdepth_T3 >= 5),]
DHL$ID <- "P6.DHL"

savi.table <- rbind(FL,DHL)

savi.Sel <- savi.table[which(savi.table$ID %in% case),]

savi <- data.frame(savi.Sel$chr,savi.Sel$ID,savi.Sel$GeneID,savi.Sel$Effect_Impact,savi.Sel$Effect)
colnames(savi) <- c("chr","case","GeneID","Impact","Effect")

#analysis for selected genes
selGene <- interestedGene

gene.Matrix <- rep('N',length(case))
for( i in 2:length(selGene)){
  gene.Matrix <- cbind(gene.Matrix,rep('N',length(case)))
}

for(i in 1:length(case)){
  for(j in 1:length(selGene)){
    temp.P <- savi$Effect[which(savi$case == as.character(case[i]) & savi$GeneID == selGene[j] & savi$Impact!='HIGH' & savi$Impact!='MODERATE')]
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
mmm[which(xxx >= 4) ] <-  "c5"
mmm[which(xxx >= 8) ] <- "d10"
mmm[which(xxx >= 16) ] <- "e20"
mmm[which(xxx >= 32) ] <- "f30"
mmm[which(xxx >= 64) ] <- "g40"
plot_1b.data <- cbind(plot_1b.data,mmm)
new_table <- merge(plot_1b.data,sample,1,1)

#orderID<-c(1:nrow(plot_1b.data))
F1b.PLoT<-ggplot()+theme_classic()
F1b.PLoT<-F1b.PLoT+geom_tile(data = new_table,aes(x=reorder(ID,-order),y=gene,fill=mmm),color='black',width=1,height=1,size=0.5,stat='identity')

F1b.PLoT<-F1b.PLoT+scale_fill_manual(name="Mutation Number",values =c(a0='white',b1='#fee0d2',c2='#fcbba1',c5='#fb6a4a',d10='#ef3b2c',e20='#cb181d',f30='#a50f15',g40='#67000d'),labels=c(a0='0',b1='1',c2='>=2',c5='>=5',d10='>=10',e20='>=20',f30='>=30',g40='>=40'))
F1b.PLoT<-F1b.PLoT+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(0,1,0,1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.68,'cm'),legend.position='none',axis.ticks.x = element_blank(),axis.line.x = element_blank(),axis.line.y = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=16,face='bold.italic'),axis.ticks.y = element_blank(), axis.text.y=element_text(size=12,vjust=0.5,hjust=1,face='plain',color='black'),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F1b.PLoT<-F1b.PLoT+scale_x_discrete(position = "bottom")+xlab(NULL)+ylab(NULL)+coord_flip()
F1b.PLoT
gene.scale1b = 1 + (ncol(mutation_gene_table) - 14)/14



yyy <- rep("Grade",length(sample$case))
data4<-data.frame(sample$case,yyy,sample$order,sample$Grade)
colnames(data4)<-c("xxx","yyy","zzz","his")

F4b.PLoT<-ggplot()+theme_classic()
F4b.PLoT<-F4b.PLoT+geom_tile(data = data4,aes(x=reorder(xxx,zzz),y=yyy,fill=his),color='black',width=1,height=1,size=0.5,stat='identity')
F4b.PLoT<-F4b.PLoT+scale_fill_manual(name=NULL,values=c(High=gg_color_hue(6)[1] ,Low=gg_color_hue(6)[3]),guide = guide_legend(override.aes=list(size=1),nrow=1),na.translate = F)
F4b.PLoT<-F4b.PLoT+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(1,0.1,0.1,0.1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.4,'cm'),legend.key.height=unit(0.4,'cm'),legend.position='none',axis.ticks.x = element_blank(),axis.line.x = element_blank(),axis.line.y = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=12,face='bold.italic'),axis.ticks.y = element_blank(), axis.text.y=element_blank(),legend.title=element_blank(),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F4b.PLoT<-F4b.PLoT+scale_x_discrete(position = "top")+xlab(NULL)+ylab(NULL)+guides(size=FALSE)
gene.scale4 = 1/14


###P7
sample <- data.frame(Sample_features)
sample <- sample[which(sample$patient=='P7'),]
pat <- unique(sample$patient)
case <- unique(sample$case)
#raw data preprocessing

somatic <- read.delim('./Mut/annotated.frequent.txt',header = T)
somatic <- somatic[which(somatic$CaseID %in% pat),]
FL <- somatic[which(somatic$T1_freq >= 20 & somatic$altdepth_T1 >= 5),]
FL$ID <- "P7.FL"
DHL <- somatic[which(somatic$T3_freq >= 20 & somatic$altdepth_T3 >= 5),]
DHL$ID <- "P7.DHL"

savi.table <- rbind(FL,DHL)

savi.Sel <- savi.table[which(savi.table$ID %in% case),]

savi <- data.frame(savi.Sel$chr,savi.Sel$ID,savi.Sel$GeneID,savi.Sel$Effect_Impact,savi.Sel$Effect)
colnames(savi) <- c("chr","case","GeneID","Impact","Effect")

#analysis for selected genes
selGene <- interestedGene

gene.Matrix <- rep('N',length(case))
for( i in 2:length(selGene)){
  gene.Matrix <- cbind(gene.Matrix,rep('N',length(case)))
}

for(i in 1:length(case)){
  for(j in 1:length(selGene)){
    temp.P <- savi$Effect[which(savi$case == as.character(case[i]) & savi$GeneID == selGene[j] & savi$Impact!='HIGH' & savi$Impact!='MODERATE')]
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
mmm[which(xxx >= 4) ] <-  "c5"
mmm[which(xxx >= 8) ] <- "d10"
mmm[which(xxx >= 16) ] <- "e20"
mmm[which(xxx >= 32) ] <- "f30"
mmm[which(xxx >= 64) ] <- "g40"
plot_1b.data <- cbind(plot_1b.data,mmm)
new_table <- merge(plot_1b.data,sample,1,1)

#orderID<-c(1:nrow(plot_1b.data))
F1b.PlOT<-ggplot()+theme_classic()
F1b.PlOT<-F1b.PlOT+geom_tile(data = new_table,aes(x=reorder(ID,-order),y=gene,fill=mmm),color='black',width=1,height=1,size=0.5,stat='identity')

F1b.PlOT<-F1b.PlOT+scale_fill_manual(name="Mutation Number",values =c(a0='white',b1='#fee0d2',c2='#fcbba1',c5='#fb6a4a',d10='#ef3b2c',e20='#cb181d',f30='#a50f15',g40='#67000d'),labels=c(a0='0',b1='1',c2='>=2',c5='>=5',d10='>=10',e20='>=20',f30='>=30',g40='>=40'))
F1b.PlOT<-F1b.PlOT+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(0,1,0,1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.68,'cm'),legend.position='none',axis.ticks.x = element_blank(),axis.line.x = element_blank(),axis.line.y = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=16,face='bold.italic'),axis.ticks.y = element_blank(), axis.text.y=element_text(size=12,vjust=0.5,hjust=1,face='plain',color='black'),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F1b.PlOT<-F1b.PlOT+scale_x_discrete(position = "bottom")+xlab(NULL)+ylab(NULL)+coord_flip()
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
F4b.PlOT<-F4b.PlOT+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(1,0.1,0.1,0.1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.4,'cm'),legend.key.height=unit(0.4,'cm'),legend.position='none',axis.ticks.x = element_blank(),axis.line.x = element_blank(),axis.line.y = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=12,face='bold.italic'),axis.ticks.y = element_blank(), axis.text.y=element_blank(),legend.title=element_blank(),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F4b.PlOT<-F4b.PlOT+scale_x_discrete(position = "top")+xlab(NULL)+ylab(NULL)+guides(size=FALSE)
gene.scale4 = 1/14



#P9
sample <- data.frame(Sample_features)
sample <- sample[which(sample$patient=='P9'),]
pat <- unique(sample$patient)
case <- unique(sample$case)
#raw data preprocessing

somatic <- read.delim('./Mut/annotated.frequent.txt',header = T)
somatic <- somatic[which(somatic$CaseID %in% pat),]
FL <- somatic[which(somatic$T1_freq >= 20 & somatic$altdepth_T1 >= 5),]
FL$ID <- "P9.FL"
DHL <- somatic[which(somatic$T3_freq >= 20 & somatic$altdepth_T3 >= 5),]
DHL$ID <- "P9.DHL"

savi.table <- rbind(FL,DHL)

savi.Sel <- savi.table[which(savi.table$ID %in% case),]

savi <- data.frame(savi.Sel$chr,savi.Sel$ID,savi.Sel$GeneID,savi.Sel$Effect_Impact,savi.Sel$Effect)
colnames(savi) <- c("chr","case","GeneID","Impact","Effect")

#analysis for selected genes
selGene <- interestedGene

gene.Matrix <- rep('N',length(case))
for( i in 2:length(selGene)){
  gene.Matrix <- cbind(gene.Matrix,rep('N',length(case)))
}

for(i in 1:length(case)){
  for(j in 1:length(selGene)){
    temp.P <- savi$Effect[which(savi$case == as.character(case[i]) & savi$GeneID == selGene[j] & savi$Impact!='HIGH' & savi$Impact!='MODERATE')]
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
mmm[which(xxx >= 4) ] <-  "c5"
mmm[which(xxx >= 8) ] <- "d10"
mmm[which(xxx >= 16) ] <- "e20"
mmm[which(xxx >= 32) ] <- "f30"
mmm[which(xxx >= 64) ] <- "g40"
plot_1b.data <- cbind(plot_1b.data,mmm)
new_table <- merge(plot_1b.data,sample,1,1)

#orderID<-c(1:nrow(plot_1b.data))
F1b.PloT<-ggplot()+theme_classic()
F1b.PloT<-F1b.PloT+geom_tile(data = new_table,aes(x=reorder(ID,-order),y=gene,fill=mmm),color='black',width=1,height=1,size=0.5,stat='identity')

F1b.PloT<-F1b.PloT+scale_fill_manual(name="Mutation Number",values =c(a0='white',b1='#fee0d2',c2='#fcbba1',c5='#fb6a4a',d10='#ef3b2c',e20='#cb181d',f30='#a50f15',g40='#67000d'),labels=c(a0='0',b1='1',c2='>=2',c5='>=5',d10='>=10',e20='>=20',f30='>=30',g40='>=40'))
F1b.PloT<-F1b.PloT+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(0,1,1,1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.68,'cm'),legend.position='none',axis.line.x = element_blank(),axis.line.y = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=16,face='bold.italic'),axis.ticks.y = element_blank(), axis.text.y=element_text(size=12,vjust=0.5,hjust=1,face='plain',color='black'),
                         axis.text.x=element_text(size=12,angle=90,vjust=0.5,hjust=1,face='italic',color='black'),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F1b.PloT<-F1b.PloT+scale_x_discrete(position = "bottom")+xlab(NULL)+ylab(NULL)+coord_flip()
F1b.PloT
gene.scale1b = 1 + (ncol(mutation_gene_table) - 14)/14



yyy <- rep("Grade",length(sample$case))
data4<-data.frame(sample$case,yyy,sample$order,sample$Grade)
colnames(data4)<-c("xxx","yyy","zzz","his")

F4b.PloT<-ggplot()+theme_classic()
F4b.PloT<-F4b.PloT+geom_tile(data = data4,aes(x=reorder(xxx,zzz),y=yyy,fill=his),color='black',width=1,height=1,size=0.5,stat='identity')
F4b.PloT<-F4b.PloT+scale_fill_manual(name=NULL,values=c(High=gg_color_hue(6)[1] ,Low=gg_color_hue(6)[3]),guide = guide_legend(override.aes=list(size=1),nrow=1),na.translate = F)
F4b.PloT<-F4b.PloT+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(0.1,1,1,1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.4,'cm'),legend.key.height=unit(0.4,'cm'),legend.position='none',axis.ticks.x = element_blank(),axis.line.x = element_blank(),axis.line.y = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=12,face='bold.italic'),axis.ticks.y = element_blank(), axis.text.y=element_blank(),legend.title=element_blank(),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F4b.PloT<-F4b.PloT+scale_x_discrete(position = "top")+xlab(NULL)+ylab(NULL)+guides(size=FALSE)
gene.scale4 = 1/14







gene.table<-read.table('./Mut/number_summary.txt',header=T,sep="\t")

data_plot  <- data.frame(c(as.character(gene.table$gene),as.character(gene.table$gene)),c(gene.table$rank,gene.table$rank),c(gene.table$FLNum,gene.table$DHLNum),c(gene.table$TotalNum,gene.table$TotalNum),c(rep("aFL",nrow(gene.table)),rep("DHL",nrow(gene.table))))
colnames(data_plot) <- c("gene","ranks","Num","Total","Label")
ffff_plot<-ggplot()+theme_classic()
ffff_plot<-ffff_plot+geom_bar(data=data_plot ,aes(x=reorder(gene,ranks),y=Num,fill=Label),width=0.7 ,color = "black",stat='identity',position=position_dodge())
ffff_plot<-ffff_plot+scale_y_continuous('Number of mutations',expand=c(0,0),limits = c(0, 350),breaks=seq(0,400,100))
ffff_plot<-ffff_plot+theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(1,1,0.1,1),'lines'),
                           plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),text=element_text(size=15,vjust=0.5,hjust=0.5,face='bold'),
                           legend.key.width=unit(0.7,'cm'),legend.key.height=unit(0.7,'cm'),legend.position="none",axis.text.x=element_blank(),
                           legend.text=element_text(size=15,face='bold.italic'),axis.text.y=element_text(size=12,face='plain',color='black'),axis.title.y=element_blank())
ffff_plot<-ffff_plot+ggtitle(NULL)+xlab(NULL)+scale_x_discrete(position = "bottom")+ylab(NULL)
ffff_plot<-ffff_plot+scale_fill_manual(name=NULL,values=c(DHL=gg_color_hue(6)[1] ,aFL=gg_color_hue(6)[3]))

ffff_plot




figure_1<-rbind(ggplotGrob(ffff_plot),ggplotGrob(F1b.ploT),ggplotGrob(F1b.plOT),ggplotGrob(F1b.PLoT),ggplotGrob(F1b.PlOT),ggplotGrob(F1b.pLot),ggplotGrob(F1b.PloT), size="last")

panels <- figure_1$layout$t[grep("panel", figure_1$layout$name)]


figure_1$heights[panels][1] <- unit(10/5,'null')
figure_1$heights[panels][2] <- unit(2/5,'null')
figure_1$heights[panels][3] <- unit(2/5,'null')
figure_1$heights[panels][4] <- unit(2/5,'null')
figure_1$heights[panels][5] <- unit(2/5,'null')
figure_1$heights[panels][6] <- unit(3/5,'null')
figure_1$heights[panels][7] <- unit(2/5,'null')




#grid.draw(figure_1)
ggsave(file="ExtFig2d_landscape_noncoding_v4.pdf", plot=figure_1,bg = 'white', width = 30, height = 13, units = 'cm', dpi = 600)



