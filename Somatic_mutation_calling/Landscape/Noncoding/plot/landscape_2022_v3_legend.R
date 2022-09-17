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


interestedGene <- c("CD83","TMSB4X","TCL1A","VMP1","ZFP36L1","LTB","UHRF2","FAM186A","SOCS1","BIRC3",
                    "CFAP251","BMP7","ETS1","MIR142","BCR","SGK1","BLK","SEL1L3","FOXO1","FAM102A","OSBPL10",
                    "MYC","FCHSD2","PATL2","PIM1","CXCR4","NCOA3","PALM2AKAP2","S1PR2","LHPP","BTG2",
                    "MIR3681HG","MYO16","PAX5","ZCCHC7","RHOH","IRF8","lnc-BCL6-7","DTX1","ST6GAL1",
                    "CIITA","LPP","DMD","BACH2","BCL6","DIAPH2","IGK","BCL7A","FHIT","DISC1FP1","IMMP2L","IGL","BCL2","IGH")

#UST


Sample_features <- read.table("Genetic_profiles_new.txt",sep = "\t",header = T)
Sample_features[is.na(Sample_features)] <- "f_NA"
sample <- data.frame(Sample_features)
sample <- sample[which(sample$patient=='P8'),]
pat <- unique(sample$patient)
case <- unique(sample$case)
#raw data preprocessing

somatic <- read.delim('./Mut/final/Annotated_20.frequent.txt',header = T)
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

F1b.pLot<-F1b.pLot+scale_fill_manual(name="Mutation Number",values =c(a0='white',b1='#fee0d2',c2='#fcbba1',c5='#fb6a4a',d10='#ef3b2c',e20='#cb181d',f30='#a50f15',g40='#67000d'),
                                     labels=c(a0='0',b1='1',c2='>=2',c5='>=5',d10='>=10',e20='>=20',f30='>=30',g40='>=40'), 
                                     guide = guide_legend(override.aes=list(size=0.5),nrow=8),na.translate = F)
#scale_fill_gradientn(name="Number of mutation",colours=c('white','#fe9929','#de2d26','#67001f'),limits=c(0,1.01),breaks=c(0,0.5,1))
#scale_fill_manual(name="Mutation Number",values =c(a0='white',b1='#fee0d2',c2='#fcbba1',c5='#fb6a4a',d10='#ef3b2c',e20='#cb181d',f30='#a50f15',g40='#67000d'),labels=c(a0='0',b1='1',c2='>=2',c5='>=5',d10='>=10',e20='>=20',f30='>=30',g40='>=40'))
F1b.pLot<-F1b.pLot+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1.3),plot.margin=unit(c(0.1,0.1,0.1,0.1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.7,'cm'),legend.key.height=unit(0.7,'cm'),legend.position='right',axis.ticks.x = element_blank(),axis.line.x = element_blank(),axis.line.y = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=16,face='bold.italic'),axis.ticks.y = element_blank(), axis.text.y=element_blank(),
                         axis.text.x=element_text(size=12,angle=45,vjust=1,hjust=1,face='bold',color='black'),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F1b.pLot<-F1b.pLot+scale_x_discrete(position = "bottom")+xlab(NULL)+ylab(NULL)
F1b.pLot
gene.scale1b = 1 + (ncol(mutation_gene_table) - 14)/14



plot1<-cbind(ggplotGrob(F1b.pLot),size="last")
grid.draw(plot1)

#grid.draw(figure_1)
ggsave(file="legend.pdf", plot=plot1,bg = 'white', width = 20, height = 38, units = 'cm', dpi = 600)



