# WatchDHL
library(rstudioapi)

# Getting the path of your current open file
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))
print( getwd() )

library(ggplot2)
library(plyr)
library("ggpubr")
library(ggbeeswarm)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

gene.table<-read.table('../gene_pairs_results_table.txt',header=T,sep="\t")
New <- gene.table[which(gene.table$Mutation=="Hypermutated"),] 
candidates <- New[which((New$MeanExpA>=9 & New$MeanExpB>=9 & New$SpearmanR < -0.5 ) ),]

F1a.plot<-ggplot()+theme_classic()
F1a.plot<-F1a.plot+ geom_histogram(data=gene.table,aes(x=SpearmanR,y=..density..),fill="#c69c72",binwidth=0.1,alpha=0.3,position='identity')
F1a.plot<-F1a.plot+ geom_density(data=gene.table ,aes(x=SpearmanR),color="#c69c72",size=1,alpha=1)
F1a.plot<-F1a.plot+ theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(1,1,1,1),'lines'),legend.title=element_text(size=12,face='plain',color='black'),
                          plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=12,face='bold'),
                          legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position="bottom",legend.text=element_text(size=10,hjust=0,face='plain'),
                          axis.text.x=element_text(size=12,angle=30,vjust=1,hjust=1,face='bold',color='black'),axis.text.y=element_text(size=12,face='bold',color='black'),
                          axis.title.x=element_text(size=12,vjust=0,hjust=0.5,face='plain',color='black'),axis.title.y=element_text(size=12,face='plain',color='black'))
F1a.plot<-F1a.plot+ggtitle(NULL)+xlab("Spearman's correlation coefficient")+ylab('Density of gene pairs')
cdat <- ddply(gene.table, "genePairs.direction", summarise, AF.mean=mean(SpearmanR))
cdat
F1a.plot<-F1a.plot+geom_vline(data=cdat, aes(xintercept=mean(gene.table$SpearmanR)),colour="#c69c72",linetype=2, size=0.5) + 
  scale_color_manual(name=NULL,values=c("#f25c54","#01baef"))
F1a.plot<-F1a.plot+geom_vline(data=New, aes(xintercept=SpearmanR),colour="grey40",linetype=2, size=0.2) 
F1a.plot<-F1a.plot+geom_vline(data=candidates , aes(xintercept=SpearmanR,color=genePairs.V1),linetype=1, size=0.4) 
F1a.plot<-F1a.plot+scale_y_continuous(expand=c(0,0),limits=c(0,1.5))+scale_x_continuous(expand=c(0,0),limits=c(-1,1),breaks=seq(-1,1,0.5))#
F1a.plot<-F1a.plot + scale_color_manual(name=NULL,values=c("#88419d","#06d6a0","#ef476f","#118ab2","#fc8d62"))
F1a.plot
figure_1<-rbind(ggplotGrob(F1a.plot),size="first")
ggsave(file="Fig6b_distri_SCC.pdf", plot=figure_1,bg = 'white', width = 15, height = 9, units = 'cm', dpi = 600)
the2 <- compare_means(SpearmanR ~ genePairs.direction, data = gene.table)
the2



mydata <- gene.table
F1a.plot<-ggplot(data= mydata ,aes(x=genePairs.direction,y=SpearmanR,fill=genePairs.direction))+theme_classic()
F1a.plot<-F1a.plot+
 geom_quasirandom(data= mydata ,aes(x=genePairs.direction,y=SpearmanR),width = 0.3,shape=21,size=2,alpha=0.1,stroke=0, varwidth = T)+
  geom_boxplot(data= mydata ,aes(x=genePairs.direction,y=SpearmanR),width=0.2,size=0.8,alpha=0.1,outlier.size = 1,outlier.shape = NA)
F1a.plot<-F1a.plot+theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,1,1,1),'lines'),
                         plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
                         legend.key.width=unit(1,'cm'),legend.key.height=unit(0.6,'cm'),legend.position="right",legend.text=element_text(size=14,hjust=0,face='plain'),
                         axis.text.x=element_text(size=12,angle=30,vjust=1,hjust=1,face='bold',color='black'),axis.text.y=element_text(size=12,face='bold',color='black'),
                         axis.title.x=element_text(size=12,vjust=0,hjust=0.5,face='plain',color='black'),axis.title.y=element_text(size=12,face='plain',color='black'))
F1a.plot<-F1a.plot+ggtitle(NULL)+ylab("Spearman's correlation coefficient")+xlab(NULL)

the <- compare_means(SpearmanR ~ genePairs.direction,  data =mydata  )
my_comparisons <- list( c("opposite_direction", "same_direction"))

F1a.plot<-F1a.plot+scale_y_continuous(expand=c(0,0),limits=c(-1,1.2))
F1a.plot<-F1a.plot + scale_fill_manual(name=NULL,values=c("#f4a582","#92c5de"),guide = F)

plotxxx<-cbind(ggplotGrob(F1a.plot),size="first")
ggsave(file="./ExtFig4d_box_SCC.pdf", plot=plotxxx,bg = 'white', width = 8, height = 12, units = 'cm', dpi = 600)
the



library("WGCNA")


test <-  data.frame(gene.table$SpearmanR,-gene.table$MeanExpA,-gene.table$MeanExpB)
rownames(test) <-  paste0(gene.table$genePairs.V1,"---",gene.table$genePairs.V2)
resultsAmp <- rankPvalue(test, columnweights = c(8,1,1), 
           na.last = "keep", ties.method = "average", 
           calculateQvalue = TRUE, pValueMethod = "rank")

gene.table$RankingPvalue <-  resultsAmp$pValueLowRank
write.table(gene.table , file = "Final_results.txt",sep = "\t",quote=F,row.names = F)




#