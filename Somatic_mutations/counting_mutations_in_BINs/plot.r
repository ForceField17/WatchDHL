# WatchDHL
library(rstudioapi)
library("ggplot2")
library("gridExtra")
library(grid)
# Getting the path of your current open file
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))
print( getwd() )



gistic_results <- read.table("coordinates.tsv", header=T, sep="")

amp<- gistic_results[gistic_results[,1]=="Amp",]
del<- gistic_results[gistic_results[,1]=="Del",]
chr<- gistic_results[gistic_results[,1]=="Chr",]
gene_amp<- gistic_results[gistic_results[,1]=="Gene_amp",]
gene_del<- gistic_results[gistic_results[,1]=="Gene_del",]

Mut <- read.csv("ranking_SE_Histone_3kb.txt",header = T,sep = "\t")
Mut <- Mut[which(Mut$NumRecordsAll >= 5),]
Mut$diffP <- Mut$NumP_DHL - Mut$NumP_FL
Mut$diffR <- abs(Mut$NumR_FL - Mut$NumR_DHL) + 1

MutLOW <- Mut[which(abs(Mut$diffP)<3),]
MutHIGH <- Mut[which(abs(Mut$diffP)>=3 ),]

sig <- Mut[which( (abs(Mut$diffP)>=3 & Mut$NumRecordsAll >= 5)),]

sigEnhan <- sig[which(grepl(pattern = 'H3K27',x = sig$Enhancer,fixed=T) ),]
sigSE <- sig[which(grepl(pattern = 'SE',x = sig$Enhancer,fixed=T) ),]
sigSE <- sigSE[order(sigSE$diffP,decreasing = T),]

sigEnhan <- sigEnhan[which(!(sigEnhan$BIN_id %in% sigSE$BIN_id)),]

sigH3K27 <- sig[which(grepl(pattern = 'H3K27',x = sig$Histone,fixed=T) ),]
sigH3K4 <-  sig[which(grepl(pattern = 'H3K4',x = sig$Histone,fixed=T) ),]
sigBoth <-  sig[which(grepl(pattern = 'Both',x = sig$Histone,fixed=T) ),]
sigBoth <- sigBoth[order(sigBoth$diffP,decreasing = F),]

#chr$Chromosome

ffff_plot<-ggplot()+theme_classic()
ffff_plot<-ffff_plot+  #geom_hline(yintercept=0, color="black", lty=1,size=0.5) +
  geom_point(data=MutLOW ,aes(x=as.numeric(as.character(start)),y=diffP),size=1 ,color = "grey")+
  geom_point(data=sig ,aes(x=as.numeric(as.character(start)),y=diffP,size=diffR) ,shape=21,fill="transparent",stroke=0.5,color = "black")+
  geom_point(data=sigSE ,aes(x=as.numeric(as.character(start)),y=diffP,size=diffR) ,shape=21,fill = "#e41a1c",alpha=0.8)+  
  geom_point(data=sigEnhan ,aes(x=as.numeric(as.character(start)),y=diffP,size=diffR) ,shape=21,fill = "blue",alpha=0.7)+
  scale_y_continuous("Difference in number of samples harboring mutation",expand=c(0,0),limits = c(-6, 6),breaks=seq(-8,8,2))+
  scale_x_continuous(expand=c(0,0),limits=c(0,3030851600),breaks=c(chr$X.log10.q.value.,2952009536),labels = c("1" ,2,  3  ,4,  5  ,6, 7 ,8 ,  9 ,10, 11,12 , 13,14 , 15,16 , 17,18 , 19,20 , 21,22,"X" ))+
  geom_vline( xintercept =c(chr$Position_offset),size=0.3,lty=5,color="grey")+ geom_hline(yintercept=2.5, color="#1a9850", lty=2)+
  geom_hline(yintercept=-2.5, color="#1a9850", lty=2)
 ffff_plot<-ffff_plot+geom_text_repel(data=sig ,aes(x=as.numeric(as.character(start)),y=diffP,label=gene),color='black',stat='identity',size=3.6,hjust = 1, nudge_y = 0.2)
#ffff_plot<-ffff_plot+geom_text(data=data_plot,aes(x=reorder(se,-orderID),y=0.5,label=Pat),color='black',stat='identity',size=3,hjust = 0, nudge_y = 0)
ffff_plot<-ffff_plot+theme(panel.background=element_rect(fill='transparent',color='black',size=1),plot.margin=unit(c(1,1,1,1),'lines'),
                           plot.title=element_text(size=16,vjust=0.5,hjust=0.5,face='bold',color = "black"),text=element_text(size=15,vjust=0.5,hjust=0.5,face='bold'),
                           legend.key.width=unit(0.7,'cm'),legend.key.height=unit(0.7,'cm'),legend.position="bottom",axis.line = element_blank(),axis.ticks.y.left = element_line(),axis.ticks.x = element_blank(),
                           legend.text=element_text(size=14,face='plain'),axis.text.y=element_text(size=12,hjust=1,face='bold',color='black'),axis.text.x=element_text(size=11,face='bold',color='black'),
                           axis.title.x=element_text(size=14,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=10,hjust=0.5,vjust=4,face='plain',color='black'))
ffff_plot<-ffff_plot+ggtitle(NULL)+xlab(NULL)#+coord_flip()#+scale_x_discrete(position = "bottom")
ffff_plot<-ffff_plot+scale_size(name=NULL,range=c(1,6),breaks = c(5,10,20,40,80),labels = expression(5,10,20,40,80) ,guide = guide_legend(override.aes=list(fill="white"),nrow=1))
#ffff_plot


figure_1<-rbind(ggplotGrob(ffff_plot))
ggsave(file="Figure4b.pdf", plot=figure_1,bg = 'white', width = 27, height =  11.5, units = 'cm', dpi = 600)


