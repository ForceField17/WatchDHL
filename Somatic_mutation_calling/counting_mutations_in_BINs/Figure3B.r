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
Mut <- Mut[which(Mut$NumRecordsAll>=5),]

Mut$diffP <- Mut$NumP_DHL - Mut$NumP_FL
Mut$diffR <- abs(Mut$NumR_FL - Mut$NumR_DHL) + 1


sig <- Mut#[which(abs((Mut$diffP)>=3 ) | Mut$NumRecordsAll >= 10),]

sigEnhan <- sig[which(grepl(pattern = 'H3K27',x = sig$Enhancer,fixed=T) ),]
sigSE    <- sig[which(grepl(pattern = 'SE',x = sig$Enhancer,fixed=T) ),]
sigSE    <- sigSE[order(sigSE$diffR,decreasing = F),]


#chr$Chromosome

ffff_plot<-ggplot()+theme_classic()
ffff_plot<-ffff_plot+
  #geom_point(data=sig ,aes(x=start,y=diffP,size=diffR) ,shape=21,stroke=0.5,color = "black")+
  geom_point(data=sigEnhan ,aes(x=start,y=NumP_DHL,size=NumR_DHL) ,shape=21,fill = "blue",alpha=0.8)+
  geom_point(data=sigEnhan ,aes(x=start,y=-NumP_FL,size=NumR_FL) ,shape=21,fill = "blue",alpha=0.8)+
  geom_point(data=sigSE ,aes(x=start,y=NumP_DHL,size=NumR_DHL) ,shape=21,fill = "#d73027",alpha=0.8)+
  geom_point(data=sigSE ,aes(x=start,y=-NumP_FL,size=NumR_FL) ,shape=21,fill = "#d73027",alpha=0.8)+
  scale_y_continuous(expression(paste("log" ["2"], bold(" [mutation number]"))),expand=c(0,0),limits = c(-8.5, 8.5),breaks=seq(-8,8,2))+
  scale_x_continuous(expand=c(0,0),limits=c(0,3030851600),breaks=c(chr$X.log10.q.value.,2952009536),labels = c(1 ,2,  3  ,4,  5  ,6, 7 ,8 ,  9 ,10, 11,12 , 13,14 , 15,16 , 17,18 , 19,20 , 21,22,"X" ))+
  geom_vline( xintercept =c(0,chr$Position_offset),size=0.3,lty=5,color="grey")+ geom_hline(yintercept=2.5, color="#1a9850", lty=2)+
  geom_hline(yintercept=-2.5, color="#1a9850", lty=2)
ffff_plot<-ffff_plot+theme(panel.background=element_rect(fill='transparent',color='black',size=1),plot.margin=unit(c(1,1,1,1),'lines'),
                           plot.title=element_text(size=16,vjust=0.5,hjust=0.5,face='bold',color = "black"),text=element_text(size=15,vjust=0.5,hjust=0.5,face='bold'),
                           legend.key.width=unit(0.7,'cm'),legend.key.height=unit(0.7,'cm'),legend.position="top",axis.line = element_blank(),axis.ticks.y.left = element_line(),
                           legend.text=element_text(size=14,face='plain'),axis.text.y=element_text(size=12,vjust=0.5,hjust=0.5,face='bold',color='black'),axis.text.x=element_text(size=11,face='bold',color='black'),
                           axis.title.x=element_text(size=14,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=14,hjust=0.5,vjust=4,face='plain',color='black'))
ffff_plot<-ffff_plot+ggtitle(NULL)+xlab(NULL)#+coord_flip()#+scale_x_discrete(position = "bottom")
ffff_plot<-ffff_plot+scale_size(name=NULL,breaks = c(5,10,20,40,80),labels = expression(5,10,20,40,80) ,guide = guide_legend(override.aes=list(fill="white"),nrow=1))
#ffff_plot


figure_1<-rbind(ggplotGrob(ffff_plot))
ggsave(file="Diff_points.pdf", plot=figure_1,bg = 'white', width = 27, height =  11.5, units = 'cm', dpi = 600)


