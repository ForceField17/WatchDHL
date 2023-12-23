# WatchDHL
library(rstudioapi)
library("ggplot2")
library("gridExtra")
library(grid)
# Getting the path of your current open file
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))
print( getwd() )

library(ggplot2)
library(gridExtra)
library(grid)


gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 120)[1:n]
}

i=30

########################3
gistic_results <- read.table("./coordinates.tsv", header=T, sep="")

amp<- gistic_results[gistic_results[,1]=="Amp",]
del<- gistic_results[gistic_results[,1]=="Del",]
chr<- gistic_results[gistic_results[,1]=="Chr",]
gene_amp<- gistic_results[gistic_results[,1]=="Gene_amp",]
gene_del<- gistic_results[gistic_results[,1]=="Gene_del",]

Mut <- read.csv(paste0("./kateagis_DHL_",i,".final.txt"),header = T,sep = "\t")


MutHIGH <- Mut[which(Mut$label != "zNone"),]
MutLOW <- Mut[which(Mut$label == "zNone"),]


ffff_plot<-ggplot()+theme_classic()
ffff_plot<-ffff_plot+#geom_point(data=MutHIGH ,aes(x=Position,y=log2(value)),size=1 ,color = "#b2182b")+
  geom_point(data=MutLOW ,aes(x=pos,y=(dis+1),color=label),size=0.5 )+
  geom_point(data=MutHIGH ,aes(x=pos,y=(dis+1),color=label),size=0.5)+
  scale_y_continuous(expression(paste("Intermutation distance (bp)")),expand=c(0,0),limits = c(0.99, 9000000),breaks=c(1,100,10000,1000000,100000000),trans="log10")+
  scale_x_continuous(expand=c(0,0),limits=c(0,3030851600),breaks=c(chr$X.log10.q.value.,2952009536),labels = c(1 ,2,  3  ,4,  5  ,6, 7 ,8 ,  9 ,10, 11,12 , 13,14 , 15,16 , 17,18 , 19,20 , 21,22,"X" ))+
  geom_vline( xintercept =c(chr$Position_offset),size=0.3,lty=5,color="grey")+ geom_hline(yintercept=3000, color="#1a9850", lty=2)
ffff_plot<-ffff_plot+theme(panel.background=element_rect(fill='transparent',color='black',size=1),plot.margin=unit(c(0,2,1,1),'lines'),
                           plot.title=element_text(size=16,vjust=0.5,hjust=0.5,face='bold',color = "black"),text=element_text(size=15,vjust=0.5,hjust=0.5,face='bold'),
                           legend.key.width=unit(0.7,'cm'),legend.key.height=unit(0.7,'cm'),legend.position="none",axis.ticks.x = element_blank(),axis.line = element_blank(),axis.ticks.y.left = element_line(),
                           legend.text=element_text(size=14,face='bold'),axis.text.y=element_text(size=10,vjust=0.5,hjust=0.5,face='bold',color='black'),axis.text.x=element_blank(),
                           axis.title.x=element_text(size=14,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=10,hjust=0.5,vjust=4,face='plain',color='black'))
ffff_plot<-ffff_plot+ggtitle(NULL)+xlab(NULL)#+coord_flip()#+scale_x_discrete(position = "bottom")
ffff_plot<-ffff_plot+scale_color_manual(name=NULL,values=c(H3K27="blue",SE="#e41a1c",zNone="grey"),
                                        labels=c(H3K27="Enhancer",SE="Super enhancer",zNone="others"),
                                        guide = guide_legend(override.aes=list(size=3),nrow=1),na.translate = F)

########################3
gistic_results <- read.table("./coordinates.tsv", header=T, sep="")

amp<- gistic_results[gistic_results[,1]=="Amp",]
del<- gistic_results[gistic_results[,1]=="Del",]
chr<- gistic_results[gistic_results[,1]=="Chr",]
gene_amp<- gistic_results[gistic_results[,1]=="Gene_amp",]
gene_del<- gistic_results[gistic_results[,1]=="Gene_del",]

Mut_ <- read.csv(paste0("./kateagis_FL_",i,".final.txt"),header = T,sep = "\t")
Mut_HIGH <- Mut_[which(Mut_$label!="zNone"),]
Mut_LOW <- Mut_[which(Mut_$label=="zNone"),]


gggg_plot<-ggplot()+theme_classic()
gggg_plot<-gggg_plot+#geom_point(data=MutHIGH ,aes(x=Position,y=log2(value)),size=1 ,color = "#b2182b")+
  geom_point(data=Mut_LOW ,aes(x=pos,y=(dis+1),color=label),size=0.5 )+
  geom_point(data=Mut_HIGH ,aes(x=pos,y=(dis+1),color=label),size=0.5)+
  scale_y_continuous(expression(paste("Intermutation distance (bp)")),expand=c(0,0),limits = c(0.99, 9000000),breaks=c(1,100,10000,1000000,100000000),trans="log10")+
  scale_x_continuous(expand=c(0,0),limits=c(0,3030851600),breaks=c(chr$X.log10.q.value.,2952009536),labels = c(1 ,2,  3  ,4,  5  ,6, 7 ,8 ,  9 ,10, 11,12 , 13,14 , 15,16 , 17,18 , 19,20 , 21,22,"X" ))+
  geom_vline( xintercept =c(chr$Position_offset),size=0.3,lty=5,color="grey")+ geom_hline(yintercept=3000, color="#1a9850", lty=2)
gggg_plot<-gggg_plot+theme(panel.background=element_rect(fill='transparent',color='black',size=1),plot.margin=unit(c(1,2,0.1,1),'lines'),
                           plot.title=element_text(size=16,vjust=0.5,hjust=0.5,face='bold',color = "black"),text=element_text(size=15,vjust=0.5,hjust=0.5,face='bold'),
                           legend.key.width=unit(0.7,'cm'),legend.key.height=unit(0.7,'cm'),legend.position="none",axis.ticks.x = element_blank(),axis.line = element_blank(),axis.ticks.y.left = element_line(),
                           legend.text=element_text(size=14,face='bold'),axis.text.y=element_text(size=10,vjust=0.5,hjust=0.5,face='bold',color='black'),axis.text.x=element_text(size=11,vjust=1,hjust=0.5,face='bold',color='black'),
                           axis.title.x=element_text(size=14,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=10,hjust=0.5,vjust=4,face='plain',color='black'))
gggg_plot<-gggg_plot+ggtitle(NULL)+xlab(NULL)#+coord_flip()#+scale_x_discrete(position = "bottom")
gggg_plot<-gggg_plot+scale_color_manual(name=NULL,values=c(H3K27="blue",SE="#e41a1c",zNone="grey"),
                                        labels=c(H3K27="Enhancer",SE="Super enhancer",zNone="others"),
                                        guide = guide_legend(override.aes=list(size=3),nrow=1),na.translate = F)

figure_1<-rbind(ggplotGrob(gggg_plot),ggplotGrob(ffff_plot))
panels <- figure_1$layout$t[grep("panel", figure_1$layout$name)]

figure_1$heights[panels][1] <- unit(1,'null')
figure_1$heights[panels][2] <- unit(1,'null')
ggsave(file=paste0("./Fig4a_merge.",i,".SE.pdf"), plot=figure_1,bg = 'white', width = 27, height =  9, units = 'cm', dpi = 600)




##############################################3
#chr$Chromosome
Mut <- read.csv(paste0("./kateagis_DHL_",i,".final.txt"),header = T,sep = "\t")

Mut <- Mut[which(Mut$Chr==16),]
MutHIGH <- Mut[which(Mut$label != "zNone"),]
MutLOW <- Mut[which(Mut$label == "zNone"),]


ffff_plot<-ggplot()+theme_classic()
ffff_plot<-ffff_plot+#geom_point(data=MutHIGH ,aes(x=Position,y=log2(value)),size=1 ,color = "#b2182b")+
  geom_point(data=MutLOW ,aes(x=start/1000000,y=(dis+1),color=label),size=0.5 )+
  geom_point(data=MutHIGH ,aes(x=start/1000000,y=(dis+1),color=label),size=0.5)+
  scale_y_continuous(expression(paste("Intermutation distance (bp)")),expand=c(0,0),limits = c(0.99, 9000000),breaks=c(1,100,10000,1000000,100000000),trans="log10")+
  scale_x_continuous(expand=c(0,0),limits=c(0,90.338345),breaks=c(0,30,60,90),labels = c(0 ,"30Mbp","60Mbp","90Mbp" ),position = "top")+
  geom_vline( xintercept =c(chr$Position_offset),size=0.3,lty=5,color="grey")+ geom_hline(yintercept=3000, color="#1a9850", lty=2)
ffff_plot<-ffff_plot+theme(panel.background=element_rect(fill='transparent',color='black',size=1),plot.margin=unit(c(0,2,1,1),'lines'),
                           plot.title=element_text(size=16,vjust=0.5,hjust=0.5,face='bold',color = "black"),text=element_text(size=15,vjust=0.5,hjust=0.5,face='bold'),
                           legend.key.width=unit(0.7,'cm'),legend.key.height=unit(0.7,'cm'),legend.position="none",axis.ticks.x = element_blank(),axis.line = element_blank(),axis.ticks.y.left = element_line(),
                           legend.text=element_text(size=14,face='bold'),axis.text.y=element_text(size=10,vjust=0.5,hjust=0.5,face='bold',color='black'),axis.text.x=element_blank(),
                           axis.title.x=element_text(size=14,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=10,hjust=0.5,vjust=4,face='plain',color='black'))
ffff_plot<-ffff_plot+ggtitle(NULL)+xlab(NULL)
ffff_plot<-ffff_plot+scale_color_manual(name=NULL,values=c(H3K27="blue",SE="#e41a1c",zNone="grey"),
                                        labels=c(H3K27="Enhancer",SE="Super enhancer",zNone="others"),
                                        guide = guide_legend(override.aes=list(size=3),nrow=1),na.translate = F)
#ffff_plot

########################
gistic_results <- read.table("./coordinates.tsv", header=T, sep="")

amp<- gistic_results[gistic_results[,1]=="Amp",]
del<- gistic_results[gistic_results[,1]=="Del",]
chr<- gistic_results[gistic_results[,1]=="Chr",]
gene_amp<- gistic_results[gistic_results[,1]=="Gene_amp",]
gene_del<- gistic_results[gistic_results[,1]=="Gene_del",]

Mut <- read.csv(paste0("./kateagis_FL_",i,".final.txt"),header = T,sep = "\t")


#chr$Chromosome
Mut <- Mut[which(Mut$Chr==16),]
MutHIGH <- Mut[which(Mut$label != "zNone"),]
MutLOW <- Mut[which(Mut$label == "zNone"),]


gggg_plot<-ggplot()+theme_classic()
gggg_plot<-gggg_plot+#geom_point(data=MutHIGH ,aes(x=Position,y=log2(value)),size=1 ,color = "#b2182b")+
  geom_point(data=MutLOW ,aes(x=start/1000000,y=(dis+1),color=label),size=0.5 )+
  geom_point(data=MutHIGH ,aes(x=start/1000000,y=(dis+1),color=label),size=0.5)+
  scale_y_continuous(expression(paste("Intermutation distance (bp)")),expand=c(0,0),limits = c(0.99, 9000000),breaks=c(1,100,10000,1000000,100000000),trans="log10")+
  scale_x_continuous(expand=c(0,0),limits=c(0,90.338345),breaks=c(0,30,60,90),labels = c(0 ,"30Mbp","60Mbp","90Mbp" ))+
  geom_vline( xintercept =c(chr$Position_offset),size=0.3,lty=5,color="grey")+ geom_hline(yintercept=3000, color="#1a9850", lty=2)
gggg_plot<-gggg_plot+theme(panel.background=element_rect(fill='transparent',color='black',size=1),plot.margin=unit(c(1,2,0.1,1),'lines'),
                           plot.title=element_text(size=16,vjust=0.5,hjust=0.5,face='bold',color = "black"),text=element_text(size=15,vjust=0.5,hjust=0.5,face='bold'),
                           legend.key.width=unit(0.7,'cm'),legend.key.height=unit(0.7,'cm'),legend.position="none",axis.ticks.x = element_blank(),axis.line = element_blank(),axis.ticks.y.left = element_line(),
                           legend.text=element_text(size=14,face='bold'),axis.text.y=element_text(size=10,vjust=0.5,hjust=0.5,face='bold',color='black'),axis.text.x=element_text(size=10.35,vjust=1,hjust=0.5,face='bold',color='black'),
                           axis.title.x=element_text(size=14,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=10,hjust=0.5,vjust=4,face='plain',color='black'))
gggg_plot<-gggg_plot+ggtitle(NULL)+xlab(NULL)#+coord_flip()#+scale_x_discrete(position = "bottom")
gggg_plot<-gggg_plot+scale_color_manual(name=NULL,values=c(H3K27="blue",SE="#e41a1c",zNone="grey"),
                                        labels=c(H3K27="Enhancer",SE="Super enhancer",zNone="others"),
                                        guide = guide_legend(override.aes=list(size=3),nrow=1),na.translate = F)
#gggg_plot


figure_1<-rbind(ggplotGrob(gggg_plot),ggplotGrob(ffff_plot))

panels <- figure_1$layout$t[grep("panel", figure_1$layout$name)]

figure_1$heights[panels][1] <- unit(1,'null')
figure_1$heights[panels][2] <- unit(1,'null')
ggsave(file=paste0("./merge.chr16.",i,".SE.pdf"), plot=figure_1,bg = 'white', width = 27, height =  7, units = 'cm', dpi = 600)




gistic_results <- read.table("./coordinates.tsv", header=T, sep="")

amp<- gistic_results[gistic_results[,1]=="Amp",]
del<- gistic_results[gistic_results[,1]=="Del",]
chr<- gistic_results[gistic_results[,1]=="Chr",]
gene_amp<- gistic_results[gistic_results[,1]=="Gene_amp",]
gene_del<- gistic_results[gistic_results[,1]=="Gene_del",]



Mut <- read.csv(paste0("./kateagis_DHL_",i,".final.txt"),header = T,sep = "\t")

Mut <- Mut[which(Mut$Chr==2),]
MutHIGH <- Mut[which(Mut$label != "zNone"),]
MutLOW <- Mut[which(Mut$label == "zNone"),]


ffff_plot<-ggplot()+theme_classic()
ffff_plot<-ffff_plot+#geom_point(data=MutHIGH ,aes(x=Position,y=log2(value)),size=1 ,color = "#b2182b")+
  geom_point(data=MutLOW ,aes(x=start/1000000,y=(dis+1),color=label),size=0.5 )+
  geom_point(data=MutHIGH ,aes(x=start/1000000,y=(dis+1),color=label),size=0.5)+
  scale_y_continuous(expression(paste("Intermutation distance (bp)")),expand=c(0,0),limits = c(0.99, 9000000),breaks=c(1,100,10000,1000000,100000000),trans="log10")+
  scale_x_continuous(expand=c(0,0),limits=c(0,242.193529),breaks=c(0,60,120,180,240),labels = c(0 ,"60Mbp","120Mbp","180Mbp","240Mbp" ),position = "top")+
  geom_vline( xintercept =c(0,chr$Position_offset),size=0.3,lty=5,color="grey")+ geom_hline(yintercept=3000, color="#1a9850", lty=2)
ffff_plot<-ffff_plot+theme(panel.background=element_rect(fill='transparent',color='black',size=1),plot.margin=unit(c(0,2,1,1),'lines'),
                           plot.title=element_text(size=16,vjust=0.5,hjust=0.5,face='bold',color = "black"),text=element_text(size=15,vjust=0.5,hjust=0.5,face='bold'),
                           legend.key.width=unit(0.7,'cm'),legend.key.height=unit(0.7,'cm'),legend.position="none",axis.ticks.x = element_blank(),axis.line = element_blank(),axis.ticks.y.left = element_line(),
                           legend.text=element_text(size=14,face='bold'),axis.text.y=element_text(size=10,vjust=0.5,hjust=0.5,face='bold',color='black'),axis.text.x=element_blank(),
                           axis.title.x=element_text(size=14,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=10,hjust=0.5,vjust=4,face='plain',color='black'))
ffff_plot<-ffff_plot+ggtitle(NULL)+xlab(NULL)
ffff_plot<-ffff_plot+scale_color_manual(name=NULL,values=c(H3K27="blue",SE="#e41a1c",zNone="grey"),
                                        labels=c(H3K27="Enhancer",SE="Super enhancer",zNone="others"),
                                        guide = guide_legend(override.aes=list(size=3),nrow=1),na.translate = F)
#ffff_plot


#chr$Chromosome
gistic_results <- read.table("./coordinates.tsv", header=T, sep="")

amp<- gistic_results[gistic_results[,1]=="Amp",]
del<- gistic_results[gistic_results[,1]=="Del",]
chr<- gistic_results[gistic_results[,1]=="Chr",]
gene_amp<- gistic_results[gistic_results[,1]=="Gene_amp",]
gene_del<- gistic_results[gistic_results[,1]=="Gene_del",]

Mut <- read.csv(paste0("./kateagis_FL_",i,".final.txt"),header = T,sep = "\t")

Mut <- Mut[which(Mut$Chr==2),]
MutHIGH <- Mut[which(Mut$label != "zNone"),]
MutLOW <- Mut[which(Mut$label == "zNone"),]


gggg_plot<-ggplot()+theme_classic()
gggg_plot<-gggg_plot+#geom_point(data=MutHIGH ,aes(x=Position,y=log2(value)),size=1 ,color = "#b2182b")+
  geom_point(data=MutLOW ,aes(x=start/1000000,y=(dis+1),color=label),size=0.5 )+
  geom_point(data=MutHIGH ,aes(x=start/1000000,y=(dis+1),color=label),size=0.5)+
  scale_y_continuous(expression(paste("Intermutation distance (bp)")),expand=c(0,0),limits = c(0.99, 9000000),breaks=c(1,100,10000,1000000,100000000),trans="log10")+
  scale_x_continuous(expand=c(0,0),limits=c(0,242.193529),breaks=c(0,60,120,180,240),labels = c(0 ,"60Mbp","120Mbp","180Mbp","240Mbp" ))+
  geom_vline( xintercept =c(0,chr$Position_offset),size=0.3,lty=5,color="grey")+ geom_hline(yintercept=3000, color="#1a9850", lty=2)
gggg_plot<-gggg_plot+theme(panel.background=element_rect(fill='transparent',color='black',size=1),plot.margin=unit(c(1,2,0.1,1),'lines'),
                           plot.title=element_text(size=16,vjust=0.5,hjust=0.5,face='bold',color = "black"),text=element_text(size=15,vjust=0.5,hjust=0.5,face='bold'),
                           legend.key.width=unit(0.7,'cm'),legend.key.height=unit(0.7,'cm'),legend.position="none",axis.ticks.x = element_blank(),axis.line = element_blank(),axis.ticks.y.left = element_line(),
                           legend.text=element_text(size=14,face='bold'),axis.text.y=element_text(size=10,vjust=0.5,hjust=0.5,face='bold',color='black'),axis.text.x=element_text(size=10.35,vjust=1,hjust=0.5,face='bold',color='black'),
                           axis.title.x=element_text(size=14,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=10,hjust=0.5,vjust=4,face='plain',color='black'))
gggg_plot<-gggg_plot+ggtitle(NULL)+xlab(NULL)#+coord_flip()#+scale_x_discrete(position = "bottom")
gggg_plot<-gggg_plot+scale_color_manual(name=NULL,values=c(H3K27="blue",SE="#e41a1c",zNone="grey"),
                                        labels=c(H3K27="Enhancer",SE="Super enhancer",zNone="others"),
                                        guide = guide_legend(override.aes=list(size=3),nrow=1),na.translate = F)
#gggg_plot




figure_1<-rbind(ggplotGrob(gggg_plot),ggplotGrob(ffff_plot))

panels <- figure_1$layout$t[grep("panel", figure_1$layout$name)]

figure_1$heights[panels][1] <- unit(1,'null')
figure_1$heights[panels][2] <- unit(1,'null')
ggsave(file=paste0("./ExtFig2b_merge.chr2.",i,".SE.pdf"), plot=figure_1,bg = 'white', width = 27, height =  7, units = 'cm', dpi = 600)

