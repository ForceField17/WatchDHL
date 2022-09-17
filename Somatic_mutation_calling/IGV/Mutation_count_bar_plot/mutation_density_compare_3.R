# CELLOR
setwd("/Users/songdong/Dropbox/Dropbox/DLBCL/main_figures/ZCCHC7/")
library(ggplot2)
library(gridExtra)
library(grid)
library(plyr)
library(rstatix)
library(ggpubr)


gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}



#####
FL<- read.csv("FL.specific.seg",header = F,sep = "\t")
DHL<- read.csv("DHL.specific.seg",header = F,sep = "\t")
FL <- FL[which(FL$V7!=0 & FL$V2=="chr9"),]
DHL <- DHL[which(DHL$V7!=0 & DHL$V2=="chr9"),]

gene.table <- data.frame(c(FL$V3/1000000,DHL$V3/1000000),c(rep("FL",nrow(FL)),rep("DHL",nrow(DHL))))
colnames(gene.table) <- c("Poition","Type")
F1a.plot<-ggplot()+theme_classic()
F1a.plot<-F1a.plot+geom_histogram(data=gene.table,aes(x=Poition,y=..count..,fill=Type),binwidth=0.003,alpha=0.5,position='identity',pad=T)
#F1a.plot<-F1a.plot+ geom_density(data=gene.table ,aes(x=Poition,color=Type),size=1,alpha=0.2)
F1a.plot<-F1a.plot+theme(panel.background=element_rect(fill='transparent',color='black',size=1),plot.margin=unit(c(0,1,1,1),'lines'),
                         plot.title=element_text(size=18,vjust=0.5,hjust=0.5,face='bold.italic',color='black'),text=element_text(size=14,face='bold'),
                         legend.key.width=unit(0.5,'cm'),legend.key.height=unit(0.5,'cm'),legend.position=c(0.9,0.5),legend.text=element_text(size=12,face='bold',color='black'),
                         axis.text.x=element_text(size=12,face='bold',color='black'),axis.text.y=element_text(size=12,face='bold',color='black'),axis.line = element_blank(),axis.ticks.y.left = element_line(),
                         axis.title.x=element_blank(),axis.title.y=element_blank())
F1a.plot<-F1a.plot+ggtitle(NULL)#+xlab("Age of onset")+ylab('Density')
cdat <- ddply(gene.table, "Type", summarise, AF.mean=mean(Poition))
cdat
#c(36.520268,37.655536)
#F1a.plot<-F1a.plot+geom_vline(data=cdat, aes(xintercept=AF.mean,colour=Type),linetype="dashed", size=1) #+ scale_color_manual(name=NULL,values=c(Male='#41b6c4',Female='#e78ac3'),labels=c(Male='Male',Female='Female'))
F1a.plot<-F1a.plot+scale_x_continuous(expand=c(0,0),limits=c(36.520268,37.655536),breaks=seq(0,38,0.2))+scale_y_continuous(expand=c(0,0),limits = c(0,30),breaks = c(0,15,30))#,limits=c(0,0.18))
#+scale_y_continuous(expand=c(0,0),limits=c(0,0.15))+scale_x_continuous(expand=c(0,0),limits=c(0,26),breaks=seq(0,100,5))#
F1a.plot<-F1a.plot + scale_fill_manual(name=NULL,values=c(FL=gg_color_hue(6)[3],DHL=gg_color_hue(6)[1]))+
                     scale_color_manual(name=NULL,values=c(FL=gg_color_hue(6)[3],DHL=gg_color_hue(6)[1]))
F1a.plot


#################3
CLL<- read.csv("WGS_nature.txt",header = F,sep = "\t")
colnames(CLL) <- c("Chrom",    "start",      "end",       "ref", "alt", "sig","method",   "Type",    "ID") 

gene.table2 <- data.frame( CLL$start/1000000,CLL$Type)
colnames(gene.table2) <- c("Poition","Type")
F1b.plot<-ggplot()+theme_classic()
F1b.plot<-F1b.plot+geom_histogram(data=gene.table2,aes(x=Poition,y=..count..,fill=Type),binwidth=0.003,alpha=0.5,position='identity',pad=T)

F1b.plot<-F1b.plot+theme(panel.background=element_rect(fill='transparent',color='black',size=1),plot.margin=unit(c(1,1,0,1),'lines'),
                         plot.title=element_text(size=18,vjust=0.5,hjust=0.5,face='bold.italic',color='black'),text=element_text(size=14,face='bold'),
                         legend.key.width=unit(0.5,'cm'),legend.key.height=unit(0.5,'cm'),legend.position=c(0.9,0.5),legend.text=element_text(size=12,face='bold',color='black'),
                         axis.text.x=element_blank(),axis.text.y=element_text(size=12,face='bold',color='black'),axis.line = element_blank(),axis.ticks.y.left = element_line(),
                         axis.title.x=element_blank(),axis.title.y=element_blank())
F1b.plot<-F1b.plot+ggtitle(NULL)#+xlab("Age of onset")+ylab('Density')

cdat <- ddply(gene.table2, "Type", summarise, AF.mean=mean(Poition))
cdat
F1b.plot<-F1b.plot+scale_x_continuous(expand=c(0,0),limits=c(36.520268,37.655536),breaks=seq(0,38,0.2),position = "bottom")+scale_y_continuous(expand=c(0,0),limits = c(0,30),breaks = c(0,15,30))
F1b.plot<-F1b.plot + scale_fill_manual(name=NULL,values=c(FL=gg_color_hue(6)[3],DLBCL=gg_color_hue(6)[1],CLL="blue",MBL=gg_color_hue(6)[4],MCL=gg_color_hue(6)[6],SLL=gg_color_hue(6)[5]))
F1b.plot


CLL<- read.csv("ZCCHC7.Morin.txt",header = F,sep = "\t")
colnames(CLL) <- c("Chrom",    "start",      "end",       "ref", "alt", "sig","method",   "Type",    "ID") 

gene.table2 <- data.frame( CLL$start/1000000,CLL$Type)
colnames(gene.table2) <- c("Poition","Type")
F1c.plot<-ggplot()+theme_classic()
F1c.plot<-F1c.plot+geom_histogram(data=gene.table2,aes(x=Poition,y=..count..,fill=Type),binwidth=0.003,alpha=0.5,position='identity',pad=T)

F1c.plot<-F1c.plot+theme(panel.background=element_rect(fill='transparent',color='black',size=1),plot.margin=unit(c(1,1,1,1),'lines'),
                         plot.title=element_text(size=18,vjust=0.5,hjust=0.5,face='bold.italic',color='black'),text=element_text(size=14,face='bold'),
                         legend.key.width=unit(0.5,'cm'),legend.key.height=unit(0.5,'cm'),legend.position=c(0.9,0.5),legend.text=element_text(size=12,face='bold',color='black'),
                         axis.text.x=element_blank(),axis.text.y=element_text(size=12,face='bold',color='black'),axis.line = element_blank(),axis.ticks.y.left = element_line(),
                         axis.title.x=element_blank(),axis.title.y=element_blank())
F1c.plot<-F1c.plot+ggtitle(NULL)#+xlab("Age of onset")+ylab('Density')
#F1c.plot<-F1c.plot+geom_vline(xintercept =37.033197,linetype=2,size=0.5)
#F1c.plot<-F1c.plot+geom_vline(xintercept =37.373097,linetype=2,size=0.5)
cdat <- ddply(gene.table2, "Type", summarise, AF.mean=mean(Poition))
cdat
F1c.plot<-F1c.plot+scale_x_continuous(expand=c(0,0),limits=c(36.520268,37.655536),breaks=seq(0,38,0.2))+scale_y_continuous(expand=c(0,0),limits = c(0,30),breaks = c(0,15,30,45))
F1c.plot<-F1c.plot + scale_fill_manual(name=NULL,values=c(FL=gg_color_hue(6)[3],DLBCL=gg_color_hue(6)[1],CLL="blue",MBL=gg_color_hue(6)[4],MCL=gg_color_hue(6)[6],SLL=gg_color_hue(6)[5]))
F1c.plot



figure_1<-rbind(ggplotGrob(F1b.plot),ggplotGrob(F1c.plot),ggplotGrob(F1a.plot))
ggsave(file="comparison_new.pdf", plot=figure_1,bg = 'white', width = 38, height =  16, units = 'cm', dpi = 600)



