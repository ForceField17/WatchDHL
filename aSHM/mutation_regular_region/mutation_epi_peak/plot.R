# WatchDHL
library(rstudioapi)

# Getting the path of your current open file
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))
print( getwd() )

library(ggplot2)
library(ggrepel)
library(reshape2)
library(plyr)
library(ggpubr)
library(rstatix)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
#raw data preprocessing

tableA <- read.table('Table_H3K27ac.PeakMut.txt',header=T)


xx <- (tableA$MutNumber[which(tableA$Tumor=="DHL")]  / tableA$size[which(tableA$Tumor=="DHL")])/27350*1000000
yy <- (tableA$MutNumber[which(tableA$Tumor=="FL")]  /  tableA$size[which(tableA$Tumor=="FL")])/41715*1000000
tableA$Density <- c(xx,yy)

the <- wilcox.test(tableA$Density[which(tableA$Tumor=="DHL")], tableA$Density[which(tableA$Tumor=="FL")], paired = T)
tableA1 <- tableA[which(tableA$Density==0),]
tableA2 <- tableA[which(tableA$Density>0),]

NAN_plot <- ggplot() + theme_classic() 
NAN_plot <- NAN_plot +
  geom_point(data=tableA1,aes(x=Rank,y=Density),color="grey",alpha=.6,shape=20,size=1) +
  geom_point(data=tableA2,aes(x=Rank,y=Density,color=Tumor),alpha=.6,shape=20,size=3) +
  geom_text(aes(x=max(tableA$Rank)/2,y=0.9,label=paste0("wilcox sign test pvalue = ",the$p.value)),color="black",size=3)
NAN_plot <-NAN_plot + theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,2,2,3),'lines'),
                           plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
                           legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='none',legend.text=element_text(size=14,hjust=0,face='bold'),
                           axis.text.x=element_text(size=12,face='bold',color='black'),axis.text.y=element_text(size=12,face='bold',color='black'),
                           axis.title.x=element_text(size=14,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))
NAN_plot <- NAN_plot + scale_y_continuous(expand=c(0,0),limits=c(-0.01,1),breaks = seq(0,100,0.1)) 
NAN_plot <- NAN_plot + scale_x_continuous(expand=c(0,0),limits=c(-500,max(tableA$Rank)+1),breaks = seq(0,max(tableA$Rank),5000)) 

NAN_plot<- NAN_plot +ylab("Normalized mutation density") +xlab("Rank of H3K27ac peaks")
NAN_plot <- NAN_plot + scale_fill_manual(name=NULL,values =  c(gg_color_hue(6)[1],gg_color_hue(6)[3]))
NAN_plot <- NAN_plot + scale_color_manual(name=NULL,values =  c(gg_color_hue(6)[1],gg_color_hue(6)[3]))


figure_2<-rbind(ggplotGrob(NAN_plot),size="last")

ggsave(file="H3K27ac.pdf", plot=figure_2,bg = 'white', width =28, height = 14, units = 'cm', dpi = 600)



##################################33
tableB <- read.table('Table_H3K4me3.PeakMut.txt',header=T)

xx <- (tableB$MutNumber[which(tableB$Tumor=="DHL")]  / tableB$size[which(tableB$Tumor=="DHL")])/27350*1000000
yy <- (tableB$MutNumber[which(tableB$Tumor=="FL")]  /  tableB$size[which(tableB$Tumor=="FL")])/41715*1000000
tableB$Density <- c(xx,yy)

the <- wilcox.test(tableB$Density[which(tableB$Tumor=="DHL")], tableB$Density[which(tableB$Tumor=="FL")], paired = T)
tableB1 <- tableB[which(tableB$Density==0),]
tableB2 <- tableB[which(tableB$Density>0),]

NAN_plot <- ggplot() + theme_classic() 
NAN_plot <- NAN_plot +
  geom_point(data=tableB1,aes(x=Rank,y=Density),color="grey",alpha=.6,shape=20,size=1) +
  geom_point(data=tableB2,aes(x=Rank,y=Density,color=Tumor),alpha=.6,shape=20,size=3) +
  geom_text(aes(x=max(tableB$Rank)/2,y=0.9,label=paste0("wilcox sign test pvalue = ",the$p.value)),color="black",size=3)
NAN_plot <-NAN_plot + theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,2,2,3),'lines'),
                            plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
                            legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='none',legend.text=element_text(size=14,hjust=0,face='bold'),
                            axis.text.x=element_text(size=12,face='bold',color='black'),axis.text.y=element_text(size=12,face='bold',color='black'),
                            axis.title.x=element_text(size=14,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))
NAN_plot <- NAN_plot + scale_y_continuous(expand=c(0,0),limits=c(-0.01,1),breaks = seq(0,100,0.1)) 
NAN_plot <- NAN_plot + scale_x_continuous(expand=c(0,0),limits=c(-500,max(tableB$Rank)+1),breaks = seq(0,max(tableB$Rank),5000)) 

NAN_plot<- NAN_plot +ylab("Normalized mutation density") +xlab("Rank of H3K4me3 peaks")
NAN_plot <- NAN_plot + scale_fill_manual(name=NULL,values =  c(gg_color_hue(6)[1],gg_color_hue(6)[3]))
NAN_plot <- NAN_plot + scale_color_manual(name=NULL,values =  c(gg_color_hue(6)[1],gg_color_hue(6)[3]))


figure_2<-rbind(ggplotGrob(NAN_plot),size="last")

ggsave(file="H3K4me3.pdf", plot=figure_2,bg = 'white', width =28, height = 14, units = 'cm', dpi = 600)


#################
tableC <- read.table('Table_H3K4me1.PeakMut.txt',header=T)
xx <- (tableC$MutNumber[which(tableC$Tumor=="DHL")]  / tableC$size[which(tableC$Tumor=="DHL")])/27350*1000000
yy <- (tableC$MutNumber[which(tableC$Tumor=="FL")]  /  tableC$size[which(tableC$Tumor=="FL")])/41715*1000000
tableC$Density <- c(xx,yy)

the <- wilcox.test(tableC$Density[which(tableC$Tumor=="DHL")], tableC$Density[which(tableC$Tumor=="FL")], paired = T)
tableC1 <- tableC[which(tableC$Density==0),]
tableC2 <- tableC[which(tableC$Density>0),]

NAN_plot <- ggplot() + theme_classic() 
NAN_plot <- NAN_plot +
  geom_point(data=tableC1,aes(x=Rank,y=Density),color="grey",alpha=.6,shape=20,size=1) +
  geom_point(data=tableC2,aes(x=Rank,y=Density,color=Tumor),alpha=.6,shape=20,size=3) +
  geom_text(aes(x=max(tableC$Rank)/2,y=0.9,label=paste0("wilcox sign test pvalue = ",the$p.value)),color="black",size=3)
NAN_plot <-NAN_plot + theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,2,2,3),'lines'),
                            plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
                            legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='none',legend.text=element_text(size=14,hjust=0,face='bold'),
                            axis.text.x=element_text(size=12,face='bold',color='black'),axis.text.y=element_text(size=12,face='bold',color='black'),
                            axis.title.x=element_text(size=14,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))
NAN_plot <- NAN_plot + scale_y_continuous(expand=c(0,0),limits=c(-0.02,2),breaks = seq(0,100,0.2)) 
NAN_plot <- NAN_plot + scale_x_continuous(expand=c(0,0),limits=c(-500,max(tableC$Rank)+1),breaks = seq(0,max(tableC$Rank),10000)) 

NAN_plot<- NAN_plot +ylab("Normalized mutation density") +xlab("Rank of H3K4me1 peaks")
NAN_plot <- NAN_plot + scale_fill_manual(name=NULL,values =  c(gg_color_hue(6)[1],gg_color_hue(6)[3]))
NAN_plot <- NAN_plot + scale_color_manual(name=NULL,values =  c(gg_color_hue(6)[1],gg_color_hue(6)[3]))


figure_2<-rbind(ggplotGrob(NAN_plot),size="last")

ggsave(file="H3K4me1.pdf", plot=figure_2,bg = 'white', width =28, height = 14, units = 'cm', dpi = 600)






