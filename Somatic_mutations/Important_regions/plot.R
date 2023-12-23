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

tableA <- read.table('H3K27ac/results.txt',header=F)
tableA$middle <- (tableA$V3+tableA$V5)/2
the1 <- wilcox.test(tableA$V3, tableA$V5, paired = T)
tableA1 <- data.frame(c(tableA$V1,tableA$V1),c(tableA$V2,tableA$V4),c(tableA$V3,tableA$V5),c(tableA$V7,tableA$V7),c(tableA$V11,tableA$V11),
              c(tableA$V12,tableA$V12),c(tableA$V13,tableA$V13),c(tableA$V14,tableA$V14),c(tableA$middle,tableA$middle),c(rep("FL",nrow(tableA)),rep("DHL",nrow(tableA))))
colnames(tableA1) <- c("peakID","Mutfraction","MutNumber","nearbyGene","MutDiff","Rank1","Rank2","Rank3","middle","Tumor")

#tableA1 <- tableA1[which(tableA1$Rank1 <= 10000),]
NAN_plot <- ggplot(data=tableA1,aes(x=Rank1,y=Mutfraction)) + theme_classic() 
NAN_plot <- NAN_plot +
  geom_line(aes(group = peakID),color="grey",size=0.5)+
  geom_point(data=tableA1,aes(color=Tumor),alpha=0.6,shape=20,size=2) 
NAN_plot <-NAN_plot + theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,2,0,3),'lines'),
                           plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
                           legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='none',legend.text=element_text(size=14,hjust=0,face='bold'),
                           axis.text.x=element_text(size=12,face='bold',color='black'),axis.text.y=element_text(size=12,face='bold',color='black'),
                           axis.title.x=element_text(size=14,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))
NAN_plot <- NAN_plot + scale_x_continuous(expand=c(0,0),limits=c(0.9,max(tableA1$Rank1+3000)),breaks = c(1,10,100,1000,10000),trans="log2") 
NAN_plot <- NAN_plot + scale_y_continuous(expand=c(0,0),limits=c(-0.01,0.65))
NAN_plot<- NAN_plot +ylab("Percentage of mutations loaded") +xlab("Rank of peaks")
NAN_plot <- NAN_plot + scale_fill_manual(name=NULL,values =  c(gg_color_hue(6)[1],gg_color_hue(6)[3]))
NAN_plot <- NAN_plot + scale_color_manual(name=NULL,values =  c(gg_color_hue(6)[1],gg_color_hue(6)[3]))
NAN_plot1 <- NAN_plot

figure<-rbind(ggplotGrob(NAN_plot),size="last")
ggsave(file="Fig5a_H3K27ac.pdf", plot=figure,bg = 'white', width =15, height = 8, units = 'cm', dpi = 600)
the1$p.value

tableA <- read.table('H3K4me3/results.txt',header=F)
tableA$middle <- (tableA$V3+tableA$V5)/2
the2 <- wilcox.test(tableA$V3, tableA$V5, paired = T)
tableA1 <- data.frame(c(tableA$V1,tableA$V1),c(tableA$V2,tableA$V4),c(tableA$V3,tableA$V5),c(tableA$V7,tableA$V7),c(tableA$V11,tableA$V11),
                      c(tableA$V12,tableA$V12),c(tableA$V13,tableA$V13),c(tableA$V14,tableA$V14),c(tableA$middle,tableA$middle),c(rep("FL",nrow(tableA)),rep("DHL",nrow(tableA))))
colnames(tableA1) <- c("peakID","Mutfraction","MutNumber","nearbyGene","MutDiff","Rank1","Rank2","Rank3","middle","Tumor")

#tableA1 <- tableA1[which(tableA1$Rank1 <= 10000),]
NAN_plot <- ggplot(data=tableA1,aes(x=Rank1,y=Mutfraction)) + theme_classic() 
NAN_plot <- NAN_plot +
  geom_line(aes(group = peakID),color="grey",size=0.5)+
  geom_point(data=tableA1,aes(color=Tumor),alpha=0.6,shape=20,size=2) 
NAN_plot <-NAN_plot + theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,2,0,3),'lines'),
                            plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
                            legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='none',legend.text=element_text(size=14,hjust=0,face='bold'),
                            axis.text.x=element_text(size=12,face='bold',color='black'),axis.text.y=element_text(size=12,face='bold',color='black'),
                            axis.title.x=element_text(size=14,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))
NAN_plot <- NAN_plot + scale_x_continuous(expand=c(0,0),limits=c(0.9,max(tableA1$Rank1+3000)),breaks = c(1,10,100,1000,10000),trans="log2") 
NAN_plot <- NAN_plot + scale_y_continuous(expand=c(0,0),limits=c(-0.01,0.65))
NAN_plot<- NAN_plot +ylab("Percentage of mutations loaded") +xlab("Rank of peaks")
NAN_plot <- NAN_plot + scale_fill_manual(name=NULL,values =  c(gg_color_hue(6)[1],gg_color_hue(6)[3]))
NAN_plot <- NAN_plot + scale_color_manual(name=NULL,values =  c(gg_color_hue(6)[1],gg_color_hue(6)[3]))
NAN_plot2 <- NAN_plot

figure<-rbind(ggplotGrob(NAN_plot),size="last")
ggsave(file="Fig5a_H3K4me3.pdf", plot=figure,bg = 'white', width =15, height = 8, units = 'cm', dpi = 600)
the2$p.value

tableA <- read.table('H3K4me1/results.txt',header=F)
tableA$middle <- (tableA$V3+tableA$V5)/2
the3 <- wilcox.test(tableA$V3, tableA$V5, paired = T)
tableA1 <- data.frame(c(tableA$V1,tableA$V1),c(tableA$V2,tableA$V4),c(tableA$V3,tableA$V5),c(tableA$V7,tableA$V7),c(tableA$V11,tableA$V11),
                      c(tableA$V12,tableA$V12),c(tableA$V13,tableA$V13),c(tableA$V14,tableA$V14),c(tableA$middle,tableA$middle),c(rep("FL",nrow(tableA)),rep("DHL",nrow(tableA))))
colnames(tableA1) <- c("peakID","Mutfraction","MutNumber","nearbyGene","MutDiff","Rank1","Rank2","Rank3","middle","Tumor")

#tableA1 <- tableA1[which(tableA1$Rank1 <= 10000),]
NAN_plot <- ggplot(data=tableA1,aes(x=Rank1,y=Mutfraction)) + theme_classic() 
NAN_plot <- NAN_plot +
  geom_line(aes(group = peakID),color="grey",size=0.5)+
  geom_point(data=tableA1,aes(color=Tumor),alpha=.6,shape=20,size=2) 
NAN_plot <-NAN_plot + theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,2,2,3),'lines'),
                            plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
                            legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='none',legend.text=element_text(size=14,hjust=0,face='bold'),
                            axis.text.x=element_text(size=12,face='bold',color='black'),axis.text.y=element_text(size=12,face='bold',color='black'),
                            axis.title.x=element_text(size=14,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))
NAN_plot <- NAN_plot + scale_x_continuous(expand=c(0,0),limits=c(0.9,max(tableA1$Rank1+3000)),breaks = c(1,10,100,1000,10000),trans="log2") 
NAN_plot<- NAN_plot +ylab("Percentage of mutations loaded") +xlab("Rank of peaks")
NAN_plot <- NAN_plot + scale_fill_manual(name=NULL,values =  c(gg_color_hue(6)[1],gg_color_hue(6)[3]))
NAN_plot <- NAN_plot + scale_color_manual(name=NULL,values =  c(gg_color_hue(6)[1],gg_color_hue(6)[3]))
NAN_plot3 <- NAN_plot

figure<-rbind(ggplotGrob(NAN_plot),size="last")
#ggsave(file="H3K4me1.pdf", plot=figure,bg = 'white', width =15, height = 8, units = 'cm', dpi = 600)
the3$p.value


figure_2<-rbind(ggplotGrob(NAN_plot1),ggplotGrob(NAN_plot2),size="last")
#ggsave(file="Merge_3peaks.pdf", plot=figure_2,bg = 'white', width =14, height = 17, units = 'cm', dpi = 600)
the1$p.value
the2$p.value
the3$p.value
