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

tableA <- read.table('./results.txt',header=F)
tableA$middle <- (tableA$V3+tableA$V5)/2
the1 <- wilcox.test(tableA$V3, tableA$V5, paired = T)
tableA1 <- data.frame(c(tableA$V1,tableA$V1),c(tableA$V2,tableA$V4),c(tableA$V3,tableA$V5),c(tableA$V7,tableA$V7),c(tableA$V11,tableA$V11),
              c(tableA$V12,tableA$V12),c(tableA$V13,tableA$V13),c(tableA$V14,tableA$V14),c(tableA$V15,tableA$V15),c(tableA$V16,tableA$V16),c(tableA$middle,tableA$middle),c(rep("FL",nrow(tableA)),rep("DHL",nrow(tableA))))
colnames(tableA1) <- c("peakID","Mutfraction","MutNumber","nearbyGene","MutDiff","MutTotal","Rank1","Rank2","Rank3","Rank4","middle","Tumor")


NAN_plot <- ggplot(data=tableA1,aes(x=Rank1,y=MutNumber)) + theme_classic() 
NAN_plot <- NAN_plot +
  geom_line(aes(group = peakID),color="grey",size=0.5)+
  geom_point(data=tableA1,aes(color=Tumor),alpha=.6,shape=20,size=2) 
NAN_plot <-NAN_plot + theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,2,0,3),'lines'),
                           plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
                           legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='none',legend.text=element_text(size=14,hjust=0,face='bold'),
                           axis.text.x=element_text(size=12,face='bold',color='black'),axis.text.y=element_text(size=12,face='bold',color='black'),
                           axis.title.x=element_text(size=14,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))
NAN_plot <- NAN_plot + scale_x_continuous(expand=c(0,0),limits=c(0.9,max(tableA1$Rank1+1000)),breaks = c(1,10,100,1000,10000),trans="log2") 
NAN_plot <- NAN_plot + scale_y_continuous(expand=c(0,0),limits=c(-0.2,10),breaks = seq(0,10,2))
NAN_plot<- NAN_plot +ylab("Number of Insertion") +xlab("Rank of peaks (Top 1000 mutated)")
NAN_plot <- NAN_plot + scale_fill_manual(name=NULL,values =  c(gg_color_hue(6)[1],gg_color_hue(6)[3]))
NAN_plot <- NAN_plot + scale_color_manual(name=NULL,values =  c(gg_color_hue(6)[1],gg_color_hue(6)[3]))
NAN_plot1 <- NAN_plot


figure_2<-rbind(ggplotGrob(NAN_plot1),size="last")

ggsave(file="ExtFig6b_figure_num.pdf", plot=figure_2,bg = 'white', width =16, height = 8, units = 'cm', dpi = 600)
the1$p.value

