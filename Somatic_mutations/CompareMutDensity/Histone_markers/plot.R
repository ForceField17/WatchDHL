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

# Function to produce summary statistics (mean and +/- sd)
data_summary <- function(x) {
  m <- mean(x)
  ymin <- max(m-sd(x),0)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

plot.table <- read.table('results.txt',header=T)
plot.table <- plot.table[which(plot.table$case!="P4" & plot.table$case!="P5"),]
plot.table$fraction <- plot.table$mutation/plot.table$total_mut/plot.table$size*100*1000000

plot.table2 <-plot.table[which( plot.table$tumor=="FL"),] 

NAN_plot <- ggplot(data=plot.table2,aes(x=compartment,y=fraction)) + theme_classic() 
NAN_plot <- NAN_plot +  #geom_boxplot(data=plot.table2,aes(x=compartment,y=fraction),width=0.2,alpha=0.5,size=0.6,outlier.shape = NA)+
  geom_dotplot(binaxis='y', stackdir='center', fill="grey20",binwidth=0.025,alpha=1)+
  stat_summary(fun.data=data_summary, color=gg_color_hue(6)[3],geom = "errorbar" ,width = 0)+
  stat_summary(fun=mean, aes(ymax = ..y.., ymin = ..y..), color=gg_color_hue(6)[3],geom = "errorbar" ,width = 0.3)
NAN_plot <-NAN_plot + theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,2,2,3),'lines'),
                            plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
                            legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='none',legend.text=element_text(size=14,hjust=0,face='bold'),
                            axis.text.x=element_text(size=14,angle=35,vjust=1,hjust=1,face='bold',color='black'),axis.text.y=element_text(size=14,face='bold',color='black'),
                            axis.title.x=element_text(size=14,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))
NAN_plot <- NAN_plot + scale_y_continuous(expand=c(0,0),limits=c(0,1.1),breaks = seq(0,1,0.2)) 
NAN_plot <- NAN_plot + scale_x_discrete()
NAN_plot<- NAN_plot +ylab("Normalized mutation density") +xlab(NULL)
theFL <- compare_means(fraction ~ compartment,  data = plot.table2,paired = T,method = "wilcox.test",p.adjust.method = "fdr" )

NAN_plot <- NAN_plot +stat_compare_means(label = "p.format", method = "wilcox.test",paired = T ,ref.group = "A_Random" ,label.y = 1.05,size=3)    

figure_2<-rbind(ggplotGrob(NAN_plot),size="last")


ggsave(file="ExtFig16b_FL_compartment_P126789.pdf", plot=figure_2,bg = 'white', width =16, height = 14, units = 'cm', dpi = 600)



plot.table2 <-plot.table[which( plot.table$tumor=="DHL"),] 

NAN_plot <- ggplot(data=plot.table2,aes(x=compartment,y=fraction)) + theme_classic() 
NAN_plot <- NAN_plot +  #geom_boxplot(data=plot.table2,aes(x=compartment,y=fraction),width=0.2,alpha=0.5,size=0.6,outlier.shape = NA)+
  geom_dotplot(binaxis='y', stackdir='center', fill="grey20",binwidth=0.025,alpha=1)+
  stat_summary(fun.data=data_summary, color=gg_color_hue(6)[1],geom = "errorbar" ,width = 0)+
  stat_summary(fun=mean, aes(ymax = ..y.., ymin = ..y..), color=gg_color_hue(6)[1],geom = "errorbar" ,width = 0.3)
NAN_plot <-NAN_plot + theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,2,2,3),'lines'),
                            plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
                            legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='none',legend.text=element_text(size=14,hjust=0,face='bold'),
                            axis.text.x=element_text(size=14,angle=35,vjust=1,hjust=1,face='bold',color='black'),axis.text.y=element_text(size=14,face='bold',color='black'),
                            axis.title.x=element_text(size=14,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))
NAN_plot <- NAN_plot + scale_y_continuous(expand=c(0,0),limits=c(0,1.1),breaks = seq(0,1,0.2)) 
NAN_plot <- NAN_plot + scale_x_discrete()
NAN_plot<- NAN_plot +ylab("Normalized mutation density") +xlab(NULL)
theDHL <- compare_means(fraction ~ compartment,  data = plot.table2,paired = T,method = "wilcox.test",p.adjust.method = "fdr" )

NAN_plot <- NAN_plot +stat_compare_means(label = "p.format", method = "wilcox.test",paired = T ,ref.group = "A_Random" ,label.y = 1.05,size=3)    

figure_2<-rbind(ggplotGrob(NAN_plot),size="last")

ggsave(file="ExtFig16b_DHL_compartment_P126789.pdf", plot=figure_2,bg = 'white', width =16, height = 14, units = 'cm', dpi = 600)

theFL
theDHL
