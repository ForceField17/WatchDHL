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

data_summary <- function(x) {
  m <- mean(x)
  ymin <- max(m-sd(x),0)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}


plot.table <- read.table('result.txt',header=T)
plot.table <- plot.table[which(plot.table$case!="P4" & plot.table$case!="P5"),]
plot.table$fraction <- plot.table$mutation/plot.table$total_mut*100

plot.table2 <-plot.table#[which( plot.table$tumor=="FL"),] 

NAN_plot <- ggplot(data=plot.table2,aes(x=xxx,y=fraction)) + theme_classic() 
NAN_plot <- NAN_plot + 
  geom_path(data=plot.table2,aes(x=xxx,y=fraction,group=tags),color="grey70",size=0.25) +
  geom_point(data=plot.table2,aes(x=xxx,y=as.numeric(fraction)),shape=21,fill="grey20",alpha=0.7,size=2) +
  stat_summary(fun.data=data_summary, aes(color=tumor),geom = "errorbar" ,width = 0)+
  stat_summary(fun=mean, aes(ymax = ..y.., ymin = ..y..,color=tumor),geom = "errorbar" ,width = 0.3)
NAN_plot <-NAN_plot + theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,2,2,3),'lines'),
                            plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
                            legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='none',legend.text=element_text(size=14,hjust=0,face='bold'),
                            axis.text.x=element_text(size=14,face='bold',color='black'),axis.text.y=element_text(size=14,face='bold',color='black'),
                            axis.title.x=element_text(size=14,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))
NAN_plot <- NAN_plot + scale_y_continuous(expand=c(0,0),limits=c(0,40),breaks = seq(0,100,10)) 
NAN_plot <- NAN_plot + scale_x_discrete(labels=c(B="B",A="A"))#
NAN_plot<- NAN_plot +ylab("percentage of mutations in SE") +xlab(NULL)
NAN_plot <- NAN_plot + scale_color_manual(name=NULL,values = c(FL=gg_color_hue(6)[3],DHL=gg_color_hue(6)[1]))
plot.table3 <-plot.table[which( plot.table$case!="P8.b"),] 
the<-compare_means(fraction ~ xxx,  data = plot.table3,paired = T,method = "wilcox.test" ,p.adjust.method = "fdr")
the

figure_2<-rbind(ggplotGrob(NAN_plot),size="last")

ggsave(file="ExtFig5e.pdf", plot=figure_2,bg = 'white', width =14, height = 10, units = 'cm', dpi = 600)

FL <- plot.table[which(plot.table$tumor=="FL"),]
the<-compare_means(fraction ~ xxx,  data = FL,paired = T,method = "wilcox.test" ,p.adjust.method = "fdr")
the
 
DHL <- plot.table[which(plot.table$tumor=="DHL"),]
the<-compare_means(fraction ~ xxx,  data = DHL,paired = T,method = "wilcox.test" ,p.adjust.method = "fdr")
the