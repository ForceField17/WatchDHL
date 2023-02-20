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

plot.table <- read.table('summary_INDEL.txt',header=T)
plot.table <- plot.table[which(!(plot.table$Case %in% c("P3","P4","P5"))),]
plot.table$InsP <- plot.table$Ins / plot.table$Total * 100
plot.table$DelP <- plot.table$Del / plot.table$Total * 100

NAN_plot <- ggplot() + theme_classic() 
NAN_plot <- NAN_plot + geom_violin(data=plot.table,aes(x=reorder(Tissue,-as.integer(as.factor(Tissue))),y=InsP,fill=Tissue),width=1,alpha=0.4,size=0.4) +
  geom_boxplot(data=plot.table,aes(x=reorder(Tissue,-as.integer(as.factor(Tissue))),y=InsP),width=0.2,alpha=0.5,size=0.6)+
  geom_path(data=plot.table,aes(x=reorder(Tissue,-as.integer(as.factor(Tissue))),y=InsP,group=Case),color="grey70",size=0.25) +
  geom_point(data=plot.table,aes(x=reorder(Tissue,-as.integer(as.factor(Tissue))),y=InsP,color=Tissue),alpha=.7,size=2.2) 
NAN_plot <-NAN_plot + theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,2,2,3),'lines'),
                            plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
                            legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='none',legend.text=element_text(size=14,hjust=0,face='bold'),
                            axis.text.x=element_text(size=14,face='bold',color='black'),axis.text.y=element_text(size=14,face='bold',color='black'),
                            axis.title.x=element_text(size=14,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))
NAN_plot <- NAN_plot + scale_y_continuous(expand=c(0,0),limits=c(0,10),breaks = seq(0,50,2)) 
NAN_plot<- NAN_plot +ylab("% of microinsertion among all mutations") +xlab(NULL)
NAN_plot <- NAN_plot + scale_fill_manual(name=NULL,values = c(FL=gg_color_hue(6)[3],DHL=gg_color_hue(6)[1]))
NAN_plot <- NAN_plot + scale_color_manual(name=NULL,values = c(FL=gg_color_hue(6)[3],DHL=gg_color_hue(6)[1]))
the<-compare_means(InsP ~ Tissue,  data = plot.table,paired = T,method = "wilcox.test",alternative = "greater" )
#NAN_plot<- NAN_plot +stat_compare_means(comparisons = my_comparisons)
NAN_plot <- NAN_plot +stat_pvalue_manual(
  the, label = NULL, 
  y.position = c( 9.6)
)

figure_2<-rbind(ggplotGrob(NAN_plot),size="last")
ggsave(file="ExtFig6a_Feq_Ins.pdf", plot=figure_2,bg = 'white', width =8, height = 10, units = 'cm', dpi = 600)


##########################################33
NAN_plot <- ggplot() + theme_classic() 
NAN_plot <- NAN_plot + geom_violin(data=plot.table,aes(x=reorder(Tissue,-as.integer(as.factor(Tissue))),y=DelP,fill=Tissue),width=1,alpha=0.4,size=0.4) +
  geom_boxplot(data=plot.table,aes(x=reorder(Tissue,-as.integer(as.factor(Tissue))),y=DelP),width=0.2,alpha=0.5,size=0.6)+
  geom_path(data=plot.table,aes(x=reorder(Tissue,-as.integer(as.factor(Tissue))),y=DelP,group=Case),color="grey70",size=0.25) +
  geom_point(data=plot.table,aes(x=reorder(Tissue,-as.integer(as.factor(Tissue))),y=DelP,color=Tissue),alpha=.7,size=2.2) 
NAN_plot <-NAN_plot + theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,2,2,3),'lines'),
                            plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
                            legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='none',legend.text=element_text(size=14,hjust=0,face='bold'),
                            axis.text.x=element_text(size=14,face='bold',color='black'),axis.text.y=element_text(size=14,face='bold',color='black'),
                            axis.title.x=element_text(size=14,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))
NAN_plot <- NAN_plot + scale_y_continuous(expand=c(0,0),limits=c(0,10),breaks = seq(0,50,2)) 
NAN_plot<- NAN_plot +ylab("% of microdeletion among all mutations") +xlab(NULL)
NAN_plot <- NAN_plot + scale_fill_manual(name=NULL,values = c(FL=gg_color_hue(6)[3],DHL=gg_color_hue(6)[1]))
NAN_plot <- NAN_plot + scale_color_manual(name=NULL,values = c(FL=gg_color_hue(6)[3],DHL=gg_color_hue(6)[1]))
the<-compare_means(DelP ~ Tissue,  data = plot.table,paired = T,method = "wilcox.test",alternative = "greater" )
#NAN_plot<- NAN_plot +stat_compare_means(comparisons = my_comparisons)
NAN_plot <- NAN_plot +stat_pvalue_manual(
  the, label = NULL, 
  y.position = c( 9.6)
)
#
figure_2<-rbind(ggplotGrob(NAN_plot),size="last")
ggsave(file="ExtFig5a_Feq_Del.pdf", plot=figure_2,bg = 'white', width =8, height = 10, units = 'cm', dpi = 600)
