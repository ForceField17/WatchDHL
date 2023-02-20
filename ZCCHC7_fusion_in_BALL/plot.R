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

plot.table <- read.table('data.txt',header=T)

Sample <- plot.table


#Sample <- Sample[which(Sample$Cancer=="BALL" | Sample$Cancer=="TALL"),]
#Sample <- Sample[which(Sample$Type=="Fusion"),]
Sample$Zscore_PAX5 <- log10(Sample$PAX5)
Sample$Zscore_ZCCHC7 <- log10(Sample$ZCCHC7)

a<-cor.test(Sample$Zscore_PAX5,Sample$Zscore_ZCCHC7,method = "spearman")
a$p.value
a$estimate


NAN_plot <- ggplot() + theme_classic() 
NAN_plot <- NAN_plot + #geom_violin(data=Sample,aes(x=SVgene,y=as.numeric(NumINS),fill=SVgene),width=4,alpha=0.7,size=0.4) +
  geom_boxplot(data=Sample,aes(x=reorder(Type,-as.integer(as.factor(Type))),y=Zscore_ZCCHC7,fill=Type),width=0.4,alpha=0.9,size=0.5,outlier.size = 0.3)#,outlier.shape = NA)+
  #geom_dotplot(data=Sample,aes(x=Type,y=Zscore_ZCCHC7),color="grey40",binaxis='y',binwidth=0.05, stackdir='center', dotsize=0.2)
NAN_plot <-NAN_plot + theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,2,2,2),'lines'),
                           plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
                           legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='none',legend.text=element_text(size=14,hjust=0,face='bold'),
                           axis.text.x=element_text(size=12,face='bold',color='black'),axis.text.y=element_text(size=12,face='bold',color='black'),
                           axis.title.x=element_text(size=14,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))
#NAN_plot <- NAN_plot + scale_y_continuous(expand=c(0,0),limits=c(0,3.1),breaks = seq(-10,10,0.5)) 
NAN_plot<- NAN_plot +ylab("RNA expression of ZCCHC7 (log10 FPKM)") +xlab(NULL)
NAN_plot <- NAN_plot + scale_fill_manual(name=NULL,values = c("#bc80bd","#d9d9d9"))
NAN_plot <- NAN_plot + scale_color_manual(name=NULL,values = c("#fdb462","#80b1d3"))
the<-compare_means(Zscore_ZCCHC7 ~ Type,  data = Sample,paired = F )
#NAN_plot<- NAN_plot +stat_compare_means(comparisons = my_comparisons)
NAN_plot <- NAN_plot +stat_pvalue_manual(
  the, label = NULL, 
  y.position = c( 2.9)
)

figure_2<-rbind(ggplotGrob(NAN_plot),size="last")

ggsave(file="ZCCHC7_all.pdf", plot=figure_2,bg = 'white', width =8, height = 12, units = 'cm', dpi = 600)





NAN_plot <- ggplot() + theme_classic() 
NAN_plot <- NAN_plot + #geom_violin(data=Sample,aes(x=SVgene,y=as.numeric(NumINS),fill=SVgene),width=4,alpha=0.7,size=0.4) +
  geom_boxplot(data=Sample,aes(x=reorder(Label,-as.integer(as.factor(Label))),y=Zscore_ZCCHC7,fill=Type),width=0.4,alpha=0.9,size=0.5,outlier.size = 0.3)#,outlier.shape = NA)+
#geom_dotplot(data=Sample,aes(x=Type,y=Zscore_ZCCHC7),color="grey40",binaxis='y',binwidth=0.05, stackdir='center', dotsize=0.2)
NAN_plot <-NAN_plot + theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,2,2,2),'lines'),
                            plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
                            legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='none',legend.text=element_text(size=14,hjust=0,face='bold'),
                            axis.text.x=element_text(size=12,face='bold',color='black'),axis.text.y=element_text(size=12,face='bold',color='black'),
                            axis.title.x=element_text(size=14,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))
NAN_plot <- NAN_plot + scale_y_continuous(expand=c(0,0),limits=c(0,3.1),breaks = seq(-10,10,0.5)) 
NAN_plot<- NAN_plot +ylab("RNA expression of ZCCHC7 (log10 FPKM)") +xlab(NULL)
NAN_plot <- NAN_plot + scale_fill_manual(name=NULL,values = c("#bc80bd","#d9d9d9"))
NAN_plot <- NAN_plot + scale_color_manual(name=NULL,values = c("#fdb462","#80b1d3"))
the<-compare_means(Zscore_ZCCHC7 ~ Label,  data = Sample,paired = F )
#NAN_plot<- NAN_plot +stat_compare_means(comparisons = my_comparisons)
NAN_plot <- NAN_plot +stat_pvalue_manual(
  the, label = NULL, 
  y.position = c( 2.9)
)

figure_2<-rbind(ggplotGrob(NAN_plot),size="last")

ggsave(file="ExtFig18b_ZCCHC7_Label.pdf", plot=figure_2,bg = 'white', width =16, height = 12, units = 'cm', dpi = 600)





