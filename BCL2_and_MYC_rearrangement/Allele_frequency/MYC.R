# CELLOR
setwd("/Users/songdong/Dropbox/Dropbox/DLBCL/main_figures/SV/fraction_change/")
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

plot.table2 <- read.table('BCL2_AF.txt',header=T)
plot.table2 <-plot.table2[which(plot.table2$Fusion!="A.Normal"),] 
plot.table3 <- plot.table2[which(plot.table2$Var1=="B"),]
plot.table4 <- plot.table2[which(plot.table2$Var1=="A"),]

NAN_plot <- ggplot() + theme_classic() 
NAN_plot <- NAN_plot + geom_violin(data=plot.table2,aes(x=Fusion,y=as.numeric(fraction),fill=Fusion),width=1,alpha=0.4,size=0.4) +
  geom_boxplot(data=plot.table2,aes(x=Fusion,y=as.numeric(fraction)),width=0.2,alpha=0.5,size=0.6)+
  geom_path(data=plot.table2,aes(x=Fusion,y=fraction,group=Case),color="grey70",size=0.25) +
  #geom_jitter(data=plot.table2,aes(x=Fusion,y=fraction,color=Case),width = 0.2,cex=3,alpha=.3,size=1.5)#+
  # geom_path(data=plot.table4,aes(x=Fusion,y=fraction,group=Case),color="black",size=0.3) +
  geom_point(aes(x=Fusion,y=as.numeric(fraction),color=Fusion),alpha=.7,data=plot.table2,size=2.2) 
#geom_text(aes(x=Fusion,y=as.numeric(fraction),label=Case),data=plot.table2,size=3) 
NAN_plot <-NAN_plot + theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,2,2,3),'lines'),
                           plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
                           legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='none',legend.text=element_text(size=14,hjust=0,face='bold'),
                           axis.text.x=element_text(size=14,face='bold',color='black'),axis.text.y=element_text(size=14,face='bold',color='black'),
                           axis.title.x=element_text(size=14,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))
NAN_plot <- NAN_plot + scale_y_continuous(expand=c(0,0),limits=c(-0.05,1),breaks = seq(0,1,0.25)) 
NAN_plot <- NAN_plot + scale_x_discrete(labels=c(B.FL="FL \ntumor",C.DHL="DHL \ntumor"))#
NAN_plot<- NAN_plot +ylab("Allele frequency of MYC structural variation") +xlab(NULL)
NAN_plot <- NAN_plot + scale_fill_manual(name=NULL,values = c(B.FL=gg_color_hue(6)[3],C.DHL=gg_color_hue(6)[1]))
NAN_plot <- NAN_plot + scale_color_manual(name=NULL,values = c(B.FL=gg_color_hue(6)[3],C.DHL=gg_color_hue(6)[1]))
the<-compare_means(fraction ~ Fusion,  data = plot.table2,paired = F )
#NAN_plot<- NAN_plot +stat_compare_means(comparisons = my_comparisons)
NAN_plot <- NAN_plot +stat_pvalue_manual(
  the, label = NULL, 
  y.position = c( 0.95)
)

#figure_2<-rbind(ggplotGrob(NAN_plot),size="last")

#ggsave(file="figure_rareSNPs.pdf", plot=figure_2,bg = 'white', width =13, height = 15, units = 'cm', dpi = 600)




figure_2<-rbind(ggplotGrob(NAN_plot),size="last")

ggsave(file="BCL2_b.pdf", plot=figure_2,bg = 'white', width =8, height = 12, units = 'cm', dpi = 600)







