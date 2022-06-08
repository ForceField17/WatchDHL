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

plot.table <- read.table('SV_insertion.txt',header=T)


NAN_plot <- ggplot() + theme_classic() 
NAN_plot <- NAN_plot + #geom_violin(data=plot.table,aes(x=SVgene,y=as.numeric(NumINS),fill=SVgene),width=4,alpha=0.7,size=0.4) +
  geom_boxplot(data=plot.table,aes(x=SVgene,y=as.numeric(NumINS)),width=0.2,alpha=0.5,size=0.6)+
  #geom_jitter(data=plot.table,aes(x=SVgene,y=as.numeric(NumINS),color=SVgene),width = 0.16,cex=2,alpha=1)#+
  geom_dotplot(data=plot.table,aes(x=SVgene,y=as.numeric(NumINS),fill=SVgene),color="grey40",shape=20,binaxis='y',binwidth=1, stackdir='center', dotsize=0.6)
  # geom_path(data=plot.table4,aes(x=Fusion,y=fraction,group=Case),color="black",size=0.3) +
  #geom_point(aes(x=SVgene,y=as.numeric(NumINS),color=SVgene),alpha=.7,data=plot.table,size=2.2) 
#geom_text(aes(x=Fusion,y=as.numeric(fraction),label=Case),data=plot.table2,size=3) 
NAN_plot <-NAN_plot + theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,2,2,3),'lines'),
                           plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
                           legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='none',legend.text=element_text(size=14,hjust=0,face='bold'),
                           axis.text.x=element_text(size=12,face='bold',color='black'),axis.text.y=element_text(size=12,face='bold',color='black'),
                           axis.title.x=element_text(size=14,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))
NAN_plot <- NAN_plot + scale_y_continuous(expand=c(0,0),limits=c(-0.5,28),breaks = seq(0,30,5)) 
#NAN_plot <- NAN_plot + scale_x_discrete(labels=c(B.FL="FL \ntumor",C.DHL="DHL \ntumor"))#
NAN_plot<- NAN_plot +ylab("Length of break-end insertions") +xlab(NULL)
NAN_plot <- NAN_plot + scale_fill_manual(name=NULL,values = c("#fdb462","#80b1d3"))
NAN_plot <- NAN_plot + scale_color_manual(name=NULL,values = c("#fdb462","#80b1d3"))
the<-compare_means(NumINS ~ SVgene,  data = plot.table,paired = F )
#NAN_plot<- NAN_plot +stat_compare_means(comparisons = my_comparisons)
NAN_plot <- NAN_plot +stat_pvalue_manual(
  the, label = NULL, 
  y.position = c( 26.6)
)



figure_2<-rbind(ggplotGrob(NAN_plot),size="last")

ggsave(file="breakEndInsertion.pdf", plot=figure_2,bg = 'white', width =11, height = 12, units = 'cm', dpi = 600)







