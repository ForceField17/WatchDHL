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

plot.table <- read.table('SV_read_count.txt',header=T)
plot.table <- plot.table[which(plot.table$Case!="P3"),]
plot.table2 <-plot.table[which(  (plot.table$Tissue!="B.FL.a" & (plot.table$Tissue!="A.Normal" ))),] 
plot.table2$ID <- paste0(plot.table2$Tissue," ",plot.table2$Fusion," fusion")

NAN_plot <- ggplot() + theme_classic() 
NAN_plot <- NAN_plot + geom_violin(data=plot.table2,aes(x=ID,y=as.numeric(fraction),fill=Tissue),width=2,alpha=0.4,size=0.4) +
  geom_boxplot(data=plot.table2,aes(x=ID,y=as.numeric(fraction),fill=Tissue),width=0.1,alpha=0.5,size=0.6)+
  geom_path(data=plot.table2,aes(x=ID,y=fraction,group=Case),color="grey70",size=0.25) +
  #geom_jitter(data=plot.table2,aes(x=Fusion,y=fraction,color=Case),width = 0.2,cex=3,alpha=.3,size=1.5)#+
  # geom_path(data=plot.table4,aes(x=Fusion,y=fraction,group=Case),color="black",size=0.3) +
  geom_point(aes(x=ID,y=as.numeric(fraction),color=Tissue),alpha=.7,data=plot.table2,size=2.2) 
#geom_text(aes(x=Fusion,y=as.numeric(fraction),label=Case),data=plot.table2,size=3) 
NAN_plot <-NAN_plot + theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,2,2,3),'lines'),
                            plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
                            legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='none',legend.text=element_text(size=14,hjust=0,face='bold'),
                            axis.text.x=element_text(size=14,angle=35,vjust=1,hjust=1,face='bold',color='black'),axis.text.y=element_text(size=14,face='bold',color='black'),
                            axis.title.x=element_text(size=14,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))
NAN_plot <- NAN_plot + scale_y_continuous(expand=c(0,0),limits=c(-0.05,1.1),breaks = seq(0,1,0.25)) 
NAN_plot <- NAN_plot + scale_x_discrete(labels=c(B.FL="FL \ntumor",C.DHL="DHL \ntumor"))#
NAN_plot<- NAN_plot +ylab("Allele frequency") +xlab(NULL)
NAN_plot <- NAN_plot + scale_fill_manual(name=NULL,values = c("grey",gg_color_hue(6)[3],gg_color_hue(6)[1],gg_color_hue(6)[1]))
NAN_plot <- NAN_plot + scale_color_manual(name=NULL,values = c("grey",gg_color_hue(6)[3],gg_color_hue(6)[1],gg_color_hue(6)[1]))
the <- compare_means(fraction ~ ID,  data = plot.table2,paired = T )
the <- the[c(1,5,6),]
#my_comparisons <- list( c("A.Normal BCL2 fusion", "B.FL BCL2 fusion") ,c("B.FL BCL2 fusion", "C.DHL BCL2 fusion") ,c("C.DHL BCL2 fusion",  "C.DHL MYC fusion"))
#NAN_plot <- NAN_plot  + stat_compare_means(comparisons = my_comparisons,method = "wilcox.test")
NAN_plot <- NAN_plot + stat_pvalue_manual(the, label = "p.format",  y.position = c( 0.85,0.95,1.05))



figure_2<-rbind(ggplotGrob(NAN_plot),size="last")

ggsave(file="AF_change1.pdf", plot=figure_2,bg = 'white', width =20, height = 16, units = 'cm', dpi = 600)

