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
plot.table <- plot.table[which(!(plot.table$Case %in% c("P3","DHL38","DHL42"))),]
plot.table2 <-plot.table[which( plot.table$Tissue!="A.Normal" & plot.table$Tissue!="B.FL.a" & plot.table$Fusion=="MYC" ),] 
plot.table2$ID <- paste0(plot.table2$Tissue," ",plot.table2$Fusion," fusion")

NAN_plot <- ggplot() + theme_classic() 
NAN_plot <- NAN_plot + geom_violin(data=plot.table2,aes(x=Tissue,y=as.numeric(fraction),fill=Tissue),width=1,alpha=0.4,size=0.4) +
  geom_boxplot(data=plot.table2,aes(x=Tissue,y=as.numeric(fraction)),width=0.2,alpha=0.5,size=0.6)+
  geom_path(data=plot.table2,aes(x=Tissue,y=fraction,group=Case),color="grey70",size=0.25) +
  geom_point(aes(x=Tissue,y=as.numeric(fraction),color=Tissue),alpha=.7,data=plot.table2,size=2.2) 
NAN_plot <-NAN_plot + theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,2,2,3),'lines'),
                            plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
                            legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='none',legend.text=element_text(size=14,hjust=0,face='bold'),
                            axis.text.x=element_text(size=14,face='bold',color='black'),axis.text.y=element_text(size=14,face='bold',color='black'),
                            axis.title.x=element_text(size=14,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))
NAN_plot <- NAN_plot + scale_y_continuous(expand=c(0,0),limits=c(-0.01,1),breaks = seq(0,1,0.25)) 
NAN_plot <- NAN_plot + scale_x_discrete(labels=c(B.FL="FL \ntumor",C.DHL="DHL \ntumor"))#
NAN_plot<- NAN_plot +ylab("Allele frequency of MYC structural variation") +xlab(NULL)
NAN_plot <- NAN_plot + scale_fill_manual(name=NULL,values = c(B.FL=gg_color_hue(6)[3],C.DHL=gg_color_hue(6)[1]))
NAN_plot <- NAN_plot + scale_color_manual(name=NULL,values = c(B.FL=gg_color_hue(6)[3],C.DHL=gg_color_hue(6)[1]))
the<-compare_means(fraction ~ Tissue,  data = plot.table2,paired = T,method = "wilcox.test",alternative = "greater" )
the
NAN_plot <- NAN_plot +stat_pvalue_manual(
  the, label = NULL, 
  y.position = c( 0.95)
)


figure_2<-rbind(ggplotGrob(NAN_plot),size="last")

ggsave(file="Fig1e_AF_MYC.pdf", plot=figure_2,bg = 'white', width =8, height = 12, units = 'cm', dpi = 600)


##########################################33

plot.table2 <-plot.table[which( plot.table$Tissue!="A.Normal" & plot.table$Tissue!="B.FL.a" & plot.table$Fusion=="BCL2" ),] 
plot.table2$ID <- paste0(plot.table2$Tissue," ",plot.table2$Fusion," fusion")

NAN_plot <- ggplot() + theme_classic() 
NAN_plot <- NAN_plot + geom_violin(data=plot.table2,aes(x=Tissue,y=as.numeric(fraction),fill=Tissue),width=1,alpha=0.4,size=0.4) +
  geom_boxplot(data=plot.table2,aes(x=Tissue,y=as.numeric(fraction)),width=0.2,alpha=0.5,size=0.6)+
  geom_path(data=plot.table2,aes(x=Tissue,y=fraction,group=Case),color="grey70",size=0.25) +
  geom_point(aes(x=Tissue,y=as.numeric(fraction),color=Tissue),alpha=.7,data=plot.table2,size=2.2) 
NAN_plot <-NAN_plot + theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,2,2,3),'lines'),
                            plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
                            legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='none',legend.text=element_text(size=14,hjust=0,face='bold'),
                            axis.text.x=element_text(size=14,face='bold',color='black'),axis.text.y=element_text(size=14,face='bold',color='black'),
                            axis.title.x=element_text(size=14,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))
NAN_plot <- NAN_plot + scale_y_continuous(expand=c(0,0),limits=c(-0.01,1),breaks = seq(0,1,0.25)) 
NAN_plot <- NAN_plot + scale_x_discrete(labels=c(B.FL="FL \ntumor",C.DHL="DHL \ntumor"))#
NAN_plot<- NAN_plot +ylab("Allele frequency of BCL2 structural variation") +xlab(NULL)
NAN_plot <- NAN_plot + scale_fill_manual(name=NULL,values = c(B.FL=gg_color_hue(6)[3],C.DHL=gg_color_hue(6)[1]))
NAN_plot <- NAN_plot + scale_color_manual(name=NULL,values = c(B.FL=gg_color_hue(6)[3],C.DHL=gg_color_hue(6)[1]))
the <- compare_means(fraction ~ Tissue,  data = plot.table2,paired = T ,method = "wilcox.test",alternative = "greater" )
the 
NAN_plot <- NAN_plot +stat_pvalue_manual(
  the, label = NULL, 
  y.position = c( 0.95)
)


figure_2<-rbind(ggplotGrob(NAN_plot),size="last")

ggsave(file="Fig1d_AF_BCL2.pdf", plot=figure_2,bg = 'white', width =8, height = 12, units = 'cm', dpi = 600)


