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

table <- read.table('bp_motif.txt',header=T)

tableA <- table[which(table$Locus == "BCL2" | table$Target == "BCL2"),]
tableA <- tableA[which(tableA$patient != "P3"),]

NAN_plot <- ggplot() + theme_classic() 
NAN_plot <- NAN_plot + geom_violin(data=tableA,aes(x=Locus,y=as.numeric(Min_CpG),fill=Locus),width=1.1,alpha=0.3,size=0.4) +
  geom_boxplot(data=tableA,aes(x=Locus,y=as.numeric(Min_CpG)),width=0.06,alpha=0.5,size=0.5)+
  #geom_jitter(data=tableA,aes(x=Locus,y=as.numeric(Min_CpG),color=SVgene),width = 0.16,cex=2,alpha=1)#+
  #geom_dotplot(data=tableA,aes(x=Locus,y=as.numeric(Min_CpG),fill=Tumor),color="grey40",shape=20,binaxis='y',binwidth=1, stackdir='center', dotsize=0.6)
   geom_path(data=tableA,aes(x=Locus,y=Min_CpG,group=patient ),color="grey60",size=0.3) +
  geom_point(aes(x=Locus,y=as.numeric(Min_CpG),color=Tumor),alpha=.9,shape=20,data=tableA,size=3) 
 #geom_text(aes(x=Fusion,y=as.numeric(fraction),label=Case),data=tableA2,size=3) 
NAN_plot <-NAN_plot + theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,2,2,3),'lines'),
                           plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
                           legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='none',legend.text=element_text(size=14,hjust=0,face='bold'),
                           axis.text.x=element_text(size=12,face='bold',color='black'),axis.text.y=element_text(size=12,face='bold',color='black'),
                           axis.title.x=element_text(size=14,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))
NAN_plot <- NAN_plot + scale_y_continuous(expand=c(0,0),limits=c(-0.5,17),breaks = seq(0,30,5)) 
#NAN_plot <- NAN_plot + scale_x_discrete(labels=c(B.FL="FL \ntumor",C.DHL="DHL \ntumor"))#
NAN_plot<- NAN_plot +ylab("Distance from nearest CpG to breakpoint") +xlab(NULL)
NAN_plot <- NAN_plot + scale_fill_manual(name=NULL,values = c("#ff7f00","#bc80bd"))
NAN_plot <- NAN_plot + scale_color_manual(name=NULL,values = c(gg_color_hue(6)[3],"#80b1d3"))
the<-compare_means(Min_CpG ~ Locus,  data = tableA,paired = T )
#NAN_plot<- NAN_plot +stat_compare_means(comparisons = my_comparisons)
NAN_plot <- NAN_plot +stat_pvalue_manual(
  the, label = NULL, 
  y.position = c( 16)
)


figure_2<-rbind(ggplotGrob(NAN_plot),size="last")

ggsave(file="ExtFig1d_CpG_BCL2.pdf", plot=figure_2,bg = 'white', width =11, height = 12, units = 'cm', dpi = 600)




tableB <- table[which(table$Locus != "BCL2" & table$Target != "BCL2"),]
tableB$bp <- tableB$Locus 
tableB$bp[which(tableB$Locus != "MYC")] <- "partner_gene" 
NAN_plot <- ggplot() + theme_classic() 
NAN_plot <- NAN_plot + geom_violin(data=tableB,aes(x=bp,y=as.numeric(Min_CpG)),fill="#fdb462",width=1.1,alpha=0.4,size=0.4) +
  geom_boxplot(data=tableB,aes(x=bp,y=as.numeric(Min_CpG)),width=0.06,alpha=0.5,size=0.5)+
  geom_path(data=tableB,aes(x=bp,y=Min_CpG,group=patient ),color="grey60",size=0.3) +
  geom_point(aes(x=bp,y=as.numeric(Min_CpG),color=Tumor),alpha=.6,shape=20,data=tableB,size=3) 
#geom_text(aes(x=Fusion,y=as.numeric(fraction),label=Case),data=tableB2,size=3) 
NAN_plot <-NAN_plot + theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,2,2,3),'lines'),
                            plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
                            legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='none',legend.text=element_text(size=14,hjust=0,face='bold'),
                            axis.text.x=element_text(size=12,face='bold',color='black'),axis.text.y=element_text(size=12,face='bold',color='black'),
                            axis.title.x=element_text(size=14,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))
NAN_plot <- NAN_plot + scale_y_continuous(expand=c(0,0),limits=c(-0.5,30),breaks = seq(0,30,5)) 
#NAN_plot <- NAN_plot + scale_x_discrete(labels=c(B.FL="FL \ntumor",C.DHL="DHL \ntumor"))#
NAN_plot<- NAN_plot +ylab("Distance from nearest CpG to breakpoint") +xlab(NULL)
NAN_plot <- NAN_plot + scale_fill_manual(name=NULL,values = c("#fdb462","#80b1d3"))
NAN_plot <- NAN_plot + scale_color_manual(name=NULL,values = c(gg_color_hue(6)[3],gg_color_hue(6)[1]))
the<-compare_means(Min_CpG ~ bp,  data = tableB,paired = T )
#NAN_plot<- NAN_plot +stat_compare_means(comparisons = my_comparisons)
NAN_plot <- NAN_plot +stat_pvalue_manual(
  the, label = NULL, 
  y.position = c( 28.5)
)


figure_2<-rbind(ggplotGrob(NAN_plot),size="last")

ggsave(file="CpG_MYC.pdf", plot=figure_2,bg = 'white', width =10, height = 12, units = 'cm', dpi = 600)


tableC <- table[which(table$Locus != "BCL2" & table$Target != "BCL2"),]
tableC$bp <- tableC$Locus 
tableC$bp[which(tableC$Locus != "MYC" & tableC$Locus != "IgH")] <- "others" 
NAN_plot <- ggplot() + theme_classic() 
NAN_plot <- NAN_plot + geom_violin(data=tableC,aes(x=bp,y=as.numeric(Min_CpG)),fill="#fdb462",width=1.1,alpha=0.4,size=0.4) +
  geom_boxplot(data=tableC,aes(x=bp,y=as.numeric(Min_CpG)),width=0.06,alpha=0.5,size=0.5,outlier.shape = NA)+
  geom_jitter(aes(x=bp,y=as.numeric(Min_CpG),color=Tumor),alpha=.6,shape=20,data=tableC,size=3,width = 0.2) 
#geom_text(aes(x=Fusion,y=as.numeric(fraction),label=Case),data=tableC2,size=3) 
NAN_plot <-NAN_plot + theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,2,2,3),'lines'),
                            plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
                            legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='none',legend.text=element_text(size=14,hjust=0,face='bold'),
                            axis.text.x=element_text(size=12,face='bold',color='black'),axis.text.y=element_text(size=12,face='bold',color='black'),
                            axis.title.x=element_text(size=14,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))
NAN_plot <- NAN_plot + scale_y_continuous(expand=c(0,0),limits=c(-0.5,32),breaks = seq(0,30,5)) 
#NAN_plot <- NAN_plot + scale_x_discrete(labels=c(B.FL="FL \ntumor",C.DHL="DHL \ntumor"))#
NAN_plot<- NAN_plot +ylab("Distance from nearest CpG to breakpoint") +xlab(NULL)
NAN_plot <- NAN_plot + scale_fill_manual(name=NULL,values = c("#fdb462","#80b1d3"))
NAN_plot <- NAN_plot + scale_color_manual(name=NULL,values = c(gg_color_hue(6)[3],gg_color_hue(6)[1]))
the<-compare_means(Min_CpG ~ bp,  data = tableC,paired = F )
#NAN_plot<- NAN_plot +stat_compare_means(comparisons = my_comparisons)
NAN_plot <- NAN_plot +stat_pvalue_manual(
  the, label = NULL, 
  y.position = c( 26,28,30)
)


figure_2<-rbind(ggplotGrob(NAN_plot),size="last")

#ggsave(file="CpG_MYC_2.pdf", plot=figure_2,bg = 'white', width =10, height = 12, units = 'cm', dpi = 600)



tableD <- table[which(table$Locus != "BCL2" & table$Target != "BCL2"),]
tableD <- tableD[which(tableD$Locus == "IgH" | tableD$Target == "IgH"),]
tableD <- tableD[which(tableD$patient != "P3"),]

NAN_plot <- ggplot() + theme_classic() 
NAN_plot <- NAN_plot + geom_violin(data=tableD,aes(x=Locus,y=as.numeric(Min_CpG),fill=Locus),width=1.1,alpha=0.3,size=0.4) +
  geom_boxplot(data=tableD,aes(x=Locus,y=as.numeric(Min_CpG)),width=0.06,alpha=0.5,size=0.5)+
  geom_path(data=tableD,aes(x=Locus,y=Min_CpG,group=patient ),color="grey60",size=0.3) +
  geom_point(aes(x=Locus,y=as.numeric(Min_CpG),color=Tumor),alpha=.9,shape=20,data=tableD,size=3) 
NAN_plot <-NAN_plot + theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,2,2,3),'lines'),
                            plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
                            legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='none',legend.text=element_text(size=14,hjust=0,face='bold'),
                            axis.text.x=element_text(size=12,face='bold',color='black'),axis.text.y=element_text(size=12,face='bold',color='black'),
                            axis.title.x=element_text(size=14,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))
NAN_plot <- NAN_plot + scale_y_continuous(expand=c(0,0),limits=c(-0.5,17),breaks = seq(0,30,5)) 
#NAN_plot <- NAN_plot + scale_x_discrete(labels=c(B.FL="FL \ntumor",C.DHL="DHL \ntumor"))#
NAN_plot<- NAN_plot +ylab("Distance from nearest CpG to breakpoint") +xlab(NULL)
NAN_plot <- NAN_plot + scale_fill_manual(name=NULL,values = c("#bc80bd","#00FF0040"))
NAN_plot <- NAN_plot + scale_color_manual(name=NULL,values = c(gg_color_hue(6)[1],gg_color_hue(6)[1]))
the<-compare_means(Min_CpG ~ Locus,  data = tableD,paired = T )
#NAN_plot<- NAN_plot +stat_compare_means(comparisons = my_comparisons)
NAN_plot <- NAN_plot +stat_pvalue_manual(
  the, label = NULL, 
  y.position = c( 16)
)


figure_2<-rbind(ggplotGrob(NAN_plot),size="last")

ggsave(file="ExtFig1d_CpG_IgH_MYC.pdf", plot=figure_2,bg = 'white', width =11, height = 12, units = 'cm', dpi = 600)

