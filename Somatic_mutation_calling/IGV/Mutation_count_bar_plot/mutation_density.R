# WatchDHL
library(rstudioapi)
library("ggplot2")
library("gridExtra")
library(grid)
# Getting the path of your current open file
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))
print( getwd() )



gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}



#####
FL<- read.csv("FL.specific.seg",header = F,sep = "\t")
DHL<- read.csv("DHL.specific.seg",header = F,sep = "\t")
FL <- FL[which(FL$V2=="chr9" & as.numeric(FL$V3)>36742650 & as.numeric(FL$V3)<37618428 & FL$V1!="P3.FL"),]
DHL <- DHL[which(DHL$V2=="chr9" & as.numeric(DHL$V3) > 36742650 & as.numeric(DHL$V3) < 37618428 & DHL$V1!="P3.DHL"),]

gene.table <- data.frame(c(as.numeric(FL$V3)/1000000,as.numeric(DHL$V3)/1000000),c(rep("FL",nrow(FL)),rep("DHL",nrow(DHL))))
colnames(gene.table) <- c("Poition","Type")
F1a.plot<-ggplot()+theme_classic()
F1a.plot<-F1a.plot+geom_histogram(data=gene.table,aes(x=Poition,y=..count..,fill=Type),binwidth=0.003,alpha=0.5,position='identity',pad=T)
#F1a.plot<-F1a.plot+ geom_density(data=gene.table ,aes(x=Poition,color=Type),size=1,alpha=0.2)
F1a.plot<-F1a.plot+theme(panel.background=element_rect(fill='transparent',color='transparent',size=1),plot.margin=unit(c(1,1,1,1),'lines'),
                         plot.title=element_text(size=18,vjust=0.5,hjust=0.5,face='bold.italic',color='black'),text=element_text(size=14,face='bold'),
                         legend.key.width=unit(0.5,'cm'),legend.key.height=unit(0.5,'cm'),legend.position="none",legend.text=element_blank(),
                         axis.text.x=element_text(size=12,face='bold',color='black'),axis.text.y=element_text(size=12,face='bold',color='black'),axis.ticks.y.left = element_line(),
                         axis.title.x=element_blank(),axis.title.y=element_blank())
F1a.plot<-F1a.plot+ggtitle(NULL)#+xlab("Age of onset")+ylab('Density')

F1a.plot<-F1a.plot+scale_x_continuous(expand=c(0,0),limits=c(36.742650,37.618428),breaks=seq(0,38,0.2),position = "top")+scale_y_continuous(expand=c(0,0),limits = c(30,0),breaks = c(0,10,20,30),trans = "reverse")#,limits=c(0,0.18))
#+scale_y_continuous(expand=c(0,0),limits=c(0,0.15))+scale_x_continuous(expand=c(0,0),limits=c(0,26),breaks=seq(0,100,5))#
F1a.plot<-F1a.plot + scale_fill_manual(name=NULL,values=c(FL=gg_color_hue(6)[3],DHL=gg_color_hue(6)[1]),labels=c(FL='FL tumor',DHL='DLBCL tumor'))+
                     scale_color_manual(name=NULL,values=c(FL=gg_color_hue(6)[3],DHL=gg_color_hue(6)[1]),labels=c(FL='FL tumor',DHL='DLBCL tumor'))
F1a.plot

figure_1<-cbind(rbind(ggplotGrob(F1a.plot)))
ggsave(file="zcchc7_mutation_load_3kb_new.2022.pdf", plot=figure_1,bg = 'white', width = 38, height =  5, units = 'cm', dpi = 600)






