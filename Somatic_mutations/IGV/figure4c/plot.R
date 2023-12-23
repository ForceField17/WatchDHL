# WatchDHL
library(rstudioapi)

# Getting the path of your current open file
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))
print( getwd() )

library(ggplot2)
library(gridExtra)
library(grid)
library(rstatix)
library(ggpubr)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 150)[1:n]
}


DHL <- read.table("../../SAVI/Filtering/Final_results/DHL.mutation.txt",sep = "\t",header = T)
DHL <- DHL[,c(1,2,3,4,8)]
DHL$position <- (DHL$start + DHL$end + 1)/2
DHL$Index <- 10-as.integer(as.factor(DHL$CaseID))

FL <- read.table("../../SAVI/Filtering/Final_results/FL.mutation.txt",sep = "\t",header = T)
FL <- FL[,c(1,2,3,4,8)]
FL$position <- (FL$start + FL$end + 1)/2
FL$Index <- 19-as.integer(as.factor(FL$CaseID))

Mut <- rbind(FL,DHL)

################################################## chr16:10853004-10939884
table1 <- Mut[which(Mut$chromosome=="chr16" & Mut$start >= 10853004 & Mut$end <= 10939884),]
table1$Range2 <- (table1$end - table1$start)
table1$Range2 <- -(10853004-10939884)/400
table1$Range <- rowMaxs(cbind(table1$Range1,table1$Range2))


histology.plot<-ggplot(data=table1,aes(x=position,y=Index,width=Range,fill=mutType)) +
  geom_hline(yintercept=seq(min(Mut$Index)-0.4,max(Mut$Index)+0.6,1),color="grey",size=0.25)+
  geom_hline(yintercept=c(min(Mut$Index)-0.4,min(Mut$Index)-0.4+8,max(Mut$Index)+0.6),color="black",size=0.5)+
  geom_tile(data=table1,aes(x=position,y=Index,width=Range,fill=mutType),color=NA,alpha=1,height=0.8,stat='identity')+
  xlab(NULL)+ylab(NULL) + 
  theme_classic()+
  theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(1,1,1,1),'lines'),
        plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
        text=element_text(size=2,face='bold'),legend.key.width=unit(0.5,'cm'),legend.key.height=unit(0.5,'cm'),legend.position="none",legend.title=element_blank(),axis.ticks.y = element_blank(),
        legend.margin=margin(t=0.1,r=0.1,b=0,l=0.1,unit='cm'),legend.text=element_text(size=12,face='plain',color='black'),axis.text.y=element_text(size=12,face='plain',color='black'),axis.line = element_blank(),
        axis.text.x=element_text(size=12,face='plain',color='black'),axis.title.x=element_text(size=18,vjust=-1,face='plain',color='black'),axis.title.y=element_text(size=16,hjust=0.5,vjust=2,face='plain',color='black'))

histology.plot<-histology.plot + scale_fill_manual(name="Mutation type",values=c(SNV="black",Ins="black",Del="black",IgH="red",BCL2="#7133FF"),guide = guide_legend(override.aes=list(size=3),nrow=4),na.translate = F)
histology.plot<-histology.plot + scale_x_continuous(expand=c(0,0),limits = c(10853004,10939884),breaks = seq(0,100000000,1000000))
histology.plot<-histology.plot + scale_y_continuous(expand=c(0,0),limits = c(min(Mut$Index)-0.5,max(Mut$Index)+0.7),breaks = unique(Mut$Index),labels = unique(Mut$CaseID))

figure_1<-rbind(ggplotGrob(histology.plot),size="first")
ggsave(file="Figure4c_1.pdf", plot=figure_1,bg = 'white', width = 22, height = 10, units = 'cm', dpi = 600)


################################################## chr16:85888587-85923527
table1 <- Mut[which(Mut$chromosome=="chr16" & Mut$start >= 85888587 & Mut$end <= 85923527),]
table1$Range1 <- (table1$end - table1$start)
table1$Range2 <- -(85888587-85923527)/400
table1$Range <- rowMaxs(cbind(table1$Range1,table1$Range2))


histology.plot<-ggplot(data=table1,aes(x=position,y=Index,width=Range,fill=mutType)) +
  geom_hline(yintercept=seq(min(Mut$Index)-0.4,max(Mut$Index)+0.6,1),color="grey",size=0.25)+
  geom_hline(yintercept=c(min(Mut$Index)-0.4,min(Mut$Index)-0.4+8,max(Mut$Index)+0.6),color="black",size=0.5)+
  geom_tile(data=table1,aes(x=position,y=Index,width=Range,fill=mutType),color=NA,alpha=1,height=0.8,stat='identity')+
  xlab(NULL)+ylab(NULL) + 
  theme_classic()+
  theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(1,1,1,1),'lines'),
        plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
        text=element_text(size=2,face='bold'),legend.key.width=unit(0.5,'cm'),legend.key.height=unit(0.5,'cm'),legend.position="none",legend.title=element_blank(),axis.ticks.y = element_blank(),
        legend.margin=margin(t=0.1,r=0.1,b=0,l=0.1,unit='cm'),legend.text=element_text(size=12,face='plain',color='black'),axis.text.y=element_text(size=12,face='plain',color='black'),axis.line = element_blank(),
        axis.text.x=element_text(size=12,face='plain',color='black'),axis.title.x=element_text(size=18,vjust=-1,face='plain',color='black'),axis.title.y=element_text(size=16,hjust=0.5,vjust=2,face='plain',color='black'))

histology.plot<-histology.plot + scale_fill_manual(name="Mutation type",values=c(SNV="black",Ins="black",Del="black",IgH="red",BCL2=gg_color_hue(5)[4]),guide = guide_legend(override.aes=list(size=3),nrow=4),na.translate = F)
histology.plot<-histology.plot + scale_x_continuous(expand=c(0,0),limits = c(85888587,85923527),breaks = seq(0,100000000,1000000))
histology.plot<-histology.plot + scale_y_continuous(expand=c(0,0),limits = c(min(Mut$Index)-0.5,max(Mut$Index)+0.7),breaks = unique(Mut$Index),labels = unique(Mut$CaseID))

figure_1<-rbind(ggplotGrob(histology.plot),size="first")
ggsave(file="Figure4c_2.pdf", plot=figure_1,bg = 'white',  width = 22, height = 10, units = 'cm', dpi = 600)



