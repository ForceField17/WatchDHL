# WatchDHL
library(rstudioapi)
library("ggplot2")
library("gridExtra")
library(grid)
# Getting the path of your current open file
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))
print( getwd() )

library(matrixStats)
#source("../lib/CELLO.R")

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 120)[1:n]
}

#raw data preprocessing
FL <-   read.delim("./mutDens_FL.txt",sep = c("\t"," "),header = F)
FL2 <-  read.delim("./mutDens_FL_nege.txt",sep = c("\t"," "),header = F)
DHL <-  read.delim("./mutDens_DHL.txt",sep = c("\t"," "),header = F)
DHL2 <- read.delim("./mutDens_DHL_nege.txt",sep = c("\t"," "),header = F)

colnames(FL) <-   c("V1","V2","Distance","P1","P2","P4","P5","P6","P7","P8","P8.b","P9","interval","P1t","P2t","P4t","P5t","P6t","P7t","P8t","P8.bt","P9t","V23")
colnames(FL2) <-   c("V1","V2","Distance","P1","P2","P4","P5","P6","P7","P8","P8.b","P9","interval","P1t","P2t","P4t","P5t","P6t","P7t","P8t","P8.bt","P9t","V23")
colnames(DHL) <-  c("V1","V2","Distance","P1","P2","P4","P5","P6","P7","P8","P8.b","P9","interval","P1t","P2t","P4t","P5t","P6t","P7t","P8t","P9t","V22")
colnames(DHL2) <- c("V1","V2","Distance","P1","P2","P4","P5","P6","P7","P8","P8.b","P9","interval","P1t","P2t","P4t","P5t","P6t","P7t","P8t","P9t","V22")

FL_P1t <- FL$P8t[1]
DHL_P1t <- DHL$P8t[1]

FL$Burden.avg <- rowMeans( cbind( (FL$P8.b/FL$interval*1000000)))
FL$Burden.sd  <- rowSds(  cbind( (FL$P8.b/FL$interval*1000000)))

FL2$Burden.avg <- rowMeans( cbind( (FL2$P8.b/FL2$interval*1000000)))
FL2$Burden.sd  <- rowSds(  cbind( (FL2$P8.b/FL2$interval*1000000)))

DHL$Burden.avg <- rowMeans( cbind( (DHL$P8/DHL$interval*1000000)))
DHL$Burden.sd  <- rowSds(  cbind( (DHL$P8/DHL$interval*1000000)))

DHL2$Burden.avg <- rowMeans( cbind( (DHL2$P8/DHL2$interval*1000000)))
DHL2$Burden.sd  <- rowSds(  cbind( (DHL2$P8/DHL2$interval*1000000)))

list <- c("V1","Distance","Burden.avg","Burden.sd")
FL  <-  FL[, which(colnames(FL) %in% list)]
FL2 <-  FL2[,which(colnames(FL2) %in% list)]
DHL  <-  DHL[, which(colnames(DHL) %in% list)]
DHL2 <-  DHL2[,which(colnames(DHL2) %in% list)]
samples <- rbind(FL,FL2,DHL,DHL2)

All_Mut <- samples

All_Mut$xxxx <- All_Mut$Burden.avg -All_Mut$Burden.sd
All_Mut$xxxx[which(All_Mut$xxxx<0)] <- 0
cumAge <- ggplot()+theme_classic()
cumAge <- cumAge + geom_vline(xintercept = 0,color="grey",size=0.5,linetype=2)
cumAge <- cumAge + geom_line(data=All_Mut,aes(x=Distance,y=Burden.avg,color=V1), alpha=1,size=1.5)
cumAge <- cumAge + geom_errorbar(data=All_Mut,width=0,size=0.2, aes(x=Distance,ymin=xxxx, ymax=Burden.avg+Burden.sd,color=V1),alpha=0.2)
#cumAge <- cumAge + geom_point(data=All_Mut,aes(x=Distance,y=Burden*1000000,color=V1), alpha=0.85,size=1)
cumAge<- cumAge +theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,2,2,2),'lines'),
                       plot.title=element_text(size=15,vjust=0.5,hjust=0.5,face='bold.italic',color='black'),text=element_text(size=14,face='bold'),
                       legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.5,'cm'),legend.position="right",legend.text=element_text(size=12,hjust=0,face='bold.italic'),
                       axis.text.x=element_text(size=12,angle=30,vjust=1,hjust=1,face='bold',color='black'),axis.text.y=element_text(size=12,face='bold',color='black'),
                       axis.title.x=element_text(size=14,vjust=0,hjust=0.5,face="plain",color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))
cumAge <- cumAge+ggtitle(paste0("$P8 ","FL:",FL_P1t," DHL:",DHL_P1t))+xlab("Distance to TSS (bp)")+ylab('Normalized mutation density')+
  scale_y_continuous(expand=c(0,0),limits=c(0,8),breaks = seq(0,8,2)) +
  scale_x_continuous(expand=c(0,0),limits=c(-20000,20000),breaks = seq(-200000,200000,2000))
table2 <- data.frame(c(0,0),c(0,0.0118))
colnames(table2) <- c("Xaxis","Yaxis")
cumAge <- cumAge + geom_path( data=table2, aes(x=Xaxis,y=Yaxis),alpha=1,size=1,color="black",arrow = NULL)
table3 <- data.frame(c(0,1000),c(0.0115,0.0115))
colnames(table3) <- c("Xaxis","Yaxis")
cumAge <- cumAge + geom_path( data=table3, aes(x=Xaxis,y=Yaxis),alpha=1,size=1,color="black",arrow = arrow(length = unit(0.18, "cm"),type = "closed"))
cumAge <- cumAge + scale_colour_manual(name=NULL,values=c(DHL=gg_color_hue(6)[1] ,FL=gg_color_hue(6)[3])) 

figure_1<-rbind(ggplotGrob(cumAge),size="first")
ggsave(file="new_ab_$P8b_normalized_MB.pdf", plot=figure_1,bg = 'white', width = 28, height = 14, units = 'cm', dpi = 600)

