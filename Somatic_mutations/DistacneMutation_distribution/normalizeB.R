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
  hcl(h = hues, l = 65, c = 120)[1:n]
}

#raw data preprocessing
FL <-   read.delim("Slide_Only_H3K27ac/mutDens_FL.txt",sep = c("\t"," "),header = F)
FL2 <-  read.delim("Slide_Only_H3K27ac/mutDens_FL_nege.txt",sep = c("\t"," "),header = F)
DHL <-  read.delim("Slide_Only_H3K27ac/mutDens_DHL.txt",sep = c("\t"," "),header = F)
DHL2 <- read.delim("Slide_Only_H3K27ac/mutDens_DHL_nege.txt",sep = c("\t"," "),header = F)

colnames(FL) <-   c("V1","V2","Distance","P1","P2","P4","P5","P6","P7","P8","P8.b","P9","interval","P1t","P2t","P4t","P5t","P6t","P7t","P8t","P8.bt","P9t","V23")
colnames(FL2) <-   c("V1","V2","Distance","P1","P2","P4","P5","P6","P7","P8","P8.b","P9","interval","P1t","P2t","P4t","P5t","P6t","P7t","P8t","P8.bt","P9t","V23")
colnames(DHL) <-  c("V1","V2","Distance","P1","P2","P4","P5","P6","P7","P8","P8.b","P9","interval","P1t","P2t","P4t","P5t","P6t","P7t","P8t","P9t","V22")
colnames(DHL2) <- c("V1","V2","Distance","P1","P2","P4","P5","P6","P7","P8","P8.b","P9","interval","P1t","P2t","P4t","P5t","P6t","P7t","P8t","P9t","V22")


FL$Burden.avg <- rowMeans( cbind( (FL$P1/FL$interval/FL$P1t*100*1000000), (FL$P2/FL$interval/FL$P2t*100*1000000), (FL$P6/FL$interval/FL$P6t*100*1000000), (FL$P7/FL$interval/FL$P7t*100*1000000), (FL$P8/FL$interval/FL$P8t*100*1000000), (FL$P8.b/FL$interval/FL$P8.bt*100*1000000), (FL$P9/FL$interval/FL$P9t*100*1000000)))
FL$Burden.sd  <- rowSds(  cbind( (FL$P1/FL$interval/FL$P1t*100*1000000), (FL$P2/FL$interval/FL$P2t*100*1000000), (FL$P6/FL$interval/FL$P6t*100*1000000), (FL$P7/FL$interval/FL$P7t*100*1000000), (FL$P8/FL$interval/FL$P8t*100*1000000), (FL$P8.b/FL$interval/FL$P8.bt*100*1000000), (FL$P9/FL$interval/FL$P9t*100*1000000)))

FL2$Burden.avg <- rowMeans( cbind( (FL2$P1/FL2$interval/FL2$P1t*100*1000000), (FL2$P2/FL2$interval/FL2$P2t*100*1000000), (FL2$P6/FL2$interval/FL2$P6t*100*1000000), (FL2$P7/FL2$interval/FL2$P7t*100*1000000), (FL2$P8/FL2$interval/FL2$P8t*100*1000000), (FL2$P8.b/FL2$interval/FL2$P8.bt*100*1000000), (FL2$P9/FL2$interval/FL2$P9t*100*1000000)))
FL2$Burden.sd  <- rowSds(  cbind( (FL2$P1/FL2$interval/FL2$P1t*100*1000000), (FL2$P2/FL2$interval/FL2$P2t*100*1000000), (FL2$P6/FL2$interval/FL2$P6t*100*1000000), (FL2$P7/FL2$interval/FL2$P7t*100*1000000), (FL2$P8/FL2$interval/FL2$P8t*100*1000000), (FL2$P8.b/FL2$interval/FL2$P8.bt*100*1000000), (FL2$P9/FL2$interval/FL2$P9t*100*1000000)))

DHL$Burden.avg <- rowMeans( cbind( (DHL$P1/DHL$interval/DHL$P1t*100*1000000), (DHL$P2/DHL$interval/DHL$P2t*100*1000000), (DHL$P6/DHL$interval/DHL$P6t*100*1000000), (DHL$P7/DHL$interval/DHL$P7t*100*1000000), (DHL$P8/DHL$interval/DHL$P8t*100*1000000), (DHL$P9/DHL$interval/DHL$P9t*100*1000000)))
DHL$Burden.sd  <- rowSds(  cbind( (DHL$P1/DHL$interval/DHL$P1t*100*1000000), (DHL$P2/DHL$interval/DHL$P2t*100*1000000), (DHL$P6/DHL$interval/DHL$P6t*100*1000000), (DHL$P7/DHL$interval/DHL$P7t*100*1000000), (DHL$P8/DHL$interval/DHL$P8t*100*1000000), (DHL$P9/DHL$interval/DHL$P9t*100*1000000)))

DHL2$Burden.avg <- rowMeans( cbind( (DHL2$P1/DHL2$interval/DHL2$P1t*100*1000000), (DHL2$P2/DHL2$interval/DHL2$P2t*100*1000000), (DHL2$P6/DHL2$interval/DHL2$P6t*100*1000000), (DHL2$P7/DHL2$interval/DHL2$P7t*100*1000000), (DHL2$P8/DHL2$interval/DHL2$P8t*100*1000000), (DHL2$P9/DHL2$interval/DHL2$P9t*100*1000000)))
DHL2$Burden.sd  <- rowSds(  cbind( (DHL2$P1/DHL2$interval/DHL2$P1t*100*1000000), (DHL2$P2/DHL2$interval/DHL2$P2t*100*1000000), (DHL2$P6/DHL2$interval/DHL2$P6t*100*1000000), (DHL2$P7/DHL2$interval/DHL2$P7t*100*1000000), (DHL2$P8/DHL2$interval/DHL2$P8t*100*1000000), (DHL2$P9/DHL2$interval/DHL2$P9t*100*1000000)))

list <- c("V1","Distance","Burden.avg","Burden.sd")
FL  <-  FL[, which(colnames(FL) %in% list)]
FL2 <-  FL2[,which(colnames(FL2) %in% list)]
DHL  <-  DHL[, which(colnames(DHL) %in% list)]
DHL2 <-  DHL2[,which(colnames(DHL2) %in% list)]
samples <- rbind(FL,FL2,DHL,DHL2)

myOnly1 <- cbind(rbind(FL,FL2),rbind(DHL,DHL2))
write.table(myOnly1,file ="./Source_Only1.txt", sep = "\t",quote = F,row.names = F)


All_Mut <- samples

All_Mut$xxxx <- All_Mut$Burden.avg -All_Mut$Burden.sd
All_Mut$xxxx[which(All_Mut$xxxx<0)] <- 0
cumAge <- ggplot()+theme_classic()
cumAge <- cumAge + geom_vline(xintercept = 0,color="grey",size=0.4,linetype=2)
cumAge <- cumAge + geom_line(data=All_Mut,aes(x=Distance,y=Burden.avg,color=V1), alpha=1,size=0.7)
cumAge <- cumAge + geom_errorbar(data=All_Mut,width=0,size=0.2, aes(x=Distance,ymin=xxxx, ymax=Burden.avg+Burden.sd,color=V1),alpha=0.2)
#cumAge <- cumAge + geom_point(data=All_Mut,aes(x=Distance,y=Burden*1000000,color=V1), alpha=0.85,size=1)
cumAge<- cumAge +theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,0.9,2,2),'lines'),
                       plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
                       legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.5,'cm'),legend.position="top",legend.text=element_text(size=12,hjust=0,face='bold.italic'),
                       axis.text.x=element_text(size=12,angle=30,vjust=1,hjust=1,face='bold',color='black'),axis.text.y=element_text(size=10,face='bold',color='black'),
                       axis.title.x=element_text(size=14,vjust=0,hjust=0.5,face="plain",color='black'),axis.title.y=element_blank())
cumAge <- cumAge+ggtitle(NULL)+xlab("Only H3K27ac")+ylab('Average proportion of somatic mutations (%/Mbp)')+
  scale_y_continuous(expand=c(0,0),limits=c(0,0.71),breaks = seq(0,1,0.1)) +coord_cartesian(ylim = c(0, 0.61))+
  scale_x_continuous(expand=c(0,0),limits=c(-100000,100000),breaks = seq(-200000,200000,40000))
table2 <- data.frame(c(0,0),c(0,0.0088))
colnames(table2) <- c("Xaxis","Yaxis")
#cumAge <- cumAge + geom_path( data=table2, aes(x=Xaxis,y=Yaxis),alpha=1,size=1,color="black",arrow = NULL)
table3 <- data.frame(c(0,1000),c(0.0085,0.0085))
colnames(table3) <- c("Xaxis","Yaxis")
#cumAge <- cumAge + geom_path( data=table3, aes(x=Xaxis,y=Yaxis),alpha=1,size=1,color="black",arrow = arrow(length = unit(0.18, "cm"),type = "closed"))
cumAge <- cumAge + scale_colour_manual(name=NULL,values=c(DHL="#d01c8b" ,FL=gg_color_hue(6)[3])) 



#raw data preprocessing
FL <-   read.delim("Slide_Only_H3K4me3/mutDens_FL.txt",sep = c("\t"," "),header = F)
FL2 <-  read.delim("Slide_Only_H3K4me3/mutDens_FL_nege.txt",sep = c("\t"," "),header = F)
DHL <-  read.delim("Slide_Only_H3K4me3/mutDens_DHL.txt",sep = c("\t"," "),header = F)
DHL2 <- read.delim("Slide_Only_H3K4me3/mutDens_DHL_nege.txt",sep = c("\t"," "),header = F)

colnames(FL) <-   c("V1","V2","Distance","P1","P2","P4","P5","P6","P7","P8","P8.b","P9","interval","P1t","P2t","P4t","P5t","P6t","P7t","P8t","P8.bt","P9t","V23")
colnames(FL2) <-   c("V1","V2","Distance","P1","P2","P4","P5","P6","P7","P8","P8.b","P9","interval","P1t","P2t","P4t","P5t","P6t","P7t","P8t","P8.bt","P9t","V23")
colnames(DHL) <-  c("V1","V2","Distance","P1","P2","P4","P5","P6","P7","P8","P8.b","P9","interval","P1t","P2t","P4t","P5t","P6t","P7t","P8t","P9t","V22")
colnames(DHL2) <- c("V1","V2","Distance","P1","P2","P4","P5","P6","P7","P8","P8.b","P9","interval","P1t","P2t","P4t","P5t","P6t","P7t","P8t","P9t","V22")


FL$Burden.avg <- rowMeans( cbind( (FL$P1/FL$interval/FL$P1t*100*1000000), (FL$P2/FL$interval/FL$P2t*100*1000000), (FL$P6/FL$interval/FL$P6t*100*1000000), (FL$P7/FL$interval/FL$P7t*100*1000000), (FL$P8/FL$interval/FL$P8t*100*1000000), (FL$P8.b/FL$interval/FL$P8.bt*100*1000000), (FL$P9/FL$interval/FL$P9t*100*1000000)))
FL$Burden.sd  <- rowSds(  cbind( (FL$P1/FL$interval/FL$P1t*100*1000000), (FL$P2/FL$interval/FL$P2t*100*1000000), (FL$P6/FL$interval/FL$P6t*100*1000000), (FL$P7/FL$interval/FL$P7t*100*1000000), (FL$P8/FL$interval/FL$P8t*100*1000000), (FL$P8.b/FL$interval/FL$P8.bt*100*1000000), (FL$P9/FL$interval/FL$P9t*100*1000000)))

FL2$Burden.avg <- rowMeans( cbind( (FL2$P1/FL2$interval/FL2$P1t*100*1000000), (FL2$P2/FL2$interval/FL2$P2t*100*1000000), (FL2$P6/FL2$interval/FL2$P6t*100*1000000), (FL2$P7/FL2$interval/FL2$P7t*100*1000000), (FL2$P8/FL2$interval/FL2$P8t*100*1000000), (FL2$P8.b/FL2$interval/FL2$P8.bt*100*1000000), (FL2$P9/FL2$interval/FL2$P9t*100*1000000)))
FL2$Burden.sd  <- rowSds(  cbind( (FL2$P1/FL2$interval/FL2$P1t*100*1000000), (FL2$P2/FL2$interval/FL2$P2t*100*1000000), (FL2$P6/FL2$interval/FL2$P6t*100*1000000), (FL2$P7/FL2$interval/FL2$P7t*100*1000000), (FL2$P8/FL2$interval/FL2$P8t*100*1000000), (FL2$P8.b/FL2$interval/FL2$P8.bt*100*1000000), (FL2$P9/FL2$interval/FL2$P9t*100*1000000)))

DHL$Burden.avg <- rowMeans( cbind( (DHL$P1/DHL$interval/DHL$P1t*100*1000000), (DHL$P2/DHL$interval/DHL$P2t*100*1000000), (DHL$P6/DHL$interval/DHL$P6t*100*1000000), (DHL$P7/DHL$interval/DHL$P7t*100*1000000), (DHL$P8/DHL$interval/DHL$P8t*100*1000000), (DHL$P9/DHL$interval/DHL$P9t*100*1000000)))
DHL$Burden.sd  <- rowSds(  cbind( (DHL$P1/DHL$interval/DHL$P1t*100*1000000), (DHL$P2/DHL$interval/DHL$P2t*100*1000000), (DHL$P6/DHL$interval/DHL$P6t*100*1000000), (DHL$P7/DHL$interval/DHL$P7t*100*1000000), (DHL$P8/DHL$interval/DHL$P8t*100*1000000), (DHL$P9/DHL$interval/DHL$P9t*100*1000000)))

DHL2$Burden.avg <- rowMeans( cbind( (DHL2$P1/DHL2$interval/DHL2$P1t*100*1000000), (DHL2$P2/DHL2$interval/DHL2$P2t*100*1000000), (DHL2$P6/DHL2$interval/DHL2$P6t*100*1000000), (DHL2$P7/DHL2$interval/DHL2$P7t*100*1000000), (DHL2$P8/DHL2$interval/DHL2$P8t*100*1000000), (DHL2$P9/DHL2$interval/DHL2$P9t*100*1000000)))
DHL2$Burden.sd  <- rowSds(  cbind( (DHL2$P1/DHL2$interval/DHL2$P1t*100*1000000), (DHL2$P2/DHL2$interval/DHL2$P2t*100*1000000), (DHL2$P6/DHL2$interval/DHL2$P6t*100*1000000), (DHL2$P7/DHL2$interval/DHL2$P7t*100*1000000), (DHL2$P8/DHL2$interval/DHL2$P8t*100*1000000), (DHL2$P9/DHL2$interval/DHL2$P9t*100*1000000)))

list <- c("V1","Distance","Burden.avg","Burden.sd")
FL  <-  FL[, which(colnames(FL) %in% list)]
FL2 <-  FL2[,which(colnames(FL2) %in% list)]
DHL  <-  DHL[, which(colnames(DHL) %in% list)]
DHL2 <-  DHL2[,which(colnames(DHL2) %in% list)]
samples <- rbind(FL,FL2,DHL,DHL2)

myOnly2 <- cbind(rbind(FL,FL2),rbind(DHL,DHL2))
write.table(myOnly2,file ="./Source_Only2.txt", sep = "\t",quote = F,row.names = F)


All_Mut <- samples

All_Mut$xxxx <- All_Mut$Burden.avg -All_Mut$Burden.sd
All_Mut$xxxx[which(All_Mut$xxxx<0)] <- 0
cumAge2 <- ggplot()+theme_classic()
cumAge2 <- cumAge2 + geom_vline(xintercept = 0,color="grey",size=0.4,linetype=2)
cumAge2 <- cumAge2 + geom_line(data=All_Mut,aes(x=Distance,y=Burden.avg,color=V1), alpha=1,size=0.7)
cumAge2 <- cumAge2 + geom_errorbar(data=All_Mut,width=0,size=0.2, aes(x=Distance,ymin=xxxx, ymax=Burden.avg+Burden.sd,color=V1),alpha=0.2)
#cumAge2 <- cumAge2 + geom_point(data=All_Mut,aes(x=Distance,y=Burden*1000000,color=V1), alpha=0.85,size=1)
cumAge2<- cumAge2 +theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,0.9,2,1),'lines'),
                         plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
                         legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.5,'cm'),legend.position="top",legend.text=element_text(size=12,hjust=0,face='bold.italic'),
                         axis.text.x=element_text(size=12,angle=30,vjust=1,hjust=1,face='bold',color='black'),axis.text.y=element_blank(),
                         axis.title.x=element_text(size=14,vjust=0,hjust=0.5,face="plain",color='black'),axis.title.y=element_blank())
cumAge2 <- cumAge2+ggtitle(NULL)+xlab("Only_H3K4me3")+ylab('Average proportion of somatic mutations (%/Mbp)')+
  scale_y_continuous(expand=c(0,0),limits=c(0,0.71),breaks = seq(0,1,0.1)) +coord_cartesian(ylim = c(0, 0.61))+
  scale_x_continuous(expand=c(0,0),limits=c(-100000,100000),breaks = seq(-200000,200000,40000))
table2 <- data.frame(c(0,0),c(0,0.0088))
colnames(table2) <- c("Xaxis","Yaxis")
#cumAge2 <- cumAge2 + geom_path( data=table2, aes(x=Xaxis,y=Yaxis),alpha=1,size=1,color="black",arrow = NULL)
table3 <- data.frame(c(0,1000),c(0.0085,0.0085))
colnames(table3) <- c("Xaxis","Yaxis")
#cumAge2 <- cumAge2 + geom_path( data=table3, aes(x=Xaxis,y=Yaxis),alpha=1,size=1,color="black",arrow = arrow(length = unit(0.18, "cm"),type = "closed"))
cumAge2 <- cumAge2 + scale_colour_manual(name=NULL,values=c(DHL="#d01c8b" ,FL=gg_color_hue(6)[3])) 


#raw data preprocessing
FL <-   read.delim("Slide_Both_H3K4me3_and_H3K27ac/mutDens_FL.txt",sep = c("\t"," "),header = F)
FL2 <-  read.delim("Slide_Both_H3K4me3_and_H3K27ac/mutDens_FL_nege.txt",sep = c("\t"," "),header = F)
DHL <-  read.delim("Slide_Both_H3K4me3_and_H3K27ac/mutDens_DHL.txt",sep = c("\t"," "),header = F)
DHL2 <- read.delim("Slide_Both_H3K4me3_and_H3K27ac/mutDens_DHL_nege.txt",sep = c("\t"," "),header = F)

colnames(FL) <-   c("V1","V2","Distance","P1","P2","P4","P5","P6","P7","P8","P8.b","P9","interval","P1t","P2t","P4t","P5t","P6t","P7t","P8t","P8.bt","P9t","V23")
colnames(FL2) <-   c("V1","V2","Distance","P1","P2","P4","P5","P6","P7","P8","P8.b","P9","interval","P1t","P2t","P4t","P5t","P6t","P7t","P8t","P8.bt","P9t","V23")
colnames(DHL) <-  c("V1","V2","Distance","P1","P2","P4","P5","P6","P7","P8","P8.b","P9","interval","P1t","P2t","P4t","P5t","P6t","P7t","P8t","P9t","V22")
colnames(DHL2) <- c("V1","V2","Distance","P1","P2","P4","P5","P6","P7","P8","P8.b","P9","interval","P1t","P2t","P4t","P5t","P6t","P7t","P8t","P9t","V22")


FL$Burden.avg <- rowMeans( cbind( (FL$P1/FL$interval/FL$P1t*100*1000000), (FL$P2/FL$interval/FL$P2t*100*1000000), (FL$P6/FL$interval/FL$P6t*100*1000000), (FL$P7/FL$interval/FL$P7t*100*1000000), (FL$P8/FL$interval/FL$P8t*100*1000000), (FL$P8.b/FL$interval/FL$P8.bt*100*1000000), (FL$P9/FL$interval/FL$P9t*100*1000000)))
FL$Burden.sd  <- rowSds(  cbind( (FL$P1/FL$interval/FL$P1t*100*1000000), (FL$P2/FL$interval/FL$P2t*100*1000000), (FL$P6/FL$interval/FL$P6t*100*1000000), (FL$P7/FL$interval/FL$P7t*100*1000000), (FL$P8/FL$interval/FL$P8t*100*1000000), (FL$P8.b/FL$interval/FL$P8.bt*100*1000000), (FL$P9/FL$interval/FL$P9t*100*1000000)))

FL2$Burden.avg <- rowMeans( cbind( (FL2$P1/FL2$interval/FL2$P1t*100*1000000), (FL2$P2/FL2$interval/FL2$P2t*100*1000000), (FL2$P6/FL2$interval/FL2$P6t*100*1000000), (FL2$P7/FL2$interval/FL2$P7t*100*1000000), (FL2$P8/FL2$interval/FL2$P8t*100*1000000), (FL2$P8.b/FL2$interval/FL2$P8.bt*100*1000000), (FL2$P9/FL2$interval/FL2$P9t*100*1000000)))
FL2$Burden.sd  <- rowSds(  cbind( (FL2$P1/FL2$interval/FL2$P1t*100*1000000), (FL2$P2/FL2$interval/FL2$P2t*100*1000000), (FL2$P6/FL2$interval/FL2$P6t*100*1000000), (FL2$P7/FL2$interval/FL2$P7t*100*1000000), (FL2$P8/FL2$interval/FL2$P8t*100*1000000), (FL2$P8.b/FL2$interval/FL2$P8.bt*100*1000000), (FL2$P9/FL2$interval/FL2$P9t*100*1000000)))

DHL$Burden.avg <- rowMeans( cbind( (DHL$P1/DHL$interval/DHL$P1t*100*1000000), (DHL$P2/DHL$interval/DHL$P2t*100*1000000), (DHL$P6/DHL$interval/DHL$P6t*100*1000000), (DHL$P7/DHL$interval/DHL$P7t*100*1000000), (DHL$P8/DHL$interval/DHL$P8t*100*1000000), (DHL$P9/DHL$interval/DHL$P9t*100*1000000)))
DHL$Burden.sd  <- rowSds(  cbind( (DHL$P1/DHL$interval/DHL$P1t*100*1000000), (DHL$P2/DHL$interval/DHL$P2t*100*1000000), (DHL$P6/DHL$interval/DHL$P6t*100*1000000), (DHL$P7/DHL$interval/DHL$P7t*100*1000000), (DHL$P8/DHL$interval/DHL$P8t*100*1000000), (DHL$P9/DHL$interval/DHL$P9t*100*1000000)))

DHL2$Burden.avg <- rowMeans( cbind( (DHL2$P1/DHL2$interval/DHL2$P1t*100*1000000), (DHL2$P2/DHL2$interval/DHL2$P2t*100*1000000), (DHL2$P6/DHL2$interval/DHL2$P6t*100*1000000), (DHL2$P7/DHL2$interval/DHL2$P7t*100*1000000), (DHL2$P8/DHL2$interval/DHL2$P8t*100*1000000), (DHL2$P9/DHL2$interval/DHL2$P9t*100*1000000)))
DHL2$Burden.sd  <- rowSds(  cbind( (DHL2$P1/DHL2$interval/DHL2$P1t*100*1000000), (DHL2$P2/DHL2$interval/DHL2$P2t*100*1000000), (DHL2$P6/DHL2$interval/DHL2$P6t*100*1000000), (DHL2$P7/DHL2$interval/DHL2$P7t*100*1000000), (DHL2$P8/DHL2$interval/DHL2$P8t*100*1000000), (DHL2$P9/DHL2$interval/DHL2$P9t*100*1000000)))

list <- c("V1","Distance","Burden.avg","Burden.sd")
FL  <-  FL[, which(colnames(FL) %in% list)]
FL2 <-  FL2[,which(colnames(FL2) %in% list)]
DHL  <-  DHL[, which(colnames(DHL) %in% list)]
DHL2 <-  DHL2[,which(colnames(DHL2) %in% list)]

samples <- rbind(FL,FL2,DHL,DHL2)

myxx <- cbind(rbind(FL,FL2),rbind(DHL,DHL2))
write.table(myxx,file ="./Source_xx.txt", sep = "\t",quote = F,row.names = F)


All_Mut <- samples

All_Mut$xxxx <- All_Mut$Burden.avg -All_Mut$Burden.sd
All_Mut$xxxx[which(All_Mut$xxxx<0)] <- 0
All_Mut$yyyy <- All_Mut$Burden.avg +All_Mut$Burden.sd
All_Mut$yyyy[which(All_Mut$yyyy>0.61)] <- 0.61
cumAge3 <- ggplot()+theme_classic()
cumAge3 <- cumAge3 + geom_vline(xintercept = 0,color="grey",size=0.4,linetype=2)
cumAge3 <- cumAge3 + geom_line(data=All_Mut,aes(x=Distance,y=Burden.avg,color=V1), alpha=1,size=0.7)
cumAge3 <- cumAge3 + geom_errorbar(data=All_Mut,width=0,size=0.2, aes(x=Distance,ymin=xxxx, ymax=yyyy,color=V1),alpha=0.2)
#cumAge3 <- cumAge3 + geom_point(data=All_Mut,aes(x=Distance,y=Burden*1000000,color=V1), alpha=0.85,size=1)
cumAge3<- cumAge3 +theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,2,2,1),'lines'),
                         plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
                         legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.5,'cm'),legend.position="top",legend.text=element_text(size=12,hjust=0,face='bold.italic'),
                         axis.text.x=element_text(size=12,angle=30,vjust=1,hjust=1,face='bold',color='black'),axis.text.y=element_blank(),
                         axis.title.x=element_text(size=14,vjust=0,hjust=0.5,face="plain",color='black'),axis.title.y=element_blank())
cumAge3 <- cumAge3+ggtitle(NULL)+xlab("Both H3K4me3 H3K27ac")+ylab('Average proportion of somatic mutations (%/Mbp)')+
  scale_y_continuous(expand=c(0,0),limits=c(0,0.71),breaks = seq(0,1,0.1)) +coord_cartesian(ylim = c(0, 0.61))+
  scale_x_continuous(expand=c(0,0),limits=c(-100000,100000),breaks = seq(-200000,200000,40000))
table2 <- data.frame(c(0,0),c(0,0.0088))
colnames(table2) <- c("Xaxis","Yaxis")
#cumAge3 <- cumAge3 + geom_path( data=table2, aes(x=Xaxis,y=Yaxis),alpha=1,size=1,color="black",arrow = NULL)
table3 <- data.frame(c(0,1000),c(0.0085,0.0085))
colnames(table3) <- c("Xaxis","Yaxis")
#cumAge3 <- cumAge3 + geom_path( data=table3, aes(x=Xaxis,y=Yaxis),alpha=1,size=1,color="black",arrow = arrow(length = unit(0.18, "cm"),type = "closed"))
cumAge3 <- cumAge3 + scale_colour_manual(name=NULL,values=c(DHL="#d01c8b" ,FL=gg_color_hue(6)[3])) 

figure_1<-cbind(ggplotGrob(cumAge),ggplotGrob(cumAge2),ggplotGrob(cumAge3),size="first")
ggsave(file="Fig4d_Slide_B_normalized_MB.pdf", plot=figure_1,bg = 'white', width = 20, height = 10, units = 'cm', dpi = 600)

