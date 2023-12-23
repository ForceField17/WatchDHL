# WatchDHL
library(rstudioapi)

# Getting the path of your current open file
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))
print( getwd() )

library(ggplot2)
library(DESeq2)
library(stringr)
library(rstatix)
library(ggpubr)
library(ggrepel)


gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 150)[1:n]
}

scale_rows = function(x){
  m = apply(x, 1, mean, na.rm = T)
  s = apply(x, 1, sd, na.rm = T)
  return((x - m) / s)
}

data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

#get CNV cases

#raw data preprocessing
Sample_features <- read.delim("Data/WT.txt",sep = "\t",header = T)
Sample_features <- Sample_features[which(Sample_features$Fraction.Number<13),]
theCoord <- Sample_features[which(Sample_features$Fraction.Number>=0 & Sample_features$Fraction.Number<12),]

mySub2 <- ggplot()+theme_classic()
mySub2 <- mySub2 + geom_line(data = Sample_features,   aes(x = Position, y = Absorbance), color="#6666F6",alpha=0.9,size=1) 
#mySub2 <- mySub2 + geom_point(data = theCoord,   aes(x = Position, y = Absorbance), fill="grey40",alpha=0.9,size=1, shape = 21,stroke=0.2) 
mySub2 <- mySub2 + xlab("Sucrose, 10~50% density gradient") + ylab("Absorbance, 260 nm") +
  scale_y_continuous(expand=c(0,0),limits=c(0,5.8),breaks=seq(-6,6,1)) + coord_cartesian(ylim = c(0,3))+
  scale_x_continuous(expand=c(0,0),limits=c(-1.5,78.2),breaks=theCoord$Position,labels = theCoord$Fraction.Number) 
mySub2 <- mySub2 + theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,1,0.5,1),'lines'),#axis.line.y = element_blank(),
                         legend.key.width=unit(0.3,'cm'),legend.key.height=unit(1,'cm'),legend.position='right',legend.text=element_text(size=12),
                         axis.text.y=element_text(size=10,face='plain',color='black'),axis.text.x=element_text(size=10,face='plain',color='black'),
                         axis.title.x=element_text(size=12,face='plain',color='black'),axis.title.y=element_text(size=12,hjust=0.5,vjust=2,face='plain',color='black'))
figure_4<-rbind(ggplotGrob(mySub2),size="last")
ggsave(file="./Fig7f_WT_ploysome.pdf", plot=figure_4,bg = 'white', width = 12, height = 9, units = 'cm', dpi = 600)


#raw data preprocessing
Sample_features <- read.delim("Data/HDR.txt",sep = "\t",header = T)
Sample_features <- Sample_features[which(Sample_features$Fraction.Number<13),]
theCoord <- Sample_features[which(Sample_features$Fraction.Number>=0 & Sample_features$Fraction.Number<12),]

mySub2 <- ggplot()+theme_classic()
mySub2 <- mySub2 + geom_line(data = Sample_features,   aes(x = Position, y = Absorbance), color="#fc8d62",alpha=0.9,size=1) 
#mySub2 <- mySub2 + geom_point(data = theCoord,   aes(x = Position, y = Absorbance), fill="#de77ae",alpha=0.9,size=1, shape = 21,stroke=0.2) 
mySub2 <- mySub2 + xlab("Sucrose, 10~50% density gradient") + ylab("Absorbance, 260 nm") +
  scale_y_continuous(expand=c(0,0),limits=c(0,5.8),breaks=seq(-6,6,1)) + coord_cartesian(ylim = c(0,3))+
  scale_x_continuous(expand=c(0,0),limits=c(-1.5,78.2),breaks=theCoord$Position,labels = theCoord$Fraction.Number) 
mySub2 <- mySub2 + theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,1,0.5,1),'lines'),#axis.line.y = element_blank(),
                         legend.key.width=unit(0.3,'cm'),legend.key.height=unit(1,'cm'),legend.position='right',legend.text=element_text(size=12),
                         axis.text.y=element_text(size=10,face='plain',color='black'),axis.text.x=element_text(size=10,face='plain',color='black'),
                         axis.title.x=element_text(size=12,face='plain',color='black'),axis.title.y=element_text(size=12,hjust=0.5,vjust=2,face='plain',color='black'))
figure_4<-rbind(ggplotGrob(mySub2),size="last")
ggsave(file="./Fig7f_HDR_ploysome.pdf", plot=figure_4,bg = 'white', width = 12, height = 9, units = 'cm', dpi = 600)

