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
  hcl(h = hues, l = 65, c = 100)[1:n]
}


table1 <- read.table("sample_table.txt",sep = "\t",header = T)
table1 <- table1[which(table1$patient!="P3"),]


xx1 <- c("P4","P5","P7","P9")
xx2 <- c("P1","P2","P3","P6")  # patients with IgH-MYC translocation
table2.1 <- table1[which(table1$patient %in% xx1),]
table2.2 <- table1[which(table1$patient %in% xx2),]
table2.3 <- table1[which(table1$patient == "P8" & table1$Grade == "a_FL"),]
table2.4 <- table1[which(table1$patient == "P8" & table1$TimeToTrans != 0),]
table2.4$TimeToTrans[which(table2.4$Grade=="a_FL")] <- 0.05 + table2.4$TimeToTrans[which(table2.4$Grade=="a_FL")]

histology.plot<-ggplot() +
  geom_point(data=table1,aes(x=TimeToTrans,y=reorder(patient,-Index),color=Grade),size=3.5,alpha=1,shape=1,stroke=1.7)+
  geom_path(data=table2.1,aes(x=TimeToTrans,y=reorder(patient,-Index)),alpha=0.85,size=0.7,
            arrow = arrow(length = unit(0.2, "cm")))+
  geom_path(data=table2.2,aes(x=TimeToTrans,y=reorder(patient,-Index)),alpha=0.85,size=0.7,color="blue",
            arrow = arrow(length = unit(0.2, "cm")))+
  geom_path( data=table2.3, aes(x=TimeToTrans,y=reorder(patient,-Index)),alpha=0.85,size=0.7,
            arrow = arrow(length = unit(0.2, "cm")))+
  geom_path( data=table2.4, aes(x=TimeToTrans,y=reorder(patient,-Index)),alpha=0.85,size=0.7,
             arrow = arrow(length = unit(0.2, "cm")))+
  xlab("Time after initial diagnosis (year)")+ylab(NULL) + 

  theme_classic()+
  theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(1,1,1,1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
        text=element_text(size=2,face='bold'),legend.key.width=unit(1,'cm'),legend.key.height=unit(0.8,'cm'),legend.position="none",legend.title=element_blank(),axis.ticks.y = element_blank(),
        legend.margin=margin(t=0.1,r=0.1,b=0,l=0.1,unit='cm'),legend.text=element_text(size=14,face='plain',color='black'),axis.text.y=element_text(size=14,face='bold',color='black'),axis.line.y = element_blank(),
        axis.text.x=element_text(size=16,face='bold',color='black'),axis.title.x=element_text(size=18,vjust=-1,face='plain',color='black'),axis.title.y=element_text(size=16,hjust=0.5,vjust=2,face='plain',color='black'))

histology.plot<-histology.plot + scale_color_manual(name="Grade",values=c(a_FL=gg_color_hue(6)[3],b_DHL=gg_color_hue(6)[1]),labels=c(a_FL="FL",b_DHL="DHL"),guide = guide_legend(override.aes=list(size=3.5),nrow=2),na.translate = F)
histology.plot<-histology.plot + scale_x_continuous(expand=c(0,0),limits = c(-0.3,14),breaks = seq(0,14,2))

figure_1<-rbind(ggplotGrob(histology.plot),size="first")
ggsave(file="Figure1A.pdf", plot=figure_1,bg = 'white', width = 16, height = 10, units = 'cm', dpi = 600)


