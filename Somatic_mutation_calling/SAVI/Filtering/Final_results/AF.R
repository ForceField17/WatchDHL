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


library("dplyr") 

table1 <- read.delim("./DHL.specific.bed",header = T,sep = "\t")
#xxx <- c("P3","P4","P5")
table1 <- table1[which(!(table1$CaseID %in% xxx)),]
revp_plot<-ggplot()+theme_classic()
revp_plot<-revp_plot+geom_density(data=table1 ,aes(x=T3_freq,color=CaseID),size=1,alpha=0.8)
#  geom_smooth(data = tmpsyn2, aes(x = Age, y = MB),method = "lm", se = TRUE,formula= y ~ x,color="black")
revp_plot<-revp_plot+ theme(panel.background=element_rect(fill='transparent',color='transparent',size=1),plot.margin=unit(c(1,1,2,1),'lines'),
                            plot.title=element_text(size=16,vjust=0.5,hjust=0.5,face='bold.italic'),text=element_text(size=15,vjust=0.5,hjust=0.5,face='bold'),
                            legend.key.width=unit(0.5,'cm'),legend.key.height=unit(0.5,'cm'),legend.position=c(0.86,0.78),
                            legend.text=element_text(size=12,face='plain'),axis.text.y=element_text(size=12,angle=0,vjust=0.5,hjust=0.5,face='bold',color='black'),
                            axis.text.x=element_text(size=12,angle=25,vjust=1,hjust=1,face='bold',color='black'),axis.title.x=element_text(size=16,vjust=-4,hjust=0.5,face='bold',color='black'),
                            axis.title.y=element_text(size=14,hjust=0.5,vjust=4,face='plain',color='black'))
revp_plot<-revp_plot+ggtitle("FL")+xlab(NULL)+ylab("Num")+scale_x_discrete(position = "bottom")
#revp_plot<-revp_plot+scale_y_continuous(expand=c(0,0),limits=c(0,0.1))
revp_plot<-revp_plot+scale_x_continuous(expand=c(0,0),limits=c(0,100),breaks=seq(0,100,20))#
revp_plot

figure_3<-rbind(ggplotGrob(revp_plot),size="first")
ggsave(file="DHL.pdf", plot=figure_3,bg = 'white', width = 36, height = 15, units = 'cm', dpi = 600)






