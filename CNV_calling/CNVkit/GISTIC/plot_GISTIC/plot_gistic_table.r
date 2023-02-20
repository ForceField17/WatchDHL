# WatchDHL
library(rstudioapi)
library("ggplot2")
library("gridExtra")
library(grid)
# Getting the path of your current open file
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))
print( getwd() )


gistic_results <- read.table("DHL_gistic_table.tsv", header=T, sep="")

amp<- gistic_results[gistic_results[,1]=="Amp",]
del<- gistic_results[gistic_results[,1]=="Del",]
chr<- gistic_results[gistic_results[,1]=="Chr",]
gene_amp<- gistic_results[gistic_results[,1]=="Gene_amp",]
gene_del<- gistic_results[gistic_results[,1]=="Gene_del",]


#chr$Chromosome

ffff_plot<-ggplot()+theme_classic()
ffff_plot<-ffff_plot+geom_line(data=amp ,aes(x=-Position_offset,y=X.log10.q.value.),size=1 ,color = "#b2182b")+
  scale_y_continuous(expression(paste("-log" ["10"], italic(" q-value"))),expand=c(0,0),limits = c(0, 8),breaks=seq(0,8,2))+
  scale_x_continuous(expand=c(0,0),breaks=-chr$X.log10.q.value.[c(1 ,  3  ,  5  , 7  ,  9 , 11 , 13 , 15 , 17 , 19 , 21 )],labels = c(1 ,  3  ,  5  , 7  ,  9 , 11 , 13 , 15 , 17 , 19 , 21 ))+
  geom_vline( xintercept =c(-chr$Position_offset),size=0.3,lty=5,color="grey")+ geom_hline(yintercept=2, color="#1a9850", lty=2)

#ffff_plot<-ffff_plot+geom_text(data=data_plot,aes(x=reorder(se,-orderID),y=Npat,label=Npat),color='black',stat='identity',size=3.6,hjust = 1, nudge_y = 0.2)
#ffff_plot<-ffff_plot+geom_text(data=data_plot,aes(x=reorder(se,-orderID),y=0.5,label=Pat),color='black',stat='identity',size=3,hjust = 0, nudge_y = 0)
ffff_plot<-ffff_plot+theme(panel.background=element_rect(fill='transparent',color='black',size=1.3),plot.margin=unit(c(1,1,1,0.1),'lines'),
                           plot.title=element_text(size=14,vjust=0.5,hjust=0.5,face='bold',color = "#b2182b"),text=element_text(size=15,vjust=0.5,hjust=0.5,face='bold'),
                           legend.key.width=unit(0.7,'cm'),legend.key.height=unit(0.7,'cm'),legend.position=c(0.7,0.2),axis.ticks.y = element_blank(),
                           legend.text=element_text(size=15,face='bold.italic'),axis.text.y=element_text(size=14,vjust=0.5,hjust=0.5,face='bold',color='black'),axis.text.x=element_text(size=14,face='plain',color='black'),
                           axis.title.x=element_text(size=14,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=20,hjust=0.5,vjust=4,face='bold',color='black'))
ffff_plot<-ffff_plot+ggtitle("Amplifications")+xlab(NULL)+coord_flip()#+scale_x_discrete(position = "bottom")

ffff_plot


dddd_plot<-ggplot()+theme_classic()
dddd_plot<-dddd_plot+geom_line(data=del ,aes(x=-Position_offset,y=X.log10.q.value.),size=1 ,color = "#2166ac")+
  scale_y_continuous(expression(paste("-log" ["10"], italic(" q-value"))),expand=c(0,0),limits = c(8, 0),breaks=seq(0,8,2),trans = 'reverse')+
  scale_x_continuous(expand=c(0,0),breaks=-chr$X.log10.q.value.[c( 2 , 4 ,  6 ,  8  ,10 , 12,14 , 16, 18 , 20 , 22)],labels = c( 2 , 4 ,  6 ,  8  ,10 , 12,14 , 16, 18 , 20 , 22),position = "top")+
  geom_vline( xintercept =c(-chr$Position_offset),size=0.3,lty=5,color="grey")+ geom_hline(yintercept=2, color="#1a9850", lty=2)

#ffff_plot<-ffff_plot+geom_text(data=data_plot,aes(x=reorder(se,-orderID),y=Npat,label=Npat),color='black',stat='identity',size=3.6,hjust = 1, nudge_y = 0.2)
#ffff_plot<-ffff_plot+geom_text(data=data_plot,aes(x=reorder(se,-orderID),y=0.5,label=Pat),color='black',stat='identity',size=3,hjust = 0, nudge_y = 0)
dddd_plot<-dddd_plot+theme(panel.background=element_rect(fill='transparent',color='black',size=1.3),plot.margin=unit(c(1,0.05,1,1),'lines'),
                           plot.title=element_text(size=14,vjust=0.5,hjust=0.5,face='bold',color ="#2166ac" ),text=element_text(size=15,vjust=0.5,hjust=0.5,face='bold'),
                           legend.key.width=unit(0.7,'cm'),legend.key.height=unit(0.7,'cm'),legend.position=c(0.7,0.2),axis.ticks.y = element_blank(),
                           legend.text=element_text(size=15,face='bold.italic'),axis.text.y=element_text(size=14,vjust=0.5,hjust=0.5,face='bold',color='black'),axis.text.x=element_text(size=14,face='plain',color='black'),
                           axis.title.x=element_text(size=14,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=20,hjust=0.5,vjust=4,face='bold',color='black'))
dddd_plot<-dddd_plot+ggtitle("Deletions")+xlab(NULL)+coord_flip()#+scale_x_discrete(position = "bottom")

dddd_plot


plot1<-cbind(ggplotGrob(dddd_plot),ggplotGrob(ffff_plot),size="last")
grid.draw(plot1)
ggsave(file="Standard_GisticPlot_DHL.pdf", plot=plot1,bg = 'white', width = 14, height = 20, units = 'cm', dpi = 300)


##########################
gistic_results <- read.table("FL_gistic_table.tsv", header=T, sep="")

amp<- gistic_results[gistic_results[,1]=="Amp",]
del<- gistic_results[gistic_results[,1]=="Del",]
chr<- gistic_results[gistic_results[,1]=="Chr",]
gene_amp<- gistic_results[gistic_results[,1]=="Gene_amp",]
gene_del<- gistic_results[gistic_results[,1]=="Gene_del",]


#chr$Chromosome

ffff_plot<-ggplot()+theme_classic()
ffff_plot<-ffff_plot+geom_line(data=amp ,aes(x=-Position_offset,y=X.log10.q.value.),size=1 ,color = "#b2182b")+
  scale_y_continuous(expression(paste("-log" ["10"], italic(" q-value"))),expand=c(0,0),limits = c(0, 8),breaks=seq(0,8,2))+
  scale_x_continuous(expand=c(0,0),breaks=-chr$X.log10.q.value.[c(1 ,  3  ,  5  , 7  ,  9 , 11 , 13 , 15 , 17 , 19 , 21 )],labels = c(1 ,  3  ,  5  , 7  ,  9 , 11 , 13 , 15 , 17 , 19 , 21 ))+
  geom_vline( xintercept =c(-chr$Position_offset),size=0.3,lty=5,color="grey")+ geom_hline(yintercept=2, color="#1a9850", lty=2)

#ffff_plot<-ffff_plot+geom_text(data=data_plot,aes(x=reorder(se,-orderID),y=Npat,label=Npat),color='black',stat='identity',size=3.6,hjust = 1, nudge_y = 0.2)
#ffff_plot<-ffff_plot+geom_text(data=data_plot,aes(x=reorder(se,-orderID),y=0.5,label=Pat),color='black',stat='identity',size=3,hjust = 0, nudge_y = 0)
ffff_plot<-ffff_plot+theme(panel.background=element_rect(fill='transparent',color='black',size=1.3),plot.margin=unit(c(1,1,1,0.1),'lines'),
                           plot.title=element_text(size=14,vjust=0.5,hjust=0.5,face='bold',color = "#b2182b"),text=element_text(size=15,vjust=0.5,hjust=0.5,face='bold'),
                           legend.key.width=unit(0.7,'cm'),legend.key.height=unit(0.7,'cm'),legend.position=c(0.7,0.2),axis.ticks.y = element_blank(),
                           legend.text=element_text(size=15,face='bold.italic'),axis.text.y=element_text(size=14,vjust=0.5,hjust=0.5,face='bold',color='black'),axis.text.x=element_text(size=14,face='plain',color='black'),
                           axis.title.x=element_text(size=14,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=20,hjust=0.5,vjust=4,face='bold',color='black'))
ffff_plot<-ffff_plot+ggtitle("Amplifications")+xlab(NULL)+coord_flip()#+scale_x_discrete(position = "bottom")

ffff_plot


dddd_plot<-ggplot()+theme_classic()
dddd_plot<-dddd_plot+geom_line(data=del ,aes(x=-Position_offset,y=X.log10.q.value.),size=1 ,color = "#2166ac")+
  scale_y_continuous(expression(paste("-log" ["10"], italic(" q-value"))),expand=c(0,0),limits = c(8, 0),breaks=seq(0,8,2),trans = 'reverse')+
  scale_x_continuous(expand=c(0,0),breaks=-chr$X.log10.q.value.[c( 2 , 4 ,  6 ,  8  ,10 , 12,14 , 16, 18 , 20 , 22)],labels = c( 2 , 4 ,  6 ,  8  ,10 , 12,14 , 16, 18 , 20 , 22),position = "top")+
  geom_vline( xintercept =c(-chr$Position_offset),size=0.3,lty=5,color="grey")+ geom_hline(yintercept=2, color="#1a9850", lty=2)

#ffff_plot<-ffff_plot+geom_text(data=data_plot,aes(x=reorder(se,-orderID),y=Npat,label=Npat),color='black',stat='identity',size=3.6,hjust = 1, nudge_y = 0.2)
#ffff_plot<-ffff_plot+geom_text(data=data_plot,aes(x=reorder(se,-orderID),y=0.5,label=Pat),color='black',stat='identity',size=3,hjust = 0, nudge_y = 0)
dddd_plot<-dddd_plot+theme(panel.background=element_rect(fill='transparent',color='black',size=1.3),plot.margin=unit(c(1,0.05,1,1),'lines'),
                           plot.title=element_text(size=14,vjust=0.5,hjust=0.5,face='bold',color ="#2166ac" ),text=element_text(size=15,vjust=0.5,hjust=0.5,face='bold'),
                           legend.key.width=unit(0.7,'cm'),legend.key.height=unit(0.7,'cm'),legend.position=c(0.7,0.2),axis.ticks.y = element_blank(),
                           legend.text=element_text(size=15,face='bold.italic'),axis.text.y=element_text(size=14,vjust=0.5,hjust=0.5,face='bold',color='black'),axis.text.x=element_text(size=14,face='plain',color='black'),
                           axis.title.x=element_text(size=14,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=20,hjust=0.5,vjust=4,face='bold',color='black'))
dddd_plot<-dddd_plot+ggtitle("Deletions")+xlab(NULL)+coord_flip()#+scale_x_discrete(position = "bottom")

dddd_plot


plot1<-cbind(ggplotGrob(dddd_plot),ggplotGrob(ffff_plot),size="last")
grid.draw(plot1)
  ggsave(file="Standard_GisticPlot_FL.pdf", plot=plot1,bg = 'white', width = 14, height = 20, units = 'cm', dpi = 300)

