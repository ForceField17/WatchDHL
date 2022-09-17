setwd("/Users/songdong/Dropbox/Dropbox/DLBCL/main_figures/landscape_noncoding/Mut/")
library("ggplot2")
library("gridExtra")
library(grid)
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 120)[1:n]
}

ZZ=c('1_exp.','2_ESFF1')

gene.table<-read.table('./final/result_rank_mut_number_rmP3P4.txt',header=T,sep="\t")


data_plot  <- data.frame(c(as.character(gene.table$gene),as.character(gene.table$gene)),c(gene.table$rank,gene.table$rank),c(gene.table$FL,gene.table$DHL),c(gene.table$Total,gene.table$Total),c(rep("FL",nrow(gene.table)),rep("DHL",nrow(gene.table))))
colnames(data_plot) <- c("gene","ranks","Num","Total","Label")

orderID<-c(1:nrow(data_plot))

ffff_plot<-ggplot()+theme_classic()
ffff_plot<-ffff_plot+geom_bar(data=data_plot ,aes(x=reorder(gene,Total),y=Num,fill=Label),width=0.7 ,color = "black",stat='identity',position=position_dodge())
 ffff_plot<-ffff_plot+scale_y_continuous('Number of mutations',expand=c(0,0),limits = c(0, 400),breaks=seq(0,400,100))
ffff_plot<-ffff_plot+theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(1,1,1,1),'lines'),
                           plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),text=element_text(size=15,vjust=0.5,hjust=0.5,face='bold'),
                           legend.key.width=unit(0.7,'cm'),legend.key.height=unit(0.7,'cm'),legend.position="none",
                           legend.text=element_text(size=15,face='bold.italic'),axis.text.y=element_text(size=12,face='bold.italic',color='black'),axis.text.x=element_text(size=12,angle=45,vjust=1,hjust=1,face='bold',color='black'),
                           axis.title.x=element_text(size=20,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=20,hjust=0.5,vjust=4,face='bold',color='black'))
ffff_plot<-ffff_plot+ggtitle(NULL)+xlab(NULL)+scale_x_discrete(position = "bottom")+coord_flip()+ylab("Cases")
#ffff_plot<-ffff_plot+scale_y_continuous(trans="log10",expand=c(0,0.1))
#ffff_plot<-ffff_plot+annotation_logticks(base = 10, sides = "b", scaled = TRUE,short = unit(0.1, "cm"), mid = unit(0.2, "cm"), 
#                                     long = unit(0.3, "cm"), colour = "black", size = 0.5, linetype = 1, alpha = 1)
ffff_plot<-ffff_plot+scale_fill_manual(name=NULL,values=c(DHL=gg_color_hue(6)[1] ,FL=gg_color_hue(6)[3]))

ffff_plot

plot1<-cbind(ggplotGrob(ffff_plot),size="last")
	grid.draw(plot1)
ggsave(file="freqent_rank_rmP3P4.pdf", plot=plot1,bg = 'white', width = 18, height = 36, units = 'cm', dpi = 300)

