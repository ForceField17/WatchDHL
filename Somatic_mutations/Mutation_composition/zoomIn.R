# WatchDHL
library(rstudioapi)

# Getting the path of your current open file
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))
print( getwd() )

library(ggplot2)
library(gridExtra)
library(grid)



gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

#raw data preprocessing
table <- read.table('Coding.txt',header = T)

table$Fraction <- table$Num/table$Total

table2 <- table
revp_plot<-ggplot()+theme_classic()
revp_plot<-revp_plot+geom_bar(data=table2,aes(x=reorder(His,Fraction),y=Fraction,fill=Mut),width=0.65, color = "black",stat='identity',position=position_stack())+
     scale_y_continuous(expand=c(0,0),limits=c(0,0.03),breaks = seq(0,0.05,0.01),labels = c("0%","1%","2%","3%","4%","5%"))
#revp_plot<-revp_plot+geom_text(data = table2,aes(x=His,y=Pos,label=xxx))
revp_plot<-revp_plot+ theme(panel.background=element_rect(fill='transparent',color='transparent',size=1),plot.margin=unit(c(2,4,2,2),'lines'),
                            plot.title=element_text(size=16,vjust=0.5,hjust=0.5,face='bold.italic'),text=element_text(size=15,vjust=0.5,hjust=0.5,face='bold'),
                            legend.key.width=unit(0.5,'cm'),legend.key.height=unit(0.5,'cm'),legend.position="bottom",
                            legend.text=element_text(size=11,face='bold'),axis.text.y=element_text(size=12,angle=0,vjust=0.5,hjust=0.5,face='bold',color='black'),
                            axis.text.x=element_text(size=12,angle=25,vjust=1,hjust=1,face='bold',color='black'),axis.title.x=element_text(size=14,vjust=-4,hjust=0.5,face='plain',color='black'),
                            axis.title.y=element_text(size=14,hjust=0.5,vjust=4,face='plain',color='black'))
revp_plot<-revp_plot+ggtitle(NULL)+xlab(NULL)+ylab("Fraction of somatic mutations")+scale_x_discrete(position = "bottom")#+coord_flip()
#revp_plot<-revp_plot+annotation_logticks(base = 10, sides = "b", scaled = TRUE,short = unit(0.1, "cm"), mid = unit(0.2, "cm"), 
#                                     long = unit(0.3, "cm"), colour = "black", size = 0.5, linetype = 1, alpha = 1)
revp_plot<-revp_plot+scale_fill_manual(name=NULL,values=c(BBB='#b2df8a',CCC='#cab2d6',AAA='#fdbf6f'),labels=c(BBB="5' UTR",CCC='Coding',AAA="3' UTR"),guide=guide_legend(override.aes=list(size=0.5),nrow=2),na.translate = F)
revp_plot
figure_2<-rbind(ggplotGrob(revp_plot),size="first")
ggsave(file="Fig2a_UTRs.pdf", plot=figure_2,bg = 'white', width = 8, height = 13, units = 'cm', dpi = 600)


a <- data.frame(table$Num[1:3],table$Num[4:6])
chisq.test(a)