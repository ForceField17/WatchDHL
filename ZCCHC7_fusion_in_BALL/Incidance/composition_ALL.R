# WatchDHL
library(rstudioapi)

# Getting the path of your current open file
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))
print( getwd() )




gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

#raw data preprocessing
table <- read.table('Fraction.txt',header = T)

table$Fraction <- table$Num/table$Total
table2 <- table#[which(table$Mut=="BBB"),]
revp_plot<-ggplot()+theme_classic()
revp_plot<-revp_plot+geom_bar(data=table2,aes(x=reorder(His,Total),y=Fraction,fill=Mut),width=0.65, color = "black",stat='identity',position=position_stack())+scale_y_continuous(expand=c(0,0),limits=c(0,1.002),breaks = seq(0,1,0.2),labels = c("0%","20%","40%","60%","80%","100%"))
#revp_plot<-revp_plot+geom_text(data = table2,aes(x=His,y=Pos,label=xxx))
revp_plot<-revp_plot+ theme(panel.background=element_rect(fill='transparent',color='transparent',size=1),plot.margin=unit(c(1,2,2,1),'lines'),
                            plot.title=element_text(size=16,vjust=0.5,hjust=0.5,face='bold.italic'),text=element_text(size=15,vjust=0.5,hjust=0.5,face='bold'),
                            legend.key.width=unit(0.5,'cm'),legend.key.height=unit(0.5,'cm'),legend.position="top",
                            legend.text=element_text(size=12,face='plain'),axis.text.y=element_text(size=12,angle=0,vjust=0.5,hjust=0.5,face='bold',color='black'),
                            axis.text.x=element_text(size=12,face='plain',color='black'),axis.title.x=element_text(size=14,vjust=-4,hjust=0.5,face='plain',color='black'),
                            axis.title.y=element_text(size=12,hjust=0.5,vjust=4,face='plain',color='black'))
revp_plot<-revp_plot+ggtitle(NULL)+xlab(NULL)+ylab("Fraction of PAX5-ZCCHC7 fusion")+scale_x_discrete(position = "bottom")+coord_flip()
#revp_plot<-revp_plot+annotation_logticks(base = 10, sides = "b", scaled = TRUE,short = unit(0.1, "cm"), mid = unit(0.2, "cm"), 
#                                     long = unit(0.3, "cm"), colour = "black", size = 0.5, linetype = 1, alpha = 1)
revp_plot<-revp_plot+scale_fill_manual(name=NULL,values=c(BBB='#d9d9d9',AAA='#bc80bd'),labels=c(BBB='Wild type',AAA='PAX5-ZCCHC7 fusion'),guide=guide_legend(override.aes=list(size=0.5)),na.translate = F)
revp_plot
figure_2<-rbind(ggplotGrob(revp_plot),size="first")
ggsave(file="ExtFig18a_incidance.pdf", plot=figure_2,bg = 'white', width = 18, height = 8, units = 'cm', dpi = 600)


