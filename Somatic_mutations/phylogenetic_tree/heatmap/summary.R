# WatchDHL
library(rstudioapi)

# Getting the path of your current open file
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))
print( getwd() )

library(ggplot2)
library(gridExtra)
library(grid)
library("dplyr") 
library(reshape2)
library(oncoprint)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

#raw data preprocessing
Samples <- read.table('P8.mutation.heatmap.txt',header = T,sep="\t")
Samples$COLOR <- "G"
Samples$COLOR[which(Samples$FL.a==1 & Samples$FL.b==1 & Samples$DHL==1)] <- "A"
Samples$COLOR[which(Samples$FL.a==1 & Samples$FL.b==0 & Samples$DHL==0)] <- "B"
Samples$COLOR[which(Samples$FL.a==0 & Samples$FL.b==1 & Samples$DHL==1)] <- "C"
Samples$COLOR[which(Samples$FL.a==0 & Samples$FL.b==1 & Samples$DHL==0)] <- "D"
Samples$COLOR[which(Samples$FL.a==0 & Samples$FL.b==0 & Samples$DHL==1)] <- "E"

sample <- Samples[order(-Samples$FL.a),]
rownames(sample) <- c(1:nrow(sample))


ID<-paste0(sample$chromosome,sample$position,sample$ref,sample$alt)


mutGenes.table.somatic <- data.frame(sample$germline,sample$FL.a,sample$FL.b,sample$DHL)
rownames(mutGenes.table.somatic) <- ID
colnames(mutGenes.table.somatic) <- c("germline","FL.a","FL.b","DHL")



M<-t(mutGenes.table.somatic)
memoSort(M,sortGenes = FALSE)
#oncoPrint(M)

new_sort<-colnames(memoSort(M,sortGenes = FALSE))
sortedList <- data.frame(new_sort,c(1:length(new_sort)))
colnames(sortedList)<-c("ID","order")

xxx <- data.frame(ID,sample$germline,sample$FL.a,sample$FL.b,sample$DHL,sample$COLOR)
rownames(xxx) <- ID
colnames(xxx) <- c("ID","germline","FL.a","FL.b","DHL","COLOR")

sample <- merge(xxx,sortedList,1,1)

plot_1b.data <- melt(xxx[,c(1:5)])
colnames(plot_1b.data) <- c('ID','gene','value')

new_table <- merge(plot_1b.data,sample,1,1)
new_table$FILL <- new_table$COLOR
new_table$FILL[which(new_table$value==0)] <- "F"

F0b.plot<-ggplot()+theme_bw()
F0b.plot<-F0b.plot+geom_tile(data = new_table,aes(x=reorder(ID,order),y=reorder(gene,-as.integer(gene)),fill=as.factor(FILL)),color=NA,alpha=1,width=1.5,height=1,stat='identity')
F0b.plot<-F0b.plot+scale_fill_manual(name=NULL,values=c(A="#B8DDAE",B="#B2D867",C="#F5C342",D="#53B94C",E="#E77D72",F="transparent",G="#bc80bd"),guide = NULL)
F0b.plot<-F0b.plot+theme(panel.background=element_rect(fill='transparent',color='black',size=1),plot.margin=unit(c(1,1,1,1),'lines'),plot.title=element_text(size=14,vjust=0.5,hjust=0.5,face='plain'),
                         text=element_text(size=16,face='bold'),legend.key.width=unit(0.5,'cm'),legend.key.height=unit(0.5,'cm'),legend.position='none',axis.line.x = element_blank(),axis.ticks.x = element_blank(),axis.line.y = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=14,face='bold'),axis.text.y=element_text(size=14,vjust=0.5,hjust=1,face='plain',color='black'),legend.title=element_blank(),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F0b.plot<-F0b.plot+xlab(NULL)+ylab(NULL)+ggtitle(paste0("P8 ( ",nrow(Samples)," mutant positions)"))
F0b.plot<-F0b.plot+scale_y_discrete(expand=c(0,0))
F0b.plot<-F0b.plot+scale_x_discrete(expand=c(0,0))
F0b.plot

#ggsave(file="P8.pdf", plot=figure_1,bg = 'transparent', width = 18, height = 4, units = 'cm', dpi = 600)

length(which(Samples$FL.a==1 & Samples$FL.b==1 & Samples$DHL==1))
length(which(Samples$FL.a==1 & Samples$FL.b==0 & Samples$DHL==0))
length(which(Samples$FL.a==0 & Samples$FL.b==1 & Samples$DHL==1))
length(which(Samples$FL.a==0 & Samples$FL.b==1 & Samples$DHL==0))
length(which(Samples$FL.a==0 & Samples$FL.b==0 & Samples$DHL==1))
length(which(Samples$FL.a==1 & Samples$FL.b==0 & Samples$DHL==1))
length(which(Samples$FL.a==1 & Samples$FL.b==1 & Samples$DHL==0))
P8 <- F0b.plot


#############3
#raw data preprocessing
Samples <- read.table('P1.mutation.heatmap.txt',header = T,sep="\t")
Samples$COLOR <- "G"
Samples$COLOR[which(Samples$FL==1 & Samples$DHL==1)] <- "C"
Samples$COLOR[which(Samples$FL==1 & Samples$DHL==0)] <- "D"
Samples$COLOR[which(Samples$FL==0 & Samples$DHL==1)] <- "E"

sample <- Samples[order(-Samples$FL),]
rownames(sample) <- c(1:nrow(sample))

ID<-paste0(sample$chromosome,sample$position,sample$ref,sample$alt)

mutGenes.table.somatic <- data.frame(sample$germline,sample$FL,sample$DHL)
rownames(mutGenes.table.somatic) <- ID
colnames(mutGenes.table.somatic) <- c("germline","FL","DHL")


M<-t(mutGenes.table.somatic)

new_sort<-colnames(memoSort(M,sortGenes = FALSE))
sortedList <- data.frame(new_sort,c(1:length(new_sort)))
colnames(sortedList)<-c("ID","order")

xxx <- data.frame(ID,sample$germline,sample$FL,sample$DHL,sample$COLOR)
rownames(xxx) <- ID
colnames(xxx) <- c("ID","germline","FL","DHL","COLOR")

sample <- merge(xxx,sortedList,1,1)

plot_1b.data <- melt(xxx[,c(1:4)])
colnames(plot_1b.data) <- c('ID','gene','value')

new_table <- merge(plot_1b.data,sample,1,1)
new_table$FILL <- new_table$COLOR
new_table$FILL[which(new_table$value==0)] <- "F"

F0b.plot<-ggplot()+theme_bw()
F0b.plot<-F0b.plot+geom_tile(data = new_table,aes(x=reorder(ID,order),y=reorder(gene,-as.integer(gene)),fill=as.factor(FILL)),color=NA,alpha=1,width=1.5,height=1,stat='identity')
F0b.plot<-F0b.plot+scale_fill_manual(name=NULL,values=c(A="#B8DDAE",B="#B2D867",C="#F5C342",D="#53B94C",E="#E77D72",F="transparent",G="#bc80bd"),guide = NULL)
F0b.plot<-F0b.plot+theme(panel.background=element_rect(fill='transparent',color='black',size=1),plot.margin=unit(c(1,1,1,1),'lines'),plot.title=element_text(size=14,vjust=0.5,hjust=0.5,face='plain'),
                         text=element_text(size=16,face='bold'),legend.key.width=unit(0.5,'cm'),legend.key.height=unit(0.5,'cm'),legend.position='none',axis.line.x = element_blank(),axis.ticks.x = element_blank(),axis.line.y = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=14,face='bold'),axis.text.y=element_text(size=14,vjust=0.5,hjust=1,face='plain',color='black'),legend.title=element_blank(),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F0b.plot<-F0b.plot+xlab(NULL)+ylab(NULL)+ggtitle(paste0("P1 ( ",nrow(Samples)," mutant positions)"))
F0b.plot<-F0b.plot+scale_y_discrete(expand=c(0,0))
F0b.plot<-F0b.plot+scale_x_discrete(expand=c(0,0))
F0b.plot

length(which(Samples$FL==1 & Samples$DHL==1))
length(which(Samples$FL==1 & Samples$DHL==0))
length(which(Samples$FL==0 & Samples$DHL==1))

P1 <- F0b.plot


#############3
#raw data preprocessing
Samples <- read.table('P2.mutation.heatmap.txt',header = T,sep="\t")
Samples$COLOR <- "G"
Samples$COLOR[which(Samples$FL==1 & Samples$DHL==1)] <- "C"
Samples$COLOR[which(Samples$FL==1 & Samples$DHL==0)] <- "D"
Samples$COLOR[which(Samples$FL==0 & Samples$DHL==1)] <- "E"

sample <- Samples[order(-Samples$FL),]
rownames(sample) <- c(1:nrow(sample))

ID<-paste0(sample$chromosome,sample$position,sample$ref,sample$alt)

mutGenes.table.somatic <- data.frame(sample$germline,sample$FL,sample$DHL)
rownames(mutGenes.table.somatic) <- ID
colnames(mutGenes.table.somatic) <- c("germline","FL","DHL")


M<-t(mutGenes.table.somatic)

new_sort<-colnames(memoSort(M,sortGenes = FALSE))
sortedList <- data.frame(new_sort,c(1:length(new_sort)))
colnames(sortedList)<-c("ID","order")

xxx <- data.frame(ID,sample$germline,sample$FL,sample$DHL,sample$COLOR)
rownames(xxx) <- ID
colnames(xxx) <- c("ID","germline","FL","DHL","COLOR")

sample <- merge(xxx,sortedList,1,1)

plot_1b.data <- melt(xxx[,c(1:4)])
colnames(plot_1b.data) <- c('ID','gene','value')

new_table <- merge(plot_1b.data,sample,1,1)
new_table$FILL <- new_table$COLOR
new_table$FILL[which(new_table$value==0)] <- "F"

F0b.plot<-ggplot()+theme_bw()
F0b.plot<-F0b.plot+geom_tile(data = new_table,aes(x=reorder(ID,order),y=reorder(gene,-as.integer(gene)),fill=as.factor(FILL)),color=NA,alpha=1,width=1.5,height=1,stat='identity')
F0b.plot<-F0b.plot+scale_fill_manual(name=NULL,values=c(A="#B8DDAE",B="#B2D867",C="#F5C342",D="#53B94C",E="#E77D72",F="transparent",G="#bc80bd"),guide = NULL)
F0b.plot<-F0b.plot+theme(panel.background=element_rect(fill='transparent',color='black',size=1),plot.margin=unit(c(1,1,1,1),'lines'),plot.title=element_text(size=14,vjust=0.5,hjust=0.5,face='plain'),
                         text=element_text(size=16,face='bold'),legend.key.width=unit(0.5,'cm'),legend.key.height=unit(0.5,'cm'),legend.position='none',axis.line.x = element_blank(),axis.ticks.x = element_blank(),axis.line.y = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=14,face='bold'),axis.text.y=element_text(size=14,vjust=0.5,hjust=1,face='plain',color='black'),legend.title=element_blank(),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F0b.plot<-F0b.plot+xlab(NULL)+ylab(NULL)+ggtitle(paste0("P2 ( ",nrow(Samples)," mutant positions)"))
F0b.plot<-F0b.plot+scale_y_discrete(expand=c(0,0))
F0b.plot<-F0b.plot+scale_x_discrete(expand=c(0,0))
F0b.plot

length(which(Samples$FL==1 & Samples$DHL==1))
length(which(Samples$FL==1 & Samples$DHL==0))
length(which(Samples$FL==0 & Samples$DHL==1))

P2 <- F0b.plot

#############3
#raw data preprocessing
Samples <- read.table('P6.mutation.heatmap.txt',header = T,sep="\t")
Samples$COLOR <- "G"
Samples$COLOR[which(Samples$FL==1 & Samples$DHL==1)] <- "C"
Samples$COLOR[which(Samples$FL==1 & Samples$DHL==0)] <- "D"
Samples$COLOR[which(Samples$FL==0 & Samples$DHL==1)] <- "E"

sample <- Samples[order(-Samples$FL),]
rownames(sample) <- c(1:nrow(sample))

ID<-paste0(sample$chromosome,sample$position,sample$ref,sample$alt)

mutGenes.table.somatic <- data.frame(sample$germline,sample$FL,sample$DHL)
rownames(mutGenes.table.somatic) <- ID
colnames(mutGenes.table.somatic) <- c("germline","FL","DHL")


M<-t(mutGenes.table.somatic)

new_sort<-colnames(memoSort(M,sortGenes = FALSE))
sortedList <- data.frame(new_sort,c(1:length(new_sort)))
colnames(sortedList)<-c("ID","order")

xxx <- data.frame(ID,sample$germline,sample$FL,sample$DHL,sample$COLOR)
rownames(xxx) <- ID
colnames(xxx) <- c("ID","germline","FL","DHL","COLOR")

sample <- merge(xxx,sortedList,1,1)

plot_1b.data <- melt(xxx[,c(1:4)])
colnames(plot_1b.data) <- c('ID','gene','value')

new_table <- merge(plot_1b.data,sample,1,1)
new_table$FILL <- new_table$COLOR
new_table$FILL[which(new_table$value==0)] <- "F"

F0b.plot<-ggplot()+theme_bw()
F0b.plot<-F0b.plot+geom_tile(data = new_table,aes(x=reorder(ID,order),y=reorder(gene,-as.integer(gene)),fill=as.factor(FILL)),color=NA,alpha=1,width=1.5,height=1,stat='identity')
F0b.plot<-F0b.plot+scale_fill_manual(name=NULL,values=c(A="#B8DDAE",B="#B2D867",C="#F5C342",D="#53B94C",E="#E77D72",F="transparent",G="#bc80bd"),guide = NULL)
F0b.plot<-F0b.plot+theme(panel.background=element_rect(fill='transparent',color='black',size=1),plot.margin=unit(c(1,1,1,1),'lines'),plot.title=element_text(size=14,vjust=0.5,hjust=0.5,face='plain'),
                         text=element_text(size=16,face='bold'),legend.key.width=unit(0.5,'cm'),legend.key.height=unit(0.5,'cm'),legend.position='none',axis.line.x = element_blank(),axis.ticks.x = element_blank(),axis.line.y = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=14,face='bold'),axis.text.y=element_text(size=14,vjust=0.5,hjust=1,face='plain',color='black'),legend.title=element_blank(),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F0b.plot<-F0b.plot+xlab(NULL)+ylab(NULL)+ggtitle(paste0("P6 ( ",nrow(Samples)," mutant positions)"))
F0b.plot<-F0b.plot+scale_y_discrete(expand=c(0,0))
F0b.plot<-F0b.plot+scale_x_discrete(expand=c(0,0))
F0b.plot

length(which(Samples$FL==1 & Samples$DHL==1))
length(which(Samples$FL==1 & Samples$DHL==0))
length(which(Samples$FL==0 & Samples$DHL==1))

P6 <- F0b.plot

#############3
#raw data preprocessing
Samples <- read.table('P7.mutation.heatmap.txt',header = T,sep="\t")
Samples$COLOR <- "G"
Samples$COLOR[which(Samples$FL==1 & Samples$DHL==1)] <- "C"
Samples$COLOR[which(Samples$FL==1 & Samples$DHL==0)] <- "D"
Samples$COLOR[which(Samples$FL==0 & Samples$DHL==1)] <- "E"

sample <- Samples[order(-Samples$FL),]
rownames(sample) <- c(1:nrow(sample))

ID<-paste0(sample$chromosome,sample$position,sample$ref,sample$alt)

mutGenes.table.somatic <- data.frame(sample$germline,sample$FL,sample$DHL)
rownames(mutGenes.table.somatic) <- ID
colnames(mutGenes.table.somatic) <- c("germline","FL","DHL")


M<-t(mutGenes.table.somatic)

new_sort<-colnames(memoSort(M,sortGenes = FALSE))
sortedList <- data.frame(new_sort,c(1:length(new_sort)))
colnames(sortedList)<-c("ID","order")

xxx <- data.frame(ID,sample$germline,sample$FL,sample$DHL,sample$COLOR)
rownames(xxx) <- ID
colnames(xxx) <- c("ID","germline","FL","DHL","COLOR")

sample <- merge(xxx,sortedList,1,1)

plot_1b.data <- melt(xxx[,c(1:4)])
colnames(plot_1b.data) <- c('ID','gene','value')

new_table <- merge(plot_1b.data,sample,1,1)
new_table$FILL <- new_table$COLOR
new_table$FILL[which(new_table$value==0)] <- "F"

F0b.plot<-ggplot()+theme_bw()
F0b.plot<-F0b.plot+geom_tile(data = new_table,aes(x=reorder(ID,order),y=reorder(gene,-as.integer(gene)),fill=as.factor(FILL)),color=NA,alpha=1,width=1.5,height=1,stat='identity')
F0b.plot<-F0b.plot+scale_fill_manual(name=NULL,values=c(A="#B8DDAE",B="#B2D867",C="#F5C342",D="#53B94C",E="#E77D72",F="transparent",G="#bc80bd"),guide = NULL)
F0b.plot<-F0b.plot+theme(panel.background=element_rect(fill='transparent',color='black',size=1),plot.margin=unit(c(1,1,1,1),'lines'),plot.title=element_text(size=14,vjust=0.5,hjust=0.5,face='plain'),
                         text=element_text(size=16,face='bold'),legend.key.width=unit(0.5,'cm'),legend.key.height=unit(0.5,'cm'),legend.position='none',axis.line.x = element_blank(),axis.ticks.x = element_blank(),axis.line.y = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=14,face='bold'),axis.text.y=element_text(size=14,vjust=0.5,hjust=1,face='plain',color='black'),legend.title=element_blank(),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F0b.plot<-F0b.plot+xlab(NULL)+ylab(NULL)+ggtitle(paste0("P7 ( ",nrow(Samples)," mutant positions)"))
F0b.plot<-F0b.plot+scale_y_discrete(expand=c(0,0))
F0b.plot<-F0b.plot+scale_x_discrete(expand=c(0,0))
F0b.plot

length(which(Samples$FL==1 & Samples$DHL==1))
length(which(Samples$FL==1 & Samples$DHL==0))
length(which(Samples$FL==0 & Samples$DHL==1))

P7 <- F0b.plot

#############3
#raw data preprocessing
Samples <- read.table('P9.mutation.heatmap.txt',header = T,sep="\t")
Samples$COLOR <- "G"
Samples$COLOR[which(Samples$FL==1 & Samples$DHL==1)] <- "C"
Samples$COLOR[which(Samples$FL==1 & Samples$DHL==0)] <- "D"
Samples$COLOR[which(Samples$FL==0 & Samples$DHL==1)] <- "E"

sample <- Samples[order(-Samples$FL),]
rownames(sample) <- c(1:nrow(sample))

ID<-paste0(sample$chromosome,sample$position,sample$ref,sample$alt)

mutGenes.table.somatic <- data.frame(sample$germline,sample$FL,sample$DHL)
rownames(mutGenes.table.somatic) <- ID
colnames(mutGenes.table.somatic) <- c("germline","FL","DHL")


M<-t(mutGenes.table.somatic)

new_sort<-colnames(memoSort(M,sortGenes = FALSE))
sortedList <- data.frame(new_sort,c(1:length(new_sort)))
colnames(sortedList)<-c("ID","order")

xxx <- data.frame(ID,sample$germline,sample$FL,sample$DHL,sample$COLOR)
rownames(xxx) <- ID
colnames(xxx) <- c("ID","germline","FL","DHL","COLOR")

sample <- merge(xxx,sortedList,1,1)

plot_1b.data <- melt(xxx[,c(1:4)])
colnames(plot_1b.data) <- c('ID','gene','value')

new_table <- merge(plot_1b.data,sample,1,1)
new_table$FILL <- new_table$COLOR
new_table$FILL[which(new_table$value==0)] <- "F"

F0b.plot<-ggplot()+theme_bw()
F0b.plot<-F0b.plot+geom_tile(data = new_table,aes(x=reorder(ID,order),y=reorder(gene,-as.integer(gene)),fill=as.factor(FILL)),color=NA,alpha=1,width=1.5,height=1,stat='identity')
F0b.plot<-F0b.plot+scale_fill_manual(name=NULL,values=c(A="#B8DDAE",B="#B2D867",C="#F5C342",D="#53B94C",E="#E77D72",F="transparent",G="#bc80bd"),guide = NULL)
F0b.plot<-F0b.plot+theme(panel.background=element_rect(fill='transparent',color='black',size=1),plot.margin=unit(c(1,1,1,1),'lines'),plot.title=element_text(size=14,vjust=0.5,hjust=0.5,face='plain'),
                         text=element_text(size=16,face='bold'),legend.key.width=unit(0.5,'cm'),legend.key.height=unit(0.5,'cm'),legend.position='none',axis.ticks.x = element_blank(),
                         legend.direction='horizontal',legend.text=element_text(size=14,face='bold'),axis.text.y=element_text(size=14,vjust=0.5,hjust=1,face='plain',color='black'),legend.title=element_blank(),
                         axis.line = element_blank(),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='transparent'))
F0b.plot<-F0b.plot+xlab(NULL)+ylab(NULL)+ggtitle(paste0("P9 ( ",nrow(Samples)," mutant positions)"))
F0b.plot<-F0b.plot+scale_y_discrete(expand=c(0,0))
F0b.plot<-F0b.plot+scale_x_discrete(expand=c(0,0))
F0b.plot

length(which(Samples$FL==1 & Samples$DHL==1))
length(which(Samples$FL==1 & Samples$DHL==0))
length(which(Samples$FL==0 & Samples$DHL==1))

P9 <- F0b.plot


figure_1<-rbind(ggplotGrob(P1),ggplotGrob(P2),ggplotGrob(P6),ggplotGrob(P7),ggplotGrob(P8),ggplotGrob(P9), size="last")

panels <- figure_1$layout$t[grep("panel", figure_1$layout$name)]


figure_1$heights[panels][1] <- unit(3/5,'null')
figure_1$heights[panels][2] <- unit(3/5,'null')
figure_1$heights[panels][3] <- unit(3/5,'null')
figure_1$heights[panels][4] <- unit(3/5,'null')
figure_1$heights[panels][5] <- unit(4/5,'null')
figure_1$heights[panels][6] <- unit(3/5,'null')

#grid.draw(figure_1)
ggsave(file="ExtFig17_Mutation_number_heatmap_all.pdf", plot=figure_1,bg = 'white', width = 18, height = 20, units = 'cm', dpi = 600)


