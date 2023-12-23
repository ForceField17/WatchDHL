# WatchDHL
library(rstudioapi)

# Getting the path of your current open file
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))
print( getwd() )

library(ggplot2)
library(DESeq2)
library(stringr)
library("pheatmap")
library("RColorBrewer")
library(fgsea)
library(ComplexHeatmap)
library(circlize)
library(rstatix)
library(ggpubr)
library(ggrepel)
library(preprocessCore)
library(clusterProfiler)
library(enrichplot)

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

#raw data preprocessing
Sample_features <- read.delim("Data/Info_Nascent.txt",sep = "\t",header = T)
samples <- data.frame(Sample_features)
rownames(samples) <- Sample_features$case
################
counts1 <- read.csv("Data/Final_Nascent_protein_matrix.txt", sep="", head=T, row.names = "Gene",quote = "\"")
counts1 <- counts1[,which((colnames(counts1 ) %in% samples$case ))]

nrow(counts)
counts <- round(counts1/10000)

xxx <- counts[,c(4:6,10:12,16:18)]
keep <- rowSums(xxx == min(counts)) <= 4
counts <- counts[keep,]
nrow(counts)
expMatrix_new <- log2(counts) 


#ANOVA
data <- expMatrix_new

AnovaP_genetic   <- rep(NA,nrow(data))
AnovaP_OPP       <- rep(NA,nrow(data))

Rank_genetic     <- rep(NA,nrow(data))
Rank_OPP         <- rep(NA,nrow(data))

nVar1 <- 0
nVar2 <- 0
for(i in 1:nrow(data)){
  test <- data.frame(as.numeric(data[i,]), samples$Protein,samples$genetic)
  colnames(test) <- c("value","protein","genetic")
  rownames(test) <- colnames(data)
  G1 <- test[which(test$genetic=="WT"),1]
  G2 <- test[which(test$genetic=="Alt"),1]
  
  G3 <- test[which(test$protein=="background"),1]
  G4 <- test[which(test$protein=="nascent"),1]
  
  res.aov2 <- aov(value ~  genetic + protein, data = test)
  myresults <- summary(res.aov2)
  
  AnovaP_genetic[i] <- myresults[[1]][[5]][1]
  AnovaP_OPP[i] <- myresults[[1]][[5]][2]
  
  if(mean(G2) > mean(G1)){
      Rank_genetic[i] <- -1   #Alt highly expressed ranking high
  }else{
      Rank_genetic[i] <- 1
  }
  
  if(mean(G4) > mean(G3)){
    Rank_OPP[i] <- -1   #nascent highly expressed ranking high
  }else{
    Rank_OPP[i] <- 1
  }
   
}

geneName <-rownames(expMatrix_new)

FDR_genetic <- p.adjust(AnovaP_genetic,method = "fdr")
FDR_OPP     <- p.adjust(AnovaP_OPP,method = "fdr")

Rank_genetic_GSEA <- Rank_genetic * log10(AnovaP_genetic)
Rank_genetic <- Rank_genetic * log10(FDR_genetic)
Rank_OPP_GSEA <- Rank_OPP * log10(AnovaP_OPP)
Rank_OPP <- Rank_OPP * log10(FDR_OPP)

final_data <- data.frame(geneName,AnovaP_genetic,FDR_genetic,Rank_genetic,AnovaP_OPP,FDR_OPP,Rank_OPP,Rank_OPP_GSEA,Rank_genetic_GSEA)
write.table(final_data , file = "./figs/Statistical_two_way_ANOVA.txt",sep = "\t",quote = F,row.names = F)


#################
TSG <- read.table("./lib/TS_list.txt",header=F)
ONCO <- read.table("./lib/OC_list.txt",header=F)
#################
plotTable <- final_data
plotTable_up <- plotTable[which( plotTable$Rank_genetic > 0  & plotTable$FDR_genetic <  0.01 ) ,]
plotTable_down <- plotTable[which( plotTable$Rank_genetic  <  0  & plotTable$FDR_genetic <  0.01 ) ,]
plotTable_other <- plotTable[which(!(plotTable$geneName %in% plotTable_up$geneName | plotTable$geneName %in% plotTable_down$geneName)) ,]

plotTable_Onco <- plotTable[(plotTable$geneName %in% ONCO$V1),]
plotTable_TSG <- plotTable[(plotTable$geneName %in% TSG$V1),]

plotTable_sig1 <- plotTable_up[(plotTable_up$geneName %in% plotTable_Onco$geneName ),]
plotTable_sig2 <- plotTable_down[(plotTable_down$geneName %in% plotTable_TSG$geneName),]
plotTable_sig <- rbind(plotTable_sig1,plotTable_sig2)

plotTable_sig3 <- plotTable_up[(plotTable_up$geneName %in% plotTable_TSG$geneName ),]
plotTable_sig4 <- plotTable_down[(plotTable_down$geneName %in% plotTable_Onco$geneName),]

write.table(plotTable_sig1,file = "./figs/ANOVA_upregulated_oncoproteins.txt",sep = "\t", quote=F, row.names = F)
write.table(plotTable_sig2,file = "./figs/ANOVA_downregulated_TSPs.txt",sep = "\t", quote=F, row.names = F)

mySub2 <- ggplot()+ theme_classic() 
mySub2 <- mySub2 + geom_hline(yintercept = 0,linetype=1,size=0.4,color="grey20")
mySub2 <- mySub2 + geom_point(data = plotTable_other,   aes(x =   (Rank_genetic), y = Rank_OPP),fill="grey",alpha=0.2,size=1.5, shape = 22,stroke=0.2) 
mySub2 <- mySub2  +
  geom_vline(xintercept = min(plotTable_up$Rank_genetic),linetype=2,size=0.4,color="grey50") +
  geom_vline(xintercept = max(plotTable_down$Rank_genetic),linetype=2,size=0.4,color="grey50")
mySub2 <- mySub2 + geom_point(data = plotTable_up,      aes(x =   (Rank_genetic), y = Rank_OPP),fill=gg_color_hue(3)[1],alpha=0.2,size=1.5, shape = 21,stroke=0.2) 
mySub2 <- mySub2 + geom_point(data = plotTable_down,    aes(x =   (Rank_genetic), y = Rank_OPP),fill=gg_color_hue(3)[3],alpha=0.2,size=1.5, shape = 21,stroke=0.2) 
mySub2 <- mySub2 + geom_point(data = plotTable_sig1,    aes(x =   (Rank_genetic), y = Rank_OPP),fill=gg_color_hue(3)[1],alpha=0.9,size=1.5, shape = 21,stroke=0.2) 
mySub2 <- mySub2 + geom_point(data = plotTable_sig2,    aes(x =   (Rank_genetic), y = Rank_OPP),fill=gg_color_hue(3)[3],alpha=0.9,size=1.5, shape = 21,stroke=0.2) 

mySub2 <- mySub2 +geom_text(data =plotTable_sig,   aes(x = Rank_genetic, y = Rank_OPP,label = geneName), size =0.7,color ='black')

mySub2 <- mySub2 + xlab("Log10(FDR of genetic factor)") + ylab("Log10(FDR of OPP label factor)") +
  scale_y_continuous(expand=c(0,0),limits=c(-17,17),breaks=seq(-20,20,5)) + 
  scale_x_continuous(expand=c(0,0),limits=c(-8.5,8.5),breaks=seq(-8,8,4))
mySub2 <- mySub2 + theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,1,0.5,1),'lines'),plot.title=element_text(size=14,vjust=0.5,hjust=0.5,face='bold.italic'),
                         legend.key.width=unit(1.5,'cm'),legend.key.height=unit(0.5,'cm'),legend.position='top',
                         legend.margin=margin(t=0.1,r=0.1,b=0,l=0.1,unit='cm'),legend.text=element_text(size=12),axis.text.y=element_text(size=10,face='plain',color='black'),
                         axis.text.x=element_text(size=10,face='plain',color='black'),axis.title.x=element_text(size=10,face='plain',color='black'),axis.title.y=element_text(size=10,hjust=0.5,vjust=2,face='plain',color='black'))
figure_4<-rbind(ggplotGrob(mySub2),size="last")
ggsave(file="./figs/Fig7e_TwoWay_ANOVA.pdf", plot=figure_4,bg = 'white', width = 10, height = 7, units = 'cm', dpi = 600)

a <- data.frame( c( nrow(plotTable_sig2), nrow(plotTable_TSG[which(!(plotTable_TSG$geneName %in% plotTable_down$geneName)),]) ),
                 c( nrow(plotTable_down[which(!(plotTable_down$geneName %in% plotTable_TSG$geneName)),]), nrow(plotTable[which(!(plotTable$geneName %in% plotTable_down$geneName | plotTable$geneName %in% plotTable_TSG$geneName)),])  ))
colnames(a) <- c("TSG","nonTSG")
rownames(a) <- c("Down","notDown")
a
fisher.test(a,alternative = "two.sided")

b <- data.frame( c( nrow(plotTable_sig1), nrow(plotTable_Onco[which(!(plotTable_Onco$geneName %in% plotTable_up$geneName)),]) ),
                 c( nrow(plotTable_up[which(!(plotTable_up$geneName %in% plotTable_Onco$geneName)),]), nrow(plotTable[which(!(plotTable$geneName %in% plotTable_up$geneName | plotTable$geneName %in% plotTable_Onco$geneName)),])  ))
colnames(b) <- c("ONCO","nonONCO")
rownames(b) <- c("Up","notUp")
b
fisher.test(b,alternative = "two.sided")

a <- data.frame( c( nrow(plotTable_sig3), nrow(plotTable_TSG[which(!(plotTable_TSG$geneName %in% plotTable_up$geneName)),]) ),
                 c( nrow(plotTable_up[which(!(plotTable_up$geneName %in% plotTable_TSG$geneName)),]), nrow(plotTable[which(!(plotTable$geneName %in% plotTable_up$geneName | plotTable$geneName %in% plotTable_TSG$geneName)),])  ))
colnames(a) <- c("TSG","nonTSG")
rownames(a) <- c("Down","notDown")
a
fisher.test(a,alternative = "two.sided")

b <- data.frame( c( nrow(plotTable_sig4), nrow(plotTable_Onco[which(!(plotTable_Onco$geneName %in% plotTable_down$geneName)),]) ),
                 c( nrow(plotTable_down[which(!(plotTable_down$geneName %in% plotTable_Onco$geneName)),]), nrow(plotTable[which(!(plotTable$geneName %in% plotTable_down$geneName | plotTable$geneName %in% plotTable_Onco$geneName)),])  ))
colnames(b) <- c("ONCO","nonONCO")
rownames(b) <- c("Up","notUp")
b
fisher.test(b,alternative = "two.sided")


######################
ranks <-  final_data$Rank_genetic_GSEA
names(ranks) <- final_data$geneName
geneList <- na.omit(ranks)
geneList = sort(geneList, decreasing = T)
str(ranks)
kegg <- read.gmt("./lib/c2.cp.kegg.v2022.1.Hs.symbols.gmt")
.Random.seed = 420
gse <- GSEA(geneList=geneList, TERM2GENE= kegg, pvalueCutoff = 0.1, minGSSize = 9,maxGSSize = 400,
            eps = 1e-50, pAdjustMethod = "BH", verbose = FALSE, seed = T, by = "fgsea")
gse_order<-gse[order(gse$enrichmentScore,decreasing=T)]
head(gse_order)
write.table(gse_order,file="figs/KEGG/results_GSEA_kegg.txt",sep="\t",quote=F,row.names = F)

fgseaResTidy <- gse_order[,c(1:8)]
fgseaResTidy$adj <- ifelse(fgseaResTidy$pvalue <= 0.05, "sig", "c_ns")
fgseaResTidy$adj[fgseaResTidy$NES>0 & fgseaResTidy$pvalue <= 0.05] <- "a_up"
fgseaResTidy$adj[fgseaResTidy$NES<0 & fgseaResTidy$pvalue <= 0.05] <- "b_down"
fgseaResTidy$ID <- gsub("_"," ",tolower(fgseaResTidy$ID))
fgseaResTidy$ID <- gsub("kegg ","",tolower(fgseaResTidy$ID))
TheTable <- fgseaResTidy[which(fgseaResTidy$pvalue <= 0.01 & fgseaResTidy$p.adjust <= 0.05),]
TheTable$dir <- ifelse(TheTable$NES >= 0, 1, -1)
TheTable$ADP <- TheTable$NES

plot2 <- ggplot() + theme_classic()
plot2 <- plot2 + geom_bar(data=TheTable, aes(x=reorder(ID,-ADP), y=ADP,fill = adj),alpha=0.7,stat='identity',width=0.5) + 
  geom_hline(yintercept = 0,linetype=1,color="black")+ scale_y_continuous(expand=c(0,0),limits=c(-2.5,2.5),breaks=seq(-4,4,1),position = "left") +
  coord_flip() +  labs(x="Pathway", y="NES", title="kegg enrichment") +scale_x_discrete(position = "top")+
  scale_fill_manual(NULL,values = c(c_ns = "grey", a_up = "#f4a582",b_down = gg_color_hue(6)[5]),labels=c(c_ns = "NS", a_up = "Up",b_down = "Down") ,guide = NULL) +
  theme(panel.background=element_rect(fill='transparent',color='transparent',size=1),plot.margin=unit(c(1,1,1,1),'lines'),plot.title=element_text(size=14,vjust=0.5,hjust=0.5,face='plain'),
        text=element_text(size=14,face='bold'),legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='right',legend.title=element_text(size=12,face='plain'),
        legend.margin=margin(t=0.1,r=0.1,b=0,l=0.1,unit='cm'),legend.text=element_text(size=12,face='plain'),axis.text.y=element_text(size=12,face='plain',color='black'),axis.line.y = element_blank(),
        axis.text.x=element_text(size=12,face='plain',color='black'),axis.title.x=element_text(size=14,face='plain',color='black'),axis.title.y=element_blank())
xxa <- plot2
figure_2<-rbind(ggplotGrob(plot2),size="last")
ggsave(file="./figs/ExtFig10a_KEGG_GSEA_summary.pdf", plot=figure_2,bg = 'white', width = 16, height = 12, units = 'cm', dpi = 600)

###########################################333
gse$ID
for(x in c(1:12)){
  part1 <- gseaplot2(gse, geneSetID = c(x) ,subplots=c(1),pvalue_table = F,color = "purple",title = sub("gobp","GO:BP",gsub("_"," ",tolower(gse$ID[x]))))
  part1 <- part1 + theme_classic()
  part1 <- part1 + scale_y_continuous(breaks=seq(-1,1,0.2)) 
  part1 <- part1 + ylab("Enrichment score")
  part1 <- part1 + geom_hline(yintercept = 0,linetype=2,color="grey")
  part1 <- part1   + theme(panel.background=element_rect(fill='white',color='black',size=1),plot.margin=unit(c(1,1,-0.11,1),'lines'),plot.title=element_text(size=14,vjust=0.5,hjust=0.5,face='plain'),
                           text=element_text(size=12,face='plain'),legend.key.width=unit(1,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='none',axis.line = element_blank(),axis.ticks.x = element_blank(),
                           legend.margin=margin(t=0.1,r=0.1,b=0,l=0.1,unit='cm'),legend.text=element_text(size=12,face='plain'),axis.text.y=element_text(size=12,face='plain',color='black'),
                           axis.text.x=element_blank(),axis.title.x=element_blank(),axis.title.y=element_text(size=12,hjust=0.5,vjust=2,face='plain',color='black'))
  part1
  
  part2 <- gseaplot2(gse, geneSetID = c(x) ,subplots=c(2))
  part2 <- part2 + theme_classic()
  part2 <- part2 + xlab("Protein ranking by genetic factor")
  part2 <- part2 +  coord_cartesian(ylim = c(0, 0.6))
  part2 <- part2   + theme(panel.background=element_rect(fill='transparent',color='black',size=1),plot.margin=unit(c(-0.14,1,1,3),'lines'),plot.title=element_text(size=16,vjust=0.5,hjust=0.5,face='bold'),
                           text=element_text(size=12,face='plain'),legend.key.width=unit(1,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='none',axis.line = element_blank(),axis.ticks.x = element_blank(),
                           legend.margin=margin(t=0.1,r=0.1,b=0,l=0.1,unit='cm'),legend.text=element_text(size=12,face='plain'),axis.text.y=element_blank(),axis.ticks.y = element_blank(),
                           axis.text.x=element_text(size=10,face='plain',color='black'),axis.title.x=element_text(size=12,face='plain',color='black'),axis.title.y=element_blank())
  part2
  
  figure <- rbind(ggplotGrob(part1),ggplotGrob(part2),size="first")
  panels <- figure$layout$t[grep("panel", figure$layout$name)]
  figure$heights[panels][1] <- unit(3,'null')
  figure$heights[panels][2] <- unit(0.5,'null')
  ggsave(file=paste0("./figs/KEGG/KEGGpathways_",x,".pdf"), plot=figure,bg = 'white', width = 9, height = 7, units = 'cm', dpi = 600)
}

