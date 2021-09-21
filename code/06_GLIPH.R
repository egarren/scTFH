rm(list=ls())
library(Seurat)
library(cowplot)
library(ggplot2)
library(plyr)
library(dplyr)
library(biomaRt)
library(plotly)
library(scales)
library(EnhancedVolcano)
library(data.table)
library(ggpubr)
library(limma)
library(VennDiagram)
library(viridis)
library(pheatmap)
library(phylotools)
library(ggforce)
library(ggseqlogo)
library(immunarch)
library(cowplot)
library(ggrepel)
library(yingtools2)
library(qgraph)
library(GGally)
library(network)
library(sna)
library(tools)
library(ggalluvial)
library(tidyverse)
library(RColorBrewer)
library(gridExtra)
library(circlize)
library(useful)
library(patchwork)
library(vegan)
my.ttest <- function(...) {
  obj<-try(t.test(...), silent=TRUE)
  if (is(obj, "try-error")) return(NA) else return(obj$p.value)
}
data_concater <- function(x){
  x<- levels(factor(x))
  paste(x, collapse = "+")
}
data_concater2 <- function(x){
  x<- levels(factor(x))
  paste(x[1])
}
gm_mean = function(x, na.rm=TRUE){  #geometric mean functions for graphing
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
} 
gm_ymax<- function(x, na.rm = TRUE)
{
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))+exp(sd(log(x), na.rm = na.rm))
}

load("clone.ab.data.RData") 

###GLIPH export
dir.create("./GLIPH2")
setwd("./GLIPH2")
df<-clone.data
df.TRA<-df[df$chain=="TRA",]
df.TRB<-df[df$chain=="TRB",]
cells<-inner_join(df.TRB[,c("barcode","cdr3","v_gene","j_gene")],
                  df.TRA[,c("barcode","cdr3","v_gene","j_gene")],by="barcode")
cells<-left_join(cells,unique(clone.data[,c("barcode","mouse_ID","condition","my.clusters")]),by="barcode")
cells$cdr3_ID<-paste0(cells$cdr3.x,"+",cells$cdr3.y)
cdr3s<-as.data.table(cells)[, lapply(.SD, data_concater), by=cdr3_ID]
cells$clone_ID<-paste0(cells$cdr3.x,cells$v_gene.x,cells$j_gene.x,
                       cells$cdr3.y,cells$v_gene.y,cells$j_gene.y,cells$mouse_ID)
cells.cloneID<-transform(cells,freq=ave(seq(nrow(cells)),clone_ID,FUN=length))
cells.cloneID<-as.data.table(cells.cloneID)[, lapply(.SD, data_concater), by=clone_ID]
cells.cloneID$freq<-as.numeric(as.character(cells.cloneID$freq))
cells.GLIPH2<-data.frame(CDR3b=cells.cloneID$cdr3.x,TRBV=cells.cloneID$v_gene.x,TRBJ=cells.cloneID$j_gene.x,
                         CDR3a=cells.cloneID$cdr3.y,subject.condition=paste0(cells.cloneID$mouse_ID,".",cells.cloneID$condition),
                         count=cells.cloneID$freq)
write.table(cells.GLIPH2,file="cdr3.paired.GLIPH2.txt",sep="\t",row.names=F,col.names=F,quote=F)
save(list=c("cells","cdr3s","cells.GLIPH2","cells.cloneID"),file="meta.for.GLIPH.RData")
setwd("../")
#submit file to: http://50.255.35.37:8080/

#load GLIPH data
GLIPH.df<-read.csv("./P1749.csv",header=T) #input path to GLIPH results here
GLIPH.df$index<-paste0("GLIPH_",GLIPH.df$index)
GLIPH.groups<-unique(GLIPH.df[,1:12])
#add metadata
GLIPH.df$cdr3_ID<-paste0(GLIPH.df$TcRb,"+",GLIPH.df$TcRa)
GLIPH.df.filter<-GLIPH.df[GLIPH.df$number_subject>3&
                            GLIPH.df$vb_score<0.05&GLIPH.df$final_score<1e-5,]
GLIPH.filter.groups<-unique(GLIPH.df.filter[,1:12])
head(GLIPH.filter.groups[order(GLIPH.filter.groups$final_score),],n=9)

#Add Seurat to GLIPH data
#adding barcode to GLIPH df
df<-clone.data
cells<-inner_join(df[df$chain=="TRB",c("barcode","cdr3")],df[df$chain=="TRA",c("barcode","cdr3")],by="barcode")
cells$cdr3_ID<-paste0(cells$cdr3.x,"+",cells$cdr3.y)
cells.meta<-left_join_replace(cells,clone.data.ab,by="barcode")
GLIPH.barcode.ref<-inner_join(cells[,c("barcode","cdr3_ID")],unique(GLIPH.df.filter[,c(1:19,31)]),by="cdr3_ID") #17
#adding Seurat using barcode
GLIPH.data<-left_join_replace(x = GLIPH.barcode.ref,y = cells.meta, by = "barcode")
GLIPH.data<-left_join_replace(GLIPH.data,metadata,by="mouse_ID")
GLIPH.data$condition2<-"WT"
GLIPH.data$condition2[GLIPH.data$condition=="m564"]<-"564Igi"
for(i in 1:nrow(GLIPH.filter.groups)){
  index<-GLIPH.filter.groups$index[i]
  select.index<-GLIPH.data[GLIPH.data$index==index,]
  select.index.seurat<-select.index[!is.na(select.index$my.clusters),]
  GLIPH.filter.groups$total.size[i]<-nrow(select.index)
  GLIPH.filter.groups$tfr.score[i]<-100*nrow(select.index.seurat[select.index.seurat$my.clusters=="1",])/nrow(select.index.seurat)
  GLIPH.filter.groups$tcm.score[i]<-100*nrow(select.index.seurat[select.index.seurat$my.clusters=="3",])/nrow(select.index.seurat)
  GLIPH.filter.groups$top.cluster[i]<-names(sort(table(select.index.seurat$my.clusters2),decreasing=T))[1]
  GLIPH.filter.groups$gm.clone.size[i]<-gm_mean(table(select.index.seurat$cdr3_ID))
  GLIPH.filter.groups$clone.div.shannon[i]<-diversity(table(select.index.seurat$cdr3_ID),index="shannon")
  GLIPH.filter.groups$clone.div.simpson[i]<-diversity(table(select.index.seurat$cdr3_ID),index="simpson")
  GLIPH.filter.groups$clone.div.invsimp[i]<-diversity(table(select.index.seurat$cdr3_ID),index="invsimp")
}
GLIPH.df.filter<-left_join(GLIPH.df.filter,GLIPH.filter.groups[,c("index","total.size","tfr.score","tcm.score","top.cluster","gm.clone.size",
                                                                  "clone.div.shannon","clone.div.simpson","clone.div.invsimp")],by="index")
GLIPH.data<-left_join_replace(GLIPH.data,GLIPH.filter.groups[,c("index","total.size","tfr.score","tcm.score","top.cluster","gm.clone.size",
                                                                "clone.div.shannon","clone.div.simpson","clone.div.invsimp")],by="index")
GLIPH.data.seurat<-GLIPH.data[!is.na(GLIPH.data$my.clusters),]
GLIPH.barcode.ref<-left_join(GLIPH.barcode.ref,GLIPH.filter.groups[,c("index","total.size","tfr.score","tcm.score","top.cluster","gm.clone.size",
                                                                      "clone.div.shannon","clone.div.simpson","clone.div.invsimp")],by="index")
clusters<-names(table(GLIPH.data.seurat$my.clusters))
mice<-names(table(metadata$mouse_ID))
AID.mice<-names(table(metadata[metadata$condition=="AID",]$mouse_ID))
m564.mice<-names(table(metadata[metadata$condition=="m564",]$mouse_ID))

##Venn diagram
#between conditions
df<-GLIPH.data
dz.list<-list()
for(i in unique(df$condition)){
  dz.list[[i]]<-unique(df$index[df$condition==i])
}
venn.diagram(
  x = dz.list,category.names = c("WT","564Igi"),filename = "condition.GLIPH.venn.png",
  imagetype="png",height=700,width=700,margin=0.1,
  lwd=1,col=c("black", 'red'),fill=c(alpha("black",0.3), alpha('red',0.3)),
  # lty="blank",fill=sample(brewer.pal(9, "Set1"),size=length(dz.list)),
  cex=0.3,fontface="bold",fontfamily="sans",
  cat.cex=0.6,cat.fontface="bold",cat.default.pos="outer",cat.fontfamily="sans",
  cat.col = c("black", 'red'),cat.dist = c(0.1, 0.1),
  ext.pos=180,ext.line.lwd=0.25,ext.dist=-0.2
)
#between clusters
df<-GLIPH.data.seurat[GLIPH.data.seurat$my.clusters %in% c("0","1","2","3"),]
dz.list<-list()
for(i in c("Tfh-activated","Tfr","Sostdc1","Tfh-CM")){
  dz.list[[i]]<-unique(df$index[df$my.clusters2==i])
}
venn.colors<-hue_pal()(6)[1:4]
venn.diagram(
  x = dz.list,category.names = names(dz.list),filename = "Tfh.Tfr.GLIPH.venn.png",
  imagetype="png",height=1200,width=1200,margin=0.1,
  lwd=1,col=venn.colors,fill=alpha(venn.colors,0.3),
  # lty="blank",fill=sample(brewer.pal(9, "Set1"),size=length(dz.list)),
  cex=0.3,fontface="bold",fontfamily="sans",
  cat.cex=0.5,cat.fontface="bold",cat.default.pos="outer",cat.fontfamily="sans",
  cat.col = venn.colors#,cat.dist = c(0.1, 0.1,0.1,0.1)
)


#Rarefaction plot
for(i in c("0","1","2","3","4","5","all")){
  df.cluster<-GLIPH.data.seurat[GLIPH.data.seurat$my.clusters==i,] #gating on given cluster
  if(i=="all"){df.cluster<-GLIPH.data}
  raref.df<-as.data.frame.matrix(table(df.cluster$mouse_ID,df.cluster$index))
  colors<-rep("black",nrow(raref.df))
  colors[rownames(raref.df) %in% m564.mice]<-"red"
  test <- try(rarecurve(raref.df, step = 20, sample = min(rowSums(raref.df)), col = alpha(colors,0.6), cex = 0.2,
                        ylab="Clones",label=F,lwd=4))
  if(!is(test,"try-error")){
    png(paste0("clust.",i,".gliph.raref.png"),width=4,height=4,units="in",res=200)
    rarecurve(raref.df, step = 20, sample = min(rowSums(raref.df)), col = alpha(colors,0.6), cex = 0.2,
              ylab="Clones",label=F,lwd=4)
    dev.off()
  }
}


#NMDS
df<-GLIPH.data
community_matrix<-as.data.frame.matrix(table(df$mouse_ID,df$index))
example_NMDS=metaMDS(community_matrix, k=2,trymax=100)
png("gliph.NMDS.stress.png",width=4,height=4,units="in",res=200)
stressplot(example_NMDS)
dev.off()
plot(example_NMDS)
treat=rep("WT",nrow(community_matrix))
treat[rownames(community_matrix) %in% m564.mice]<-"564Igi"
colors=rep("black",nrow(community_matrix))
colors[rownames(community_matrix) %in% m564.mice]<-"red"
elevation=runif(10,0.5,1.5)
png("gliph.NMDS.png",width=6,height=6,units="in",res=200)
ordisurf(example_NMDS,elevation,main="",col="forestgreen")
ordiellipse(example_NMDS$point[grep("WT",treat),],draw="polygon",groups=treat[treat=="WT"],col="grey",alpha=50)
ordiellipse(example_NMDS$point[grep("564Igi",treat),],draw="polygon",groups=treat[treat=="564Igi"],col="red",alpha=50)
orditorp(example_NMDS,display="species",col="black",pch=21,bg=adjustcolor("grey",0.5),pcex=1,air=50000,cex=0.1)
orditorp(example_NMDS,display="sites",col=colors,air=0.0001,cex=1.25)
dev.off()



#public GLIPH repeteroire, by condition
df<-GLIPH.data
imm.list<-list()
for(i in names(table(df$mouse_ID))){ #creating immunarch list from collapsed data
  df2<-df[df$mouse_ID==i&!is.na(df$index),]
  df2<-transform(df2,Clones=ave(seq(nrow(df2)),index,FUN=length))
  df2<-unique(df2[,c("index","Clones")])
  df2$Proportion<-df2$Clones/sum(df2$Clones)
  imm.list[[i]]<-tibble(Clones=df2$Clones, Proportion=df2$Proportion, CDR3.aa=df2$index )
}
immdata = repLoad("./rep")
immdata$meta$Sample<-immdata$meta$mouse_ID
pr = pubRep(imm.list, "aa", .coding = T, .verbose = F)
pr.AID = pubRepFilter(pr, immdata$meta, c(condition = "AID"))
pr.564 = pubRepFilter(pr, immdata$meta, c(condition = "m564"))
pr.AID[is.na(pr.AID)]<-0
pr.564[is.na(pr.564)]<-0
pr.AID[["avgfreq.AID"]] = rowMeans(public_matrix(pr.AID), na.rm = T)
pr.564[["avgfreq.564"]] = rowMeans(public_matrix(pr.564), na.rm = T)
pr.AID[["sum.AID"]] = rowSums(public_matrix(pr.AID)[,1:5], na.rm = T)+1
pr.564[["sum.564"]] = rowSums(public_matrix(pr.564)[,1:5], na.rm = T)+1
pr.res.GLIPH = dplyr::full_join(pr.AID, pr.564, by = c("CDR3.aa")) #,"J.name"
pr.res.GLIPH[is.na(pr.res.GLIPH)]<-0
pr.res.GLIPH$sum.AID[pr.res.GLIPH$sum.AID==0]<-1
pr.res.GLIPH$sum.564[pr.res.GLIPH$sum.564==0]<-1
pr.res.GLIPH[["Samples.sum"]] = pr.res.GLIPH[["Samples.x"]] + pr.res.GLIPH[["Samples.y"]]
pr.res.GLIPH[["freq.ratio"]] = apply(pr.res.GLIPH[, c("avgfreq.AID", "avgfreq.564")],1, function(x) log10(x[1])/log10(x[2]))
pr.res.GLIPH[["log2FC"]]= apply(pr.res.GLIPH[, c("sum.AID", "sum.564")],1, function(x) log2((x[2])/(x[1])))
pr.res.GLIPH<-left_join(pr.res.GLIPH,GLIPH.filter.groups[,c("index","final_score","number_subject",
                                                            "number_unique_cdr3","total.size","tfr.score","tcm.score",
                                                            "top.cluster","gm.clone.size","clone.div.shannon",
                                                            "clone.div.simpson","clone.div.invsimp")],
                        by=c("CDR3.aa"="index"))
labels<-pr.res.GLIPH[abs(pr.res.GLIPH$log2FC)>4|pr.res.GLIPH$sum.564>100,]
ggplot(pr.res.GLIPH,aes(x = sum.564, y =  sum.AID,size=number_subject))+theme_classic()+
  scale_x_continuous(trans='log10') +scale_y_continuous(trans='log10')+ scale_size(range = c(4, 10))+
  geom_point(alpha=0.25,aes(colour=number_unique_cdr3))+ 
  guides(fill = guide_legend(override.aes = list(size = 7)))+
  labs(x = "564Igi", y = "WT", size="Samples",color="Clones")+
  theme(legend.direction = "vertical", legend.box = "horizontal")+
  geom_label_repel(data = labels,aes(label = CDR3.aa),size = 3,alpha=0.7,
                   color="royalblue",min.segment.length=unit(0,'lines'),nudge_x=0.1,nudge_y=0.1,segment.size=0.25)+
  scale_color_gradient(low="grey",high="darkgreen",trans="log",breaks=c(1,5,10,15,25),labels=format(c(1,5,10,15,25)))
ggsave2("overlap.AID.564.GLIPH.scatter.all.labeled.png",width=12, height=8,device="png")
labels<-pr.res.GLIPH[pr.res.GLIPH$CDR3.aa %in% c("GLIPH_330","GLIPH_4117","GLIPH_4614","GLIPH_140","GLIPH_4695","GLIPH_31",
                                                 "GLIPH_506","GLIPH_4448"),]
ggplot(pr.res.GLIPH,aes(x = sum.564, y =  sum.AID,size=number_subject))+theme_classic()+
  scale_x_continuous(trans='log10') +scale_y_continuous(trans='log10')+ scale_size(range = c(4, 10))+
  geom_point(alpha=0.25,aes(colour=number_unique_cdr3))+ 
  guides(fill = guide_legend(override.aes = list(size = 7)))+
  labs(x = "564Igi", y = "WT", size="Samples",color="Clones")+
  theme(legend.direction = "vertical", legend.box = "horizontal")+
  geom_label_repel(data = labels,aes(label = CDR3.aa),size = 3,alpha=0.7,
                   color="royalblue",min.segment.length=unit(0,'lines'),nudge_x=0.1,nudge_y=0.1,segment.size=0.25)+
  scale_color_gradient(low="grey",high="darkgreen",trans="log",breaks=c(1,5,10,15,25),labels=format(c(1,5,10,15,25)))
ggsave2("overlap.AID.564.GLIPH.scatter.png",width=6, height=4,device="png")
for(i in c("tfr.score","tcm.score","gm.clone.size","clone.div.shannon","clone.div.simpson","clone.div.invsimp")){
  ggplot(pr.res.GLIPH,aes(x = sum.564, y =  sum.AID,size=number_subject))+theme_classic()+
    scale_x_continuous(trans='log10') +scale_y_continuous(trans='log10')+ scale_size(range = c(4, 10))+
    geom_point(alpha=0.25,aes(colour=.data[[i]]))+ #tfr.score, most common cluster, avg clone size, final_score
    guides(fill = guide_legend(override.aes = list(size = 7)),size=F)+#,color=F)+ #,size=F
    labs(x = "564Igi", y = "WT", size="Samples",color=i)+scale_color_viridis()+ #,color="Clones"
    theme(legend.direction = "vertical", legend.box = "horizontal",legend.title=element_blank())#+
  ggsave2(paste0("overlap.AID.564.GLIPH.",i,".scatter.png"),width=5, height=4,device="png")
}
show_col(hue_pal()(6))
colors<-hue_pal()(6)[1:5]
pr.res.GLIPH$top.cluster <- factor(pr.res.GLIPH$top.cluster, levels = c("Tfh-activated","Tfr","Sostdc1","Tfh-CM","Tfh-effector"))
ggplot(pr.res.GLIPH,aes(x = sum.564, y =  sum.AID,size=number_subject))+theme_classic()+
  scale_x_continuous(trans='log10') +scale_y_continuous(trans='log10')+ scale_size(range = c(4, 10))+
  geom_point(aes(colour=factor(top.cluster),fill = factor(top.cluster)), shape=21) + 
  scale_color_manual(values = colors)+
  scale_fill_manual(values = alpha(colors,0.5))+
  guides(fill = guide_legend(override.aes = list(size = 7)),size=F)+#,color=F)+ #,size=F
  labs(x = "564Igi", y = "WT")+ #,color="Clones"
  theme(legend.direction = "vertical", legend.box = "horizontal",legend.title=element_blank())#+
ggsave2(paste0("overlap.AID.564.GLIPH.clusters.scatter.png"),width=5.5, height=4,device="png")


### WT vs 564 specific expanded GLIPHs (log2FC, 2,3)
AID.public.expand.GLIPH<-pr.res.GLIPH[pr.res.GLIPH$log2FC< -3 &pr.res.GLIPH$sum.AID>10,]
m564.public.expand.GLIPH<-pr.res.GLIPH[pr.res.GLIPH$log2FC>3&pr.res.GLIPH$sum.564>10,]
pr.res.GLIPH$annotate<-ifelse(pr.res.GLIPH$CDR3.aa %in% AID.public.expand.GLIPH$CDR3.aa, "WT-expanded","Public")
pr.res.GLIPH$annotate[pr.res.GLIPH$CDR3.aa %in% m564.public.expand.GLIPH$CDR3.aa]<-"564Igi-expanded"
pr.res.GLIPH$annotate <- factor(pr.res.GLIPH$annotate, levels = c("Public","WT-expanded","564Igi-expanded"))
ggplot(pr.res.GLIPH,aes(x = sum.564, y =  sum.AID,size=Samples.sum))+theme_classic()+
  scale_x_continuous(trans='log10') +scale_y_continuous(trans='log10')+ scale_size(range = c(4, 10))+
  geom_point(aes(colour=factor(annotate),fill = factor(annotate)), shape=21) + 
  scale_color_manual(values = c(alpha("black",0.5),"black","red"))+
  scale_fill_manual(values = c("#1C00ff00",alpha("black",0.2),alpha("red",0.2)))+
  guides(fill = guide_legend(override.aes = list(size = 7)),color=F,size=F)+ #,size=F
  labs(x = "564Igi", y = "WT", size="Samples",fill="Expansion")+
  theme(legend.direction = "vertical", legend.box = "horizontal")#+#theme(legend.title = element_blank())+
ggsave2("GLIPH.expansion.scatter.png",width=6, height=4,device="png")
pr.res.GLIPH$annotate <- factor(pr.res.GLIPH$annotate, levels = c("Public","WT-expanded","564Igi-expanded"))
plot.list<-list()
for(i in c("tfr.score","gm.clone.size","clone.div.shannon")){ 
  plot.list[[i]]<-ggbarplot(pr.res.GLIPH, x = "annotate", y = i, color="annotate", add = c("mean_se","jitter"),add.params = list(size = 0.5),
            position = position_dodge(0.8), legend="right",palette=c("gray50","black","red"),xlab=F)+
    theme(legend.title = element_blank(),axis.text.x=element_blank())+NoLegend()#+
}
plot.list[[1]]+plot.list[[2]]+plot.list[[3]]+plot_layout(guides="collect")&theme(legend.position="right")
ggsave2("bar.dot.cluster.GLIPH.metrics.png",width=7.5, height=2.5,device="png")



#studying specific GLIPH
sig.indexes<-head(GLIPH.filter.groups[order(GLIPH.filter.groups$final_score),],n=9)$index
sig.indexes<-c("GLIPH_506","GLIPH_4117","GLIPH_140","GLIPH_2359","GLIPH_1908","GLIPH_3500","GLIPH_4448","GLIPH_4695")
for(i in sig.indexes){
  pattern.df<-GLIPH.data[GLIPH.data$index==i,]
  df<-unique(pattern.df[,c("index","TcRa","TcRb","Freq","mouse_ID","condition2")])
  df<- df[order(df$Freq,decreasing=T),]
  write.csv(df,file=paste0(i,".csv"))
  #comparing clone size in specific GLIPH
  pattern.df<-GLIPH.data[GLIPH.data$index==i,]
  df.list<-list()
  for(j in c("AID","m564")){
    cond.gliph<-pattern.df[pattern.df$condition==j,]
    freq.tab<-as.data.frame(table(cond.gliph$cdr3_ID))
    if(j=="AID"){freq.tab$condition<-"WT"}else{freq.tab$condition<-"564Igi"}
    df.list[[j]]<-freq.tab #add to list for each condition
  }
  df<-rbindlist(df.list)
  df$condition <- factor(df$condition, levels = c("WT","564Igi"))
  ggbarplot(df, x = "condition", y = "Freq", color="condition", add = c("mean","jitter"),
            position = position_dodge(0.8), legend="none",palette=c("black","red"),xlab=F)+
    labs(x="", y = "Clone Size") #+scale_y_log10()
  ggsave2(paste0("bar.dot.condition.cdr3.freq.",i,".png"),width=2, height=3,device="png")
  mu <- ddply(df, "condition", summarise, grp.mean=gm_mean(Freq))
  df$condition <- factor(df$condition, levels = c("WT","564Igi"))
  ggplot(df, aes(x=Freq,color=condition,fill=condition)) + 
    geom_histogram(aes(y=..density..),bins=10,position="identity",alpha=0.2) +geom_density(fill=NA,adjust=5,size=1)+
    geom_vline(data=mu, aes(xintercept=grp.mean, color=condition),linetype="dashed")+theme_classic()+
    labs(x = "Clone Size", y = "Density")+scale_x_log10()+theme(legend.title = element_blank(),legend.position="bottom")+
    scale_color_manual(values = c("black", "red"))+
    scale_fill_manual(values = c("black", "red"))
  ggsave2(paste0("histogram.condition.cdr3.freq.",i,".png"),width=5, height=3,device="png")
  
}


GLIPH.data$GLIPH.exp<-ifelse(GLIPH.data$index %in% AID.public.expand.GLIPH$CDR3.aa, "WT-expanded","Public")
GLIPH.data$GLIPH.exp[GLIPH.data$index %in% m564.public.expand.GLIPH$CDR3.aa]<-"564Igi-expanded"
#studying expanded GLIPH
for(i in c("WT-expanded","564Igi-expanded","Public")){
  pattern.df<-GLIPH.data[GLIPH.data$GLIPH.exp==i,]
  df<-unique(pattern.df[,c("condition2","mouse_ID","TcRa","TcRb","cdr3.freq")])
  df<- df[order(df$cdr3.freq,decreasing=T),]
  write.csv(df,file=paste0("GLIPH.exp.",i,".csv"))
}

#comparing clone size in expanded GLIPH
df.list<-list()
for(j in c("WT-expanded","564Igi-expanded","Public")){
  cond.gliph<-GLIPH.data[GLIPH.data$GLIPH.exp==j,]
  freq.tab<-as.data.frame(table(cond.gliph$cdr3_ID))
  freq.tab$condition<-j
  df.list[[j]]<-freq.tab
}
df<-rbindlist(df.list)
df$condition <- factor(df$condition, levels = c("Public","WT-expanded","564Igi-expanded"))
ggbarplot(df, x = "condition", y = "Freq", color="condition", add = c("mean","jitter"),
          position = position_dodge(0.8), legend="none",palette=c("gray","black","red"),xlab=F)+
  labs(x="", y = "Clone Size")
ggsave2("bar.dot.condition.cdr3.freq.gliph.exp.png",width=3, height=3,device="png")
mu <- ddply(df, "condition", summarise, grp.mean=gm_mean(Freq))
ggplot(df, aes(x=Freq,color=condition,fill=condition)) + 
  geom_histogram(aes(y=..density..),bins=10,position="identity",alpha=0.2) +geom_density(fill=NA,adjust=5,size=1)+
  geom_vline(data=mu, aes(xintercept=grp.mean, color=condition),linetype="dashed")+theme_classic()+
  labs(x = "Clone Size", y = "Density")+scale_x_log10()+theme(legend.title = element_blank(),legend.position="bottom")+
  scale_color_manual(values = c("gray30","black", "red"))+
  scale_fill_manual(values = c("gray30","black", "red"))
ggsave2("histogram.condition.cdr3.freq.gliph.exp.png",width=5, height=3,device="png")


### GLIPH-tcr network
diff.GLIPH<-c("GLIPH_1442","GLIPH_397","GLIPH_606","GLIPH_206","GLIPH_47","GLIPH_115","GLIPH_10","GLIPH_4","GLIPH_97","GLIPH_566")
df<-GLIPH.data[GLIPH.data$index %in% c(diff.GLIPH) &GLIPH.data$cdr3.freq>1,]
# extract link between clonotype_id.long and GLIPH for network
tmp3=data.frame(df$index,df$cdr3,stringsAsFactors = F)
tmp3=tmp3[which(duplicated(paste(tmp3[,1],tmp3[,2],sep="_"))==F),]
colnames(tmp3)=c("x","y")
# create table for qgraph
toQgraph=rbind(tmp3)
# color of Dx
l=unique(c(toQgraph[,1],toQgraph[,2]))
col.tmp=rep(adjustcolor("blue",.4),length(l)) 
col.tmp[which(l %in% unique(df$cdr3[which(df$condition=="AID")]))]<-adjustcolor("black",.2)
col.tmp[which(l %in% unique(df$cdr3[which(df$condition=="m564")]))]<-adjustcolor("red",.2)
col.tmp[which(l %in% unique(intersect(df$cdr3[which(df$condition=="AID")],
                                      df$cdr3[which(df$condition=="m564")])))]<-adjustcolor("forestgreen",.9)
# size of nodes
size.tmp=rep(1,length(l))
size.tmp[which(l %in% unique(df$index))]<-10
for(i in which(size.tmp==1)){
  size.tmp[i]=log(gm_mean(df$cdr3.freq[which(df$cdr3==l[i])]))*1.5
}
# color of edges
line.tmp=rep("grey",nrow(toQgraph))
#labels
labels.cex.tmp=rep(0.0001,length(l))
labels.cex.tmp[which(l %in% c(diff.GLIPH))]<-1 #,exp.GLIPH
#graph
dim(toQgraph)
png("diff.GLIPH.network.by.cdr3s.png",width=10,height=10,units="in",res=400)
qgraph(toQgraph,color=col.tmp,vsize=size.tmp,edge.color=line.tmp,edge.width=0.15,
       directed=F,repulsion=1,borders=F,labels=TRUE,label.cex=labels.cex.tmp)
dev.off()


#public GLIPH repeteroire, by Tfh vs Tfr
df<-GLIPH.data.seurat
imm.list<-list()
for(i in names(table(df$mouse_ID))){ #creating immunarch list from collapsed data
  for(j in c("0","1")){
    df2<-df[df$mouse_ID==i&!is.na(df$index)&df$my.clusters==j,]
    df2<-transform(df2,Clones=ave(seq(nrow(df2)),index,FUN=length))
    df2<-unique(df2[,c("index","Clones")])
    df2$Proportion<-df2$Clones/sum(df2$Clones)
    imm.list[[paste0(i,j)]]<-tibble(Clones=df2$Clones, Proportion=df2$Proportion, CDR3.aa=df2$index)
  }
}
immdata = repLoad("./clust.rep")
immdata$meta$Sample<-paste0(immdata$meta$mouse_ID,immdata$meta$Cluster)
pr = pubRep(imm.list, "aa", .coding = T, .verbose = F)
pr.Tfh = pubRepFilter(pr, immdata$meta, c(Cluster=0))
pr.Tfr = pubRepFilter(pr, immdata$meta, c(Cluster=1))
pr.Tfh[is.na(pr.Tfh)]<-0
pr.Tfr[is.na(pr.Tfr)]<-0
pr.Tfh[["avgfreq.Tfh"]] = rowMeans(public_matrix(pr.Tfh), na.rm = T)
pr.Tfr[["avgfreq.Tfr"]] = rowMeans(public_matrix(pr.Tfr), na.rm = T)
pr.Tfh[["sum.Tfh"]] = rowSums(public_matrix(pr.Tfh)[,1:10], na.rm = T)+1
pr.Tfr[["sum.Tfr"]] = rowSums(public_matrix(pr.Tfr)[,1:10], na.rm = T)+1
pr.res.GLIPH = dplyr::full_join(pr.Tfh, pr.Tfr, by = c("CDR3.aa")) #,"J.name"
pr.res.GLIPH[is.na(pr.res.GLIPH)]<-0
pr.res.GLIPH$sum.Tfh[pr.res.GLIPH$sum.Tfh==0]<-1
pr.res.GLIPH$sum.Tfr[pr.res.GLIPH$sum.Tfr==0]<-1
pr.res.GLIPH[["Samples.sum"]] = pr.res.GLIPH[["Samples.x"]] + pr.res.GLIPH[["Samples.y"]]
pr.res.GLIPH[["freq.ratio"]] = apply(pr.res.GLIPH[, c("avgfreq.Tfh", "avgfreq.Tfr")],1, function(x) log10(x[1])/log10(x[2]))
pr.res.GLIPH[["log2FC"]]= apply(pr.res.GLIPH[, c("sum.Tfh", "sum.Tfr")],1, function(x) log2((x[2])/(x[1])))
pr.res.GLIPH<-left_join(pr.res.GLIPH,GLIPH.filter.groups[,c("index","final_score","number_subject",
                                                            "number_unique_cdr3","total.size","tfr.score")],
                        by=c("CDR3.aa"="index"))
GLIPH.public.FC<-pr.res.GLIPH[,c("CDR3.aa","Samples.sum","log2FC","final_score","number_subject",
                                 "number_unique_cdr3","total.size","tfr.score")] 
labels<-pr.res.GLIPH[abs(pr.res.GLIPH$log2FC)>3|pr.res.GLIPH$sum.Tfr>100,]
ggplot(pr.res.GLIPH,aes(x = sum.Tfr, y =  sum.Tfh,size=number_subject))+theme_classic()+
  scale_x_continuous(trans='log10') +scale_y_continuous(trans='log10')+ scale_size(range = c(4, 10))+
  geom_point(alpha=0.25,aes(colour=number_unique_cdr3))+ 
  guides(fill = guide_legend(override.aes = list(size = 7)))+
  labs(x = "Tfr", y = "Tfh", size="Samples",color="Clones")+
  theme(legend.direction = "vertical", legend.box = "horizontal")+
  geom_label_repel(data = labels,aes(label = CDR3.aa),size = 3,
                   color="royalblue",min.segment.length=unit(0,'lines'),nudge_x=0.1,nudge_y=0.1,segment.size=0.25)+
  scale_color_gradient(low="grey",high="darkgreen",trans="log",breaks=c(1,5,10,15,25),labels=format(c(1,5,10,15,25)))
ggsave2("overlap.Tfh.Tfr.GLIPH.scatter.png",width=6, height=4,device="png")


### GLIPH-tcr network (Tfh. Tfr)
diff.GLIPH<-unique(pr.res.GLIPH[abs(pr.res.GLIPH$log2FC)>3 | pr.res.GLIPH$sum.Tfr>100,]$CDR3.aa) 
df<-GLIPH.data.seurat[GLIPH.data.seurat$index %in% c(diff.GLIPH) &GLIPH.data.seurat$cdr3.freq>1,]
df<-GLIPH.data.seurat[GLIPH.data.seurat$cdr3.freq>1,]
# extract link between clonotype_id.long and GLIPH for network
tmp3=data.frame(df$index,df$cdr3,stringsAsFactors = F)
tmp3=tmp3[which(duplicated(paste(tmp3[,1],tmp3[,2],sep="_"))==F),]
colnames(tmp3)=c("x","y")
# create table for qgraph
toQgraph=rbind(tmp3)
# color of Dx
l=unique(c(toQgraph[,1],toQgraph[,2]))
col.tmp=rep(adjustcolor("black",.4),length(l)) #adjustcolor("purple",.3)
col.tmp[which(l %in% unique(df$index))]<-adjustcolor("blue",.4)
for(i in 6:1){
  col.tmp[which(l %in% unique(df$cdr3[which(df$my.clusters==(i-1))]))]<-adjustcolor(hue_pal()(6)[i],.4)
}
# size of nodes
size.tmp=rep(1,length(l))
size.tmp[which(l %in% unique(df$index))]<-10
for(i in which(size.tmp==1)){
  size.tmp[i]=log(gm_mean(df$cdr3.freq[which(df$cdr3==l[i])]))
}
# color of edges
line.tmp=rep("grey",nrow(toQgraph))
#labels
labels.cex.tmp=rep(0.0001,length(l))
labels.cex.tmp[which(l %in% c(diff.GLIPH))]<-1 
#graph
dim(toQgraph)
png("diff.GLIPH.network.by.cdr3s.Tfh.Tfr.png",width=10,height=10,units="in",res=400)
qgraph(toQgraph,color=col.tmp,vsize=size.tmp,edge.color=line.tmp,edge.width=0.15,
       directed=F,repulsion=1,borders=F,labels=TRUE,label.cex=labels.cex.tmp)
dev.off()


load("graphed.RData")

#Add GLIPH to Seurat data
GLIPH.barcode.ref.ordered<-GLIPH.barcode.ref[order(GLIPH.barcode.ref$final_score),]
meta<-GLIPH.barcode.ref.ordered[!duplicated(GLIPH.barcode.ref.ordered$barcode),]
rownames(meta)<-meta$barcode
T_all.combined@meta.data$orig.barcode<-rownames(T_all.combined@meta.data)
T_all.combined@meta.data$condition2<-"WT"
T_all.combined@meta.data$condition2[T_all.combined@meta.data$condition=="m564"]<-"564Igi"
T_all.combined@meta.data$condition2 <- factor(x = T_all.combined@meta.data$condition2, levels = c("WT", "564Igi"))
T_all.combined<-AddMetaData(T_all.combined,cbind(meta[,1:14],dplyr::select(meta, total.size)))
Idents(T_all.combined) <- "condition"
m564.combined<-subset(T_all.combined, idents='m564')
AID.combined<-subset(T_all.combined, idents='AID')
save.image("temp.GLIPH.seurat.RData")


#mapping GLIPH group size 
notNA<-T_all.combined@meta.data$barcode[!is.na(T_all.combined@meta.data$index)]
Idents(T_all.combined)<-"barcode"
Tall.GLIPH.groups<-subset(x = T_all.combined, idents=notNA) 
p<-FeaturePlot(Tall.GLIPH.groups, features= "total.size",split.by = "condition2",pt.size=0.1, 
               order=T,cols=c("grey","royalblue"),combine=F) #NoLegend()+
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] +NoAxes()+#NoLegend()+
    theme(panel.border = element_rect(colour = "black", size=1),
          plot.title=element_blank(),axis.title.y.right=element_blank(),
          axis.line=element_blank())
}
cowplot::plot_grid(p[[1]],p[[2]],ncol=2)
ggsave2("GLIPH.group.size.umap.bycondition.png",width=9, height=3.25,device="png")

##Mapping individual GLIPH groups within condition
pat.list<-c("GLIPH_506","GLIPH_4117","GLIPH_4448")
Idents(AID.combined) <- "index"
Idents(m564.combined) <- "index"
plot.list<-list()
for(i in pat.list){
  plot.list[[paste0(i,"AID")]]<-DimPlot(
    object = AID.combined, cols.highlight="royalblue",
    cells.highlight = Cells(subset(AID.combined, idents=i)),sizes.highlight=0.3) + 
    NoLegend() + NoAxes()+#ggtitle("WT")+
    theme(plot.title = element_text(size = 10,hjust = 0.5),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5))
  plot.list[[paste0(i,"564")]]<-DimPlot(
    object = m564.combined, cols.highlight="royalblue",
    cells.highlight = Cells(subset(m564.combined, idents=i)),sizes.highlight=0.3) + 
    NoLegend() + NoAxes()+#ggtitle("564Igi")+
    theme(plot.title = element_text(size = 10,hjust = 0.5),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5))
}
CombinePlots(plots = plot.list, ncol = 2)
ggsave2("GLIPH.oi.umap.png",width=4.5, height=6,device="png")

#Weblogo for GLIPH groups
plot.list<-list()
for(j in pat.list){
  df<-GLIPH.data[GLIPH.data$index ==j,]
  freq.tab<-as.data.frame(table(df$TcRa))  #ms_cdr3 frequency table for given cluster
  freq.tab$aa_length<-nchar(as.character(freq.tab$Var1))
  len<-freq.tab$aa_length[1]
  df2<-freq.tab[freq.tab$aa_length==len,]
  plot.list[[paste0("TRA",j)]]<-ggplot()+geom_logo(as.character(df2$Var1),method="probability")+
    theme_logo()+ggtitle(paste0("TRA_",j))+
    theme(plot.title = element_text(hjust = 0.5)) 
  freq.tab<-as.data.frame(table(df$TcRb))  #ms_cdr3 frequency table for given cluster
  freq.tab$aa_length<-nchar(as.character(freq.tab$Var1))
  len<-freq.tab$aa_length[1]
  df2<-freq.tab[freq.tab$aa_length==len,]
  plot.list[[paste0("TRB",j)]]<-ggplot()+geom_logo(as.character(df2$Var1),method="probability")+
    theme_logo()+ggtitle(paste0("TRB_",j))+
    theme(plot.title = element_text(hjust = 0.5))
}
CombinePlots(plot.list,ncol=2,legend="bottom")
ggsave2("cdr3.weblogo.topGLIPHs.png",width=6,height=6,device="png")

#######GLIPH DE
# # DE b/w condition-specific expanded clones
#DE between AID vs 564 specific expanded clones
AID.pats<-unique(AID.public.expand.GLIPH$CDR3.aa)
m564.pats<-unique(m564.public.expand.GLIPH$CDR3.aa)
T_all.combined@meta.data$GLIPH.exp<-ifelse(T_all.combined@meta.data$index %in% AID.pats, "WT-expanded","Public")
T_all.combined@meta.data$GLIPH.exp[T_all.combined@meta.data$index %in% m564.pats]<-"564Igi-expanded"
T_all.combined@meta.data$GLIPH.exp<- factor(T_all.combined@meta.data$GLIPH.exp, levels = c("Public","WT-expanded","564Igi-expanded"))
DimPlot(T_all.combined, reduction = "umap", group.by="GLIPH.exp",split.by = "condition2",
        cols=c("grey","black","red"),order=T)+ #NoLegend()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text=element_blank(),axis.ticks=element_blank(),axis.line=element_blank(),axis.title=element_blank())
ggsave2("GLIPH.expansion.umap.png",width=7, height=2.75,device="png")
Idents(T_all.combined)<-"GLIPH.exp"
GLIPH.expanded.autoimmune.response <- try(FindMarkers(T_all.combined, ident.1 = "564Igi-expanded", ident.2 = "WT-expanded",
                                                      min.pct=0,logfc.threshold = -Inf,test.use = "MAST"))
if(!is(GLIPH.expanded.autoimmune.response,"try-error")){
  GLIPH.expanded.autoimmune.response$log2FC<-log2(exp(GLIPH.expanded.autoimmune.response$avg_logFC))
  avg.GLIPH.expanded <- log1p(AverageExpression(T_all.combined)$RNA)
  avg.GLIPH.expanded$gene <- rownames(T_all.combined)
}

# DE b/w conditions in individual clones
pat<-pr.res.GLIPH.seurat.condition[pr.res.GLIPH.seurat.condition$sum.AID>50,]$CDR3.aa
DefaultAssay(T_all.combined) <- "RNA"
Idents(T_all.combined) <- "cellID"
for(i in pat){
  barcodes<-T_all.combined@meta.data$cellID[T_all.combined@meta.data$index %in% i]
  df<-subset(T_all.combined, idents=barcodes)
  Idents(df) <- "condition" #setting idents to condition metadata
  df.autoimmune.response <- try(FindMarkers(df, ident.1 = "m564", ident.2 = "AID",
                                        min.pct=0,logfc.threshold = -Inf,test.use = "MAST"))
  if(!is(df.autoimmune.response,"try-error")){
    df.autoimmune.response$log2FC<-log2(exp(df.autoimmune.response$avg_logFC))
    avg.df <- log1p(AverageExpression(df)$RNA)
    avg.df$gene <- rownames(df)
    assign(paste0(i,".GLIPH.autoimmune.response"),df.autoimmune.response)
    assign(paste0(i,".GLIPH.avg.exp"),avg.df)
  }
}

# 
##GLIPH expansion anlaysis
# identifying expanded clones
m564.expanded <- names(table(m564.combined@meta.data$index))[table(m564.combined@meta.data$index) > 50]
m564.unexpanded <- names(table(m564.combined@meta.data$index))[table(m564.combined@meta.data$index) >= 1&
                                                                   table(m564.combined@meta.data$index) <=5]
AID.expanded <- names(table(AID.combined@meta.data$index))[table(AID.combined@meta.data$index) > 50]
AID.unexpanded <- names(table(AID.combined@meta.data$index))[table(AID.combined@meta.data$index) >= 1&
                                                                  table(AID.combined@meta.data$index) <=5]
#expanded DE (by condition)
Idents(m564.combined) <- "index"
expandedDE.m564<-FindMarkers(m564.combined, ident.1 = m564.expanded, ident.2 = m564.unexpanded,
                             min.pct=0,logfc.threshold = -Inf,test.use = "MAST")
expandedDE.m564$gene<-rownames(expandedDE.m564)
expandedDE.m564$log2FC<-log2(exp(expandedDE.m564$avg_logFC))
Idents(AID.combined) <- "index"
expandedDE.AID<-FindMarkers(AID.combined, ident.1 = AID.expanded, ident.2 = AID.unexpanded,
                            min.pct=0,logfc.threshold = -Inf,test.use = "MAST")
expandedDE.AID$gene<-rownames(expandedDE.AID)
expandedDE.AID$log2FC<-log2(exp(expandedDE.AID$avg_logFC))
exp.vs.expandedDE.GLIPH <-merge(expandedDE.m564[,c("gene","log2FC")], expandedDE.AID[,c("gene","log2FC")], by="gene")
rownames(exp.vs.expandedDE.GLIPH)<-exp.vs.expandedDE.GLIPH$gene


###Expansion within GLIPH group
m<-"GLIPH_506"
for(m in pat){
  DefaultAssay(T_all.combined) <- "RNA"
  Idents(T_all.combined) <- "index"
  seurat.temp<-subset(T_all.combined, idents=m)
  Idents(seurat.temp) <- "condition"
  m564.temp<-subset(seurat.temp, idents='m564')
  AID.temp<-subset(seurat.temp, idents='AID')
  cor.df<-list()
  for(k in c("AID","m564")){
    obj<-get(paste0(k,".temp"))
    cond.freqs<-as.data.frame(table(obj@meta.data$cdr3_ID))
    if(nrow(cond.freqs)>5){
      rownames(cond.freqs)<-cond.freqs$Var1
      colnames(cond.freqs)<-c("cdr3_ID","cdr3_ID.freq")
      obj@meta.data<-left_join(obj@meta.data,cond.freqs,by="cdr3_ID")
      df<-as.data.frame(t(as.matrix(obj[["RNA"]]@data)))
      mod_score<-as.numeric(obj@meta.data$cdr3_ID.freq)
      mat<-cbind(mod_score,df)
      # ggscatter(mat,x="Cd74",y="Foxp3",add="reg.line",conf.int=T,cor.coef=T,cor.method="spearman")
      # correlations<-apply(matrix,1,function(x){cor(mod_score,x)})
      # cors<-apply(mat,2,cor.test,mod_score,method="spearman") #non-parametric (pearson for parametric)
      cor.scores<-data.frame(gene=rep(NA,ncol(mat)),pval=rep(NA,ncol(mat)),rho=rep(NA,ncol(mat)))
      for(j in 1:ncol(mat)){
        cor<-cor.test(mat[[j]],mat$mod_score,method="spearman")
        cor.scores$pval[j]<-cor$p.value
        cor.scores$rho[j]<-cor$estimate
        cor.scores$gene[j]<-colnames(mat)[j]
      }
      cor.scores[is.na(cor.scores)]<-0
      assign(paste0(k,".",m,".size.correl"),cor.scores[cor.scores$gene!="mod_score",])
      colnames(cor.scores)<-c("gene",paste0(k,".pval"),paste0(k,".rho"))
      cor.df[[k]]<-cor.scores
    }
  }
  if(length(cor.df)==2){
    cor.compare<-merge(cor.df[[1]], cor.df[[2]], by="gene")
    assign(paste0(m,".size.cor.compare"),cor.compare[cor.compare$gene!="mod_score",])
  }
 
  # identifying expanded clones
  m564.exp.in.GLIPH <- names(table(m564.temp@meta.data$cdr3_ID))[table(m564.temp@meta.data$cdr3_ID) > 2]
  m564.unexp.in.GLIPH <- names(table(m564.temp@meta.data$cdr3_ID))[table(m564.temp@meta.data$cdr3_ID) ==1]
  AID.exp.in.GLIPH <- names(table(AID.temp@meta.data$cdr3_ID))[table(AID.temp@meta.data$cdr3_ID) > 2]
  AID.unexp.in.GLIPH <- names(table(AID.temp@meta.data$cdr3_ID))[table(AID.temp@meta.data$cdr3_ID) ==1]
  #exp.in.GLIPH DE (by condition)
  if(length(m564.unexp.in.GLIPH)>2 & length(AID.unexp.in.GLIPH)>2 & length(m564.exp.in.GLIPH)>2 & length(AID.exp.in.GLIPH)>2){
    Idents(m564.temp) <- "cdr3_ID"
    exp.in.GLIPHDE.m564<-FindMarkers(m564.temp, ident.1 = m564.exp.in.GLIPH, ident.2 = m564.unexp.in.GLIPH, 
                                     min.pct=0,logfc.threshold = -Inf,test.use = "MAST")
    exp.in.GLIPHDE.m564$gene<-rownames(exp.in.GLIPHDE.m564)
    exp.in.GLIPHDE.m564$log2FC<-log2(exp(exp.in.GLIPHDE.m564$avg_logFC))
    Idents(AID.temp) <- "cdr3_ID"
    exp.in.GLIPHDE.AID<-FindMarkers(AID.temp, ident.1 = AID.exp.in.GLIPH, ident.2 = AID.unexp.in.GLIPH, 
                                    min.pct=0,logfc.threshold = -Inf,test.use = "MAST")
    exp.in.GLIPHDE.AID$gene<-rownames(exp.in.GLIPHDE.AID)
    exp.in.GLIPHDE.AID$log2FC<-log2(exp(exp.in.GLIPHDE.AID$avg_logFC))
    exp.vs.exp.in.GLIPHDE <-merge(exp.in.GLIPHDE.m564[,c("gene","log2FC")], exp.in.GLIPHDE.AID[,c("gene","log2FC")], by="gene")
    rownames(exp.vs.exp.in.GLIPHDE)<-exp.vs.exp.in.GLIPHDE$gene
    assign(paste0(m,".exp.vs.exp.in.GLIPHDE"),exp.vs.exp.in.GLIPHDE)
  }
}


save.image("GLIPH.analyzed.RData")
save(T_all.combined,file="T_all.RData")

clusters<-names(table(GLIPH.data.seurat$my.clusters))
save(list=c(ls(pattern="DE"),ls(pattern=".markers"),ls(pattern=".response"),"clusters"),file="DE.dfs.RData")


#Volcano Plots
for(i in c(ls(pattern=".GLIPH.autoimmune.response"),"GLIPH.expanded.autoimmune.response")){
  df<-get(i)
  # rownames(df)<-df$gene
  if(!is(df,"try-error")){
    if(!("p_val_adj" %in% colnames(df))){
      if("ttest" %in% colnames(df) ){df$p_val_adj<-p.adjust(df$ttest,method="BH")}
      if("p_val" %in% colnames(df) ){df$p_val_adj<-p.adjust(df$p_val,method="BH")}
    }
    if(!("log2FC" %in% colnames(df))){
      if("delta" %in% colnames(df)){df$log2FC<-df$delta}
      if("avg_logFC" %in% colnames(df)){df$log2FC<-log2(exp(df$avg_logFC))}
    }
    if(is.numeric(df$log2FC)&is.numeric(df$p_val_adj)){
      gene.list<-rownames(df[abs(df$log2FC)>0.2&df$p_val_adj<0.01,])
      if(length(gene.list)>15){gene.list<-rownames(df[abs(df$log2FC)>0.2&df$p_val_adj<0.005,])}
      if(length(gene.list)>15){gene.list<-rownames(df[abs(df$log2FC)>0.21&df$p_val_adj<0.0035,])}
      if(length(gene.list)>15){gene.list<-rownames(df[abs(df$log2FC)>0.21&df$p_val_adj<0.001,])}
      if(length(gene.list)>15){gene.list<-rownames(df[abs(df$log2FC)>0.22&df$p_val_adj<10e-5,])}
      if(length(gene.list)>15){gene.list<-rownames(df[abs(df$log2FC)>0.35&df$p_val_adj<10e-6,])}
      if(length(gene.list)>15){gene.list<-rownames(df[abs(df$log2FC)>0.35&df$p_val_adj<10e-8,])}
      if(length(gene.list)>15){gene.list<-rownames(df[abs(df$log2FC)>0.35&df$p_val_adj<10e-10,])}
      if(length(gene.list)>15){gene.list<-rownames(df[abs(df$log2FC)>0.35&df$p_val_adj<10e-15,])}
      if(length(gene.list)>15){gene.list<-rownames(df[abs(df$log2FC)>0.5&df$p_val_adj<10e-15,])}
      EnhancedVolcano(df,lab = rownames(df),
                      x = 'log2FC',y = 'p_val_adj',title = i,col=c("black","black","black","red3"),
                      selectLab=gene.list,xlab=bquote(~Log[2]~ (frac("564Igi expanded","WT expanded"))),
                      pCutoff = 0.01,FCcutoff = 0.2,pointSize = 0.5,labSize = 3,axisLabSize=10,colAlpha = 1, #transparencyxlim = c(-1.5, 1.5),
                      legendVisible=F,drawConnectors = TRUE,widthConnectors = 0.2,colConnectors = 'grey30',
                      subtitle="", caption="",border="full",cutoffLineWidth=0,
                      gridlines.major=F,gridlines.minor=F,titleLabSize=10
      )
      ggsave2(paste0(i,".volcano.png"),width=3.5, height=4,device="png")
    }
  }
}

###Antigen Prediction

## export masterdb to GLIPH
load("masterdb3.RData")
df<-unique(masterdb[!is.na(masterdb$chain)&masterdb$chain=="TRB"&masterdb$TRBV!=""&
                      masterdb$Species=="Mouse"&masterdb$Tcell=="CD4",]) #
df<-unique(masterdb[!is.na(masterdb$chain)&masterdb$chain=="TRB"&masterdb$TRBV!=""&
                      masterdb$Tcell=="CD4",]) #masterdb$Species=="Mouse"&
df<-unique(masterdb)
df[df==""]<-NA
masterdb.GLIPH<-unique(data.frame(CDR3b=df$cdr3.b,TRBV=df$TRBV,TRBJ=df$TRBJ,CDR3a="NA",
                                  subject.condition=paste0(df$study2,":",df$disease2),
                                  count=df$freq))
write.table(masterdb.GLIPH,file="./GLIPH2/ref.GLIPH.txt",sep="\t",row.names=F,col.names=F,quote=F)
#submit file to: http://50.255.35.37:8080/

#load GLIPH ref data
GLIPH.ref.df<-read.csv("./GLIPH2/P829.csv",header=T) #input path to GLIPH results here
GLIPH.ref.df$index<-paste0("GLIPH.ref_",GLIPH.ref.df$index)
GLIPH.ref.groups<-unique(GLIPH.ref.df[,1:12])
GLIPH.ref.df.filter<-GLIPH.ref.df
GLIPH.ref.filter.groups<-unique(GLIPH.ref.df.filter[,1:12])
head(GLIPH.ref.filter.groups[order(GLIPH.ref.filter.groups$final_score),],n=9)

##add back masterdb metadata
load("masterdb3.RData")
df<-unique(masterdb[!is.na(masterdb$chain)&masterdb$chain=="TRB"&masterdb$TRBV!=""& #P654
                      masterdb$Species=="Mouse"&masterdb$Tcell=="CD4",]) #
df<-unique(masterdb[!is.na(masterdb$chain)&masterdb$chain=="TRB"&masterdb$TRBV!=""& #P733
                      masterdb$Tcell=="CD4",]) #masterdb$Species=="Mouse"&
df<-unique(masterdb) #P735
df[df==""]<-NA
df$Sample_cdr3<-paste0(df$study2,":",df$disease2,df$cdr3.b,df$TRBV)
GLIPH.ref.df.filter$Sample_cdr3<-paste0(GLIPH.ref.df.filter$Sample,GLIPH.ref.df.filter$TcRb,GLIPH.ref.df.filter$V)
GLIPH.ref.df.filter.meta<-left_join(GLIPH.ref.df.filter[,c(1:19,31)],df,by="Sample_cdr3")
GLIPH.ref<-unique(GLIPH.ref.df.filter.meta[,c("pattern","type","TcRb","disease","antigen",
                                              "epitope","disease2","study2","final_score","Species")])

load("GLIPH.data.RData")
##making public repertoire df
df<-GLIPH.data
length(unique(intersect(GLIPH.ref.df.filter$pattern,GLIPH.data$pattern)))
length(unique(GLIPH.data$pattern))
imm.list<-list()
for(i in names(table(df$mouse_ID))){ #creating immunarch list from collapsed data
  df2<-df[df$mouse_ID==i&!is.na(df$index),]
  df2<-transform(df2,Clones=ave(seq(nrow(df2)),index,FUN=length))
  df2<-unique(df2[,c("index","Clones")])
  df2$Proportion<-df2$Clones/sum(df2$Clones)
  imm.list[[i]]<-tibble(Clones=df2$Clones, Proportion=df2$Proportion, CDR3.aa=df2$index )
}
immdata = repLoad("./rep")
immdata$meta$Sample<-immdata$meta$mouse_ID
pr = pubRep(imm.list, "aa", .coding = T, .verbose = F)
pr.AID = pubRepFilter(pr, immdata$meta, c(condition = "AID"))
pr.564 = pubRepFilter(pr, immdata$meta, c(condition = "m564"))
pr.AID[is.na(pr.AID)]<-0
pr.564[is.na(pr.564)]<-0
pr.AID[["avgfreq.AID"]] = rowMeans(public_matrix(pr.AID), na.rm = T)
pr.564[["avgfreq.564"]] = rowMeans(public_matrix(pr.564), na.rm = T)
pr.AID[["sum.AID"]] = rowSums(public_matrix(pr.AID)[,1:5], na.rm = T)+1
pr.564[["sum.564"]] = rowSums(public_matrix(pr.564)[,1:5], na.rm = T)+1
pr.res.GLIPH = dplyr::full_join(pr.AID, pr.564, by = c("CDR3.aa")) #,"J.name"
pr.res.GLIPH[is.na(pr.res.GLIPH)]<-0
pr.res.GLIPH$sum.AID[pr.res.GLIPH$sum.AID==0]<-1
pr.res.GLIPH$sum.564[pr.res.GLIPH$sum.564==0]<-1
pr.res.GLIPH[["Samples.sum"]] = pr.res.GLIPH[["Samples.x"]] + pr.res.GLIPH[["Samples.y"]]
pr.res.GLIPH[["freq.ratio"]] = apply(pr.res.GLIPH[, c("avgfreq.AID", "avgfreq.564")],1, function(x) log10(x[1])/log10(x[2]))
pr.res.GLIPH[["log2FC"]]= apply(pr.res.GLIPH[, c("sum.AID", "sum.564")],1, function(x) log2((x[2])/(x[1])))
pr.res.GLIPH<-left_join(pr.res.GLIPH,GLIPH.filter.groups[,c("index","final_score","number_subject","pattern",
                                                            "number_unique_cdr3","total.size","tfr.score")],
                        by=c("CDR3.aa"="index"))
labels<-pr.res.GLIPH[abs(pr.res.GLIPH$log2FC)>4|pr.res.GLIPH$sum.564>1700|pr.res.GLIPH$CDR3.aa %in% c("GLIPH_47"),]#,]
ggplot(pr.res.GLIPH,aes(x = sum.564, y =  sum.AID,size=number_subject))+theme_classic()+
  scale_x_continuous(trans='log10') +scale_y_continuous(trans='log10')+ scale_size(range = c(4, 10))+
  geom_point(alpha=0.25,aes(colour=number_unique_cdr3))+ #tfr.score
  guides(fill = guide_legend(override.aes = list(size = 7)))+#,color=F)+ #,size=F
  labs(x = "564Igi", y = "WT", size="Samples",color="Clones")+
  theme(legend.direction = "vertical", legend.box = "horizontal")+
  geom_label_repel(data = labels,aes(label = pattern),size = 3,alpha=0.7,
                   color="royalblue",min.segment.length=unit(0,'lines'),nudge_x=0.1,nudge_y=0.1,segment.size=0.25)+
  scale_color_gradient(low="grey",high="darkgreen",trans="log",breaks=c(5,25,100,400),labels=format(c(5,25,100,400)))
ggsave2("overlap.AID.564.GLIPH.scatter.indexlabel.png",width=6, height=4,device="png")
pr.res.GLIPH.db<-pr.res.GLIPH
#annotating with antigen prediction
df.db<-unique(GLIPH.ref[GLIPH.ref$antigen!="unknown"&!is.na(GLIPH.ref$antigen),]) #&GLIPH.ref$disease2!="species"
df.db<-df.db[order(df.db$final_score),]
pr.res.db<-unique(join(pr.res.GLIPH, df.db[,!names(df.db)%in%c("final_score")], by="pattern", type="left", match="first"))
match.df<-unique(pr.res.db[,c("pattern","disease","antigen","epitope","disease2")])
match.freq<-1-sum(is.na(match.df$disease))/nrow(match.df)
pr.res.db$disease2[is.na(pr.res.db$disease2)]<-"unknown"
labels<-unique(pr.res.db[(abs(pr.res.db$log2FC)>4&pr.res.db$Samples.sum>1&!is.na(pr.res.db$antigen)&
                            (pr.res.db$sum.AID>10|pr.res.db$sum.564>10))|(pr.res.db$sum.564>100&!is.na(pr.res.db$antigen)),])
labels2<-labels[labels$CDR3.aa %in% c("GLIPH_3662"),] 
brewer.pal(7,"Dark2") #display.brewer.pal
ggplot(pr.res.db,aes(x = sum.564, y =  sum.AID,size=Samples.sum))+theme_classic()+
  scale_x_continuous(trans='log10') +scale_y_continuous(trans='log10')+ scale_size(range = c(4, 10))+
  geom_point(aes(colour=factor(disease2),fill = factor(disease2)), shape=21) + 
  scale_color_manual(values = c("#1B9E77","#7570B3","#66A61E",alpha("black",0.5),"#A6761D"))+
  scale_fill_manual(values = c(alpha("#1B9E77",0.5),alpha("#7570B3",0.5),alpha("#66A61E",0.5),"#1C00ff00",alpha("#A6761D",0.5)))+
  guides(fill = guide_legend(override.aes = list(size = 5)),color=F)+ #,size=F
  labs(x = "564Igi", y = "WT", size="Samples",fill="Disease")+
  theme(legend.direction = "vertical", legend.box = "horizontal")
ggsave2("annot.index.GLIPH.scatter.bycondition.png",width=6.5, height=4,device="png")


#UMAP for specific GLIPH
load("GLIPH.analyzed.RData")
df<-T_all.combined
df@meta.data<-unique(join(df@meta.data, df.db[,!names(df.db)%in%c("final_score")], by="pattern", type="left", match="first"))
df@meta.data$disease2[is.na(df@meta.data$disease2)]<-"unknown"
Idents(df)<-"disease2"
DimPlot(df, reduction = "umap", split.by = "condition2",order=T, #"pattern","type","TcRb","disease","antigen","epitope","disease2","study2","final_score","Species"
        cells.highlight=list(autoimmune=WhichCells(df, idents = c("autoimmune")),
                             virus=WhichCells(df, idents = c( "virus"))),
        cols.highlight=c("#A6761D","#1B9E77"),cols="grey",sizes.highlight=0.3)+ #NoLegend()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text=element_blank(),axis.ticks=element_blank(),axis.line=element_blank(),axis.title=element_blank())
ggsave2("GLIPH.dz2.umap.png",width=7, height=3,device="png")
for(j in c("disease","antigen","epitope","pattern")){
  Idents(df)<-j
  cell.list<-list()
  for(i in names(head(sort(table(df[[j]]), decreasing = T), n=9))){
    cell.list[[i]]<-WhichCells(df, idents = i)
  }
  DimPlot(df, reduction = "umap", split.by = "condition2",order=T, #"pattern","type","TcRb","disease","antigen","epitope","disease2","study2","final_score","Species"
          cells.highlight=cell.list,cols.highlight=brewer.pal(n = length(cell.list), name = "Paired"),
          sizes.highlight=0.3)+ #NoLegend()+cols="grey",
    theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
          axis.text=element_blank(),axis.ticks=element_blank(),axis.line=element_blank(),axis.title=element_blank())
  ggsave2(paste0("GLIPH.",j,".umap.png"),width=7.5, height=2.75,device="png")
}

##GLIPH characteristic frequency graphs
seurat.meta.first<-unique(join(T_all.combined@meta.data, unique(df.db[,c("disease","antigen","epitope","pattern","disease2")]), 
                               by="pattern", type="left", match="first"))
seurat.meta.all<-unique(left_join(x = T_all.combined@meta.data, y = unique(df.db[,c("disease","antigen","epitope","pattern","disease2")]),
                                  by = "pattern",keep=F))
meta<-seurat.meta.first
for(j in c("disease","antigen","epitope","disease2")){
  df<-meta[!is.na(meta[[j]])&meta[[j]]!="unknown",] 
  prediction.Counts <- as.data.frame(table(df[[j]],df$condition2))
  top9<-names(sort(table(df[[j]]),decreasing=T))[1:9]
  prediction.Counts<-prediction.Counts[prediction.Counts$Var1 %in% top9,]
  color.order<-c("#CAB2D6","#FF7F00","#FDBF6F","#E31A1C","#FB9A99","#33A02C","#B2DF8A","#1F78B4","#A6CEE3")
  ggplot(prediction.Counts, aes(fill=Var1,y=Freq, x=Var2,alluvium=Var1,stratum=Var1)) +
    geom_lode()+geom_flow()+geom_stratum(alpha=0) +theme_classic()+
    theme(legend.title = element_blank(),axis.title=element_blank())+ 
    scale_fill_manual(values=color.order)
  ggsave2(paste0(j,".prediction.prop.flow.png"),width=4, height=2.5,device="png")
}



