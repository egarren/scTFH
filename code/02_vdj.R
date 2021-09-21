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
library(tools)
library(dplyr)
library(alakazam)
library(ggalluvial)
library(tidyverse)
library(RColorBrewer)
library(gridExtra)
library(circlize)
library(useful)
library(ggseqlogo)
library(yingtools2)
library(qgraph)
library(vegan)
gm_mean = function(x, na.rm=TRUE){  #geometric mean functions for graphing
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
} 
gm_ymax<- function(x, na.rm = TRUE)
{
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))+exp(sd(log(x), na.rm = na.rm))
}
barcoder<-function(df, prefix){
  df$ms.barcoder<-prefix
  df$orig.barcode<-df$barcode
  df$barcode <- gsub("\\-1", "", df$barcode)
  df$barcode <- paste0(prefix, df$barcode)
  df$ms_clone <- paste0(prefix, df$raw_clonotype_id) 
  df$ms_v <- paste0(prefix, df$v_gene) 
  df$ms_j<- paste0(prefix, df$j_gene) 
  df$ms_cdr3 <- paste0(prefix, df$cdr3) 
  df
}
data_concater <- function(x){
  x<- levels(factor(x))
  paste(x, collapse = "+")
}
data_concater2 <- function(x){
  x<- levels(factor(x))
  paste(x[1])
}
data_concater3 <- function(x){
  x<- levels(factor(x))
  if(length(x)>1){paste(x[1:2], collapse = "+")}else{paste(x[1])}
}
quantile_breaks <- function(xs, n = 100) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n),na.rm=T)
  breaks[!duplicated(breaks)]
}

load("graphed.RData")

metadata<-read.csv("Tcell.metadata.csv", header=T) 
metadata[] <- lapply(metadata, as.character)

#load cdr3 data
mice<-metadata$T.sampleID
for (i in mice){
  csv.path<-paste0("./vdj/",i,"/outs/filtered_contig_annotations.csv") #path to 10X output
  fasta.path<-paste0("./vdj/",i,"/outs/filtered_contig.fasta") #path to 10X output
  df<-merge(read.csv(csv.path, header=T),read.fasta(fasta.path),by.x="contig_id",by.y="seq.name")
  colnames(df)[colnames(df)=="length"]<-"vdj_length"
  assign(paste0(i,".vdj"),df)
}

#rename barcodes
for (i in mice){
  assign(paste0(i,".vdj"),barcoder(get(paste0(i,".vdj")),prefix=paste0(i,"_")))
}
Tall.vdj<-rbind(T1.vdj,T2.vdj,T3.vdj,T4.vdj,T7.vdj,T8.vdj,T11.vdj,T12.vdj,T34.vdj,T36.vdj)

#extract and add Seurat metadata
clone.data<-Tall.vdj[Tall.vdj$productive=="True"&!is.na(Tall.vdj$productive)&Tall.vdj$chain %in% c("TRA","TRB"),] 
clone.data<-transform(clone.data, vab.freq = ave(seq(nrow(clone.data)), v_gene, FUN=length))
T_all.combined@meta.data$orig.barcode<-rownames(T_all.combined@meta.data)
clone.data<-left_join(x = clone.data, y = subset(T_all.combined@meta.data, select=-c(ms.barcoder)),by = c("barcode"="orig.barcode"),keep=F)
clone.data<-left_join_replace(clone.data,metadata,by="ms.barcoder")
clone.data.seurat<-clone.data[!is.na(clone.data$my.clusters),] 
clone.meta<- as.data.table(clone.data)[, lapply(.SD, data_concater), by=cdr3]
clone.meta$cdr3.freq <- as.numeric(as.character(clone.meta$cdr3.freq))
clone.meta$ms_cdr3.freq <- as.numeric(as.character(clone.meta$ms_cdr3.freq))
clone.meta$aa_length<-nchar(clone.meta$cdr3)
clusters<-unique(clone.data.seurat$my.clusters)
mice<-unique(clone.data$mouse_ID)
AID.mice<-unique(clone.data[clone.data$condition=="AID",]$mouse_ID)
m564.mice<-unique(clone.data[clone.data$condition=="m564",]$mouse_ID)
save(list=c(ls(pattern="clone."),"metadata","Tall.vdj","clusters",ls(pattern="mice")),file="clone.data.RData")

#TCR venn diagram
df<-clone.data
dz.list<-list()
for(i in unique(df$chain)){
  dz.list[[i]]<-unique(df$barcode[df$chain==i])
}
venn.diagram(
  x = dz.list,category.names = names(dz.list),filename = "TR.venn.png",
  imagetype="png",height=700,width=700,margin=0.1,
  lwd=1,col=c("#440154ff", '#21908dff'),fill=c(alpha("#440154ff",0.3), alpha('#21908dff',0.3)),
  # lty="blank",fill=sample(brewer.pal(9, "Set1"),size=length(dz.list)),
  cex=0.3,fontface="bold",fontfamily="sans",
  cat.cex=0.6,cat.fontface="bold",cat.default.pos="outer",cat.fontfamily="sans",
  cat.col = c("#440154ff", '#21908dff'),cat.dist = c(0.1, 0.1)
)

##CDR3 property analysis by cluster
#WebLogo
#all clones w/n all mice, comparing conditions
plot.list<-list()
for(j in c("AID","m564")){
  for(i in c("TRA","TRB")){
    df<-clone.data[clone.data$chain==i & clone.data$condition ==j,]
    freq.tab<-as.data.frame(table(df$cdr3))  #ms_cdr3 frequency table for given cluster
    freq.tab$aa_length<-nchar(as.character(freq.tab$Var1))
    for(k in 14){ #10:18
      if(j=="AID"){name="WT"}else{name="564Igi"}
      df2<-freq.tab[freq.tab$aa_length==k,]
      plot.list[[paste0(i,j,k)]]<-ggplot()+geom_logo(as.character(df2$Var1),method="probability")+
        theme_logo()+ggtitle(paste0(i,"_",name))+
        theme(plot.title = element_text(hjust = 0.5)) #,"_",k, axis.text.x=element_blank(),
    }
  }
}
CombinePlots(plot.list,ncol=2,legend="bottom")
ggsave2("cdr3.weblogo.bycondition.png",width=6,height=4,device="png")

#aa length
#all clones w/n all mice, comparing conditions
plot.list<-list()
plot.list2<-list()
for(i in c("TRA","TRB")){
  df.list<-list()
  for(j in c("AID","m564")){
    df<-clone.data[clone.data$chain==i & clone.data$condition ==j,]
    freq.tab<-as.data.frame(table(df$cdr3))  #ms_cdr3 frequency table for given cluster
    freq.tab$aa_length<-nchar(as.character(freq.tab$Var1))
    if(j=="AID"){freq.tab$condition<-"WT"}else{freq.tab$condition<-"564Igi"}
    df.list[[j]]<-freq.tab #add to list for each condition
  }
  df<-rbindlist(df.list)
  df$condition <- factor(df$condition, levels = c("WT","564Igi"))
  mu <- ddply(df, "condition", summarise, grp.mean=mean(aa_length))
  plot.list[[i]]<-ggplot(df, aes(x=aa_length,color=condition,fill=condition)) + 
    geom_histogram(aes(y=..density..),binwidth=1,position="identity",alpha=0.2) +#geom_density(fill=NA)+
    geom_vline(data=mu, aes(xintercept=grp.mean, color=condition),linetype="dashed")+ggtitle(i)+theme_classic()+
    labs(x = "CDR3 Length", y = "") +theme(plot.title = element_text(hjust = 0.5),legend.title = element_blank())+
    scale_color_manual(values = c("black", "red"))+
    scale_fill_manual(values = c("black", "red"))
  plot.list2[[i]]<-ggbarplot(df, x = "condition", y = "aa_length", add = c("mean","jitter"),
                             position = position_dodge(0.8), legend="right")+
    stat_compare_means(aes(group = condition))+
    labs(x="condition", y = "CDR3 Length")+ggtitle(i)
}
CombinePlots(plot.list,ncol=2,legend="right")
ggsave2("cdr3.aa_length.histo.png",width=7,height=3,device="png")
CombinePlots(plot.list2,ncol=2,legend="right")
ggsave2("cdr3.aa_length.bar.dot.png",width=8,height=4,device="png")

#Heatmap by mouse v usage
plot.list<-list()
for(i in c("TRA","TRB")){
  df<-clone.data[clone.data$chain==i,]
  tab<-table(df$v_gene,df$mouse_ID)
  Proportions <- scale(tab,scale=colSums(tab),center=FALSE)*100 
  annot.col<-data.frame(BMChimera=metadata$condition)
  annot.col$BMChimera[annot.col$BMChimera %in% "AID"]<-"WT"
  annot.col$BMChimera[annot.col$BMChimera %in% "m564"]<-"564Igi"
  rownames(annot.col)<-metadata$mouse_ID
  plot.list[[i]]<-pheatmap(as.data.frame.matrix(Proportions),annotation_col=annot.col,main=i,annotation_names_col=F,
                           annotation_colors=list(BMChimera=c("WT"="grey","564Igi"="red")),fontsize_row=5)
}
save_pheatmap_png <- function(x, filename, width=700, height=1200, res = 150) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_png(plot.list[["TRA"]], "TRA.v_gene.heatmap.png")
save_pheatmap_png <- function(x, filename, width=700, height=600, res = 150) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_png(plot.list[["TRB"]], "TRB.v_gene.heatmap.png")

###VDJtools export
dir.create("./vdjtools")
setwd("./vdjtools")
mice<-metadata$T.sampleID
#by condition 
file.list=list()
for (i in mice) { 
  db<-clone.data[clone.data$ms.barcoder==paste0(i,"_"),]
  db<-transform(db, ab.cdr3.freq = ave(seq(nrow(db)), cdr3, FUN=length)) #calculate cdr3 frequencies
  db<- as.data.table(db)[, lapply(.SD, data_concater2), by=cdr3]   #collapse by cdr3
  db$ab.cdr3.freq <- as.numeric(as.character(db$ab.cdr3.freq))
  vdj <- db %>% transmute(count=ab.cdr3.freq,   #make vdjtools table
                          freq=ab.cdr3.freq/sum(ab.cdr3.freq, na.rm=TRUE),
                          cdr3nt=cdr3_nt,
                          cdr3aa=cdr3,
                          v=v_gene,
                          d=d_gene,
                          j=j_gene)
  write.table(vdj, file=paste0(i,"_vdjtools.tab"), quote=FALSE, sep="\t", row.names=FALSE)
  file.list[i]<-paste0(i,"_vdjtools.tab")
}
vdjtools.metadata<-data.frame(file.name=unlist(file.list,use.names=F),sample.id=names(file.list))
vdjtools.metadata<-left_join(x = vdjtools.metadata, y = metadata, by = c("sample.id"="T.sampleID"),keep=F)
write.table(vdjtools.metadata, file = "vdjtools.metadata.txt", sep = "\t",quote=F,row.names = FALSE)
setwd("../")
#by cluster
dir.create("./vdjtools.cluster")
setwd("./vdjtools.cluster")
vdjtools.metadata2<-data.frame()
for (i in mice) { 
  db<-clone.data[clone.data$ms.barcoder==paste0(i,"_"),]
  for(j in clusters){
    db<-db[!is.na(db$my.clusters),]
    db2<-db[db$my.clusters==j,]
    db2<-transform(db2, ab.cdr3.freq = ave(seq(nrow(db2)), cdr3, FUN=length)) #calculate cdr3 frequencies
    db2<- as.data.table(db2)[, lapply(.SD, data_concater2), by=cdr3]   #collapse by cdr3
    db2$ab.cdr3.freq <- as.numeric(as.character(db2$ab.cdr3.freq))
    vdj <- db2 %>% transmute(count=ab.cdr3.freq,   #make vdjtools table
                             freq=ab.cdr3.freq/sum(ab.cdr3.freq, na.rm=TRUE),
                             cdr3nt=cdr3_nt,
                             cdr3aa=cdr3,
                             v=v_gene,
                             d=d_gene,
                             j=j_gene)
    write.table(vdj, file=paste0(i,".clust",j,"_vdjtools.tab"), quote=FALSE, sep="\t", row.names=FALSE)
    vdjtools.metadata2[paste0(i,".clust",j),1]<-paste0(i,".clust",j,"_vdjtools.tab")
    vdjtools.metadata2[paste0(i,".clust",j),2]<-paste0(i,".clust",j)
    vdjtools.metadata2[paste0(i,".clust",j),3]<-i
    vdjtools.metadata2[paste0(i,".clust",j),4]<-j
  }
}
colnames(vdjtools.metadata2)<-c("file.name","sample.id","mouse","Cluster")
vdjtools.metadata2<-left_join(x = vdjtools.metadata2, y = metadata, by = c("mouse"="T.sampleID"),keep=F)
write.table(vdjtools.metadata2, file = "vdjtools.clust.metadata.txt", sep = "\t",quote=F,row.names = FALSE)
setwd("../")

##DeepTCR export
#scTfh data
dir.create("./deepTCR")
setwd("./deepTCR")
dir.create("./scTfh.data")
setwd("./scTfh.data")
for(i in names(table(clone.data$condition))){
  df<-clone.data[clone.data$condition==i,]
  dir.create(paste0("./",i))
  setwd(paste0("./",i))
  for(j in names(table(df$mouse_ID))){
    df2<-df[df$mouse_ID==j,]
    df.TRA<-df2[df2$chain=="TRA",]
    df.TRB<-df2[df2$chain=="TRB",]
    cells<-inner_join(df.TRB[,c("barcode","cdr3","v_gene","j_gene","d_gene")],
                      df.TRA[,c("barcode","cdr3","v_gene","j_gene")],by="barcode")
    cells<-left_join(cells,unique(clone.data[,c("barcode","mouse_ID","condition","my.clusters")]),by="barcode")
    cells$clone_ID<-paste0(cells$cdr3.x,cells$v_gene.x,cells$j_gene.x,
                           cells$cdr3.y,cells$v_gene.y,cells$j_gene.y,cells$mouse_ID)
    cells.cloneID<-transform(cells,freq=ave(seq(nrow(cells)),clone_ID,FUN=length))
    cells.cloneID<-as.data.table(cells.cloneID)[, lapply(.SD, data_concater), by=clone_ID]
    cells.cloneID$freq<-as.numeric(as.character(cells.cloneID$freq))
    write.table(cells.cloneID,file=paste0(j,".deepTCR.tsv"),sep="\t",row.names=F,col.names=F,quote=F)
  }
  setwd("../")
}
setwd("../")
#making count matrix (for ML regression)
df<-clone.data
df.TRA<-df[df$chain=="TRA",]
df.TRB<-df[df$chain=="TRB",]
cells<-inner_join(df.TRB[,c("barcode","cdr3","v_gene","j_gene")],
                  df.TRA[,c("barcode","cdr3","v_gene","j_gene")],by="barcode")
cells<-left_join(cells,unique(clone.data[,c("barcode","mouse_ID","condition","my.clusters")]),by="barcode")
cells$cdr3_ID<-paste0(cells$cdr3.x,"+",cells$cdr3.y)
df<-as.data.frame.matrix(table(cells$cdr3_ID,cells$condition))
df$cdr3_ID<-rownames(df)
df2<-left_join(df,unique(cells[,c("cdr3_ID","cdr3.x","cdr3.y")]),by="cdr3_ID")
df3<-data.frame(alpha=df2$cdr3.y,beta=df2$cdr3.x,AID=df2$AID,m564=df2$m564)
df3$alpha<-gsub("[^GALMFWKQESPVICYHRNDT]+", "", df3$alpha)
df3$beta<-gsub("[^GALMFWKQESPVICYHRNDT]+", "", df3$beta)
write.csv(df3, "BMchim.counts_regression.csv",row.names=F,col.names=T,quote=F)
#training data
load("masterdb4.RData")
dir.create("./training.data")
setwd("./training.data")
# write.table(df,file="masterdb.tsv",sep="\t",row.names=F,col.names=F,quote=F)
for(i in c("disease","antigen","epitope","disease2")){
  df<-get(paste0("masterdb.",i))
  df<-unique(df[!is.na(df$chain)&df$chain=="TRB"&df$TRBV!="",]) #&masterdb$Species=="Mouse"
  dir.create(paste0("./",i))
  setwd(paste0("./",i))
  for(k in unique(df[[i]])){
    df2<-df[df[[i]]==k&!is.na(df[[i]]),]
    df3<-unique(df2[,c("cdr3.b","TRBV","TRBJ")])
    if(nrow(df3)>1){
      dir.create(paste0("./",k))
      setwd(paste0("./",k))
      write.table(df3,file=paste0(k,".deepTCR.tsv"),sep="\t",row.names=F,col.names=F,quote=F)
      setwd("../")
    }
  }
  setwd("../")
}
setwd("../")
setwd("../")

#paired data
clone.data.ab<-Tall.vdj[Tall.vdj$productive=="True"&!is.na(Tall.vdj$productive)&Tall.vdj$chain %in% c("TRA","TRB"),] 
clone.data.ab<-as.data.table(clone.data.ab)[, lapply(.SD, data_concater3), by=barcode]
clone.data.ab<-clone.data.ab[clone.data.ab$chain=="TRA+TRB",]
clone.data.ab<-transform(clone.data.ab, cdr3.freq = ave(seq(nrow(clone.data.ab)), cdr3, FUN=length))
clone.data.ab<-transform(clone.data.ab, ms_cdr3.freq = ave(seq(nrow(clone.data.ab)), ms_cdr3, FUN=length))
clone.data.ab<-transform(clone.data.ab, vab.freq = ave(seq(nrow(clone.data.ab)), v_gene, FUN=length))
T_all.combined@meta.data$orig.barcode<-rownames(T_all.combined@meta.data)
clone.data.ab<-left_join(x = clone.data.ab, y = subset(T_all.combined@meta.data, select=-c(ms.barcoder)), by = c("barcode"="orig.barcode"),keep=F)
clone.data.ab<-left_join_replace(clone.data.ab,metadata,by="ms.barcoder")
clone.data.ab.seurat<-clone.data.ab[!is.na(clone.data.ab$my.clusters),] 
clone.meta.ab<-as.data.table(clone.data.ab)[, lapply(.SD, data_concater), by=cdr3]
clone.meta.ab$cdr3.freq<-as.numeric(as.character(clone.meta.ab$cdr3.freq))
# clone.meta.ab$aa_length<-nchar(clone.meta.ab$cdr3)
clone.data.ab.ms_cdr3<-as.data.table(clone.data.ab)[, lapply(.SD, data_concater), by=ms_cdr3]
clone.data.ab.ms_cdr3$cdr3.freq <- as.numeric(as.character(clone.data.ab.ms_cdr3$cdr3.freq))
clone.data.ab.ms_cdr3$ms_cdr3.freq <- as.numeric(as.character(clone.data.ab.ms_cdr3$ms_cdr3.freq))
clone.data.ab$clust_cdr3<-paste0(clone.data.ab$my.clusters2,"_",clone.data.ab$cdr3)
clone.data.ab.clust_cdr3<-as.data.table(clone.data.ab)[, lapply(.SD, data_concater), by=clust_cdr3]
clone.data.ab.clust_cdr3$cdr3.freq <- as.numeric(as.character(clone.data.ab.clust_cdr3$cdr3.freq))
clusters<-unique(clone.data.seurat$my.clusters)
mice<-unique(clone.data$mouse_ID)
AID.mice<-unique(clone.data[clone.data$condition=="AID",]$mouse_ID)
m564.mice<-unique(clone.data[clone.data$condition=="m564",]$mouse_ID)
save(list=c(ls(pattern="clone."),"metadata","Tall.vdj","clusters",ls(pattern="mice")),file="clone.ab.data.RData")

##Venn diagram
#between conditions
df<-clone.data.ab
dz.list<-list()
for(i in unique(df$condition)){
  dz.list[[i]]<-unique(df$cdr3[df$condition==i])
}
venn.diagram(
  x = dz.list,category.names = c("WT","564Igi"),filename = "condition.cdr3.venn.png",
  imagetype="png",height=700,width=700,margin=0.1,
  lwd=1,col=c("black", 'red'),fill=c(alpha("black",0.3), alpha('red',0.3)),
  # lty="blank",fill=sample(brewer.pal(9, "Set1"),size=length(dz.list)),
  cex=0.3,fontface="bold",fontfamily="sans",
  cat.cex=0.6,cat.fontface="bold",cat.default.pos="outer",cat.fontfamily="sans",
  cat.col = c("black", 'red'),cat.dist = c(0.1, 0.1),
  ext.pos=180,ext.line.lwd=0.25,ext.dist=-0.2
)
#between clusters
df<-clone.data.ab[clone.data.ab$my.clusters %in% c("0","1","2","3"),]
dz.list<-list()
for(i in c("Tfh-activated","Tfr","Sostdc1","Tfh-CM")){
  dz.list[[i]]<-unique(df$cdr3[df$my.clusters2==i])
}
venn.colors<-hue_pal()(7)[1:4]
venn.diagram(
  x = dz.list,category.names = names(dz.list),filename = "Tfh.Tfr.cdr3.venn.png",
  imagetype="png",height=1200,width=1200,margin=0.1,
  lwd=1,col=venn.colors,fill=alpha(venn.colors,0.3),
  # lty="blank",fill=sample(brewer.pal(9, "Set1"),size=length(dz.list)),
  cex=0.5,fontface="bold",fontfamily="sans",
  cat.cex=0.5,cat.fontface="bold",cat.default.pos="outer",cat.fontfamily="sans",
  cat.col = venn.colors#,cat.dist = c(0.1, 0.1,0.1,0.1)
)

##clone size
# all clones w/n all mice, comparing conditions
df.list<-list()
for(j in c("AID","m564")){
  df<-clone.data.ab[clone.data.ab$condition ==j,]
  freq.tab<-as.data.frame(table(df$cdr3))  #ms_cdr3 frequency table for given cluster
  if(j=="AID"){freq.tab$condition<-"WT"}else{freq.tab$condition<-"564Igi"}
  df.list[[j]]<-freq.tab #add to list for each cluster
}
df<-rbindlist(df.list)
ggbarplot(df, x = "condition", y = "Freq", color="condition", add = c("mean","jitter"),
          position = position_dodge(0.8), legend="none",palette=c("black","red"),xlab=F)+
  labs(x="", y = "Clone Size") #+scale_y_log10()
ggsave2("bar.dot.condition.cdr3.freq.png",width=2, height=3,device="png")
mu <- ddply(df, "condition", summarise, grp.mean=gm_mean(Freq))
df$condition <- factor(df$condition, levels = c("WT","564Igi"))
ggplot(df, aes(x=Freq,color=condition,fill=condition)) + 
  geom_histogram(aes(y=..density..),bins=10,position="identity",alpha=0.2) +#geom_density(fill=NA)+
  geom_vline(data=mu, aes(xintercept=grp.mean, color=condition),linetype="dashed")+theme_classic()+
  labs(x = "Clone Size", y = "Density")+scale_x_log10()+theme(legend.title = element_blank())+
  scale_color_manual(values = c("black", "red"))+
  scale_fill_manual(values = c("black", "red"))
ggsave2("histogram.condition.cdr3.freq.png",width=4, height=3,device="png")
#all clones w/n all mice w/n cluster
df.list <- list()
plot.list<-list()
mycluster.names<-names(sort(table(clone.data.ab.seurat$my.clusters2),decreasing=T))
for (i in mycluster.names){
  for(j in c("AID","m564")){
    df<-clone.data.ab.seurat[clone.data.ab.seurat$my.clusters2==i&clone.data.ab.seurat$condition==j,] #gating on given cluster
    freq.tab<-as.data.frame(table(df$cdr3))  #ms_cdr3 frequency table for given cluster
    if(j=="AID"){freq.tab$condition<-"WT"}else{freq.tab$condition<-"564Igi"}
    freq.tab$my.clusters2<-i
    df.list[[paste0(i,j)]]<-freq.tab #add to list for each cluster
  }
  df2<-rbindlist(df.list)
}
df<-rbindlist(df.list)
df$condition <- factor(df$condition, levels = c("WT","564Igi"))
ggbarplot(df, x = "my.clusters2", y = "Freq", color="condition", add = c("mean","jitter"),
          position = position_dodge(0.8), legend="right",palette=c("black","red"),xlab=F)+
  labs(y = "Clone Size")+theme(legend.title = element_blank(),axis.text.x = element_text(angle = 45, hjust = 1))+ 
  scale_x_discrete(limits = mycluster.names) #+scale_y_log10()
ggsave2("bar.dot.cluster.cdr3.freq.png",width=4, height=3,device="png")

##top 10 cdr3 percent (diversity), within cluster, by mouse
df.list <- list()
for (i in unique(clone.data.ab.seurat$my.clusters2)){
  df.cluster<-clone.data.ab.seurat[clone.data.ab.seurat$my.clusters2==i,] #gating on given cluster
  ms.clone.freq<-as.data.frame(table(df.cluster$ms_cdr3)) #ms_cdr3 frequency table for given cluster
  ms.clone.freq$ms_cdr3<-as.character(ms.clone.freq$Var1) 
  for(k in c("m564","AID")){
    meta2<-metadata[metadata$condition==k,]
    meta<-data.frame(mice=names(table(meta2$ms.barcoder))) #blank dataframe of mice (to fill in upcoming loop)
    if(k=="AID"){meta$condition<-"WT"}else{meta$condition<-"564Igi"}
    meta$top10.pct<-0
    meta$my.clusters2<-i
    for (j in 1:nrow(meta)) {
      df<-ms.clone.freq[grepl(meta$mice[j],ms.clone.freq$ms_cdr3), ] #frequency info for cdrs3 from given mouse
      df <- df[order(df$Freq,decreasing=T),] #calculating what % of all clones the top 10 clones represent
      pct <- 100*(df$Freq/sum(df$Freq))
      meta$top10.pct[j]<-sum(pct[1:10]) #add to list for each mouse
    }
    df.list[[paste0(i,k)]]<-meta
  }
}
df<-rbindlist(df.list)
df$condition <- factor(df$condition, levels = c("WT","564Igi"))
ggbarplot(df, x = "my.clusters2", y = "top10.pct", color="condition", add = c("mean_se","jitter"),
          position = position_dodge(0.8), legend="right",palette=c("black","red"),xlab=F)+
  stat_compare_means(aes(group = condition),label = "p.format", method="t.test",size=2.5,label.y=110)+
  labs(y = "% of repertoire in \n 10 larges clones")+
  theme(legend.title = element_blank(),axis.text.x = element_text(angle = 45, hjust = 1))+ 
  scale_x_discrete(limits = mycluster.names) #+scale_y_log10()
ggsave2("bar.dot.cluster.cdr3.top10freq.png",width=5, height=3,device="png")


#Rarefaction plot (https://peat-clark.github.io/BIO381/veganTutorial.html)
for(i in c("0","1","2","3","4","5","all")){
  df.cluster<-clone.data.ab.seurat[clone.data.ab.seurat$my.clusters==i,] #gating on given cluster
  if(i=="all"){df.cluster<-clone.data.ab}
  raref.df<-as.data.frame.matrix(table(df.cluster$mouse_ID,df.cluster$cdr3))
  colors<-rep("black",nrow(raref.df))
  colors[rownames(raref.df) %in% m564.mice]<-"red"
  # rarecurve(raref.df)
  png(paste0("clust.",i,".raref.png"),width=3.5,height=3.5,units="in",res=200)
  rarecurve(raref.df, step = 20, sample = min(rowSums(raref.df)), col = alpha(colors,0.6), cex = 0.2,
            ylab="Clones",label=F,lwd=4)
  dev.off()
}


#non-metric multidimensional scaling (by condition)
df<-clone.data.ab
community_matrix<-as.data.frame.matrix(table(df$mouse_ID,df$cdr3))
example_NMDS=metaMDS(community_matrix, k=2,trymax=100,threshold=0.9999)
png("tcr.NMDS.stres.png",width=4,height=4,units="in",res=200)
stressplot(example_NMDS,p.col="forestgreen",l.col="black")
dev.off()
plot(example_NMDS)
treat=rep("WT",nrow(community_matrix))
treat[rownames(community_matrix) %in% m564.mice]<-"564Igi"
colors=rep("black",nrow(community_matrix))
colors[rownames(community_matrix) %in% m564.mice]<-"red"
elevation=runif(10,0.5,1.5)
png("tcr.NMDS.png",width=6,height=6,units="in",res=200)
ordisurf(example_NMDS,elevation,main="",col="forestgreen")
ordiellipse(example_NMDS$point[grep("WT",treat),],draw="polygon",groups=treat[treat=="WT"],col="grey",alpha=50)
ordiellipse(example_NMDS$point[grep("564Igi",treat),],draw="polygon",groups=treat[treat=="564Igi"],col="red",alpha=50)
orditorp(example_NMDS,display="species",col="black",pch=21,bg=adjustcolor("grey",0.5),pcex=1,air=50000,cex=0.1)
orditorp(example_NMDS,display="sites",col=colors,air=0.0001,cex=1.25)
dev.off()


#non-metric multidimensional scaling (by co=luster)
df<-clone.data.ab.seurat
community_matrix<-as.data.frame.matrix(table(df$my.clusters2,df$cdr3))
example_NMDS=metaMDS(community_matrix, k=5,trymax=1000,threshold=0.99)
png("tcr.NMDS.cluster.stres.png",width=4,height=4,units="in",res=200)
stressplot(example_NMDS,p.col="forestgreen",l.col="black")
dev.off()
# plot(example_NMDS)
colors=hue_pal()(6)
elevation=runif(14,0.5,1.5)
png("tcr.NMDS.cluster.png",width=6,height=6,units="in",res=200)
ordiplot(example_NMDS,type="n")
orditorp(example_NMDS,display="species",col="black",pch=21,bg=adjustcolor("grey",0.5),pcex=1,air=50000,cex=0.1)
orditorp(example_NMDS,display="sites",col=colors,air=0.0001,cex=1.25)
dev.off()



##Clonotype Distribution Heatmaps
#cluster propotions by CDR3
clone.data.ab.top<-clone.data.ab.seurat[clone.data.ab.seurat$cdr3.freq>30,]
clone.data.ab.top<-clone.data.ab.top[clone.data.ab.top$cdr3!="None",]
Counts <- table(as.character(clone.data.ab.top$my.clusters2),as.character(clone.data.ab.top$cdr3))
Proportions <- scale(Counts,scale=colSums(Counts),center=FALSE)*100 
pheatmap(as.data.frame.matrix(Proportions),cutree_cols = 6)
#add heatmap metadata
clone.meta.top<- as.data.table(clone.data.ab.top)[, lapply(.SD, data_concater), by=cdr3]
clone.meta.top$cdr3.freq <- as.numeric(as.character(clone.meta.top$cdr3.freq))
annot.col<-data.frame(cdr3=colnames(Proportions))
annot.col<-merge(annot.col,clone.meta.top[,c("cdr3","cdr3.freq","condition")],by="cdr3")
names(annot.col)[names(annot.col) == "condition"] <- "BMChimera"
names(annot.col)[names(annot.col) == "cdr3.freq"] <- "Clone Size"
annot.col$BMChimera[annot.col$BMChimera %in% "AID"]<-"WT"
annot.col$BMChimera[annot.col$BMChimera %in% "m564"]<-"564Igi"
annot.col$BMChimera[annot.col$BMChimera %in% "AID+m564"]<-"Both"
row.names(annot.col) <- annot.col$cdr3
annot.col$cdr3 <- NULL
p1<-pheatmap(as.data.frame.matrix(Proportions),annotation_col=annot.col,fontsize_col=3,cutree_cols = 5,annotation_names_col=F,
             annotation_colors=list(BMChimera=c("WT"="grey","564Igi"="red","Both"=alpha("red",0.5))))
save_pheatmap_png <- function(x, filename, width=1400, height=700, res = 150) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_png(p1, "clone.cluster.heatmap.png")


#public clones
public.clones<-unique(unlist(ms.overlap.list,use.names=F))
clone.data.ab.public<-clone.data.ab[clone.data.ab$cdr3 %in% public.clones,]
Counts <- table(as.character(clone.data.ab.public$mouse_ID),as.character(clone.data.ab.public$cdr3))
Proportions <- scale(Counts,scale=colSums(Counts),center=FALSE)*100 
pheatmap(as.data.frame.matrix(Proportions))
#add heatmap metadata
clone.meta.top<- as.data.table(clone.data.ab.public)[, lapply(.SD, data_concater), by=cdr3]
clone.meta.top$cdr3.freq <- as.numeric(as.character(clone.meta.top$cdr3.freq))
annot.col<-data.frame(cdr3=colnames(Proportions))
annot.col<-merge(annot.col,clone.meta.top[,c("cdr3","cdr3.freq","condition")],by="cdr3")
names(annot.col)[names(annot.col) == "condition"] <- "BMChimera"
names(annot.col)[names(annot.col) == "cdr3.freq"] <- "Clone Size"
annot.col$BMChimera[annot.col$BMChimera %in% "AID"]<-"WT"
annot.col$BMChimera[annot.col$BMChimera %in% "m564"]<-"564Igi"
annot.col$BMChimera[annot.col$BMChimera %in% "AID+m564"]<-"Both"
row.names(annot.col) <- annot.col$cdr3
annot.col$cdr3 <- NULL
mat_breaks <- quantile_breaks(as.data.frame.matrix(Proportions), n = 101)
colfunc <- colorRampPalette(c("#EDF8E9", "#006D2C"))
p1<-pheatmap(as.data.frame.matrix(Proportions),annotation_col=annot.col,annotation_names_col=F,fontsize_col=5,
             # color=cividis(length(mat_breaks) - 1),breaks=mat_breaks,
             annotation_colors=list(BMChimera=c("WT"="grey","564Igi"="red","Both"=alpha("red",0.5)))) #
save_pheatmap_png <- function(x, filename, width=1200, height=800, res = 150) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_png(p1, "clone.public.cluster.heatmap.png")

###CDR3 Network Analysis by condition
df<-clone.data.ab.ms_cdr3[clone.data.ab.ms_cdr3$cdr3.freq>10,]
# extract link between diag and patients ID for network
tmp1=data.frame(Diagnosis=df$condition,ID=df$mouse_ID,stringsAsFactors = F)
tmp1=tmp1[which(duplicated(paste(tmp1[,1],tmp1[,2],sep="_"))==F),]
colnames(tmp1)=c("x","y")
# extract link between clonotype_id.long and patients ID for network
tmp2=data.frame(df$mouse_ID,df$cdr3,stringsAsFactors = F)
tmp2=tmp2[which(duplicated(paste(tmp2[,1],tmp2[,2],sep="_"))==F),]
colnames(tmp2)=c("x","y")
# create table for qgraph
toQgraph=rbind(tmp1,tmp2)
# color of Dx
l=unique(c(toQgraph[,1],toQgraph[,2]))
col.tmp=rep(adjustcolor("grey85",.3),length(l))
col.tmp[which(l %in% c("AID"))]<-adjustcolor("black",.8)
col.tmp[which(l %in% c("m564"))]<-adjustcolor("red",.8)
# color of subjects
col.tmp[which(l %in% unique(df$mouse_ID[which(df$condition=="AID")]))]<-adjustcolor("black",.3)
col.tmp[which(l %in% unique(df$mouse_ID[which(df$condition=="m564")]))]<-adjustcolor("red",.3)
# size of nodes
size.tmp=rep(1,length(l))
size.tmp[which(l %in% c("AID","m564"))]<-30
size.tmp[which(l %in% c(names(table(df$mouse_ID))))]<-15
for(i in which(size.tmp==1)){
  size.tmp[i]=gm_mean(df$cdr3.freq[which(df$cdr3==l[i])])/10
}
# color of edges
line.tmp=rep("grey",nrow(toQgraph))
line.tmp[which(c(toQgraph[,1],toQgraph[,2]) %in% c("AID"))]<-"black"
line.tmp[which(c(toQgraph[,1],toQgraph[,2]) %in% c("m564"))]<-"red"
labels.cex.tmp=l
labels.cex.tmp=rep(0.00000000001,length(l))
labels.cex.tmp[which(l %in% c("AID","m564"))]<-1.5
labels.cex.tmp[which(l %in% c(names(table(df$mouse_ID))))]<-1
#graph
dim(toQgraph)
toQgraph$x[toQgraph$x %in% "AID"]<-"WT"
toQgraph$x[toQgraph$x %in% "m564"]<-"564Igi"
png("cdr3.network.by.mouse.png",width=5,height=5,units="in",res=400)
qgraph(toQgraph,color=col.tmp,vsize=size.tmp/2,edge.color=line.tmp,edge.width=0.25,
       labels=TRUE,label.cex=labels.cex.tmp,directed=F,repulsion=20)
dev.off()


###Stacked bar plots and pie charts
##making reference dfs
#clone frequency within condition
cdr3.list<-list()
v_gene.list<-list()
for(i in names(table(clone.data.ab$condition))){
  df<-clone.data.ab[clone.data.ab$condition==i,]
  cdr3.df<-as.data.frame(table(df$cdr3))
  colnames(cdr3.df)<-c("cdr3","cdr3.freq")
  cdr3.df$condition<-i
  cdr3.list[[i]]<-cdr3.df
  v_gene.df<-as.data.frame(table(df$v_gene))
  colnames(v_gene.df)<-c("v_gene","v_gene.freq")
  v_gene.df$condition<-i
  v_gene.df<-left_join(x = v_gene.df %>% group_by(v_gene) %>% mutate(id = row_number()),
                       y = clone.data.ab %>% select(v_gene,chain) %>% group_by(v_gene) %>% mutate(id = row_number()), 
                       by = c("v_gene","id"),keep=F)
  v_gene.list[[i]]<-v_gene.df
}
cdr3.data.bycondition<-bind_rows(cdr3.list)
cdr3.data.bycondition<-cdr3.data.bycondition%>% arrange(cdr3.freq)
v_gene.data.bycondition<-bind_rows(v_gene.list)
v_gene.data.bycondition<-v_gene.data.bycondition%>% arrange(v_gene.freq)
#clone frequency within mice
cdr3.list<-list()
v_gene.list<-list()
for(i in names(table(clone.data.ab$mouse_ID))){
  df<-clone.data.ab[clone.data.ab$mouse_ID==i,]
  cdr3.df<-as.data.frame(table(df$cdr3))
  colnames(cdr3.df)<-c("cdr3","cdr3.freq")
  cdr3.df$mouse_ID<-i
  cdr3.list[[i]]<-cdr3.df
  v_gene.df<-as.data.frame(table(df$v_gene))
  colnames(v_gene.df)<-c("v_gene","v_gene.freq")
  v_gene.df$mouse_ID<-i
  v_gene.df<-left_join(x = v_gene.df %>% group_by(v_gene) %>% mutate(id = row_number()),
                       y = clone.data.ab %>% select(v_gene,chain) %>% group_by(v_gene) %>% mutate(id = row_number()), 
                       by = c("v_gene","id"),keep=F)
  v_gene.list[[i]]<-v_gene.df
}
cdr3.data.bymouse<-bind_rows(cdr3.list)
cdr3.data.bymouse<-cdr3.data.bymouse%>% arrange(cdr3.freq)
cdr3.data.bymouse<-left_join(x = cdr3.data.bymouse, y = metadata, by = "mouse_ID",keep=F)
v_gene.data.bymouse<-bind_rows(v_gene.list)
v_gene.data.bymouse<-v_gene.data.bymouse%>% arrange(v_gene.freq)
v_gene.data.bymouse<-left_join(x = v_gene.data.bymouse, y = metadata, by = "mouse_ID",keep=F)
#clone frequency within condition + cluster
cdr3.list<-list()
v_gene.list<-list()
for(i in names(table(clone.data.ab.seurat$condition))){
  df2<-clone.data.ab.seurat[clone.data.ab.seurat$condition==i,]
  for(j in names(table(df2$my.clusters))){
    df3<-df2[df2$my.clusters==j,]
    cdr3.df<-as.data.frame(table(df3$cdr3))
    colnames(cdr3.df)<-c("cdr3","cdr3.freq")
    cdr3.df$condition<-i
    cdr3.df$my.clusters<-j
    cdr3.list[[paste0(i,".",j)]]<-cdr3.df
    v_gene.df<-as.data.frame(table(df3$v_gene))
    colnames(v_gene.df)<-c("v_gene","v_gene.freq")
    v_gene.df$condition<-i
    v_gene.df$my.clusters<-j
    v_gene.df<-left_join(x = v_gene.df %>% group_by(v_gene) %>% mutate(id = row_number()),
                         y = clone.data.ab %>% select(v_gene,chain) %>% group_by(v_gene) %>% mutate(id = row_number()), 
                         by = c("v_gene","id"),keep=F)
    v_gene.list[[paste0(i,".",j)]]<-v_gene.df
  }
}
cdr3.data.bycluster<-bind_rows(cdr3.list)
cdr3.data.bycluster<-cdr3.data.bycluster%>% arrange(cdr3.freq)
v_gene.data.bycluster<-bind_rows(v_gene.list)
v_gene.data.bycluster<-v_gene.data.bycluster%>% arrange(v_gene.freq)
save(list=c(ls(pattern="data")),file="clone.freq.prop.RData")

##pie charts
#by cluster
df<-cdr3.data.bycluster
png("cdr3.bycluster.pie.png",width=15,height=5.5,units="in",res=400)
# par(mfrow=c(2,length(names(table(df$my.clusters)))),mai=c(1,.1,.25,.1) )
layout(matrix(c(1:12,rep(13,6)),ncol=6,byrow=T),heights=c(4,4,1))
par(mai=c(.1,.1,.1,.1) )
for(i in names(table(df$condition))){
  for(k in names(table(df$my.clusters))){
    df2<-df[df$condition==i&df$my.clusters==k,]
    # df2<-df2[df2$cdr3%in% names(head(sort(table(df2$cdr3),decreasing=T),n=100)),] #select only top 100 clones
    df2$cdr3.prop2<- df2$cdr3.freq/sum(df2$cdr3.freq)
    df2.lbl<-df2[df2$cdr3.prop2>0.02,] #which clones to label
    if(length(df2.lbl$cdr3)!=0){
      df2.lbl$lbl<-paste0(df2.lbl$cdr3, "\n", round(100*df2.lbl$cdr3.prop2, 1),"%")
      df2<-left_join(df2,df2.lbl)
    } else {df2$lbl<-NA}
    # pie(df2$cdr3.freq,labels=df2$lbl,col=sample(brewer.pal(12,"Paired")),border=NA, xlab=j, cex=.4,line=-1)
    df2$color<-"darkgrey"
    for(j in 1:length(df2$cdr3.freq)){
      if(between(df2$cdr3.freq[j],2,4)){df2$color[j]<-"grey30" }
      if(between(df2$cdr3.freq[j],5,9)){df2$color[j]<-"tan1" }
      if(between(df2$cdr3.freq[j],10,19)){df2$color[j]<-"blue4"}
      if(between(df2$cdr3.freq[j],20,49)){df2$color[j]<-"red2"}
      if(df2$cdr3.freq[j]>=50){df2$color[j]<-"forestgreen"}
    }
    # pie(df2$cdr3.freq,labels=df2$lbl,col=sample(brewer.pal(12,"Paired")),border=NA, xlab=j, cex=0.4,line=-1)
    pie(df2$cdr3.freq,col=df2$color,border=NA,labels=NA)#,labels=df2$lbl,cex=0.3)
    # mtext(k,side=1,line=-3,cex=0.75)
    mtext(paste0("n = ",sum(df2$cdr3.freq)),side=1,line=-2,cex=1)
  }
  # mtext(i,side=3,cex=1.5,line=-1)
}
par(mai=c(0,0,0,0))
plot.new()
legend(x="center", ncol=6,legend=c("Unique","\u2265 2","\u2265 5","\u2265 10","\u2265 20","\u2265 50"),
       fill=c("darkgrey","grey30","tan1","blue4","red2","forestgreen"),bty="n",cex=2)
dev.off()



