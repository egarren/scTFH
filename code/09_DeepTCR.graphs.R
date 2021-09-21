rm(list=ls())
library(pheatmap)
library(reticulate)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(cowplot)
library(data.table)
library(RColorBrewer)
library(Seurat)


meta.list<-list()
for(i in c("AID","m564")){
  file.names<-list.files(path=paste0("../scTfh.data/",i))
  if(i=="AID"){cond<-"WT"}else{cond<-"564Igi"}
  meta.list[[i]]<-data.frame(condition2=cond,file.name=file.names)
}
meta<-do.call("rbind",meta.list)

#ML feature heatmap
df<-read.csv("unsupervised.rep.features.csv",header=T,row.names=1)
deep.feats<-t(df)
pheatmap(deep.feats)
annot.col<-data.frame(file.name=colnames(deep.feats))
annot.col<-merge(annot.col,meta,by="file.name")
names(annot.col)[names(annot.col) == "condition2"] <- "BMChimera"
row.names(annot.col) <- annot.col$file.name
annot.col$file.name <- NULL
p<-pheatmap(deep.feats,annotation_col=annot.col,fontsize_col=3,annotation_names_col=F,show_colnames=F,show_rownames=F,
         annotation_colors=list(BMChimera=c("WT"="grey","564Igi"="red")),main="DeepTCR")
save_pheatmap_png <- function(x, filename, width=700, height=700, res = 150) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_png(p,"unsupervised.rep.sample.heatmap.png")

#KNN scores
df<-read.csv("unsupervised.rep.KNN.csv",header=T,row.names=1)
df<-df[df$Metric %in% c("AUC","Recall"),]
df$condition2<-"WT"
df$condition2[df$Classes=="m564"]<-"564Igi"
df$condition2<-factor(df$condition2,levels=c("WT","564Igi"))
ggplot(df,aes(x=condition2,y=Value,fill=condition2))+geom_violin(trim=F)+theme_classic()+
  scale_fill_manual(values=c("grey","red"))+stat_summary(fun.data=mean_sdl, mult=1,geom="pointrange", color="black")+
  theme(legend.position="none",axis.title.x=element_blank())+labs(y="AUC")+ylim(0,1)
ggsave2("unsupervisred.rep.KNN.violin.png",width=2,height=3,device="png")
ggbarplot(df, x = "Metric", y = "Value", color="condition2", add = c("mean_se","jitter"),add.params = list(size=0.5),
          position = position_dodge(0.8), legend="right",palette=c("black","red"),xlab=F,size=0.5)+
  stat_compare_means(aes(group = condition2),label = "p.format", method="t.test",label.y=1.5,size=4)+
  theme(legend.title=element_blank())
ggsave2("unsupervisred.rep.KNN.bar.png",width=4,height=4,device="png")

#Structural Diversity
df<-read.csv("unsupervised.rep.structural.diversity.csv",header=T,row.names=1)
df$condition2<-"WT"
df$condition2[df$Class=="m564"]<-"564Igi"
df$condition2<-factor(df$condition2,levels=c("WT","564Igi"))
ggplot(df,aes(x=condition2,y=Entropy,fill=condition2))+geom_violin()+theme_classic()+ 
  stat_summary(fun.data="mean_sdl", mult=1, geom="pointrange",color="black")+
  scale_fill_manual(values=c("grey","red"))#+geom_boxplot(width=0.1,fill="white",position=position_dodge(1))
ggbarplot(df, x = "condition2", y = "Entropy", color="condition2", add = c("mean_se","jitter"),add.params = list(size=0.5),
          position = position_dodge(0.8), legend="none",palette=c("black","red"),xlab=F,size=0.5)+
  stat_compare_means(aes(group = condition2),label = "p.format", method="t.test",label.y=2.5,size=4)+
  theme(legend.title=element_blank())
ggsave2("unsupervisred.rep.entropy.png",width=2.5,height=4,device="png")

pd <- import("pandas")
df<- pd$read_pickle("./unsupervised.rep/KNN_sample.pkl")

##UMAP plots
#create metadata info
#add metadata
load("../../clone.data.RData")
data_concater <- function(x){
  x<- levels(factor(x))
  paste(x, collapse = "+")
}
df.TRA<-clone.data[clone.data$chain=="TRA",]
df.TRB<-clone.data[clone.data$chain=="TRB",]
cells<-inner_join(df.TRB[,c("barcode","cdr3","v_gene","j_gene","d_gene")],
                  df.TRA[,c("barcode","cdr3","v_gene","j_gene")],by="barcode")
cells<-left_join(cells,unique(clone.data[,c("barcode","mouse_ID","condition","my.clusters")]),by="barcode")
cells$clone_ID<-paste0(cells$cdr3.x,"+",cells$cdr3.y,"+",cells$mouse_ID)
cells.cloneID<-transform(cells,cloneID.freq=ave(seq(nrow(cells)),clone_ID,FUN=length))
cells.cloneID<-as.data.table(cells.cloneID)[, lapply(.SD, data_concater), by=clone_ID]
cells.cloneID$cloneID.freq<-as.numeric(as.character(cells.cloneID$cloneID.freq))
cells.cloneID$top.clust<-NA
for(i in 1:nrow(cells.cloneID)){
  ID<-cells.cloneID$clone_ID[i]
  temp.cells<-cells[cells$clone_ID==ID,]
  cells.cloneID$top.clust[i]<-names(sort(table(temp.cells$my.clusters), decreasing = T))[1]
}
cells.cloneID$top.clust[cells.cloneID$my.clusters==""]<-NA
#load unsupervised umap from DeepTCR
df<-read.csv("unsup.rep.umap.csv",header=F)
colnames(df)<-c("x","y","condition","mouse","cdr3b","cdr3a","freq","count")
df$mouse_ID<-gsub("\\..*","",df$mouse)
df$clone_ID<-paste0(df$cdr3b,"+",df$cdr3a,"+",df$mouse_ID)
df$condition2<-"WT"
df$condition2[df$condition=="m564"]<-"564Igi"
df$condition2 <- factor(x = df$condition2, levels = c("WT", "564Igi"))
df<-left_join(df,unique(cells.cloneID[,c("clone_ID","cloneID.freq","top.clust")]),by="clone_ID")
for(i in c("condition2","mouse","top.clust")){
  p<-ggplot(df,aes(x, y,size=cloneID.freq,color=.data[[i]]))+theme_classic()+
    geom_point(alpha=0.05)+scale_size(trans="log10",range=c(1,10))+#+
    guides(color = guide_legend(override.aes = list(size = 7,alpha=0.5)))+
    labs(x="UMAP_1",y="UMAP_2",size="Clone Size",color="")+
    theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
          axis.text=element_blank(),axis.ticks=element_blank(),axis.line=element_blank())#axis.text=element_blank(),#,color=F)+ #,size=F, legend.position="none",
  if(i=="condition2"){p+scale_color_manual(values=alpha(c("black","red"),0.05))}#+
     # scale_color_discrete(labels=c("WT","564Igi"))}
  if(i=="top.clust"){p+scale_color_discrete(labels=c("Tfr","Tfh-activated","Sostdc1","Gm42418","Tfh-CM", "Tfh-effector","Tfh-ISG"))}
  ggsave2(paste0("unsup.rep.umap.",i,".png"),width=6.5,height=5,device="png")
}

#load supervised AUC from DeepTCR
AUC<-read.csv("supervised.seq.AUC.csv",header=T)
hist(AUC$AUC)
mean(AUC$AUC)
#loading training UMAP
df2<-read.csv("sup.seq.train.csv",header=F)
colnames(df2)<-c("x","y","condition","mouse","cdr3b")
df2$ML<-"train"
df2$epitope<-df2$condition
df2<-left_join(df2,AUC,by=c("epitope"="Class"))
df2<-transform(df2,epitope.freq=ave(seq(nrow(df2)),epitope,FUN=length))
load("masterdb4.RData")
df2.meta<-left_join(df2,unique(masterdb.epitope[,c("epitope","disease","disease2",
                                                   "antigen","chain","Species","Tcell")]),by="epitope")
df2.select<-df2.meta[df2.meta$AUC>0.9,]
#Training UMAP
set.seed(42)
for(i in c("epitope","disease","disease2","antigen","chain","Species","Tcell")){
  label<-names(sort(table(df2.select[[i]]), decreasing = T))[1:10]
  ggplot(df2.select[sample(nrow(df2.select)),],aes(x, y,color=.data[[i]]))+theme_classic()+
    geom_point(size=0.1)+scale_color_manual(values=colorRampPalette(brewer.pal(8, "Paired"))(
      length(unique(df2.select[[i]]))),breaks=label)+
    guides(color = guide_legend(override.aes = list(size=5)))+ #size = 7,
    labs(x="UMAP_1",y="UMAP_2",color="")+
    theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
          axis.text=element_blank(),axis.ticks=element_blank(),axis.line=element_blank())
  ggsave2(paste0("sup.seq.train.umap.",i,".png"),width=6.5,height=4,device="png")
}
rownames(df2.select)<-make.unique(df2.select$antigen)
p<-ggplot(df2.select[sample(nrow(df2.select)),],aes(x, y,color=epitope))+theme_classic()+
  geom_point(size=0.1)+scale_color_manual(values=colorRampPalette(brewer.pal(8, "Paired"))(
    length(unique(df2.select$epitope))),breaks=label)
HoverLocator(p)
#load supervised umap from DeepTCR
df<-read.csv("sup.seq.test.csv",header=F)
colnames(df)<-c("x","y","condition","mouse","cdr3b")
df$ML<-"test"
df$mouse_ID<-gsub("\\..*","",df$mouse)
df$clone_ID3<-df$cdr3b
#Test UMAP w/ scTfh annotation
cells$clone_ID3<-cells$cdr3.x
cells.cloneID3<-transform(cells,cloneID3.freq=ave(seq(nrow(cells)),clone_ID3,FUN=length))
cells.cloneID3<-as.data.table(cells.cloneID3)[, lapply(.SD, data_concater), by=clone_ID3]
cells.cloneID3$cloneID3.freq<-as.numeric(as.character(cells.cloneID3$cloneID3.freq))
cells.cloneID3$top.clust<-NA
for(i in 1:nrow(cells.cloneID3)){
  ID3<-cells.cloneID3$clone_ID3[i]
  temp.cells<-cells[cells$clone_ID3==ID3,]
  cells.cloneID3$top.clust[i]<-names(sort(table(temp.cells$my.clusters), decreasing = T))[1]
}
cells.cloneID3$top.clust[cells.cloneID3$my.clusters==""]<-NA
df<-left_join(df,unique(cells.cloneID3[,c("clone_ID3","cloneID3.freq","top.clust","condition")]),by="clone_ID3")
df$condition3<-"Public"
df$condition3[df$condition.y=="AID"]<-"WT"
df$condition3[df$condition.y=="m564"]<-"564Igi"
#Training + Testing UMAP
df3<-bind_rows(df,df2.select)
df3$condition3[is.na(df3$condition3)]<-"Training"
df3$condition3<-factor(df3$condition3,levels=c("Training","Public","WT","564Igi"))
ggplot(df3[sample(nrow(df3)),],aes(x, y,color=condition3))+theme_classic()+
  geom_point(size=0.1)+scale_color_manual(values=c("grey","forestgreen","black","red"))+
  guides(color = guide_legend(override.aes = list(size=5)))+ #size = 7,
  labs(x="UMAP_1",y="UMAP_2",color="")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text=element_blank(),axis.ticks=element_blank(),axis.line=element_blank())
ggsave2("sup.seq.test.umap.png",width=5,height=4,device="png")


#ML feature heatmap
df<-read.csv("supervised.seq.features.csv",header=T,row.names=1)
rownames(df)<-gsub("\\..*","",rownames(df))
select.epitopes<-unique(df2.select$epitope)
df.select<-df[select.epitopes,]
pheatmap(df.select)
load("masterdb4.RData")
masterdb2<-masterdb.epitope
annot.col<-data.frame(epitope=rownames(df.select))
annot.col<-merge(annot.col,unique(masterdb2[,c("epitope","disease","disease2")]),by="epitope")
annot.col<-annot.col[!duplicated(annot.col$epitope),]
names(annot.col)[names(annot.col) == "disease2"] <- "Class"
names(annot.col)[names(annot.col) == "disease"] <- "Disease"
row.names(annot.col) <- annot.col$epitope
annot.col$epitope <- NULL
disease.colors<-colorRampPalette(brewer.pal(8, "Paired"))(length(unique(annot.col$Disease)))
names(disease.colors)<-unique(annot.col$Disease)
disease2.colors<-colorRampPalette(brewer.pal(8, "Set1"))(length(unique(annot.col$Class)))
names(disease2.colors)<-unique(annot.col$Class)
p<-pheatmap(df.select,annotation_row=annot.col,annotation_names_row=F,show_rownames=T,show_colnames=F,
            annotation_colors=list(Disease=disease.colors,Class=disease2.colors))#,main="DeepTCR"
save_pheatmap_png <- function(x, filename, width=2500, height=700, res = 150) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_png(p,"supervised.epitopes.heatmap.png")


