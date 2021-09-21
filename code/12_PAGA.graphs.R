rm(list=ls())
library(plyr)
library(Seurat)
library(dplyr)
library(Matrix)
library(knitr)
library(viridis)
library(ggplot2)
library(cowplot)
library(org.Mm.eg.db)
library(tidyverse)
library(RDAVIDWebService)
library(Seurat)
library(cowplot)
library(fgsea)
library(clusterProfiler)
library(pathview)
library(topGO)
library(scde)
library(biomaRt)
library(GO.db)
library(DBI)
library(VISION)
library(msigdbr)
library(KEGGREST)
library(DOSE)
library(ggpubr)
library(cogena)
library(viridisLite)
library(viridis)
library(gridExtra)


load("graphed.RData")

#import scanpy diffmap 
scanpy.diffmap <- read.table("scanpy.diffmap.txt", sep = "\t")
scanpy.diffmap <- scanpy.diffmap[,c(2,3)]
colnames(scanpy.diffmap) <- c("scanpy.DC1", "scanpy.DC2")
rownames(scanpy.diffmap) <- T_all.combined@meta.data$cellID
scanpy.diffmap$scanpy.DC1<-scanpy.diffmap$scanpy.DC1*1000
scanpy.diffmap$scanpy.DC2<-scanpy.diffmap$scanpy.DC2*1000
scanpy.diffmap <- as.matrix(scanpy.diffmap)
T_all.combined[["scanpy.diffmap"]] <- CreateDimReducObject(
  embeddings = scanpy.diffmap, key = "scanpyDC_", assay = DefaultAssay(T_all.combined))
#import scanpy scanpy.pseudotime
scanpy.dpt <- read.table("scanpy.dpt.txt", sep = "\t")
rownames(scanpy.dpt) <- scanpy.dpt[,1]
scanpy.dpt <- scanpy.dpt[,2, drop = FALSE]
colnames(scanpy.dpt) <- "scanpy.pseudotime"
T_all.combined <- AddMetaData(T_all.combined, metadata = scanpy.dpt)
T_all.combined@meta.data$scanpy.pseudo.rank <- rank(T_all.combined@meta.data$scanpy.pseudotime) 
save(T_all.combined,file="scanpy.pseudotime.RData")

#plot
Idents(T_all.combined) <- "my.clusters2"
DimPlot(T_all.combined, reduction= "scanpy.diffmap")+ labs(x="DC_1",y="DC_2")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.ticks=element_blank(),axis.line=element_blank(),axis.text=element_blank()) #
ggsave2("scanpy.diffmap.png",width=6, height=4,device="png")
FeaturePlot(T_all.combined, features= "scanpy.pseudotime", cols= viridis(100, begin = 0), reduction = "scanpy.diffmap")+ 
  labs(x="DC_1",y="DC_2",color="Pseudotime")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),plot.title = element_blank(),
        axis.ticks=element_blank(),axis.line=element_blank(),axis.text=element_blank())
ggsave2("scanpy.diffmap.pseudo.png",width=5.5, height=4,device="png")
FeaturePlot(T_all.combined, features= "scanpy.pseudotime", cols= viridis(100, begin = 0))+labs(color="Pseudotime")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),plot.title = element_blank(),
        axis.ticks=element_blank(),axis.line=element_blank(),axis.text=element_blank())
ggsave2("umap.scanpy.pseudo.png",width=5.5, height=4,device="png")
VlnPlot(T_all.combined, features = "scanpyDC_1",pt.size=0)+ NoLegend()+labs(y="DC_1")+
  theme(axis.title.x=element_blank(),plot.title=element_blank())
ggsave2("cluster.by.scanpyDC1.png",width=4, height=4,device="png")
VlnPlot(T_all.combined, features = "scanpy.pseudotime",pt.size=0)+ NoLegend()+labs(y="Pseudotime")+
  theme(axis.title.x=element_blank(),plot.title=element_blank())
ggsave2("cluster.by.scanpypseudo.png",width=4, height=4,device="png")
VlnPlot(T_all.combined, features = "scanpy.pseudo.rank",pt.size=0)+ NoLegend()+labs(y="Pseudotime Rank")+
  theme(axis.title.x=element_blank(),plot.title=element_blank())
ggsave2("cluster.by.scanpypseudo.rank.png",width=4, height=4,device="png")
FeatureScatter(T_all.combined, feature1 = "scanpyDC_1", feature2 = "PC_1")+labs(color="")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),plot.title = element_blank(),
        axis.ticks=element_blank(),axis.line=element_blank(),axis.text=element_blank())
ggsave2("scanpy.DC.PC.png",width=6, height=4,device="png")


##Pseudotime correlation

dir.create("./pseudo.correl")
setwd("./pseudo.correl")

load("../graphed.RData")
T_all.combined@meta.data$condition2<-"WT"
T_all.combined@meta.data$condition2[T_all.combined@meta.data$condition=="m564"]<-"564Igi"
temp.obj<-T_all.combined

load("../scanpy.pseudotime.RData")
temp.obj@meta.data$scanpy.pseudotime<-T_all.combined@meta.data$scanpy.pseudotime
temp.obj@meta.data$scanpy.pseudo.rank <- rank(T_all.combined@meta.data$scanpy.pseudotime) 

T_all.combined<-temp.obj


#Plots
Idents(T_all.combined)<-"my.clusters2"
for(h in c("scanpy")){ 
  FeaturePlot(T_all.combined, features= paste0(h,".pseudotime"), cols= viridis(100, begin = 0))+labs(color="Pseudotime")+
    theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),plot.title = element_blank(),
          axis.ticks=element_blank(),axis.line=element_blank(),axis.text=element_blank())
  ggsave2(paste0(h,".umap.pseudo.png"),width=5.5, height=4,device="png")
  VlnPlot(T_all.combined, features = paste0(h,".pseudotime"),pt.size=0)+ NoLegend()+labs(y="Pseudotime")+
    theme(axis.title.x=element_blank(),plot.title=element_blank())
  ggsave2(paste0(h,".vln.pseudo.png"),width=4, height=4,device="png")
  VlnPlot(T_all.combined, features = paste0(h,".pseudo.rank"),pt.size=0)+ NoLegend()+labs(y="Pseudotime Rank")+
    theme(axis.title.x=element_blank(),plot.title=element_blank())
  ggsave2(paste0(h,".vln.pseudo.rank.png"),width=4, height=4,device="png")
  FeaturePlot(T_all.combined, features =paste0(h,".pseudotime"), split.by = "condition2",
              cols=viridis(100, begin = 0),order=T,pt.size=0.001,combine=T)
  ggsave2(paste0(h,".umap2.pseudo.png"),width=6.5, height=3,device="png")
  p<-FeaturePlot(T_all.combined, features = paste0(h,".pseudotime"), split.by = "condition2",cols=viridis(100, begin = 0),
                 order=T,pt.size=0.001,combine=F)
  for(j in 1:length(p)) {
    p[[j]] <- p[[j]] + NoLegend()+NoAxes()+
      theme(panel.border = element_rect(colour = "black", size=1),
            plot.title=element_blank(),axis.title.y.right=element_blank(),
            axis.line=element_blank())
  }
  cowplot::plot_grid(p[[1]],p[[2]],ncol=2)
  ggsave2(paste0(h,".umap2.pseudo.clean.png"),width=6, height=2.75,device="png")
  VlnPlot(T_all.combined, features = paste0(h,".pseudotime"), split.by = "condition2",
          group.by = "my.clusters2",cols=c("red","grey"),pt.size = 0,combine=F)
  ggsave2(paste0(h,".vln2.pseudo.png"),width=5, height=4,device="png")
  VlnPlot(T_all.combined, features = paste0(h,".pseudo.rank"), split.by = "condition2",
          group.by = "my.clusters2",cols=c("red","grey"),pt.size = 0,combine=F)
  ggsave2(paste0(h,".vln2.pseudo.rank.png"),width=5, height=4,device="png")
}

#Pseudotime correlations
for(h in c("scanpy.pseudotime")){ #"slingshot.pseudotime",,"destiny.pseudotime","monocle.pseudotime"
  cor.df<-list()
  for(k in c("AID","m564","T_all")){
    Idents(T_all.combined)<-"condition"
    if(k=="T_all"){obj<-T_all.combined}else{obj<-subset(T_all.combined,idents=k)}
    df<-as.data.frame(t(as.matrix(obj[["RNA"]]@data)))
    pseudo_score<-as.numeric(obj@meta.data[[h]])
    mat<-cbind(pseudo_score,df)
    cor.scores<-data.frame(gene=rep(NA,ncol(mat)),pval=rep(NA,ncol(mat)),rho=rep(NA,ncol(mat)))
    for(j in 1:ncol(mat)){
      cor<-cor.test(mat[[j]],mat$pseudo_score,method="spearman")
      cor.scores$pval[j]<-cor$p.value
      cor.scores$rho[j]<-cor$estimate
      cor.scores$gene[j]<-colnames(mat)[j]
    }
    cor.scores[is.na(cor.scores)]<-0
    assign(paste0(h,".",k,".pseudo.correl"),cor.scores[cor.scores$gene!="pseudo_score",])
    colnames(cor.scores)<-c("gene",paste0(k,".pval"),paste0(k,".rho"))
    cor.df[[k]]<-cor.scores
  }
  cor.compare<-merge(cor.df[[1]], cor.df[[2]], by="gene")
  assign(paste0(h,".pseudo.cor.compare"),cor.compare[cor.compare$gene!="pseudo_score",])
  
  #by cluster
  for(l in unique(T_all.combined@meta.data$my.clusters)){
    Idents(T_all.combined)<-"my.clusters"
    clust<-subset(T_all.combined,idents=l)
    cor.df<-list()
    for(k in c("AID","m564","T_all")){
      Idents(clust)<-"condition"
      if(k %in% c("AID","m564")){obj<-subset(clust,idents=k)}else{obj<-clust}
      df<-as.data.frame(t(as.matrix(obj[["RNA"]]@data)))
      pseudo_score<-as.numeric(obj@meta.data[[h]])
      mat<-cbind(pseudo_score,df)
      # ggscatter(mat,x="Cd74",y="Foxp3",add="reg.line",conf.int=T,cor.coef=T,cor.method="spearman")
      # correlations<-apply(matrix,1,function(x){cor(pseudo_score,x)})
      # cors<-apply(mat,2,cor.test,pseudo_score,method="spearman") #non-parametric (pearson for parametric)
      cor.scores<-data.frame(gene=rep(NA,ncol(mat)),pval=rep(NA,ncol(mat)),rho=rep(NA,ncol(mat)))
      for(j in 1:ncol(mat)){
        cor<-cor.test(mat[[j]],mat$pseudo_score,method="spearman")
        cor.scores$pval[j]<-cor$p.value
        cor.scores$rho[j]<-cor$estimate
        cor.scores$gene[j]<-colnames(mat)[j]
      }
      cor.scores[is.na(cor.scores)]<-0
      assign(paste0(h,".clust",l,".",k,".pseudo.correl"),cor.scores[cor.scores$gene!="pseudo_score",])
      colnames(cor.scores)<-c("gene",paste0(k,".pval"),paste0(k,".rho"))
      cor.df[[k]]<-cor.scores
    }
    cor.compare<-merge(cor.df[[1]], cor.df[[2]], by="gene")
    assign(paste0(h,".clust",l,".pseudo.cor.compare"),cor.compare[cor.compare$gene!="pseudo_score",])
  }
}



for(i in c(ls(pattern=".pseudo.cor.compare"))){
  df<-get(i)
  df[is.na(df)]<-0
  if(!is(df,"try-error")){
    if("AID.rho" %in% colnames(df)){df<-cbind("gene"=df$gene,df[,grepl("rho",names(df))])}
    colnames(df)[apply(sapply(c("AID","log2FC.y"), function (y) sapply(colnames(df), 
                                                                       function (x) grepl(y, x))), 1, any)]<-"AID"
    colnames(df)[apply(sapply(c("m564","log2FC.x"), function (y) sapply(colnames(df), 
                                                                        function (x) grepl(y, x))), 1, any)]<-"m564"
    rownames(df)<-df$gene
    df$sig<-"unsig"
    df$sig[abs(df$m564)>0.05&abs(df$AID)>0.05]<-"co.DE"
    df$sig[abs(df$m564)>0.05&abs(df$AID)<=0.05]<-"564.DE"
    df$sig[abs(df$m564)<=0.05&abs(df$AID)>0.05]<-"AID.DE"
    up<-rownames(df[abs(df$m564)>0.05&abs(df$AID)>0.05&!is.na(df$gene),])
    if(length(up)>10){up<-rownames(df[abs(df$m564)>0.05&abs(df$AID)>0.05&!is.na(df$gene),])}
    if(length(up)>10){up<-rownames(df[abs(df$m564)>0.075&abs(df$AID)>0.075&!is.na(df$gene),])}
    if(length(up)>10){up<-rownames(df[abs(df$m564)>0.1&abs(df$AID)>0.1&!is.na(df$gene),])}
    if(length(up)>10){up<-rownames(df[abs(df$m564)>0.11&abs(df$AID)>0.11&!is.na(df$gene),])}
    if(length(up)>10){up<-rownames(df[abs(df$m564)>0.12&abs(df$AID)>0.12&!is.na(df$gene),])}
    if(length(up)>10){up<-rownames(df[abs(df$m564)>0.15&abs(df$AID)>0.15&!is.na(df$gene),])}
    if(length(up)>10){up<-rownames(df[abs(df$m564)>0.17&abs(df$AID)>0.17&!is.na(df$gene),])}
    if(length(up)>10){up<-rownames(df[abs(df$m564)>0.2&abs(df$AID)>0.2&!is.na(df$gene),])}
    if(length(up)>10){up<-rownames(df[abs(df$m564)>0.22&abs(df$AID)>0.22&!is.na(df$gene),])}
    if(length(up)>10){up<-rownames(df[abs(df$m564)>0.25&abs(df$AID)>0.25&!is.na(df$gene),])}
    if(length(up)>10){up<-rownames(df[abs(df$m564)>0.3&abs(df$AID)>0.3&!is.na(df$gene),])}
    diff<-rownames(df[abs(df$m564-df$AID)>0.1&!is.na(df$gene),])
    # gene.list<-rownames(df[abs(df$m564-df$AID)>0.1&!is.na(df$gene),])
    if(length(diff)>10){diff<-rownames(df[abs(df$m564-df$AID)>0.11&!is.na(df$gene),])}
    if(length(diff)>10){diff<-rownames(df[abs(df$m564-df$AID)>0.12&!is.na(df$gene),])}
    if(length(diff)>10){diff<-rownames(df[abs(df$m564-df$AID)>0.13&!is.na(df$gene),])}
    if(length(diff)>10){diff<-rownames(df[abs(df$m564-df$AID)>0.14&!is.na(df$gene),])}
    if(length(diff)>10){diff<-rownames(df[abs(df$m564-df$AID)>0.2&!is.na(df$gene),])}
    if(length(diff)>10){diff<-rownames(df[abs(df$m564-df$AID)>0.3&!is.na(df$gene),])}
    if(length(diff)>10){diff<-rownames(df[abs(df$m564-df$AID)>0.4&!is.na(df$gene),])}
    gene.list<-c(up,diff)
    df$sig <- factor(df$sig, levels = c("unsig","co.DE","564.DE","AID.DE"))
    p2 <- ggplot(df, aes(m564,AID)) + geom_point(aes(colour=sig),size=0.5)+#geom_point(fill=NA,colour=alpha("black",0.5),pch=21,size=3) + 
      ggtitle(i)+theme_classic()+
      labs(x=expression(paste("564Igi correlation (",rho,")")),y=expression(paste("WT correlation (",rho,")")))+
      theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),legend.position="none",
            axis.ticks=element_blank(),axis.line=element_blank())+ #axis.text=element_blank(),
      geom_hline(yintercept=0,color=alpha("black",0.5))+geom_vline(xintercept=0,color=alpha("black",0.5))+#geom_abline(intercept = 0, slope = 1)+
      scale_color_manual(values = c("grey","black","red","orange"))
    if(length(gene.list)>0 & length(gene.list)<30){p2 <- LabelPoints(plot = p2, points = gene.list, repel = TRUE,xnudge=0,ynudge=0,size=2.5,segment.size=0.1)}
    p2
    ggsave2(paste0(i,".scatter.png"),width=4, height=4,device="png")
  }
}


