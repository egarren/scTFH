rm(list=ls())
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
library(DOSE)
library(ggpubr)
library(data.table)
data_concater <- function(x){
  x<- levels(factor(x))
  paste(x, collapse = "+")
}


#Size correlations
load("vdj.analyzed.RData")
cdr3.temp<-T_all.combined@meta.data$cdr3
cdr3.temp2<-AID.combined@meta.data$cdr3
cdr3.temp3<-m564.combined@meta.data$cdr3
load("GLIPH.analyzed.RData")
T_all.combined@meta.data$cdr3<-cdr3.temp
AID.combined@meta.data$cdr3<-cdr3.temp2
m564.combined@meta.data$cdr3<-cdr3.temp3
for(i in c("cdr3")){ #,"index","pattern","cdr3"
  cor.df<-list()
  for(k in c("AID","m564","T_all")){
    obj<-get(paste0(k,".combined"))
    cond.freqs<-as.data.frame(table(obj@meta.data[[i]]))
    rownames(cond.freqs)<-cond.freqs$Var1
    colnames(cond.freqs)<-c(i,paste0(i,".cond.freq"))
    obj@meta.data<-left_join(obj@meta.data,cond.freqs,by=i)
    df<-as.data.frame(t(as.matrix(obj[["RNA"]]@data)))
    mod_score<-as.numeric(obj@meta.data[[paste0(i,".cond.freq")]])
    mat<-cbind(mod_score,df)
    cor.scores<-data.frame(gene=rep(NA,ncol(mat)),pval=rep(NA,ncol(mat)),rho=rep(NA,ncol(mat)))
    for(j in 1:ncol(mat)){
      cor<-cor.test(mat[[j]],mat$mod_score,method="spearman")
      cor.scores$pval[j]<-cor$p.value
      cor.scores$rho[j]<-cor$estimate
      cor.scores$gene[j]<-colnames(mat)[j]
    }
    cor.scores[is.na(cor.scores)]<-0
    assign(paste0(i,".",k,".size.correl"),cor.scores[cor.scores$gene!="mod_score",])
    colnames(cor.scores)<-c("gene",paste0(k,".pval"),paste0(k,".rho"))
    cor.df[[k]]<-cor.scores
  }
  cor.compare<-merge(cor.df[[1]], cor.df[[2]], by="gene")
  assign(paste0(i,".size.cor.compare"),cor.compare[cor.compare$gene!="mod_score",])
  #by cluster
  for(l in unique(T_all.combined@meta.data$my.clusters)){
    Idents(T_all.combined)<-"my.clusters"
    clust<-subset(T_all.combined,idents=l)
    cor.df<-list()
    for(k in c("AID","m564","T_all")){
      Idents(clust)<-"condition"
      if(k %in% c("AID","m564")){obj<-subset(clust,idents=k)}else{obj<-clust}
      cond.freqs<-as.data.frame(table(obj@meta.data[[i]]))
      rownames(cond.freqs)<-cond.freqs$Var1
      colnames(cond.freqs)<-c(i,paste0(i,".cond.freq"))
      obj@meta.data<-left_join(obj@meta.data,cond.freqs,by=i)
      df<-as.data.frame(t(as.matrix(obj[["RNA"]]@data)))
      mod_score<-as.numeric(obj@meta.data[[paste0(i,".cond.freq")]])
      mat<-cbind(mod_score,df)
      cor.scores<-data.frame(gene=rep(NA,ncol(mat)),pval=rep(NA,ncol(mat)),rho=rep(NA,ncol(mat)))
      for(j in 1:ncol(mat)){
        cor<-cor.test(mat[[j]],mat$mod_score,method="spearman")
        cor.scores$pval[j]<-cor$p.value
        cor.scores$rho[j]<-cor$estimate
        cor.scores$gene[j]<-colnames(mat)[j]
      }
      cor.scores[is.na(cor.scores)]<-0
      assign(paste0(i,".clust",l,".",k,".size.correl"),cor.scores[cor.scores$gene!="mod_score",])
      colnames(cor.scores)<-c("gene",paste0(k,".pval"),paste0(k,".rho"))
      cor.df[[k]]<-cor.scores
    }
    cor.compare<-merge(cor.df[[1]], cor.df[[2]], by="gene")
    assign(paste0(i,".clust",l,".size.cor.compare"),cor.compare[cor.compare$gene!="mod_score",])
  }
}

save.image("temp.size.correl.RData")

#within specific GLIPHs
l="GLIPH"
cor.df<-list()
for(k in c("WT-expanded","564Igi-expanded")){
  Idents(T_all.combined)<-"GLIPH.exp"
  obj<-subset(T_all.combined,idents=k)
  cond.freqs<-as.data.frame(table(obj@meta.data[["cdr3"]]))
  rownames(cond.freqs)<-cond.freqs$Var1
  colnames(cond.freqs)<-c("cdr3","cdr3.cond.freq")
  obj@meta.data<-left_join(obj@meta.data,cond.freqs,by="cdr3")
  df<-as.data.frame(t(as.matrix(obj[["RNA"]]@data)))
  mod_score<-as.numeric(obj@meta.data[["cdr3.cond.freq"]])
  mat<-cbind(mod_score,df)
  cor.scores<-data.frame(gene=rep(NA,ncol(mat)),pval=rep(NA,ncol(mat)),rho=rep(NA,ncol(mat)))
  for(j in 1:ncol(mat)){
    cor<-cor.test(mat[[j]],mat$mod_score,method="spearman")
    cor.scores$pval[j]<-cor$p.value
    cor.scores$rho[j]<-cor$estimate
    cor.scores$gene[j]<-colnames(mat)[j]
  }
  cor.scores[is.na(cor.scores)]<-0
  assign(paste0("cdr3.",l,".",k,".size.correl"),cor.scores[cor.scores$gene!="mod_score",])
  colnames(cor.scores)<-c("gene",paste0(k,".pval"),paste0(k,".rho"))
  cor.df[[k]]<-cor.scores
}
cor.compare<-merge(cor.df[[1]], cor.df[[2]], by="gene")
assign(paste0("cdr3.",l,".size.cor.compare"),cor.compare[cor.compare$gene!="mod_score",])

##graphing
dir.create("./size.correl")
setwd("./size.correl")
 
#ranks bar
i="cdr3.T_all.size.correl" #"cdr3.clust4.T_all.size.correl"
for(i in ls(pattern=".size.correl")){
# for(i in ls(pattern="-expanded.size.correl")){
  df<-get(i)
  head(df[order(df$rho),])
  tail(df[order(df$rho),])
  ranks <- df$rho
  names(ranks) <- df$gene
  png(paste0(i,".ranks.png"),width=5,height=4,units="in",res=200)
  barplot(sort(ranks, decreasing = T),cex.names=0.7,ylab=expression(paste("Correlation (",rho,")")))
  dev.off()
}

#correl vs correl scatter
i="cdr3_ID.sel.gliphGLIPH_506.size.cor.compare"
for(i in c(ls(pattern=".size.cor.compare"))){
  df<-get(i)
  df[is.na(df)]<-0
  if(!is(df,"try-error")){
    if("AID.rho" %in% colnames(df)){df<-cbind("gene"=df$gene,df[,grepl("rho",names(df))])}
    if("WT-expanded.rho" %in% colnames(df)){df<-cbind("gene"=df$gene,df[,grepl("rho",names(df))])}
    colnames(df)[apply(sapply(c("AID","log2FC.y","WT"), function (y) sapply(colnames(df), 
                                                                       function (x) grepl(y, x))), 1, any)]<-"AID"
    colnames(df)[apply(sapply(c("m564","log2FC.x","564Igi"), function (y) sapply(colnames(df), 
                                                                        function (x) grepl(y, x))), 1, any)]<-"m564"
    rownames(df)<-df$gene
    df$sig<-"unsig"
    df$sig[abs(df$m564)>0.05&abs(df$AID)>0.05]<-"co.DE"
    df$sig[abs(df$m564)>0.05&abs(df$AID)<=0.05]<-"564.DE"
    df$sig[abs(df$m564)<=0.05&abs(df$AID)>0.05]<-"AID.DE"
    up<-rownames(df[abs(df$m564)>0.05&abs(df$AID)>0.05&!is.na(df$gene),])
    if(length(up)>15){up<-rownames(df[abs(df$m564)>0.05&abs(df$AID)>0.05&!is.na(df$gene),])}
    if(length(up)>15){up<-rownames(df[abs(df$m564)>0.075&abs(df$AID)>0.075&!is.na(df$gene),])}
    if(length(up)>15){up<-rownames(df[abs(df$m564)>0.1&abs(df$AID)>0.1&!is.na(df$gene),])}
    if(length(up)>15){up<-rownames(df[abs(df$m564)>0.11&abs(df$AID)>0.11&!is.na(df$gene),])}
    if(length(up)>15){up<-rownames(df[abs(df$m564)>0.12&abs(df$AID)>0.12&!is.na(df$gene),])}
    if(length(up)>15){up<-rownames(df[abs(df$m564)>0.15&abs(df$AID)>0.15&!is.na(df$gene),])}
    if(length(up)>15){up<-rownames(df[abs(df$m564)>0.2&abs(df$AID)>0.2&!is.na(df$gene),])}
    if(length(up)>15){up<-rownames(df[abs(df$m564)>0.25&abs(df$AID)>0.25&!is.na(df$gene),])}
    if(length(up)>15){up<-rownames(df[abs(df$m564)>0.3&abs(df$AID)>0.3&!is.na(df$gene),])}
    diff<-rownames(df[abs(df$m564-df$AID)>0.1&!is.na(df$gene),])
    # gene.list<-rownames(df[abs(df$m564-df$AID)>0.1&!is.na(df$gene),])
    if(length(diff)>15){diff<-rownames(df[abs(df$m564-df$AID)>0.11&!is.na(df$gene),])}
    if(length(diff)>15){diff<-rownames(df[abs(df$m564-df$AID)>0.12&!is.na(df$gene),])}
    if(length(diff)>15){diff<-rownames(df[abs(df$m564-df$AID)>0.13&!is.na(df$gene),])}
    if(length(diff)>15){diff<-rownames(df[abs(df$m564-df$AID)>0.14&!is.na(df$gene),])}
    if(length(diff)>15){diff<-rownames(df[abs(df$m564-df$AID)>0.2&!is.na(df$gene),])}
    if(length(diff)>15){diff<-rownames(df[abs(df$m564-df$AID)>0.25&!is.na(df$gene),])}
    if(length(diff)>15){diff<-rownames(df[abs(df$m564-df$AID)>0.3&!is.na(df$gene),])}
    if(length(diff)>15){diff<-rownames(df[abs(df$m564-df$AID)>0.4&!is.na(df$gene),])}
    if(length(diff)>15){diff<-rownames(df[abs(df$m564-df$AID)>0.6&!is.na(df$gene),])}
    if(length(diff)>15){diff<-rownames(df[abs(df$m564-df$AID)>0.8&!is.na(df$gene),])}
    if(length(diff)>15){diff<-rownames(df[abs(df$m564-df$AID)>1&!is.na(df$gene),])}
    if(length(diff)>15){diff<-rownames(df[abs(df$m564-df$AID)>1.5&!is.na(df$gene),])}
    gene.list<-c(up,diff)
    df$sig <- factor(df$sig, levels = c("unsig","co.DE","564.DE","AID.DE"))
    p2 <- ggplot(df, aes(m564,AID)) + geom_point(aes(colour=sig),size=0.5)+#geom_point(fill=NA,colour=alpha("black",0.5),pch=21,size=3) + 
      ggtitle(i)+theme_classic()+
      labs(x=expression(paste("564Igi correlation (",rho,")")),y=expression(paste("WT correlation (",rho,")")))+
      theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),legend.position="none",
            axis.ticks=element_blank(),axis.line=element_blank())+ #axis.text=element_blank(),
      geom_hline(yintercept=0,color=alpha("black",0.5))+geom_vline(xintercept=0,color=alpha("black",0.5))+#geom_abline(intercept = 0, slope = 1)+
      scale_color_manual(values = c("grey","black","red","orange"))
    if(length(gene.list)>0 & length(gene.list)<=30){p2 <- LabelPoints(plot = p2, points = gene.list, repel = TRUE,xnudge=0,ynudge=0,size=2.5,segment.size=0.1)}
    p2
    ggsave2(paste0(i,".scatter.png"),width=4, height=4,device="png")
  }
}


setwd("../")
load("vdj.RData")
#Clone-normalization (average expression within all cells of same cdr3)
Idents(T_all.combined) <- "cdr3" #setting idents to condition metadata
avg.T_all.cdr3.raw <- AverageExpression(T_all.combined)$RNA
#create seurat object
avg.T_all.cdr3<-CreateSeuratObject(counts = avg.T_all.cdr3.raw, min.cells = 3, min.features  = 200, 
                                   project = "cdr3", assay = "RNA")
avg.T_all.cdr3<-AddMetaData(avg.T_all.cdr3, metadata=rownames(avg.T_all.cdr3@meta.data),col.name = "cdr3")
#add meta data from raw seurat
df<-T_all.combined@meta.data[!is.na(T_all.combined@meta.data$cdr3),]
colnames(df)<-paste0("raw_",colnames(df))
df<- as.data.table(df)[, lapply(.SD, data_concater), by=raw_cdr3]
df$raw_cdr3.freq <- as.numeric(as.character(df$raw_cdr3.freq))
rownames(df)<-df$raw_cdr3
avg.T_all.cdr3@meta.data<-left_join(x = avg.T_all.cdr3@meta.data, y = df, by = c("cdr3"="raw_cdr3"),keep=F)
rownames(avg.T_all.cdr3@meta.data)<-avg.T_all.cdr3@meta.data$cdr3
table(avg.T_all.cdr3@meta.data$raw_condition2)
Idents(avg.T_all.cdr3)<-"raw_condition2"
p<-list()
for(i in c("Ifngr1","Selplg","Eomes","Sell","Tcf7","Pdcd1","Tnfsf8","Ctsb","S100a6","Ly6a")){ #Nr4a1
  p[[i]]<-FeatureScatter(object = avg.T_all.cdr3, feature1 = "raw_cdr3.freq", feature2 = i,
                 col=alpha(c("black","grey","red"),0.6),smooth=T)+ #,pt.size=0.01
    scale_x_continuous(trans='log10')+scale_y_continuous(trans='log10')+#geom_smooth(method=lm)+ ,limits=c(10,NA)
    xlab("Clone Size")+ylab("Average expression")+ NoLegend()+ggtitle(i)+
    theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
          axis.ticks=element_blank(),axis.line=element_blank(),axis.title=element_blank())
}
cowplot::plot_grid(plotlist=p,ncol=5)
ggsave2("cdr3.size.correl.scatter.png",width=18, height=7,device="png")

load("GLIPH.analyzed.RData")
T_all.combined@meta.data$cdr3<-cdr3.temp
AID.combined@meta.data$cdr3<-cdr3.temp2
m564.combined@meta.data$cdr3<-cdr3.temp3
#Clone-normalization within GLIPH (average expression within all cells of same cdr3)
Idents(T_all.combined)<-"GLIPH.exp"
obj<-subset(T_all.combined,idents=c("WT-expanded","564Igi-expanded"))
cond.freqs<-as.data.frame(table(obj@meta.data[["cdr3"]]))
rownames(cond.freqs)<-cond.freqs$Var1
colnames(cond.freqs)<-c("cdr3","cdr3.cond.freq")
obj@meta.data<-left_join(obj@meta.data,cond.freqs,by="cdr3")
Idents(obj) <- "cdr3" #setting idents to condition metadata
avg.gliph.cdr3.raw <- AverageExpression(obj)$RNA
#create seurat object
avg.gliph.cdr3<-CreateSeuratObject(counts = avg.gliph.cdr3.raw, min.cells = 3, min.features  = 200, 
                                   project = "cdr3", assay = "RNA")
avg.gliph.cdr3<-AddMetaData(avg.gliph.cdr3, metadata=rownames(avg.gliph.cdr3@meta.data),col.name = "cdr3")
#add meta data from raw seurat
df<-obj@meta.data[!is.na(obj@meta.data$cdr3),]
colnames(df)<-paste0("raw_",colnames(df))
df<- as.data.table(df)[, lapply(.SD, data_concater), by=raw_cdr3]
df$raw_cdr3.cond.freq<- as.numeric(as.character(df$raw_cdr3.cond.freq))
rownames(df)<-df$raw_cdr3
avg.gliph.cdr3@meta.data<-left_join(x = avg.gliph.cdr3@meta.data, y = df, by = c("cdr3"="raw_cdr3"),keep=F)
rownames(avg.gliph.cdr3@meta.data)<-avg.gliph.cdr3@meta.data$cdr3
table(avg.gliph.cdr3@meta.data$raw_condition2)
#scatter graph
Idents(avg.gliph.cdr3)<-"raw_condition2"
p<-list()
for(i in c("Malat1","Hcst","Ptprc","Pkm","Ascl2","Lag3")){ #"Thy1","Cd40lg","Ikzf2","Foxp3","Il1r2""Tcf7","Tmsb10","Itgb1","Gm30211","Il1r2","Fabp5"
  p[[i]]<-FeatureScatter(object = avg.gliph.cdr3, feature1 = "raw_cdr3.cond.freq", feature2 = i,
                         col=alpha(c("black","grey","red"),0.6),smooth=T,pt.size=2)+ #,pt.size=0.01
    scale_x_continuous(trans='log10')+scale_y_continuous(trans='log10')+
    xlab("Clone Size")+ylab("Average expression")+ NoLegend()+ggtitle(i)+ 
    theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
          axis.ticks=element_blank(),axis.line=element_blank(),axis.title=element_blank())
}
cowplot::plot_grid(plotlist=p,ncol=3)
ggsave2("cdr3.gliph.exp.size.correl.scatter.png",width=9, height=6,device="png")
FeatureScatter(object = avg.gliph.cdr3, feature1 = "raw_cdr3.cond.freq", feature2 = "Ikzf2",
               col=alpha(c("black","grey","red"),0.6),smooth=T)+ #,pt.size=0.01
  scale_x_continuous(trans='log10')+scale_y_continuous(trans='log10')+
  xlab("Clone Size")+ylab("Average expression")+ NoLegend()+ggtitle("Ikzf2")+ 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.ticks=element_blank(),axis.line=element_blank())
ggsave2("cdr3.gliph.exp.size.correl.Ikzf2.scatter.png",width=5, height=5,device="png")



