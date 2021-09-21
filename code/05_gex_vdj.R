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
data_concater3 <- function(x){
  x<- levels(factor(x))
  if(length(x)>1){paste(x[1:2], collapse = "+")}else{paste(x[1])}
}

load("clone.ab.data.RData")
load("graphed.RData")

#paired data
ab.barcode<-Tall.vdj[Tall.vdj$productive=="True"&!is.na(Tall.vdj$productive)&Tall.vdj$chain %in% c("TRA","TRB"),] 
ab.barcode<-as.data.table(ab.barcode)[, lapply(.SD, data_concater3), by=barcode]
ab.barcode<-ab.barcode[ab.barcode$chain=="TRA+TRB",]
ab.barcode<-transform(ab.barcode, cdr3.freq = ave(seq(nrow(ab.barcode)), cdr3, FUN=length))
ab.barcode<-transform(ab.barcode, ms_cdr3.freq = ave(seq(nrow(ab.barcode)), ms_cdr3, FUN=length))
ab.barcode<-transform(ab.barcode, vab.freq = ave(seq(nrow(ab.barcode)), v_gene, FUN=length))
T_all.combined@meta.data$orig.barcode<-rownames(T_all.combined@meta.data)
rownames(ab.barcode)<-ab.barcode$barcode
T_all.combined<-AddMetaData(T_all.combined, metadata=ab.barcode)
T_all.combined@meta.data$condition2<-"WT"
T_all.combined@meta.data$condition2[T_all.combined@meta.data$condition=="m564"]<-"564Igi"
T_all.combined@meta.data$condition2 <- factor(x = T_all.combined@meta.data$condition2, levels = c("WT", "564Igi"))

#key objects
Idents(T_all.combined) <- "condition"
m564.combined<-subset(T_all.combined, idents='m564')
AID.combined<-subset(T_all.combined, idents='AID')
# identifying top clones
m564.top.clones<-names(head(sort(table(m564.combined@meta.data$cdr3), decreasing = T), n=9))
AID.top.clones<-names(head(sort(table(AID.combined@meta.data$cdr3), decreasing = T), n=9))
save.image("vdj.RData")

#mapping expanded clones by cdr3aa
Tall.vdj.clones<-subset(x = T_all.combined, subset = cdr3!="None")
FeaturePlot(Tall.vdj.clones, features= "cdr3.freq",pt.size=0.1, order=T,cols=c("grey","darkgreen"))+ #NoLegend()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text=element_blank(),axis.ticks=element_blank(),axis.line=element_blank(),
        axis.title=element_blank(),plot.title=element_blank(),legend.position="none")
ggsave2("umap.combined.cdr3.freq.png",width=3, height=3,device="png")
p<-FeaturePlot(Tall.vdj.clones, features= "cdr3.freq",split.by = "condition2",pt.size=0.1, 
               order=T,cols=c("grey","darkgreen"),combine=F) #NoLegend()+
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] +NoAxes()+#NoLegend()+
    theme(panel.border = element_rect(colour = "black", size=1),plot.title=element_blank(),
          axis.title.y.right=element_blank(),
          axis.line=element_blank())
}
cowplot::plot_grid(p[[1]],p[[2]],ncol=2)
ggsave2("umap.cdr3.freq.png",width=9, height=3.5,device="png")


##Mapping individual clones within condition
plot.list <- list()
#AID
#mapping top clones
Idents(AID.combined) <- "cdr3"
# plot.list <- list()
for (i in unique(x = names(head(sort(table(AID.combined@meta.data$cdr3), decreasing = T), n=4)))) {
  plot.list[[i]] <- DimPlot(
    object = AID.combined, cols.highlight="forestgreen",
    cells.highlight = Cells(subset(AID.combined, idents=i)),sizes.highlight=0.3) + 
    NoLegend() + NoAxes()+ggtitle(i)+
    theme(plot.title = element_text(size = 5,hjust = 0.5),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5))
}
#564
#mapping top clones
Idents(m564.combined) <- "cdr3"
for (i in unique(x = names(head(sort(table(m564.combined@meta.data$cdr3), decreasing = T), n=4)))) {
  plot.list[[i]] <- DimPlot(
    object = m564.combined, cols.highlight="forestgreen",
    cells.highlight = Cells(subset(m564.combined, idents=i)),
    sizes.highlight=0.3) + NoLegend() + NoAxes()+ggtitle(i)+
    theme(plot.title = element_text(size = 5,hjust = 0.5),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5))
}
CombinePlots(plots = plot.list, ncol = 4)
ggsave2("topclones.umap.png",width=10, height=5,device="png")

#top clone distribution to different clusters, by mouse
mice<-unique(T_all.combined@meta.data$mouse_ID)
freqlist=list()
for (i in mice) {
  df<-T_all.combined@meta.data[T_all.combined@meta.data$mouse_ID==i,]
  clones<-names(head(sort(table(df$cdr3), decreasing = T), n=5))
  df2<-df[df$cdr3 %in% clones, ]
  freq<-as.data.frame(table(df$my.clusters2))
  freq$pct<-100*freq$Freq/sum(freq$Freq)
  freq$mouse<-i
  freq$condition<-unique(df$condition2)
  freqlist[[i]]<-freq
}
topclone.freq.mice<-do.call(rbind,freqlist)
ggbarplot(topclone.freq.mice, x = "Var1", y = "pct",add = c("mean_se", "jitter"),
          color = "condition",position = position_dodge(0.8), legend="right",palette=c("black","red"),xlab=F)+
  stat_compare_means(aes(group = condition),label = "p.format", method="t.test",size=2.5,label.y=60)+
  # stat_compare_means(method="anova",label.y=0)+
  labs(y = "Cluster distribution of cells belonging \n to 5 largest clonotypes per mouse")+
  theme(legend.title = element_blank(),axis.text.x = element_text(angle = 45, hjust = 1))#,axis.title.y=element_text(size=10))
ggsave2("topclone.clust.distr.mice.png",width=5, height=4,device="png")

#######CloneDE
immdata = repLoad("./rep")
df<-clone.data.ab
imm.list<-list()
for(i in unique(df$T.sampleID)){ #creating immunarch list from collapsed data
  df2<-df[df$T.sampleID==i,]
  df2<-transform(df2,Clones=ave(seq(nrow(df2)),cdr3,FUN=length))
  df2<-as.data.table(df2)[, lapply(.SD, data_concater), by=cdr3]
  df2$Clones<-as.numeric(as.character(df2$Clones))
  df2$Proportion<-df2$Clones/sum(df2$Clones)
  imm.list[[i]]<-tibble(
    Clones=df2$Clones, Proportion=df2$Proportion, CDR3.nt=df2$cdr3_nt,CDR3.aa=df2$cdr3,
    V.name=df2$v_gene,D.name=df2$d_gene,J.name=df2$j_gene
  )
}
pr = pubRep(imm.list, "aa", .coding = T, .verbose = F)
pr.AID = pubRepFilter(pr, immdata$meta, c(condition = "AID"))
pr.564 = pubRepFilter(pr, immdata$meta, c(condition = "m564"))
pr.AID[is.na(pr.AID)]<-0
pr.564[is.na(pr.564)]<-0
pr.AID[["avgfreq.AID"]] = rowMeans(public_matrix(pr.AID), na.rm = T)
pr.564[["avgfreq.564"]] = rowMeans(public_matrix(pr.564), na.rm = T)
pr.AID[["sum.AID"]] = rowSums(public_matrix(pr.AID)[,1:5], na.rm = T)+1
pr.564[["sum.564"]] = rowSums(public_matrix(pr.564)[,1:5], na.rm = T)+1
pr.res.cdr3 = dplyr::full_join(pr.AID, pr.564, by = "CDR3.aa")
pr.res.cdr3.graph<-pr.res.cdr3
pr.res.cdr3.graph[is.na(pr.res.cdr3.graph)]<-0
pr.res.cdr3.graph$sum.AID[pr.res.cdr3.graph$sum.AID==0]<-1
pr.res.cdr3.graph$sum.564[pr.res.cdr3.graph$sum.564==0]<-1
pr.res.cdr3.graph[["Samples.sum"]] = pr.res.cdr3.graph[["Samples.x"]] + pr.res.cdr3.graph[["Samples.y"]]
pr.res.cdr3.graph[["freq.ratio"]] = apply(pr.res.cdr3.graph[, c("sum.AID", "sum.564")],1, function(x) log10(x[1])/log10(x[2]))
pr.res.cdr3.graph[["log2FC"]]= apply(pr.res.cdr3.graph[, c("sum.AID", "sum.564")],1, function(x) log2((x[2])/(x[1])))
rownames(pr.res.cdr3.graph)<-pr.res.cdr3.graph$CDR3.aa
###expanded clone 564 vs WT scatter
AID.public.expand.cdr3<-pr.res.cdr3.graph[pr.res.cdr3.graph$log2FC< -3 &pr.res.cdr3.graph$sum.AID>10,]#[(pr.res.cdr3$avgfreq.AID-10)>10*pr.res.cdr3$avgfreq.564,]
m564.public.expand.cdr3<-pr.res.cdr3.graph[pr.res.cdr3.graph$log2FC>3&pr.res.cdr3.graph$sum.564>10,]#pr.res.cdr3[(pr.res.cdr3$avgfreq.564-10)>10*pr.res.cdr3$avgfreq.AID,]
pr.res.cdr3.graph$annotate<-ifelse(pr.res.cdr3.graph$CDR3.aa %in% AID.public.expand.cdr3$CDR3.aa, "WT expanded","Public")
pr.res.cdr3.graph$annotate[pr.res.cdr3.graph$CDR3.aa %in% m564.public.expand.cdr3$CDR3.aa]<-"564Igi expanded"
labels<-unique(pr.res.cdr3.graph[abs(pr.res.cdr3.graph$log2FC)>3&pr.res.cdr3.graph$Samples.sum>1&
                           (pr.res.cdr3.graph$sum.AID>50|pr.res.cdr3.graph$sum.564>50),])
pr.res.cdr3.graph$annotate <- factor(pr.res.cdr3.graph$annotate, levels = c("Public","WT expanded","564Igi expanded"))
ggplot(pr.res.cdr3.graph,aes(x = sum.564, y =  sum.AID,size=Samples.sum))+theme_classic()+
  scale_x_continuous(trans='log10') +scale_y_continuous(trans='log10')+ scale_size(range = c(1, 10))+
  geom_point(aes(colour=factor(annotate),fill = factor(annotate)), shape=21) + 
  scale_color_manual(values = c(alpha("black",0.5),"black","red"))+
  scale_fill_manual(values = c("#1C00ff00",alpha("black",0.2),alpha("red",0.2)))+
  guides(fill = guide_legend(override.aes = list(size = 7)),color=F)+ #,size=F
  labs(x = "564Igi", y = "WT", size="Samples",fill="Expansion")+
  theme(legend.direction = "vertical", legend.box = "horizontal")#+#theme(legend.title = element_blank())+
ggsave2("overlap.scatter.TCRab.expansion.png",width=5.5, height=3,device="png")

# DE b/w condition-specific expanded clones
AID.public.expand.cdr3.barcode<-AID.combined@meta.data$cellID[AID.combined@meta.data$cdr3 %in% AID.public.expand.cdr3$CDR3.aa]
m564.public.expand.cdr3.barcode<-m564.combined@meta.data$cellID[m564.combined@meta.data$cdr3 %in% m564.public.expand.cdr3$CDR3.aa]
DefaultAssay(T_all.combined) <- "RNA"
Idents(T_all.combined) <- "cellID"
expandedDE.public <- FindMarkers(T_all.combined, ident.1 = m564.public.expand.cdr3.barcode, ident.2 = AID.public.expand.cdr3.barcode,
                        min.pct=0,logfc.threshold = -Inf,test.use = "MAST")
expandedDE.public$log2FC<-log2(exp(expandedDE.public$avg_logFC))
T_all.combined@meta.data$public.expand<-ifelse(T_all.combined@meta.data$cellID %in% AID.public.expand.cdr3.barcode, "AID expanded","")
T_all.combined@meta.data$public.expand[T_all.combined@meta.data$cellID %in% m564.public.expand.cdr3.barcode]<-"m564 expanded"
Idents(T_all.combined) <- "public.expand"
avg.public.expandedDE <- log1p(AverageExpression(T_all.combined)$RNA)
avg.public.expandedDE$gene<-rownames(T_all.combined)

# DE b/w conditions in individual clones
public.clones<-unique(pr.res.cdr3.graph[(pr.res.cdr3.graph$sum.AID>10|pr.res.cdr3.graph$sum.564>10)&
                                          pr.res.cdr3.graph$sum.AID>1&pr.res.cdr3.graph$sum.564>1&pr.res.cdr3.graph$Samples.sum>1,])
DefaultAssay(T_all.combined) <- "RNA"
Idents(T_all.combined) <- "cellID"
for(i in public.clones$CDR3.aa){
  barcodes<-T_all.combined@meta.data$cellID[T_all.combined@meta.data$cdr3 %in% i]
  df<-subset(T_all.combined, idents=barcodes)
  Idents(df) <- "condition" #setting idents to condition metadata
  df.autoimmune.response <- try(FindMarkers(df, ident.1 = "m564", ident.2 = "AID",min.pct=0,logfc.threshold = -Inf,test.use = "MAST"))
  if(!is(df.autoimmune.response,"try-error")){
    df.autoimmune.response$log2FC<-log2(exp(df.autoimmune.response$avg_logFC))
    assign(paste0(i,".cdr3.autoimmune.response"),df.autoimmune.response)
    avg.df <- log1p(AverageExpression(df)$RNA)
    avg.df$gene <- rownames(df)
    assign(paste0(i,".cdr3.avg.exp"),avg.df)
  }
}

# identifying expanded clones
m564.expanded <- names(table(m564.combined@meta.data$cdr3))[table(m564.combined@meta.data$cdr3) > 20]
m564.unexpanded <- names(table(m564.combined@meta.data$cdr3))[table(m564.combined@meta.data$cdr3) ==1]
AID.expanded <- names(table(AID.combined@meta.data$cdr3))[table(AID.combined@meta.data$cdr3) > 20]
AID.unexpanded <- names(table(AID.combined@meta.data$cdr3))[table(AID.combined@meta.data$cdr3) ==1]

#expanded DE (by condition)
Idents(m564.combined) <- "cdr3"
expandedDE.m564<-FindMarkers(m564.combined, ident.1 = m564.expanded, ident.2 = m564.unexpanded, 
                             min.pct=0,logfc.threshold = -Inf,test.use = "MAST")
expandedDE.m564$gene<-rownames(expandedDE.m564)
expandedDE.m564$log2FC<-log2(exp(expandedDE.m564$avg_logFC))
Idents(AID.combined) <- "cdr3"
expandedDE.AID<-FindMarkers(AID.combined, ident.1 = AID.expanded, ident.2 = AID.unexpanded, 
                            min.pct=0,logfc.threshold = -Inf,test.use = "MAST")
expandedDE.AID$gene<-rownames(expandedDE.AID)
expandedDE.AID$log2FC<-log2(exp(expandedDE.AID$avg_logFC))
exp.vs.expandedDE <-merge(expandedDE.m564[,c("gene","log2FC")], expandedDE.AID[,c("gene","log2FC")], by="gene")
rownames(exp.vs.expandedDE)<-exp.vs.expandedDE$gene


#by clusters
for (i in names(table(T_all.combined$my.clusters))){
  #cluster subset
  Idents(m564.combined) <- "my.clusters"
  m564.clust.combined<-subset(m564.combined, idents=i)
  Idents(AID.combined) <- "my.clusters"
  AID.clust.combined<-subset(AID.combined, idents=i)
  # identifying expanded clones
  if(i %in% c("0","1","2")){
    m564.clust.expanded <- names(table(m564.clust.combined@meta.data$cdr3))[table(m564.clust.combined@meta.data$cdr3) > 5]
    AID.clust.expanded <- names(table(AID.clust.combined@meta.data$cdr3))[table(AID.clust.combined@meta.data$cdr3) > 5]
  }else{
    m564.clust.expanded <- names(table(m564.clust.combined@meta.data$cdr3))[table(m564.clust.combined@meta.data$cdr3) > 1]
    AID.clust.expanded <- names(table(AID.clust.combined@meta.data$cdr3))[table(AID.clust.combined@meta.data$cdr3) > 1]
  }
  m564.clust.unexpanded <- names(table(m564.clust.combined@meta.data$cdr3))[table(m564.clust.combined@meta.data$cdr3) ==1]
  AID.clust.unexpanded <- names(table(AID.clust.combined@meta.data$cdr3))[table(AID.clust.combined@meta.data$cdr3) ==1]
  #expanded DE (by condition)
  Idents(m564.clust.combined) <- "cdr3"
  Idents(AID.clust.combined) <- "cdr3"
  expandedDE.m564.clust<-try(FindMarkers(m564.clust.combined, ident.1 = m564.clust.expanded, ident.2 = m564.clust.unexpanded, 
                               min.pct=0,logfc.threshold = -Inf,test.use = "MAST"), silent=TRUE)
  expandedDE.AID.clust<-try(FindMarkers(AID.clust.combined, ident.1 = AID.clust.expanded, ident.2 = AID.clust.unexpanded, 
                                    min.pct=0,logfc.threshold = -Inf,test.use = "MAST"), silent=TRUE)
  if(!is(expandedDE.m564.clust,"try-error")&!is(expandedDE.AID.clust,"try-error")){
    expandedDE.m564.clust$gene<-rownames(expandedDE.m564.clust)
    expandedDE.m564.clust$log2FC<-log2(exp(expandedDE.m564.clust$avg_logFC))
    expandedDE.AID.clust$gene<-rownames(expandedDE.AID.clust)
    expandedDE.AID.clust$log2FC<-log2(exp(expandedDE.AID.clust$avg_logFC))
    clust.exp.vs.exp <-merge(expandedDE.m564.clust[,c("gene","log2FC")], expandedDE.AID.clust[,c("gene","log2FC")], by="gene")
    rownames(clust.exp.vs.exp)<-clust.exp.vs.exp$gene
    assign(paste0("clust",i,".expandedDE.m564"),expandedDE.m564.clust)
    assign(paste0("clust",i,".expandedDE.AID"),expandedDE.AID.clust)
    assign(paste0("clust",i,".exp.vs.expandedDE"),clust.exp.vs.exp)
  }
  #expanded DE (by mouse)
  Idents(T_all.combined) <- "my.clusters"
  clust.combined<-subset(T_all.combined, idents=i)
  Idents(clust.combined) <- "mouse_ID"
  mice<-names(table(clust.combined$mouse_ID))
  clust.expandedDE.ms<-data.frame(gene=rownames(clust.combined))
  df.list <- list()
  for (j in mice){
    ms.combined<-subset(clust.combined, idents=j)
    if(i %in% c("0","1","2")){
      ms.expanded <- names(table(ms.combined@meta.data$cdr3))[table(ms.combined@meta.data$cdr3) > 5]
    }
    else{
      ms.expanded <- names(table(ms.combined@meta.data$cdr3))[table(ms.combined@meta.data$cdr3) > 1]
    }
    ms.unexpanded <- names(table(ms.combined@meta.data$cdr3))[table(ms.combined@meta.data$cdr3) ==1]
    Idents(ms.combined) <- "cdr3"
    df<-try(FindMarkers(ms.combined, ident.1 = ms.expanded, ident.2 = ms.unexpanded, 
                        min.pct=0,logfc.threshold = -Inf,test.use = "MAST"), silent=TRUE)
    if(!is(df,"try-error")){
      df$gene<-rownames(df)
      df$temp<-log2(exp(df$avg_logFC))
      clust.expandedDE.ms <-merge(clust.expandedDE.ms, df[,c("gene","temp")], by="gene")
      names(clust.expandedDE.ms)[names(clust.expandedDE.ms) == 'temp'] <- j
      df.list[[j]]<-df
    }else{
      clust.expandedDE.ms$temp<-rep(NA,length(clust.expandedDE.ms$gene))
      names(clust.expandedDE.ms)[names(clust.expandedDE.ms) == 'temp'] <- j
    }
  }
  clust.expandedDE.ms$AID.avgFC<-rowMeans(clust.expandedDE.ms[,AID.mouseID], na.rm=T)
  clust.expandedDE.ms$m564.avgFC<-rowMeans(clust.expandedDE.ms[,m564.mouseID], na.rm=T)
  clust.expandedDE.ms$delta<-clust.expandedDE.ms$m564.avgFC-clust.expandedDE.ms$AID.avgFC
  for(k in 1:nrow(clust.expandedDE.ms)) {
    clust.expandedDE.ms$ttest[k] <- my.ttest(clust.expandedDE.ms[k,AID.mouseID],clust.expandedDE.ms[k,m564.mouseID])
  }
  rownames(clust.expandedDE.ms)<-clust.expandedDE.ms$gene
  clust.expandedDE.ms$p_val_adj<-clust.expandedDE.ms$ttest
  clust.expandedDE.ms$log2FC<-clust.expandedDE.ms$delta
  assign(paste0("clust",i,".expandedDE.ms"),clust.expandedDE.ms)
  save.image(paste0("temp.clust.",i,".RData"))
}

save.image("vdj.analyzed.RData")
save(list=c(ls(pattern="DE"),ls(pattern=".markers"),ls(pattern=".response"),"clusters"),file="DE.dfs.RData")

#Volcano Plots
for(i in c(ls(pattern=".cdr3.autoimmune.response"))){
  df<-get(i)
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
      if(length(gene.list)>15){gene.list<-rownames(df[abs(df$log2FC)>0.2&df$p_val_adj<0.001,])}
      if(length(gene.list)>15){gene.list<-rownames(df[abs(df$log2FC)>0.25&df$p_val_adj<10e-10,])}
      if(length(gene.list)>15){gene.list<-rownames(df[abs(df$log2FC)>0.35&df$p_val_adj<10e-15,])}
      if(length(gene.list)>15){gene.list<-rownames(df[abs(df$log2FC)>0.5&df$p_val_adj<10e-15,])}
      EnhancedVolcano(df,lab = rownames(df),
                      x = 'log2FC',y = 'p_val_adj',title = i,col=c("black","black","black","red3"),
                      selectLab=gene.list,xlab=bquote(~Log[2]~ (frac("564Igi","WT"))),
                      pCutoff = 0.01,FCcutoff = 0.2,pointSize = 0.5,labSize = 3,axisLabSize=10,colAlpha = 1, #transparencyxlim = c(-1.5, 1.5),
                      legendVisible=F,drawConnectors = TRUE,widthConnectors = 0.2,colConnectors = 'grey30',
                      subtitle="", caption="",border="full",cutoffLineWidth=0,
                      gridlines.major=F,gridlines.minor=F,titleLabSize=10
      )
      ggsave2(paste0(i,".volcano.png"),width=3.5, height=4,device="png")
    }
  }
}

i="expandedDE.public"
df<-get(i)
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
    if(length(gene.list)>15){gene.list<-rownames(df[abs(df$log2FC)>0.2&df$p_val_adj<0.001,])}
    if(length(gene.list)>15){gene.list<-rownames(df[abs(df$log2FC)>0.25&df$p_val_adj<10e-10,])}
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

