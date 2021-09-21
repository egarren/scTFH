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
library(ggalluvial)
my.ttest <- function(...) {
  obj<-try(t.test(...), silent=TRUE)
  if (is(obj, "try-error")) return(NA) else return(obj$p.value)
}


#load data
metadata<-read.csv("Tcell.metadata.csv", header=T) 
metadata[] <- lapply(metadata, as.character)
mice<-metadata$T.sampleID
for (i in mice){
  path<-paste0("./gex/",i,"/outs/filtered_feature_bc_matrix") #path to 10X output
  seurat.object<-CreateSeuratObject(counts = Read10X(data.dir = path,gene.column=1) , 
                                    min.cells = 3, min.features  = 200, project = i, assay = "RNA")
  assign(i,seurat.object)
}

##Add Count Matrices (eg. all 564 samples, all AID samples)
seurat.list<-lapply(mice,get)
T_all<-merge(seurat.list[[1]], y=seurat.list[2:length(seurat.list)],add.cell.ids=mice,project="564vsAID")
# add sample metadata
T_all<-AddMetaData(T_all, metadata=rownames(T_all@meta.data),col.name = "cellID")
T_all@meta.data<-left_join(x = T_all@meta.data, y = metadata, by = c("orig.ident"="T.sampleID"),keep=F)
rownames(T_all@meta.data)<-T_all@meta.data$cellID
rm(list=setdiff(ls(), c("T_all", "metadata")))
save.image("temp.merged.RData")

##QC and Filter
#adding feature (gene) metadata
T_all@assays[["RNA"]]@meta.features$original_ensembl<-rownames(T_all@assays[["RNA"]]@meta.features)
genes.meta<-getBM(attributes=c("ensembl_gene_id", "mgi_symbol", "start_position", "end_position", "chromosome_name", 
                      "percentage_gene_gc_content", "external_gene_name", "gene_biotype","go_id","name_1006"),filters=
                    "ensembl_gene_id",values=list(rownames(T_all@assays[["RNA"]]@meta.features)),
                  mart=useMart("ensembl", dataset = "mmusculus_gene_ensembl",host="ensembl.org"),useCache=F) #useast.
m <- match(T_all@assays[["RNA"]]@meta.features$original_ensembl, genes.meta$ensembl_gene_id)
T_all@assays[["RNA"]]@meta.features<-cbind(T_all@assays[["RNA"]]@meta.features,genes.meta[m,])
#renaming features to gene symbols
rownames(T_all@assays[["RNA"]]@meta.features) <- make.names(T_all@assays[["RNA"]]@meta.features$mgi_symbol,unique=T)
rownames(T_all@assays[["RNA"]]@data) <- make.names(T_all@assays[["RNA"]]@meta.features$mgi_symbol,unique=T)
rownames(T_all@assays[["RNA"]]@counts) <- make.names(T_all@assays[["RNA"]]@meta.features$mgi_symbol,unique=T)
# #TCR and Ig identification and filtering
ig_list <- c("IG_C_gene", "IG_C_pseudogene", "IG_D_gene", "IG_D_pseudogene", "IG_J_gene", "IG_LV_gene", 
             "IG_pseudogene", "IG_V_gene", "IG_V_pseudogene")
tr_list <-c("TR_V_gene", "TR_V_pseudogene", "TR_D_gene", "TR_J_gene", "TR_J_pseudogene", "TR_C_gene")
T_all <- subset(T_all, features=rownames(T_all[!(T_all@assays[["RNA"]]@meta.features$gene_biotype %in% ig_list),]))
T_all <- subset(T_all, features=rownames(T_all[!(T_all@assays[["RNA"]]@meta.features$gene_biotype %in% tr_list),]))
T_all <- subset(T_all, features=rownames(T_all[!(is.na(T_all@assays[["RNA"]]@meta.features$mgi_symbol)),]))
#eliminate 18S rRNA contaminating genes 
rRNA.genes<-c("AY036118","Gm42418")
T_all <- subset(T_all, features=rownames(T_all[!(T_all@assays[["RNA"]]@meta.features$mgi_symbol %in% rRNA.genes),]))
T_all <- subset(T_all, features=rownames(T_all[!(rownames(T_all) %in% rRNA.genes),]))
# mitochondrial DNA percentage
T_all <- PercentageFeatureSet(T_all, pattern = "^mt", col.name = "percent.mt")
T_all <- PercentageFeatureSet(T_all, pattern = "^Rpl", col.name = "percent.Rpl")
T_all <- PercentageFeatureSet(T_all, pattern = "^Rps", col.name = "percent.Rps")
#cell cycle regression
m_cc<-readRDS("m_cc.rds") 
T_all <- CellCycleScoring(T_all, s.features = m_cc$s.genes, g2m.features = m_cc$g2m.genes, set.ident = TRUE)
T_all$CC.Difference <- T_all$S.Score - T_all$G2M.Score
#HSP regression
HSP_genes<-genes.meta[genes.meta$go_id=="GO:0034605",]$mgi_symbol
T_all <- AddModuleScore(object = T_all,features = list(HSP_genes),name = 'HSP.score')
Idents(T_all)<-"condition"
VlnPlot(object = T_all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt",
                                     "percent.Rpl","percent.Rps","HSP.score1"),pt.size=0, ncol = 3)
ggsave2("vln.QC.png",device="png")
#filter for cells with 200-2500 gene counts, <10% mito genes
save.image("temp.preQC.RData")
dim(T_all)
T_all<- subset(x = T_all, subset = nFeature_RNA > 1000 & nFeature_RNA < 3500 & nCount_RNA>2000 & 
                  percent.mt >  -Inf & percent.mt < 7 & percent.Rpl < 20 & percent.Rps < 15 &
                 S.Score <0.15 & G2M.Score<0.15) 
dim(T_all)
rm(list=setdiff(ls(), c("T_all", "metadata","genes.meta")))

##Integrating by condition
T_all.list <- SplitObject(T_all, split.by = "condition")
for (i in 1:length(T_all.list)) {
  T_all.list[[i]] <- SCTransform(T_all.list[[i]], vars.to.regress = c("nCount_RNA", "percent.mt","HSP.score1",
                                                                      "percent.Rpl","percent.Rps")) #", "CC.Difference",
}
T_all.features <- SelectIntegrationFeatures(object.list = T_all.list, nfeatures = 3000)
T_all.list <- PrepSCTIntegration(object.list = T_all.list, anchor.features = T_all.features)
T_all.anchors <- FindIntegrationAnchors(object.list = T_all.list, normalization.method = "SCT",anchor.features = T_all.features)
T_all.combined <- IntegrateData(anchorset = T_all.anchors, normalization.method = "SCT")
T_all.combined <- RunPCA(T_all.combined)
rm(list=setdiff(ls(), c("T_all.combined", "metadata","genes.meta")))

## t-SNE and Clustering
T_all.combined <- RunUMAP(T_all.combined, reduction = "pca", dims = 1:25) #change based on elbow plot
T_all.combined <- FindNeighbors(T_all.combined, reduction = "pca", dims = 1:25) #change based on elbow plot
T_all.combined <- FindClusters(T_all.combined, resolution = 0.15) #adjust resolution (bigger=more clusters), initially used 0.2
T_all.combined <- RunTSNE(object = T_all.combined, dims.use = 1:25, do.fast = TRUE) #change number of PCs to use, change perplexity (https://distill.pub/2016/misread-tsne/)
DimPlot(T_all.combined, reduction = "umap")
for (i in c("condition","Phase","mouse_ID","gender","batch","seurat_clusters")){
  Idents(T_all.combined)<-i
  DimPlot(T_all.combined, reduction = "umap",pt.size=0.1)
  ggsave2(paste0(i,".umap.png"),width=6, height=5,device="png")
}
for (i in c("nCount_RNA","nFeature_RNA","percent.mt","percent.Rps","percent.Rpl","CC.Difference","HSP.score1","S.Score","G2M.Score")){
  FeaturePlot(T_all.combined, features= i,split.by = "condition",pt.size=0.1, order=T)
  ggsave2(paste0(i,".umap.png"),width=10, height=5,device="png")
}

clusters<-as.data.frame(table(T_all.combined@meta.data$seurat_clusters))
letters<-c("A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P")
for(i in 1:nrow(clusters)){
  T_all.combined@meta.data$cluster_letter[T_all.combined@meta.data$seurat_clusters%in%clusters$Var1[i]]<-letters[i]
}
T_all.combined.sce <- as.SingleCellExperiment(T_all.combined)
save(T_all.combined.sce,file="sce.clusters.RData")


## Cluster Analysis
rm(list=setdiff(ls(), c("T_all.combined","metadata","genes.meta")))
DefaultAssay(T_all.combined) <- "RNA"
T_all.combined<- NormalizeData(T_all.combined) #if SCTransform data
#comparing cluster frequencies
T_all.combined@meta.data$my.clusters <- Idents(T_all.combined)  # Store cluster identities in object@meta.data$my.clusters
cluster.counts<-table(T_all.combined@meta.data$my.clusters, T_all.combined@meta.data$condition)
# Finding conserved markers of clusters
clusters<-names(table(T_all.combined$my.clusters))

#set metadata and subsets
mice<-names(table(T_all.combined$mouse_ID))
Idents(T_all.combined) <- "condition"
m564.combined<-subset(T_all.combined, idents='m564')
AID.combined<-subset(T_all.combined, idents='AID')
AID.mouseID<-metadata[metadata$condition=="AID",]$mouse_ID
m564.mouseID<-metadata[metadata$condition=="m564",]$mouse_ID

##Condition DE
DefaultAssay(T_all.combined) <- "RNA"
Idents(T_all.combined) <- "condition" #setting idents to condition metadata
T_all.combined.autoimmune.response <- FindMarkers(T_all.combined, ident.1 = "m564", ident.2 = "AID", 
                                                  min.pct=0,logfc.threshold = -Inf,test.use = "MAST")
avg.T_all.combined <- log1p(AverageExpression(T_all.combined)$RNA)
avg.T_all.combined$gene <- rownames(avg.T_all.combined)
top30 <- head(T_all.combined.autoimmune.response, n = 30)
T_all.combined.autoimmune.response$log2FC<-log2(exp(T_all.combined.autoimmune.response$avg_logFC))
write.csv(T_all.combined.autoimmune.response ,"564vsAID.csv")
T_all.combined$celltype.condition <- paste(T_all.combined$my.clusters, T_all.combined$condition, sep = "_") #adding metadata identifier

###Within Cluster DE
for (i in clusters){
  Idents(T_all.combined) <- "celltype.condition" #setting idents to new metadata column
  df <- FindMarkers(T_all.combined, ident.1 = paste0(i,"_m564"), ident.2 = paste0(i,"_AID"), 
                    min.pct=0,logfc.threshold = -Inf,test.use = "MAST")
  df2 <- head(df, n = 30)
  assign(paste0("cluster",i,"top30"),df2)
  df$log2FC<-log2(exp(df$avg_logFC))
  write.csv(df,paste0("cluster",i,".autoimmune.markers.csv"))
  assign(paste0("cluster",i,".autoimmune.response"),df)
  Idents(T_all.combined) <- "my.clusters"
  temp<- subset(T_all.combined, idents = i)
  Idents(temp) <- "condition"
  df <- log1p(AverageExpression(temp)$RNA)
  df$gene <- rownames(df)
  assign(paste0("avg.cluster",i),df)
  assign(paste0("cluster",i),temp)
}


Idents(T_all.combined) <- "my.clusters"
DefaultAssay(T_all.combined) <- "SCT"

#rename clusters
new.cluster.ids <- c("Tfh-activated","Tfr","Sostdc1","Tfh-CM","Tfh-effector","Tfh-ISG")
names(new.cluster.ids) <- levels(T_all.combined)
T_all.combined <- RenameIdents(T_all.combined, new.cluster.ids)
T_all.combined@meta.data$my.clusters2 <- Idents(T_all.combined)  

# # Visualization
DimPlot(T_all.combined, reduction = "umap",pt.size=0.001)+ #NoLegend()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text=element_blank(),axis.ticks=element_blank(),axis.line=element_blank())
ggsave2("umap.png",width=7, height=5,device="png")
ggsave2("umap.highres.tiff",width=8, height=6,dpi=300,device="tiff")
DimPlot(T_all.combined, reduction = "umap", split.by = "condition")+ 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text=element_blank(),axis.ticks=element_blank(),axis.line=element_blank())
ggsave2("umap2.png",width=12, height=5,device="png")
#by mouse and condition
Idents(T_all.combined) <- "condition"
m564.combined<-subset(T_all.combined, idents='m564')
AID.combined<-subset(T_all.combined, idents='AID')
Idents(m564.combined) <- "my.clusters"
DimPlot(m564.combined, reduction = "umap", split.by = "mouse_ID")+ NoLegend()+NoAxes()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5))
ggsave2("umap.m564.mouse.png",width=10, height=2.25,device="png")
Idents(AID.combined) <- "my.clusters"
DimPlot(AID.combined, reduction = "umap", split.by = "mouse_ID")+ NoLegend()+NoAxes()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5))
ggsave2("umap.AID.mouse.png",width=10, height=2.25,device="png")
gene.list<-c("Cxcr5","Pdcd1","Cd4","Cd40lg","Icos","Bcl6","Il4","Il21","Cxcr3","Maf",
             "Cxcr4","Ccr6","Ccr7","Selplg","Btla","Il2ra","Il7r","Il21r","Tnfrsf4","Slamf1")
p<-FeaturePlot(object = T_all.combined, features = gene.list, 
               cols = c("grey", "blue"), reduction = "umap",combine=F, pt.size=0.001, order=T)
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoLegend()+NoAxes()+theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5))
}
cowplot::plot_grid(plotlist = p,ncol=4)
ggsave2("umap.GOI.genes.png",width=12,height=15,device="png")

## Cluster Analysis
Idents(T_all.combined) <- "my.clusters2"
DefaultAssay(T_all.combined) <- "integrated" 
T_all.combined<-BuildClusterTree(T_all.combined)
png("cluster.tree.png",width=4,height=6,units="in",res=300)
PlotClusterTree(T_all.combined,font=1)
dev.off()
# Cluster heatmap
top8 <- T_all.markers %>% group_by(cluster) %>% top_n(8, avg_logFC)
DefaultAssay(T_all.combined) <- "integrated"
DoHeatmap(object = T_all.combined, features = top8$gene, label = TRUE,size=5,angle=0,hjust=0.5)  
ggsave2("cluster.heatmap.png",width=6, height=8,device="png")
DefaultAssay(T_all.combined) <- "RNA"
#UMAP by cluster marker
gene.list<-c("Klf2","Tnfsf8","Itgb1","Foxp3","Sostdc1","Sox4","Ccl5","Ifit1")
p<-FeaturePlot(object = T_all.combined, features = gene.list, 
               cols = c("grey", "blue"), reduction = "umap",combine=F)
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoLegend() + NoAxes()+theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5))
}
cowplot::plot_grid(plotlist=p,ncol=4)
ggsave2("umap.clustermarkers.png",width=12,height=6,device="png")
#vln cluster marker
p<-VlnPlot(T_all.combined, features = gene.list,group.by = "my.clusters2",pt.size = 0, combine=F)
for(i in 1:length(p)) {
  p[[i]] <- p[[i]]+ NoLegend()+theme(axis.title=element_blank())
}
cowplot::plot_grid(plotlist = p,ncol=4)
ggsave2("vlnplot.clustermarkers.png",width=18, height=9,device="png")

#dot plot
top6$gene<-make.unique(top6$gene,sep="--")
DotPlot(T_all.combined, features = top6$gene)+ 
  coord_flip()+theme(legend.title = element_blank(),axis.title=element_blank())+ RotatedAxis()
ggsave2("dotplot.png",width=6, height=12,device="png")

# Cluster frequency comparison
cluster.Counts <- table(T_all.combined@meta.data$my.clusters2,T_all.combined@meta.data$condition)
cluster.prop <- as.data.frame(scale(cluster.Counts,scale=colSums(cluster.Counts),center=FALSE)*100) 
ggplot(cluster.prop, aes(fill=Var1,y=Freq, x=Var2,alluvium=Var1,stratum=Var1)) + 
  geom_lode()+geom_flow()+geom_stratum(alpha=0) +theme_classic()+
  theme(legend.title = element_blank(),axis.title=element_blank())+ 
  scale_x_discrete(labels= c("WT","564Igi"))#+ theme(legend.position = "none")
ggsave2("cluster.prop.flow.png",width=4, height=4,device="png")

#UMAP comparison
gene.list<-c("Gm42031","Ly6a","Lag3","Id3")
p<-FeaturePlot(T_all.combined, features = gene.list, split.by = "condition",
               pt.size=0.001,combine=F)
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoLegend()+NoAxes()+
    theme(panel.border = element_rect(colour = "black", size=1),
          plot.title=element_blank(),axis.title.y.right=element_blank(),
          axis.line=element_blank())
}
cowplot::plot_grid(p[[1]],p[[5]],p[[2]],p[[6]],p[[3]],p[[7]],p[[4]],p[[8]],ncol=2)
ggsave2("umap.GOI.bycondition.png",width=6, height=10,device="png")

#Population VlnPlot
DefaultAssay(T_all.combined)<-"RNA"
p<-VlnPlot(object = T_all.combined, features =gene.list, pt.size=0.001,group.by ="condition",cols=c("grey","red"),combine=F)
for(i in 1:(length(p))) {
  p[[i]] <- p[[i]]+theme(axis.title=element_blank(),axis.text.x=element_blank(),
                         plot.title = element_text(size=30))+NoLegend()
}
cowplot::plot_grid(plotlist = p,ncol=2)
ggsave2("vlnpopulation.GOI.bycondition.png",width=8, height=12,device="png")
#Cluster VlnPlot
VlnPlot(T_all.combined, features = "Vsir", split.by = "condition", group.by = "my.clusters", 
        pt.size = 0.5, ncol=2)#,legend='right')
ggsave2("vlnplot.c3.cluster.png",width=8, height=5,device="png")
p<-VlnPlot(T_all.combined, features = gene.list, split.by = "condition", 
           group.by = "my.clusters2",cols=c("grey","red"),pt.size = 0,combine=F)
for(i in 1:(length(p)-1)) {
  p[[i]] <- p[[i]]+theme(axis.title=element_blank(),axis.text.x=element_blank())# NoLegend()+
}
p[[length(p)]] <- p[[length(p)]]+theme(axis.title=element_blank())
cowplot::plot_grid(plotlist = p,ncol=1,rel_heights=c(1,1,1,1.5))
ggsave2("vlnplot.GOI.cluster.png",width=5.5, height=12,device="png")

#Volcano Plots
for(i in c(ls(pattern=".autoimmune.response"))){
  df<-get(i)
  if(is.numeric(df$avg_logFC)){
    df$log2FC<-log2(exp(df$avg_logFC))
    gene.list<-rownames(df[abs(df$log2FC)>0.2&df$p_val_adj<0.01,])
    if(length(gene.list)>20){gene.list<-rownames(df[abs(df$log2FC)>0.2&df$p_val_adj<0.001,])}
    if(length(gene.list)>20){gene.list<-rownames(df[abs(df$log2FC)>0.25&df$p_val_adj<10e-10,])}
    if(length(gene.list)>20){gene.list<-rownames(df[abs(df$log2FC)>0.3&df$p_val_adj<10e-50,])}
    if(length(gene.list)>20){gene.list<-rownames(df[abs(df$log2FC)>0.75&df$p_val_adj<10e-100,])}
    if(length(gene.list)>20){gene.list<-rownames(df[abs(df$log2FC)>1.5&df$p_val_adj<10e-200,])}
    if(length(gene.list)>30){gene.list<-NA}
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


save.image("graphed.RData")
save(list=c(ls(pattern="DE"),ls(pattern=".markers"),ls(pattern=".response"),"clusters"),file="DE.dfs.seurat.RData")
