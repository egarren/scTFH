rm(list=ls())
library(immunarch)
library(cowplot)
library(ggrepel)
library(Seurat)
library(dplyr)
library(ggpubr)
library(EnhancedVolcano)
library(ggseqlogo)
library(ggplot2)
library(ggsci)
library(viridis)
library(RColorBrewer)
library(plyr)
library(patchwork)
library(gridExtra)
my.ttest <- function(...) {
  obj<-try(t.test(...), silent=TRUE)
if (is(obj, "try-error")) return(NA) else return(obj$p.value)
}
data_concater <- function(x){
  x<- levels(factor(x))
  paste(x, collapse = "+")
}

load("clone.ab.data.RData")

##Immunarch
#by condition
dir.create("./rep")
rep.tabs<-list.files("./vdjtools/out/annot",".txt$")
file.copy(file.path("./vdjtools/out/annot",rep.tabs),"./rep")
df<-read.table("./vdjtools/out/annot/metadata.txt",header=T)
df<-cbind(Sample=gsub(".txt$","",df$file_name),df)
write.table(df, file = "./rep/metadata.txt", sep = "\t",quote=F,row.names = FALSE)

# Load the data to the package
# setwd("./immunarch.results")
immdata = repLoad("./rep")
save(immdata, file="immdata.Rdata")
#Basic anlaysis (explore)
exp_vol = repExplore(immdata$data, .method = "volume")
vis(exp_vol, .by = c("condition"), .meta = immdata$meta) #.by=c("condition',"cluster")
ggsave2("number.unique.clonotypes.png",width=4, height=6,device="png") #for publication use fixVis
exp_len = repExplore(immdata$data, .method = "len", .col = "aa")
vis(exp_len,.by=c("condition"), .meta = immdata$meta)
ggsave2("cdr3.length.png",width=20, height=6,device="png") #for publication use fixVis
exp_cnt = repExplore(immdata$data, .method = "count")
vis(exp_cnt)#.by=c("condition)"), .meta = immdata$meta)
ggsave2("clonotype.abundance.png",width=6, height=6,device="png") #for publication use fixVis
#clonality
for(i in c("clonal.prop","homeo","top","rare")){
  imm_pr = repClonality(immdata$data, .method = i)
  png(paste0(i,".clonality.png"),width=10,height=6,units="in",res=200)
  grid.arrange(vis(imm_pr),vis(imm_pr,.by=c("condition"),.meta=immdata$meta),ncol=2)
  dev.off()
}

#public repertoire (identify shared clones) using CDR3b
pr = pubRep(immdata$data, "aa", .coding = T, .verbose = F)
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
pr.res.cdr3.graph$aa_length<-nchar(as.character(pr.res.cdr3.graph$CDR3.aa))
pr.res.cdr3.graph<-left_join(pr.res.cdr3.graph,clone.meta[,c("cdr3","chain")],by=c("CDR3.aa"="cdr3"),keep=F)
plot.list<-list()
plot.list2<-list()
for(h in c("TRA","TRB")){
  pr.df<-pr.res.cdr3.graph[pr.res.cdr3.graph$chain==h&!is.na(pr.res.cdr3.graph$chain),]
  cdr3.public.FC<-pr.df[,c("CDR3.aa","Samples.sum","log2FC")]
  colnames(cdr3.public.FC)<-c("CDR3.aa","Samples.sum.mice","log2FC.mice")
  most.shared.cdr3<-pr.df[abs(pr.df$log2FC)>6&
                                        pr.df$sum.564>20&pr.df$sum.AID>20,]
  ggplot(pr.df,aes(x = sum.AID, y = sum.564,size=Samples.sum,label=CDR3.aa))+geom_point(alpha=0.25)+ 
    scale_x_continuous(trans='log2') +scale_y_continuous(trans='log2')+
    geom_abline(intercept = 0, slope = 1) +
    geom_text_repel(data = most.shared.cdr3,aes(label = CDR3.aa),size = 3,
                    color="red",min.segment.length=unit(0,'lines'),nudge_x=0.1,nudge_y=0.1,segment.size=0.25)
  ggsave2(paste0(h,".overlap.AID.564.cdr3.scatter.png"),width=7, height=6,device="png")
  #WebLogo on public repertoire
  AID.xp<-pr.df[pr.df$log2FC< -5,]#&pr.df$sum.AID>0.1,]
  m564.xp<-pr.df[pr.df$log2FC>5,]#&pr.df$sum.564>0.1,]
  df.list<-list()
  for(j in c("AID","m564")){
    df<-get(paste0(j,".xp"))
    if(j=="AID"){name="WT"}else{name="564Igi"}
    for(k in 14){ #10:18
      df2<-df[df$aa_length==k,]
      plot.list[[paste0(h,j,k)]]<-ggplot()+geom_logo(as.character(df2$CDR3.aa),method="probability")+
        theme_logo()+ggtitle(paste0(h,"_",name))+
        theme(plot.title = element_text(hjust = 0.5)) #,"_",k, axis.text.x=element_blank(),
    }
    df$condition<-name
    df.list[[j]]<-df
  }
  #cdr3 length
  aa.df<-rbindlist(df.list)
  aa.df$condition <- factor(aa.df$condition, levels = c("WT","564Igi"))
  mu <- ddply(aa.df, "condition", summarise, grp.mean=mean(aa_length))
  plot.list2[[h]]<-ggplot(aa.df, aes(x=aa_length,color=condition,fill=condition)) + 
    geom_histogram(aes(y=..density..),binwidth=1,position="identity",alpha=0.2) +#geom_density(fill=NA)+
    geom_vline(data=mu, aes(xintercept=grp.mean, color=condition),linetype="dashed")+ggtitle(h)+theme_classic()+
    labs(x = "CDR3 Length", y = "") +theme(plot.title = element_text(hjust = 0.5),legend.title = element_blank())+
    scale_color_manual(values = c("black", "red"))+
    scale_fill_manual(values = c("black", "red"))
}
CombinePlots(plots=list(plot.list[[1]],plot.list[[3]],plot.list[[2]],plot.list[[4]]),ncol=2,legend="bottom")
ggsave2("cdr3.pr.weblogo.bycondition.png",width=6,height=4,device="png")
CombinePlots(plot.list2,ncol=2,legend="right")
ggsave2("cdr3.pr.aa_length.bycondition.histo.png",width=7,height=4,device="png")


#Annotating using antigen prediction databases
load("masterdb3.RData")
pr.df<-pr.res.cdr3.graph[pr.res.cdr3.graph$chain=="TRB"&!is.na(pr.res.cdr3.graph$chain),]
df.db<-unique(masterdb[masterdb$antigen!="unknown"&!is.na(masterdb$antigen)&masterdb$disease2!="species",])
# df.db<-unique(masterdb[!is.na(masterdb$chain)&masterdb$chain=="TRB"&masterdb$TRBV!=""&masterdb$Species=="Mouse"&masterdb$Tcell=="CD4",]) #
pr.res.db<-unique(left_join(x = pr.df, y = df.db[,c("cdr3.b","antigen","epitope","study2","Species","disease2","disease")], 
                        by = c("CDR3.aa"="cdr3.b"),keep=F))
pr.res.db$disease2[is.na(pr.res.db$disease2)]<-"unknown"
labels<-unique(pr.res.db[abs(pr.res.db$log2FC)>3&pr.res.db$Samples.sum>1&!is.na(pr.res.db$antigen)&
                           (pr.res.db$sum.AID>50|pr.res.db$sum.564>50),])
labels2<-labels[labels$antigen %in% c("NP177","M45","IE1","Tat","EBNA4","p65"),]
ggplot(pr.res.db,aes(x = sum.564, y =  sum.AID,size=Samples.sum))+theme_classic()+
  scale_x_continuous(trans='log10') +scale_y_continuous(trans='log10')+ scale_size(range = c(1, 10))+
  geom_point(aes(colour=factor(disease2),fill = factor(disease2)), shape=21) + 
  scale_color_manual(values = c(brewer.pal(n = 5, name = "Dark2"),alpha("black",0.5),"#A6761D"))+
  scale_fill_manual(values = c(alpha(brewer.pal(n = 5, name = "Dark2"),0.5),"#1C00ff00",alpha("#A6761D",0.5)))+
  guides(fill = guide_legend(override.aes = list(size = 5)),color=F)+ #,size=F
  labs(x = "564Igi", y = "WT", size="Samples",fill="Disease")+
  theme(legend.direction = "vertical", legend.box = "horizontal")+#theme(legend.title = element_blank())+
  geom_label_repel(data =labels2,aes(label = antigen),size = 3,min.segment.length=unit(0,'lines'),nudge_x=0.1,nudge_y=0.1,segment.size=0.25)
ggsave2("annot.db.scatter.bycondition.png",width=6, height=3.5,device="png")

#public repertoire using a/b pairing
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
pr.res.cdr3.graph[["log2FC.avg"]]= apply(pr.res.cdr3.graph[, c("avgfreq.AID", "avgfreq.564")],1, function(x) log2((x[2])/(x[1])))
rownames(pr.res.cdr3.graph)<-pr.res.cdr3.graph$CDR3.aa
pr.df<-pr.res.cdr3.graph
cdr3.public.FC<-pr.df[,c("CDR3.aa","Samples.sum","log2FC")]
colnames(cdr3.public.FC)<-c("CDR3.aa","Samples.sum.mice","log2FC.mice")
labels<-pr.df[abs(pr.df$log2FC)>2&(pr.df$sum.AID>20|pr.df$sum.564>50)&pr.df$Samples.sum>1,]
ggplot(pr.df,aes(x = sum.564, y =  sum.AID,size=Samples.sum))+theme_classic()+
  scale_x_continuous(trans='log10') +scale_y_continuous(trans='log10')+ scale_size(range = c(1, 10))+
  geom_point(color=alpha("black",0.5),fill="#1C00ff00", shape=21) + 
  labs(x = "564Igi", y = "WT", size="Samples")+
  geom_text_repel(data =labels,aes(label = CDR3.aa),size = 3,color="forestgreen",
                   min.segment.length=unit(0,'lines'),nudge_x=0.1,nudge_y=0.1,segment.size=0.25)
ggsave2("overlap.TCRab.bycondition.png",width=5, height=4,device="png")


###by clusters
#create repertoire folders
dir.create("./clust.rep")
# dir.create("./bycluster.immunarch.results")
rep.tabs<-list.files("./vdjtools.cluster/clust.out/annot",".txt$")
file.copy(file.path("./vdjtools.cluster/clust.out/annot",rep.tabs),"./clust.rep")
df<-read.table("./vdjtools.cluster/clust.out/annot/metadata.txt",header=T)
df<-cbind(Sample=gsub(".txt$","",df$file_name),df)
write.table(df, file = "./clust.rep/metadata.txt", sep = "\t",quote=F,row.names = FALSE)

# Load the data to the package
immdata = repLoad("./clust.rep")
save(immdata, file="immdata.clust.Rdata")
#Basic anlaysis (explore)
exp_vol = repExplore(immdata$data, .method = "volume")
vis(exp_vol, .by = c("Cluster","condition"), .meta = immdata$meta,.test=F) #.by=c("condition',"cluster")
ggsave2("bycluster.number.unique.clonotypes.png",width=6, height=6,device="png") #for publication use fixVis
exp_cnt = repExplore(immdata$data, .method = "count")
vis(exp_cnt,.by=c("Cluster","condition"), .meta = immdata$meta)
ggsave2("bycluster.clonotype.abundance.png",width=6, height=6,device="png") #for publication use fixVis
#clonality
imm_pr = repClonality(immdata$data, .method = "clonal.prop")
vis(imm_pr,.by=c("Cluster","condition"),.meta=immdata$meta,.test=F)
ggsave2("bycluster.clone.prop.clonality.png",width=6, height=6,device="png") #for publication use fixVis


#public TCR-beta repertoire (identify shared clones)
pr = pubRep(immdata$data, "aa", .coding = T, .verbose = F)
# vis(pr,.by = c("condition"), .meta = immdata$meta)
# ggsave2("public.clonotypes.png",width=7, height=7,device="png")
pr.Tfh = pubRepFilter(pr, immdata$meta, c(Cluster="0"))
pr.Tfr = pubRepFilter(pr, immdata$meta, c(Cluster="1"))
pr.Tfh[is.na(pr.Tfh)]<-0
pr.Tfr[is.na(pr.Tfr)]<-0
pr.Tfh[["avgfreq.Tfh"]] = rowMeans(public_matrix(pr.Tfh), na.rm = T)
pr.Tfr[["avgfreq.Tfr"]] = rowMeans(public_matrix(pr.Tfr), na.rm = T)
pr.Tfh[["sum.Tfh"]] = rowSums(public_matrix(pr.Tfh)[,1:10], na.rm = T)+1
pr.Tfr[["sum.Tfr"]] = rowSums(public_matrix(pr.Tfr)[,1:10], na.rm = T)+1
pr.res.cdr3 = dplyr::full_join(pr.Tfh, pr.Tfr, by = "CDR3.aa")
pr.res.cdr3.graph<-pr.res.cdr3
pr.res.cdr3.graph[is.na(pr.res.cdr3.graph)]<-0
pr.res.cdr3.graph$sum.Tfh[pr.res.cdr3.graph$sum.Tfh==0]<-1
pr.res.cdr3.graph$sum.Tfr[pr.res.cdr3.graph$sum.Tfr==0]<-1
pr.res.cdr3.graph[["Samples.sum"]] = pr.res.cdr3.graph[["Samples.x"]] + pr.res.cdr3.graph[["Samples.y"]]
pr.res.cdr3.graph[["freq.ratio"]] = apply(pr.res.cdr3.graph[, c("sum.Tfh", "sum.Tfr")],1, function(x) log10(x[1])/log10(x[2]))
pr.res.cdr3.graph[["log2FC"]]= apply(pr.res.cdr3.graph[, c("sum.Tfh", "sum.Tfr")],1, function(x) log2((x[2])/(x[1])))
rownames(pr.res.cdr3.graph)<-pr.res.cdr3.graph$CDR3.aa
pr.res.cdr3.graph$aa_length<-nchar(as.character(pr.res.cdr3.graph$CDR3.aa))
pr.res.cdr3.graph<-left_join(pr.res.cdr3.graph,clone.meta[,c("cdr3","chain")],by=c("CDR3.aa"="cdr3"),keep=F)
plot.list<-list()
plot.list2<-list()
for(h in c("TRA","TRB")){
  pr.df<-pr.res.cdr3.graph[pr.res.cdr3.graph$chain==h&!is.na(pr.res.cdr3.graph$chain),]
  cdr3.public.FC<-pr.df[,c("CDR3.aa","Samples.sum","log2FC")]
  colnames(cdr3.public.FC)<-c("CDR3.aa","Samples.sum.mice","log2FC.mice")
  most.shared.cdr3<-pr.df[abs(pr.df$log2FC)>6&
                            pr.df$sum.Tfr>20&pr.df$sum.Tfh>20,]
  ggplot(pr.df,aes(x = sum.Tfh, y = sum.Tfr,size=Samples.sum,label=CDR3.aa))+geom_point(alpha=0.25)+ 
    scale_x_continuous(trans='log2') +scale_y_continuous(trans='log2')+
    geom_abline(intercept = 0, slope = 1) +
    geom_text_repel(data = most.shared.cdr3,aes(label = CDR3.aa),size = 3,
                    color="#C49A00",min.segment.length=unit(0,'lines'),nudge_x=0.1,nudge_y=0.1,segment.size=0.25)
  ggsave2(paste0(h,".overlap.Tfh.Tfr.cdr3.scatter.png"),width=7, height=6,device="png")
  #WebLogo on public repertoire
  Tfh.xp<-pr.df[pr.df$log2FC< -3,]#&pr.df$sum.Tfh>0.1,]
  Tfr.xp<-pr.df[pr.df$log2FC>3,]#&pr.df$sum.Tfr>0.1,]
  df.list<-list()
  for(j in c("Tfh","Tfr")){
    df<-get(paste0(j,".xp"))
    for(k in 14){ #10:18
      df2<-df[df$aa_length==k,]
      plot.list[[paste0(h,j,k)]]<-ggplot()+geom_logo(as.character(df2$CDR3.aa),method="probability")+
        theme_logo()+ggtitle(paste0(h,"_",j))+
        theme(plot.title = element_text(hjust = 0.5)) #,"_",k, axis.text.x=element_blank(),
    }
    df$condition<-j
    df.list[[j]]<-df
  }
  #cdr3 length
  aa.df<-rbindlist(df.list)
  aa.df$condition <- factor(aa.df$condition, levels = c("Tfr","Tfh"))
  mu <- ddply(aa.df, "condition", summarise, grp.mean=mean(aa_length))
  plot.list2[[h]]<-ggplot(aa.df, aes(x=aa_length,color=condition,fill=condition)) + 
    geom_histogram(aes(y=..density..),binwidth=1,position="identity",alpha=0.2) +#geom_density(fill=NA)+
    geom_vline(data=mu, aes(xintercept=grp.mean, color=condition),linetype="dashed")+ggtitle(h)+theme_classic()+
    labs(x = "CDR3 Length", y = "") +
    theme(plot.title = element_text(hjust = 0.5),legend.title = element_blank())+
    scale_color_manual(values = c("#F8766D", "#C49A00"))+
    scale_fill_manual(values = c("#F8766D", "#C49A00"))
}
CombinePlots(plots=list(plot.list[[1]],plot.list[[3]],plot.list[[2]],plot.list[[4]]),ncol=2,legend="bottom")
ggsave2("cdr3.pr.weblogo.TfhTfr.png",width=6,height=4,device="png")
CombinePlots(plot.list2,ncol=2,legend="right")
ggsave2("cdr3.pr.aa_length.TfhTfr.histo.png",width=6,height=4,device="png")
save.image("temp2.immunarch.RData")

#Annotating using antigen prediction databases
load("masterdb3.RData")
pr.df<-pr.res.cdr3.graph[pr.res.cdr3.graph$chain=="TRB"&!is.na(pr.res.cdr3.graph$chain),]
df.db<-unique(masterdb[masterdb$antigen!="unknown"&!is.na(masterdb$antigen)&masterdb$disease2!="species",])
pr.res.db<-unique(left_join(x = pr.df, y = df.db[,c("cdr3.b","antigen","epitope","study2","Species","disease2","disease")], 
                            by = c("CDR3.aa"="cdr3.b"),keep=F))
pr.res.db$disease2[is.na(pr.res.db$disease2)]<-"unknown"
labels<-unique(pr.res.db[abs(pr.res.db$log2FC)>2&pr.res.db$Samples.sum>1&!is.na(pr.res.db$antigen)&
                           (pr.res.db$sum.Tfh>20|pr.res.db$sum.Tfr>20),])
labels2<-labels[labels$antigen %in% c("NP177","M45","IE1","Tat","EBNA4","p65","N"),]
ggplot(pr.res.db,aes(x = sum.Tfr, y =  sum.Tfh,size=Samples.sum))+theme_classic()+
  scale_x_continuous(trans='log10') +scale_y_continuous(trans='log10')+ scale_size(range = c(1, 10))+
  geom_point(aes(colour=factor(disease2),fill = factor(disease2)), shape=21) + 
  scale_color_manual(values = c(brewer.pal(n = 5, name = "Dark2"),alpha("black",0.5),"#A6761D"))+
  scale_fill_manual(values = c(alpha(brewer.pal(n = 5, name = "Dark2"),0.5),"#1C00ff00",alpha("#A6761D",0.5)))+
  guides(fill = guide_legend(override.aes = list(size = 5)),color=F)+ #,size=F
  labs(x = "Tfr", y = "Tfh", size="Samples",fill="Disease")+
  theme(legend.direction = "vertical", legend.box = "horizontal")+#theme(legend.title = element_blank())+
  geom_label_repel(data =labels2,aes(label = antigen),size = 3,min.segment.length=unit(0,'lines'),nudge_x=0.1,nudge_y=0.1,segment.size=0.25)
ggsave2("annot.db.scatter.TfhTfr.png",width=6, height=3.5,device="png")

#using a/b pairs:
#create new imm db
df<-clone.data.ab.seurat
imm.list<-list()
for(i in unique(df$T.sampleID)){ #creating immunarch list from collapsed data
  for(j in unique(df$my.clusters)){
    df2<-df[df$T.sampleID==i&df$my.clusters==j,]
    df2<-transform(df2,Clones=ave(seq(nrow(df2)),cdr3,FUN=length))
    df2<-as.data.table(df2)[, lapply(.SD, data_concater), by=cdr3]
    df2$Clones<-as.numeric(as.character(df2$Clones))
    df2$Proportion<-df2$Clones/sum(df2$Clones)
    imm.list[[paste0(i,".clust",j)]]<-tibble(
      Clones=df2$Clones, Proportion=df2$Proportion, CDR3.nt=df2$cdr3_nt,CDR3.aa=df2$cdr3,
      V.name=df2$v_gene,D.name=df2$d_gene,J.name=df2$j_gene
    )
  }
}
pr = pubRep(imm.list, "aa", .coding = T, .verbose = F)
pr.Tfh = pubRepFilter(pr, immdata$meta, c(Cluster="0"))
pr.Tfr = pubRepFilter(pr, immdata$meta, c(Cluster="1"))
pr.Tfh[is.na(pr.Tfh)]<-0
pr.Tfr[is.na(pr.Tfr)]<-0
pr.Tfh[["avgfreq.Tfh"]] = rowMeans(public_matrix(pr.Tfh), na.rm = T)
pr.Tfr[["avgfreq.Tfr"]] = rowMeans(public_matrix(pr.Tfr), na.rm = T)
pr.Tfh[["sum.Tfh"]] = rowSums(public_matrix(pr.Tfh)[,1:10], na.rm = T)+1
pr.Tfr[["sum.Tfr"]] = rowSums(public_matrix(pr.Tfr)[,1:10], na.rm = T)+1
pr.res.cdr3 = dplyr::full_join(pr.Tfh, pr.Tfr, by = "CDR3.aa")
pr.res.cdr3.graph<-pr.res.cdr3
pr.res.cdr3.graph[is.na(pr.res.cdr3.graph)]<-0
pr.res.cdr3.graph$sum.Tfh[pr.res.cdr3.graph$sum.Tfh==0]<-1
pr.res.cdr3.graph$sum.Tfr[pr.res.cdr3.graph$sum.Tfr==0]<-1
pr.res.cdr3.graph[["Samples.sum"]] = pr.res.cdr3.graph[["Samples.x"]] + pr.res.cdr3.graph[["Samples.y"]]
pr.res.cdr3.graph[["freq.ratio"]] = apply(pr.res.cdr3.graph[, c("sum.Tfh", "sum.Tfr")],1, function(x) log10(x[1])/log10(x[2]))
pr.res.cdr3.graph[["log2FC"]]= apply(pr.res.cdr3.graph[, c("sum.Tfh", "sum.Tfr")],1, function(x) log2((x[2])/(x[1])))
rownames(pr.res.cdr3.graph)<-pr.res.cdr3.graph$CDR3.aa
pr.res.cdr3.graph$aa_length<-nchar(as.character(pr.res.cdr3.graph$CDR3.aa))
pr.res.cdr3.graph<-left_join(pr.res.cdr3.graph,clone.meta[,c("cdr3","chain")],by=c("CDR3.aa"="cdr3"),keep=F)
pr.df<-pr.res.cdr3.graph
labels<-unique(pr.df[abs(pr.df$log2FC)>4.5&pr.df$Samples.sum>1&(pr.df$sum.Tfh>20|pr.df$sum.Tfr>20),])
ggplot(pr.df,aes(x = sum.Tfr, y =  sum.Tfh,size=Samples.sum))+theme_classic()+
  scale_x_continuous(trans='log10') +scale_y_continuous(trans='log10')+ scale_size(range = c(1, 10))+
  geom_point(color=alpha("black",0.5),fill = "#1C00ff00", shape=21) + 
  guides(fill = guide_legend(override.aes = list(size = 5)))+ #,size=F
  labs(x = "Tfr", y = "Tfh", size="Samples")+
  theme(legend.direction = "vertical", legend.box = "horizontal")+#theme(legend.title = element_blank())+
  geom_label_repel(data =labels,aes(label = CDR3.aa),size = 3,color="forestgreen",
                   min.segment.length=unit(0,'lines'),nudge_x=0.1,nudge_y=0.1,segment.size=0.25)
ggsave2("overlap.TCRab.Tfh.Tfr.scatter.png",width=5, height=4,device="png")








