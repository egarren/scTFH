##DeepTCR
import matplotlib
matplotlib.use("agg") 
import os, sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from DeepTCR.DeepTCR import DeepTCR_U
from DeepTCR.DeepTCR import DeepTCR_SS

###Unsupervised CNN on repertoires (disease-based) for clustering
# Instantiate training object
DTCRUr = DeepTCR_U('./unsupervised.rep')
#Load Data from directories
DTCRUr.Get_Data(directory='../scTfh.data',Load_Prev_Data=False,aggregate_by_aa=True,
               aa_column_beta=2,v_beta_column=3,j_beta_column=4, d_beta_column=5,
               aa_column_alpha=6, v_alpha_column=7,j_alpha_column=8,count_column=12)
colors={"AID":"grey","m564":"red"}
# #Train VAE
DTCRUr.Train_VAE(Load_Prev_Data=False,var_explained=0.99)
#visualize training features/latent dimensions
features = DTCRUr.features
DTCRUr.Sample_Features()
DTCRUr.sample_features
print(features.shape)
DTCRUr.sample_features.to_csv("unsupervised.rep.features.csv")
DTCRUr.HeatMap_Sequences()
plt.savefig("unsupervised.rep.heatmap.seq.png")
DTCRUr.HeatMap_Samples()
plt.savefig("unsupervised.rep.heatmap.sample.ratio.png")
#Clustering
DTCRUr.Cluster(clustering_method='phenograph',write_to_sheets=True,sample=1000) #hierarchical, dbscan; use sample to downsample
DFs = DTCRUr.Cluster_DFs
df = pd.concat(DFs,keys=list(range(0,len(DFs))),names=['cluster','rowID'])
df.to_csv("unsupervised.rep.clusters.csv")
DTCRUr.UMAP_Plot(by_cluster=True, scale=5,show_legend=False)
plt.savefig("unsupervised.rep.umap.VAE.clusters.png")
DTCRUr.UMAP_Plot(Load_Prev_Data=True,by_sample=True, scale=5,show_legend=False)
plt.savefig("unsupervised.rep.umap.VAE.samples.png")
DTCRUr.UMAP_Plot(Load_Prev_Data=True,by_class=True,freq_weight=True,scale=2500,alpha=0.5)
plt.savefig("unsupervised.rep.umap.VAE.class.png")
DTCRUr.UMAP_Plot(Load_Prev_Data=True,by_class=True,freq_weight=True,scale=2500,alpha=0.5,show_legend=False)
plt.savefig("unsupervised.rep.umap2.VAE.class.png")
DTCRUr.UMAP_Plot_Samples(scale=100)
plt.savefig("unsupervised.rep.umap.samples.png")
#export UMAP
import umap
import numpy
import pickle
umap_obj = umap.UMAP()
features_umap = umap_obj.fit_transform(DTCRUr.features)
df=numpy.c_[features_umap,DTCRUr.class_id,DTCRUr.sample_id,DTCRUr.beta_sequences,DTCRUr.alpha_sequences,DTCRUr.freq,DTCRUr.counts]
numpy.savetxt("unsup.rep.umap.csv",df,delimiter=",",fmt="%s")# save train feature numpy array to csv
pickle.dump(DTCRUr.features,open("unsup.rep.obj","wb"))
#Repertoire visualization
DTCRUr.Repertoire_Dendrogram(n_jobs=40,distance_metric='KL')
DTCRUr.Repertoire_Dendrogram(lw=6,gridsize=25,Load_Prev_Data=True,
    gaussian_sigma=0.75,dendrogram_radius=0.2,repertoire_radius=0.3,color_dict=colors) #,sample_labels=True
plt.savefig("unsupervised.rep.repertoire.dendrogram.png")
#KNN classification (testing ML predictive capability)
DTCRUr.KNN_Repertoire_Classifier(metrics=['AUC'],distance_metric='KL',plot_metrics=True,by_class=True)
# DTCRUr.KNN_Repertoire_Classifier(metrics=['AUC','Recall'], #,"Precision","F1_Score"
#   distance_metric='KL',plot_metrics=True,by_class=True,Load_Prev_Data=True) #metrics: KL, correlation, euclidean, wasserstein, JS; Sample to speed up
plt.savefig("unsupervised.rep.KNN.scores.png")
DTCRUr.KNN_Repertoire_DF.to_csv("unsupervised.rep.KNN.csv")
#Motif Identification (for WebLogo)
DTCRUr.Motif_Identification(group='m564',by_samples=True)
DTCRUr.Motif_Identification(group='AID',by_samples=True)
#Structural Diversity
DTCRUr.Structural_Diversity(sample=1000) #sample=500
print(DTCRUr.Structural_Diversity_DF)
DTCRUr.Structural_Diversity_DF.to_csv("unsupervised.rep.structural.diversity.csv")



###Supervised CDR3 sequence classification (tetramer-sorted cells)
# Instantiate training object
DTCR_SS = DeepTCR_SS('./supervised.seq.class')
#Load Data from directories
DTCR_SS.Get_Data(directory='../training.data/epitope',Load_Prev_Data=False,aggregate_by_aa=True,
               aa_column_beta=0)
#train (method #3 - Monte Carlo)
DTCR_SS.Monte_Carlo_CrossVal(test_size=0.25,folds=5)
DTCR_SS.AUC_Curve()
DTCR_SS.AUC_Curve(figsize=(10,10),legend_font_size=1,legend_loc="lower right")
plt.savefig("supervised.seq.method3.AUC.samples.png")
df=DTCR_SS.AUC_DF
df.to_csv("supervised.seq.AUC.csv")
#visualize training features/latent dimensions
features = DTCR_SS.features
DTCR_SS.Sample_Features()
DTCR_SS.sample_features
print(features.shape)
DTCR_SS.sample_features.to_csv("supervised.seq.features.csv")
DTCR_SS.HeatMap_Sequences()
plt.savefig("supervised.seq.heatmap.seq.png")
DTCR_SS.HeatMap_Samples()
plt.savefig("supervised.seq.heatmap.sample.ratio.png")
#strongest predictive sequences
DTCR_SS.Representative_Sequences()
print(DTCR_SS.Rep_Seq['NP'])
df = pd.concat(DTCR_SS.Rep_Seq,keys=list(DTCR_SS.Rep_Seq.keys()),names=['Antigen','rowID'])
df.to_csv("supervised.seq.rep_seq.csv")
DTCR_SS.Motif_Identification('NP')
#visualize learned latent space
DTCR_SS.UMAP_Plot(by_class=True,freq_weight=True,scale=1000)
DTCR_SS.UMAP_Plot(by_class=True,freq_weight=True,scale=2000,set='train') #specify if using train, valid, or test data
DTCR_SS.UMAP_Plot(by_class=True,freq_weight=True,scale=2000,set='valid') #specify if using train, valid, or test data
DTCR_SS.UMAP_Plot(by_class=True,freq_weight=True,scale=2000,set='test') #specify if using train, valid, or test data
DTCR_SS.UMAP_Plot(by_class=True,alpha=0.5,show_legend=False,scale=1,Load_Prev_Data=True)
plt.savefig("supervised.seq.umap.png")
#repertoire relatedness
DTCR_SS.Repertoire_Dendrogram(gridsize=50,gaussian_sigma=0.75,lw=6,dendrogram_radius=0.3)
DTCR_SS.Repertoire_Dendrogram(lw=6,gridsize=50,Load_Prev_Data=True,
    gaussian_sigma=0.75,dendrogram_radius=0.3,repertoire_radius=0.4) #,sample_labels=True
plt.savefig("supervised.seq.dendrogram.png")
#Inference: test new data, transform to feature space
DTCRU_test.Get_Data(directory='../scTfh.data',Load_Prev_Data=False,aggregate_by_aa=True,
               aa_column_beta=2,v_beta_column=3,j_beta_column=4)
features_test,_ = DTCR_SS.Sequence_Inference(beta_sequences=DTCRU_test.beta_sequences,
  v_beta=DTCRU_test.v_beta,j_beta=DTCRU_test.j_beta)
features_test,_ = DTCR_SS.Sequence_Inference(beta_sequences=DTCRU_test.beta_sequences)
#compare using UMAP
import umap
import numpy
import pickle
umap_obj = umap.UMAP()
features_orig = umap_obj.fit_transform(DTCR_SS.features)
df=numpy.c_[features_orig,DTCR_SS.class_id,DTCR_SS.sample_id,DTCR_SS.beta_sequences]
numpy.savetxt("sup.seq.train.csv",df,delimiter=",",fmt="%s")# save train feature numpy array to csv
plt.clf()
plt.scatter(features_orig[:,0],features_orig[:,1],c="grey",s=0.1)
plt.savefig("supervised.seq.train.umap.samples.png")
features_new = umap_obj.transform(features_test)
df=numpy.c_[features_new,DTCRU_test.class_id,DTCRU_test.sample_id,DTCRU_test.beta_sequences]
numpy.savetxt("sup.seq.test.csv",df,delimiter=",",fmt="%s")# save train feature numpy array to csv
plt.scatter(features_new[:,0],features_new[:,1],c="red",s=0.1)
plt.savefig("supervised.seq.test.umap.samples.png")
plt.clf()
DTCR_SS.UMAP_Plot(by_class=True,alpha=0.01,show_legend=False,scale=5,Load_Prev_Data=True,sample=10) #prob_plot
plt.savefig("supervised.seq.umap.png")
pickle.dump(DTCR_SS.features,open("train.obj","wb"))
pickle.dump(features_test,open("test.obj","wb"))




