###VDJtools
# #Install
# wget https://github.com/mikessh/vdjtools/releases/download/1.2.1/vdjtools-1.2.1.zip
# unzip vdjtools-1.2.1.zip
#startup
VDJTOOLS="vdjtools-1.2.1.jar" #point to jar
cd vdjtools
#by condition
$VDJTOOLS CalcBasicStats -m vdjtools.metadata.txt out/0
$VDJTOOLS CalcSpectratype -m vdjtools.metadata.txt out/1
$VDJTOOLS CalcSegmentUsage -m vdjtools.metadata.txt -p -f condition out/2
$VDJTOOLS CalcDiversityStats -m vdjtools.metadata.txt out/7
$VDJTOOLS RarefactionPlot -m vdjtools.metadata.txt -f condition -l sample_id out/8
$VDJTOOLS CalcPairwiseDistances -p -m vdjtools.metadata.txt out/10
$VDJTOOLS ClusterSamples -p -f condition -l sample_id out/10 out/10.condition
$VDJTOOLS ClusterSamples -p -e vJSD -f condition -l sample_id out/10 out/10.jsd
$VDJTOOLS JoinSamples -p -m vdjtools.metadata.txt out/12
$VDJTOOLS Annotate -m vdjtools.metadata.txt out/annot/
# $VDJTOOLS FilterNonFunctional -m vdjtools.metadata.txt -c out/nf/

cd ../vdjtools.cluster
#by cluster
$VDJTOOLS CalcBasicStats -m vdjtools.clust.metadata.txt clust.out/0
$VDJTOOLS CalcSpectratype -m vdjtools.clust.metadata.txt clust.out/1
$VDJTOOLS CalcSegmentUsage -m vdjtools.clust.metadata.txt -p -f condition clust.out/2
$VDJTOOLS CalcDiversityStats -m vdjtools.clust.metadata.txt clust.out/7
$VDJTOOLS RarefactionPlot -m vdjtools.clust.metadata.txt -f condition -l sample_id clust.out/8
$VDJTOOLS CalcPairwiseDistances -p -m vdjtools.clust.metadata.txt clust.out/10
$VDJTOOLS ClusterSamples -p -f condition -l sample_id clust.out/10 clust.out/10.condition
$VDJTOOLS ClusterSamples -p -e vJSD -f condition -l sample_id clust.out/10 clust.out/10.jsd
$VDJTOOLS JoinSamples -p -m vdjtools.clust.metadata.txt clust.out/12
$VDJTOOLS Annotate -m vdjtools.clust.metadata.txt clust.out/annot/
# $VDJTOOLS FilterNonFunctional -m vdjtools.clust.metadata.txt -c clust.out/nf/


