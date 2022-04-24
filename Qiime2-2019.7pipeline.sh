#!/bin/bash
wd=$1
echo -e "\033[43;30m Qiime2 pipeline start! \033[0m"
cd ${wd}
mkdir results
sed 's/\r//' metadata.txt >> metadata.tsv
mv metadata.tsv ./results
mkdir ./results/taxa
mkdir ./results/ancom
source activate qiime2-2019.7
#Samples_list
cd ./seqs
echo "sample-id,absolute-filepath,direction" >> Samples_list.csv
for i in `ls *fastq.gz | awk '{print $0}'`;
do
ls $i | grep 1.fastq.gz | awk -v a=`echo $i | sed s/_R1.fastq.gz//` -v b=$PWD -v c=$i '{print a","b"/"c",""forward"}' >> Samples_list.csv ;
ls $i | grep 2.fastq.gz | awk -v a=`echo $i | sed s/_R2.fastq.gz//` -v b=$PWD -v c=$i '{print a","b"/"c",""reverse"}' >> Samples_list.csv ;
done
echo -e "\033[43;30m Samples_list complete! \033[0m"

# import data
time qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path Samples_list.csv \
  --output-path demux.qza \
  --input-format PairedEndFastqManifestPhred33

# demux
time qiime demux summarize \
  --i-data demux.qza \
  --o-visualization demux.qzv

# representative-sequences
cd ../results
nohup time qiime dada2 denoise-paired \
 --i-demultiplexed-seqs ../seqs/demux.qza \
 --p-trim-left-f 17 --p-trim-left-r 20 \
 --p-trunc-len-f 260 --p-trunc-len-r 260 \
 --o-table table.qza \
 --o-representative-sequences rep-seqs.qza \
 --o-denoising-stats denoising-stats.qza --p-n-threads 30 &

# rep-seqs
time qiime feature-table summarize \
  --i-table table.qza \
  --o-visualization table.qzv \
  --m-sample-metadata-file metadata.tsv
time qiime feature-table tabulate-seqs \
  --i-data rep-seqs.qza \
  --o-visualization rep-seqs.qzv

# exported-table
qiime tools export \
--input-path table.qza \
--output-path exported-table

# feature-table.tsv
biom convert -i exported-table/feature-table.biom -o exported-table/feature-table.tsv --to-tsv

## rooted-tree.qza
time qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza
  
# alpha-rarefaction
time qiime diversity alpha-rarefaction \
  --i-table table.qza \
  --i-phylogeny rooted-tree.qza \
  --p-max-depth 50000 \
  --m-metadata-file metadata.tsv \
  --o-visualization alpha-rarefaction.qzv

# core-metrics-results
time qiime diversity core-metrics-phylogenetic \
  --i-phylogeny rooted-tree.qza \
  --i-table table.qza \
  --p-sampling-depth 10000 \
  --m-metadata-file metadata.tsv \
  --output-dir core-metrics-results

# chao1、simpson、ace、euclidean
qiime diversity alpha \
 --i-table table.qza \
 --p-metric chao1 \
 --o-alpha-diversity core-metrics-results/chao1_vector.qza

qiime diversity alpha \
 --i-table table.qza \
 --p-metric simpson \
 --o-alpha-diversity core-metrics-results/simpson_vector.qza

qiime diversity beta \
 --i-table table.qza \
 --p-metric euclidean \
 --o-distance-matrix core-metrics-results/euclidean_distance_matrix.qza

# Beta_distance_significance
distance="euclidean jaccard bray_curtis unweighted_unifrac weighted_unifrac"
for i in $distance;
do
    for j in `awk 'NR==1{for(k=2;k<=NF;k++) print $k}' metadata.tsv`;
    do
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/${i}_distance_matrix.qza \
  --m-metadata-file metadata.tsv \
  --m-metadata-column ${j} \
  --o-visualization core-metrics-results/${i}_${j}_significance.qzv \
  --p-pairwise
    done
done

# alpha_index_significance
index="observed_otus evenness shannon faith_pd chao1 simpson"
for k in $index;
do
qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/${k}_vector.qza \
  --m-metadata-file metadata.tsv \
  --o-visualization core-metrics-results/${k}_significance.qzv
done

#export data
mkdir -p core-metrics-results/alpha_diversity core-metrics-results/beta_diversity
for i in $distance;
do
qiime tools export \
  --input-path core-metrics-results/${i}_distance_matrix.qza \
  --output-path core-metrics-results/beta_diversity/${i}_distance_matrix
done
for k in $index;
do
qiime tools export \
 --input-path core-metrics-results/${k}_vector.qza \
 --output-path core-metrics-results/alpha_diversity/${k}
done
echo -e "\033[43;30m Alpha & Beta diversity metrics have been exported! \033[0m"
 
#  annotation
time qiime feature-classifier classify-sklearn \
  --i-classifier /home/foodbio/database/Qiime2/silva-132-99-nb-341F-806R_classifier.qza\
  --i-reads rep-seqs.qza \
  --o-classification taxa/taxonomy.qza --p-n-jobs 30
  
time qiime metadata tabulate \
  --m-input-file taxa/taxonomy.qza \
  --o-visualization taxa/taxonomy.qzv
  
# exported-table annotation
time qiime tools export \
 --input-path taxa/taxonomy.qza \
 --output-path ./exported-table

time qiime feature-table filter-samples \
  --i-table table.qza \
  --p-min-frequency 10000 \
  --o-filtered-table table-10k.qza

time qiime taxa barplot \
  --i-table table-10k.qza \
  --i-taxonomy taxa/taxonomy.qza \
  --m-metadata-file metadata.tsv \
  --o-visualization taxa/taxa-bar-plots.qzv

# filtered-table
time qiime feature-table filter-features \
  --i-table table-10k.qza \
  --p-min-frequency 50 \
  --p-min-samples 4 \
  --o-filtered-table filtered-table.qza
  
# rel-abundance
mkdir rel-abundance
cd rel-abundance
for i in {2..6};
do
time qiime taxa collapse \
  --i-table ../filtered-table.qza \
  --i-taxonomy ../taxa/taxonomy.qza \
  --p-level $i \
  --o-collapsed-table table-l${i}.qza
time qiime feature-table relative-frequency \
  --i-table table-l${i}.qza \
  --o-relative-frequency-table rel-abun-l${i}.qza
time qiime tools export \
  --input-path rel-abun-l${i}.qza \
  --output-path rel-abun-l${i}
biom convert -i rel-abun-l${i}/feature-table.biom -o rel-abun-l${i}/rel-abun-l${i}.tsv --to-tsv
done
