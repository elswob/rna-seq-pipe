#!/bin/bash

. /usr/share/Modules/init/bash


module load acml/gfortran64_mp/4.4.0
module load R/gcc-4.4.6/2.15.0
module load jesmaa/python/2.7.2
module load gi/samtools/0.1.19
module load gi/bedtools/2.17.0
module load jesmaa/featureCounts/1.3.6-p1

paste HTSeq/HTSeqGenes HTSeq/Combined/*HTSeqValue >HTSeq/Combined/Allcombined

species=$1
pipe=$2
design=$3

if [ -z $3 ]; then
	echo -e "No design provided\n"
	design="n"
else
	echo -e "Using $3 for design"
	design=$3
fi

R="/home/benels/software/R/R-3.0.2/bin/R"

#################################################################################################################################


#Running EdgeR if all HTSeq samples are complete
echo -e "##################### Running Differential expression with EdgeR #####################\n" 
#module unload R/gcc-4.4.6/2.15.0
module load jesmaa/R/3.0.0  
sort Conditions.txt > Conditions1.txt
mv Conditions1.txt Conditions.txt

if [ ! -s HTSeq/Combined/Allcombined_gene_names ]; then

	if [[ $species == "Human" ]]; then
		$R --vanilla --quiet --slave --args "Human" $design <$pipe/Scripts/EdgeR.R 
	elif [[ $species == "Mouse" ]]; then
		$R --vanilla --quiet --slave --args "Mouse" $design <$pipe/Scripts/EdgeR.R 
	else
		echo -e "Specify Human or Mouse!\n"
	fi

	#remove trailing integer from ensembl IDs
	perl -i -plane 's/(ENSMUSG\d+)\.\d+/$1/' HTSeq/Combined/Allcombined_gene_names
else
	echo -e "DE already done\n"
fi
#################################################################################################################################

#Conducting Gene Enrichment analysis on the differentially expressed genes between each condition and plotting the resutls
mkdir Results/GeneEnrichment/
for i in `ls Results/DE_Genes`; do gawk '{if($3>2) print $1}' Results/DE_Genes/${i} >Results/GeneEnrichment/Upreg_${i} ;done 
for i in `ls Results/DE_Genes`; do gawk '{if($3<-2) print $1}' Results/DE_Genes/${i} >Results/GeneEnrichment/Downreg_${i};done 

cd Results/GeneEnrichment/

echo -e "##################### Conducting gene enrichment analysis #####################\n"

#check if already run
ge_check=$(find GO_Analysis/ -name "*.tsv" | wc -l)
echo -e "GO_Analysis files = $ge_check\n"
if [[ $ge_check <1 ]]; then
	$R --vanilla --quiet --slave --args $pipe $species <$pipe/Scripts/GeneEnrichment.R 
else
	echo -e "Gene enrichment already done\n"
fi

echo -e "##################### Finished with Gene enrichment analysis #####################\n"

#################################################################################################################################


cd GO_Analysis 
find ./ -name *tsv -exec mv {} ./ \;
#echo "Plotting gene enrichment results"
#find ./ -size  0 -print0 |xargs -0 rm


$R --vanilla --quiet --slave --args $pipe <$pipe/Scripts/PlotGene_enrichment_FDR.R

echo -e "##################### Finished plotting gene enrichment results #####################\n"
cd ../../../

#Plotting the various genomic features that showed the most confident differentially expressed genes between conditions
mkdir Results/GenomicFeatures/
echo "Sorting the different genomic features"
if [ $species = Human ]; then
	for i in `ls $pipe/Scripts/GenomicFeatures/`;do for a in `ls Results/DE_Genes`; do awk 'NR==FNR{a[$1]++;next} (a[$2])' $pipe/Scripts/GenomicFeatures/${i} Results/DE_Genes/${a} >Results/GenomicFeatures/${a}_${i};done ;done
elif [ $species = Mouse ]; then
	for i in `ls $pipe/Scripts/Mouse_GenomicFeatures/`;do for a in `ls Results/DE_Genes`; do awk 'NR==FNR{a[$1]++;next} (a[$2])' $pipe/Scripts/Mouse_GenomicFeatures/${i} Results/DE_Genes/${a} >Results/GenomicFeatures/${a}_${i};done ;done
fi

#rm Results/GenomicFeatures/*miRNA

mkdir Results/Gene_names
for i in `ls Results/GenomicFeatures`; do gawk '{if ($3 >0) print $2}' Results/GenomicFeatures/${i} > Results/Gene_names/Up_${i};done
for i in `ls Results/GenomicFeatures`; do gawk '{if ($3 <0) print $2}' Results/GenomicFeatures/${i} > Results/Gene_names/Down_${i};done
mkdir Results/top20
for i in `ls Results/Gene_names/`; do head -20 Results/Gene_names/${i} >Results/top20/${i};done 
mkdir Results/NormalisedMatrix
for i in `ls Results/top20/`; do awk 'NR==FNR{a[$1]++;next} (a[$2])' Results/top20/${i} Results/NormalisedCounts.tsv >Results/NormalisedMatrix/${i};done 

mv Results/Gene_names Results/GenomicFeatures/
mv Results/NormalisedMatrix Results/GenomicFeatures/
mv Results/top20 Results/GenomicFeatures/
head -1 Results/NormalisedCounts.tsv| tr '\t' '\n' >Results/GenomicFeatures/Sampelnames.txt
cd Results/GenomicFeatures/NormalisedMatrix

echo "Plotting the differentially expressed genomic features" 
find ./ -size  0 -print0 | xargs -0 rm

#################################################################################################################################

$R --vanilla --quiet --slave --args $pipe <$pipe/Scripts/PlotFeatureExpression_Log10.R
echo -e "##################### Finished plotting the various differentially expressed features #####################\n"
cd ../../../

echo "Analysis finished. Have a good day"
