#!/bin/bash

#set the variables
name=$1
species=$2
tDir=$3
strand=$4

if [ $species = Human ]; then
	chromSizes="/share/ClusterShare/biodata/contrib/benels/genome_indices/star/human/hg38_gencode_20_ercc/chrNameLength.txt"
	gName="hg20"
elif [ $species = Mouse ]; then
	chromSizes="/share/ClusterShare/biodata/contrib/benels/genome_indices/star/mouse/GRCm38_M3/chrNameLength.txt"
	gName="mm10"
fi

#get the name of the sample from the filepath (not ideal)
echo "Creating bigwig files ..."

if [ $strand = "no" ]; then
	echo "Creating un-stranded bigwig data..."	
	if [ ! -f Alignments/$name.wig ]; then
		echo " - creating wig file...";
		WIG="/share/ClusterShare/biodata/contrib/benels/software/IGV_2.3.34/IGVTools/igvtools count $tDir/star/star.${name}Aligned.out.sorted.bam Alignments/${name}.wig $gName"
		echo $WIG && eval $WIG
		#echo $WIG
	fi

	if [ ! -f Alignments/$name.bw ]; then
		echo " - creating bigwig file...";
		BW="/share/ClusterShare/software/contrib/benels/ucsc/wigToBigWig Alignments/${name}.wig $chromSizes Alignments/${name}.bw"
		echo $BW && eval $BW
		#echo $WIG
	fi	
else
	echo "Creating stranded bigwig data..."	
	if [ ! -f Alignments/$name.wig ]; then
		echo " - creating wig file..."
		WIG="/share/ClusterShare/biodata/contrib/benels/software/IGV_2.3.34/IGVTools/igvtools count --strands first $tDir/star/star.${name}Aligned.out.sorted.bam Alignments/${name}.wig $gName"
		echo $WIG && eval $WIG
		#add negative value to negative strand
	fi
	
	if [ ! -f Alignments/$name.neg.wig ]; then
		awk '{if($3 ~ /^[0-9]/ && $2>0){print $1 "\t-" $2 "\t" $3}else{print $0}}' Alignments/${name}.wig > Alignments/${name}.neg.wig
	fi

	if [ ! -f Alignments/$name.1.bw ] || [ ! -f Alignments/$name.2.bw ]; then
		echo " - creating bigwig files...";
		BW1="/share/ClusterShare/software/contrib/benels/ucsc/wigToBigWig <(cut -f1,2 Alignments/${name}.neg.wig) $chromSizes Alignments/${name}.1.bw"
		BW2="/share/ClusterShare/software/contrib/benels/ucsc/wigToBigWig <(cut -f1,3 Alignments/${name}.neg.wig) $chromSizes Alignments/${name}.2.bw"
		echo $BW1 && eval $BW1
		echo $BW2 && eval $BW2
		#echo $WIG
	fi
fi
