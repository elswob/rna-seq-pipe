. /usr/share/Modules/init/bash
#!/bin/bash

module unload fabbus/python/2.7.3
module load acml/gfortran64_mp/4.4.0
module load gi/samtools/0.1.19
module load gi/bedtools/2.17.0
module load jesmaa/featureCounts/1.3.6-p1
module load hugfre/HTSeq/0.5.4p3

tDir=$3
strand=$4

if [ $2 = Human ]; then
	Anno="/share/ClusterShare/biodata//contrib/nenbar/genomes/human/hg38_ercc/gencode_ercc.v20.annotation.gtf"
elif [ $2 = Mouse ]; then
	#Anno="/home/benels/data/mouse/gencode_m2/gencode.vM2.annotation.gtf.edit" 
	Anno="/share/ClusterShare/biodata/contrib/benels/gencode/mouse/m3/gencode.vM3.annotation.gtf" 
fi

if [ ! -d HTSeq ]; then
	mkdir HTSeq
fi

echo "Counting..."
if [ ! -f HTSeq/$1HTSeq ]; then
	if [ $strand = "no" ]; then
		COUNT="samtools view $tDir/star/star.${1}Aligned.out.bam | htseq-count -s no -m union -f sam -q - $Anno > HTSeq/$1HTSeq" 
	elif [ $strand = "reverse" ]; then
		COUNT="samtools view $tDir/star/star.${1}Aligned.out.bam | htseq-count -s reverse -m union -f sam -q - $Anno > HTSeq/$1HTSeq" 
	else
		COUNT="samtools view $tDir/star/star.${1}Aligned.out.bam | htseq-count -s forward -m union -f sam -q - $Anno > HTSeq/$1HTSeq" 
	fi
	echo $COUNT && eval $COUNT
fi

if [ ! -d HTSeq/Combined ]; then
	mkdir HTSeq/Combined
fi

if [ ! -f HTSeq/HTSeqGenes ]; then
	(echo Gene_names; awk '{print $1}' HTSeq/$1HTSeq) > HTSeq/HTSeqGenes
fi

#add the sample names to the counts files
(echo $1; awk '{print $2}' HTSeq/$1HTSeq) > HTSeq/Combined/$1_HTSeqValue
