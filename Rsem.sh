#!/bin/bash

#set the variables
name=$1
cpus=$2
tDir=$3
species=$4
strand=$5

if [ ! -d rsem_data ];  then
	mkdir rsem_data
fi

if [ ! -d $tDir/rsem ];  then
	mkdir $tDir/rsem
fi


if [ $species = Human ]; then
	rsemDir="/share/ClusterShare/biodata/contrib/benels/transcript_indices/rsem/human/Hg38_gencode_20_ercc-polyA/Hg38_gencode_20_ercc_no-polyA"
elif [ $species = Mouse ]; then
	rsemDir="/share/ClusterShare/biodata/contrib/benels/transcript_indices/rsem/mouse/M3-polyA/GRCm38_m3_no-polyA"
fi

rsem="/share/ClusterShare/biodata/contrib/benels/software/rsem/rsem-1.2.18"
#rsem="/home/nenbar/local/lib/rsem-1.2.18"

#calculate expression
echo -e "\n###################### Calculating expression ######################"
if [ $strand = "no" ]; then
	COUNT="$rsem/rsem-calculate-expression \
		--no-bam-output \
		--forward-prob 0.5 \
		-p $cpus \
		--temporary-folder $tDir/rsem/rsem.$name \
		--paired-end \
		--keep-intermediate-files \
		--bam $tDir/star/star.${name}Aligned.toTranscriptome.out.sorted.filtered.bam \
		$rsemDir \
		$tDir/rsem/$name"
elif [ $strand = "reverse" ]; then
	COUNT="$rsem/rsem-calculate-expression \
		--no-bam-output \
		--forward-prob 0 \
		--keep-intermediate-files \
		-p $cpus \
		--temporary-folder $tDir/rsem/rsem.$name \
		--paired-end \
		--bam $tDir/star/star.${name}Aligned.toTranscriptome.out.sorted.filtered.bam \
		$rsemDir \
		$tDir/rsem/$name"	
else
	COUNT="$rsem/rsem-calculate-expression \
		--no-bam-output \
		--forward-prob 1 \
		-p $cpus \
		--keep-intermediate-files \
		--temporary-folder $tDir/rsem/rsem.$name \
		--paired-end \
		--bam $tDir/star/star.${name}Aligned.toTranscriptome.out.sorted.filtered.bam \
		$rsemDir \
		$tDir/rsem/$name"	
fi
echo $COUNT && eval $COUNT
#generate rsem stats
rsem-plot-model $tDir/rsem/$name rsem_data/$name.rsem.pdf 
