#!/bin/bash

. /usr/share/Modules/init/bash

#set the variables
name=$1
cpus=$2
species=$3
tDir=$4
fastq=$5
strand=$6

if [ $species = Human ]; then
	#starDir="/share/ClusterShare/biodata/contrib/benels/transcript_indices/star/"
	starDir="/share/ClusterShare/biodata/contrib/benels/genome_indices/star/human/hg38_gencode_20_ercc"

elif [ $species = Mouse ]; then
	starDir="/share/ClusterShare/biodata/contrib/benels/genome_indices/star/mouse/GRCm38_M3"
fi

if [ ! -d $tDir/rsem ];  then
	mkdir $tDir/rsem;
fi

if [ ! -d $tDir/star ];  then
	mkdir $tDir/star;
fi

if [ ! -d star ];  then
	mkdir star
fi

#clean up any old data as star fails otherwise
rm -r $tDir/star/star.$name

#align the reads
echo -e "\n###################### Aligning the reads ######################"

echo "name = $name"
echo "cpus = $cpus"
echo "species = $species"
echo "tDir = $tDir"
echo "fastq = $fastq"
echo "strand = $strand"

#ALIGN="/share/ClusterShare/software/contrib/gi/star/2.3.1z4/STAR \
# ALIGN="/home/benels/software/rna-star/STAR-STAR_2.4.0b/STARstatic \
# 	--genomeDir $starDir \
# 	--readFilesIn $fastq/$name/*R1*.gz $fastq/$name/*R2*.gz \
# 	--outSAMattributes NH   HI      \
# 	--outFilterMultimapNmax 20   \
# 	--outFilterMismatchNmax 999   \
# 	--outFilterMismatchNoverLmax 0.04   \
# 	--alignIntronMin 20   \
# 	--alignIntronMax 1000000   \
# 	--alignMatesGapMax 1000000   \
# 	--alignSJoverhangMin 8   \
# 	--alignSJDBoverhangMin 1 \
# 	--quantMode TranscriptomeSAM \
# 	--runThreadN $cpus  \
# 	--readFilesCommand zcat"
	
#ALIGN="/home/benels/software/rna-star/STAR_2.3.1z12/STAR/STARstatic --outTmpDir $tDir/star.$name --genomeDir $starDir --readFilesIn $fastq/$name/*R1*.gz $fastq/$name/*R2*.gz --outSAMattributes NH HI AS NM MD --outFilterMultimapNmax 20 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --readFilesCommand zcat --outWigType bedGraph --runThreadN $cpus --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM --outWigStrand UnStranded --outSAMheaderCommentFile commentsENCODElong.txt --outSAMheaderHD @HD VN:1.4 --outFileNamePrefix $tDir/star.$name"

#transcripts only
#ALIGN="/share/ClusterShare/software/contrib/gi/star/2.3.1z4/STAR --alignIntronMax 1 --alignIntronMin 2  --scoreDelOpen -10000  --scoreInsOpen -10000 --alignEndsType EndToEnd --genomeDir $starDir --readFilesIn $fastq/$name/*R1*.gz $fastq/$name/*R2*.gz  --readFilesCommand zcat --outFileNamePrefix $tDir/star.$name --runThreadN $cpus"

star="/share/ClusterShare/software/contrib/nenbar/star/2.4.0d/"

#transcripts and genome (https://groups.google.com/forum/#!msg/rna-star/0C3zNg70uuI/LgM407aWSCsJ)
if [ $strand = "no" ]; then 
		echo "Running non-stranded"
		ALIGN="$star/STARstatic \
		--outTmpDir $tDir/star/star.$name \
		--genomeDir $starDir \
		--readFilesIn $fastq/$name/*R1*.gz $fastq/$name/*R2*.gz \
		--readFilesCommand zcat \
		--outFilterType BySJout \
		--outSAMattributes NH HI AS NM MD \
		--outFilterMultimapNmax 20 \
		--outFilterMismatchNmax 999 \
		--outFilterMismatchNoverReadLmax 0.04 \
		--alignIntronMin 20 \
		--alignIntronMax 1000000 \
		--alignMatesGapMax 1000000 \
		--alignSJoverhangMin 8 \
		--alignSJDBoverhangMin 1  \
		--runThreadN $cpus  \
		--quantMode TranscriptomeSAM \
		--genomeLoad NoSharedMemory \
		--outSAMstrandField intronMotif \
		--outFileNamePrefix $tDir/star/star.$name"
else
		echo "Running stranded"
		ALIGN="$star/STARstatic \
		--outTmpDir $tDir/star/star.$name \
		--genomeDir $starDir \
		--readFilesIn $fastq/$name/*R1*.gz $fastq/$name/*R2*.gz \
		--readFilesCommand zcat \
		--outFilterType BySJout \		
		--outSAMattributes All \
		--outFilterMultimapNmax 20 \
		--outFilterMismatchNmax 999 \
		--outFilterMismatchNoverReadLmax 0.04 \
		--alignIntronMin 20 \
		--alignIntronMax 1000000 \
		--alignMatesGapMax 1000000 \
		--alignSJoverhangMin 8 \
		--alignSJDBoverhangMin 1  \
		--runThreadN $cpus  \
		--quantMode TranscriptomeSAM \
		--genomeLoad NoSharedMemory \
		--outFileNamePrefix $tDir/star/star.$name"
fi		
	#--outSAMattributes NH HI AS NM MD \
	#--outWigType bedGraph \
	#--outWigStrand Unstranded \
	#--outSAMtype BAM \
	#--outSAMunmapped Within \
	#--quantMode TranscriptomeSAM \
	#--genomeLoad LoadAndRemove"

if [ ! -f $tDir/star/star.${i}Aligned.out.sam ]; then
	echo $ALIGN
	eval $ALIGN
	#get the star alignment progress files from tmp
	echo "mv $tDir/star/star.${name}Log.progress.out star/"
	mv $tDir/star/star.${name}Log.progress.out star/
	mv $tDir/star/star.${name}Log.final.out star/
fi

module load gi/novosort/precompiled/1.03.01

#convert genomic sam alignment to bam
echo -e "\n###################### Converting genomic alignment from sam to bam ######################"
if [ ! -f $tDir/star/star.${i}Aligned.out.bam ]; then
	#samtools view -bS $tDir/star/star.${name}Aligned.out.sam  > $tDir/star/star.${name}Aligned.out.bam
	/share/ClusterShare/software/contrib/gi/samtools/0.1.19/bin/samtools view -@ $cpus -bS $tDir/star/star.${name}Aligned.out.sam  > $tDir/star/star.${name}Aligned.out.bam
	rm $tDir/star/star.${name}Aligned.out.sam
fi

#sort and index genome alignment
echo " - sorting and creating bam index for genome data ..."
if [ ! -f $tDir/star/star.${name}Aligned.out.sorted.bam.bai ]; then
	#SORT="samtools sort $tDir/star/star.${name}Aligned.out.bam $tDir/star/star.${name}Aligned.out.sorted"
	SORT="novosort -t $tDir/star/ -c $cpus -m 16G $tDir/star/star.${name}Aligned.out.bam > $tDir/star/star.${name}Aligned.out.sorted.bam"
	echo $SORT && eval $SORT 
	#rm $tDir/star/star.${name}Aligned.out.bam
	INDEX="samtools index $tDir/star/star.${name}Aligned.out.sorted.bam"
	echo $INDEX && eval $INDEX
	#echo $INDEX
fi

#sort and filter the transcript alignment
if [ ! -f $tDir/star/star.${name}Aligned.toTranscriptome.out.sorted.filtered.bam ]; then
	#SAM="samtools sort -n $tDir/star/star.${name}Aligned.toTranscriptome.out.bam $tDir/star/star.${name}Aligned.toTranscriptome.out.sorted"
	SAM="novosort -t $tDir/star/ -n $tDir/star/star.${name}Aligned.toTranscriptome.out.bam -c $cpus -m 16G > $tDir/star/star.${name}Aligned.toTranscriptome.out.sorted.bam"
	echo $SAM
	eval $SAM	
	rm $tDir/star/star.${name}Aligned.toTranscriptome.out.bam
	FILTER="/share/ClusterShare/software/contrib/gi/samtools/0.1.19/bin/samtools view -@ $cpus -b -f 3 $tDir/star/star.${name}Aligned.toTranscriptome.out.sorted.bam > $tDir/star/star.${name}Aligned.toTranscriptome.out.sorted.filtered.bam"
	echo $FILTER
	eval $FILTER
	rm $tDir/star/star.${name}Aligned.toTranscriptome.out.sorted.bam
fi
