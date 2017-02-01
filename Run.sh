#!/bin/bash


if [ -z $1 ]; then
	echo "Need to provide number of processors to be used e.g. 6" && exit
fi

if [ -z $2 ]; then
	echo "Write Human or Mouse" && exit
fi

if [ $2 != Human ] && [ $2 != Mouse ] ; then
	echo "Select Mouse or Human" && exit 
fi

if [ -z $3 ]; then
	echo "Provide a tmp directory" && exit
fi

if [ -z $4 ]; then
	echo "Provide the location of a directory containing the sample fastq files" && exit
fi

if [ -z $5 ]; then
	echo "Is the data strand specific? no, reverse of forward" && exit
fi

if [ $5 != no ] && [ $5 != reverse ]  && [ $5 != forward ]; then
	echo "Select no, reverse of forward for strandedness" && exit 
fi

if [ -z $6 ]; then
	echo "No DE design provided"
	design="n" 
else
	design=$6
fi

if [ ! -d logs ]; then
	mkdir logs
fi

#make name for rsem jobs using random number
r=$(( $RANDOM % 100000 ));

#set the number of cpus to use for the parallel parts, the species and the temp directory
cpus=$1
species=$2
tDir=$3
fastq=$4
strand=$5
SCRIPT=$( readlink -m $( type -p $0 ))      # Full path to script
BASE_DIR=`dirname ${SCRIPT}`                # Directory script is run in

#create Alignments directory
if [ ! -d Alignments ];  then
	mkdir Alignments/
fi

#run the aligning and counting script in parallel and multi-threaded where possible
for i in `cat Conditions.txt | cut -f1`; do
	s=$(( $RANDOM % 100000 ));
	#run star
	if [ ! -s $tDir/star/star.${i}Aligned.toTranscriptome.out.sorted.filtered.bam ]; then
		#echo "$tDir/star/star.${i}Aligned.toTranscriptome.out.sorted.bam"
		jobStats="qstat -j \`qstat | grep star.${s} | cut -d \" \" -f1\` | grep vmem"
		qsub -N star.${s} -e logs -o logs -b y -j y -cwd -pe smp $cpus -V "$BASE_DIR/Star.sh ${i} $cpus $species $tDir $fastq $strand; $jobStats" 
	fi
	#run HTSeq
	if [ ! -s HTSeq/${i}HTSeq ]; then
		#echo "Launching HTSeq job..."
		qsub -N htseq.${r} -hold_jid star.${s} -e logs -o logs -b y -j y -cwd -V bash $BASE_DIR/Htseq.sh ${i} $species $tDir $strand; 
	fi
	#create bigwigs
	if [ ! -s Alignments/${i}.bw ] && [ ! -s Alignments/${i}.1.bw ] && [ ! -s Alignments/${i}.2.bw ]; then
		jobStats="qstat -j \`qstat | grep sf.${s} | cut -d \" \" -f1\` | grep vmem" 
		qsub -hold_jid star.${s} -e logs -o logs -b y -j y -cwd -N bw.${r} -V "$BASE_DIR/bigwig.sh ${i} $species $tDir $strand; $jobStats"
	fi
	# run RSEM
	if [ ! -s $tDir/rsem/${i}.genes.results ]; then
		jobStats="qstat -j \`qstat | grep rsem.${r} | cut -d \" \" -f1\` | grep vmem"
		qsub -hold_jid star.${s} -e logs -o logs -b y -j y -cwd -pe smp $cpus -N rsem.${r} -V "$BASE_DIR/Rsem.sh ${i} $cpus $tDir $species $strand; $jobStats"
	fi
done

#Create merged bigwigs
# for i in `cat Conditions.txt | cut -f2 | sort | uniq`; do
# 	if [ $strand = "no" ]; then
# 		if [ ! -s Alignments/$i.merged.bw ]; then
# 			files=$(awk '{if ($2=="'${i}'"){print$0}}' Conditions.txt | cut -f1 | tr "\n" " " | sed 's/ /.bw /g')
# 			qsub -hold_jid bw.${r} -N bw_merge -e logs -o logs -b y -j y -cwd -V bash $BASE_DIR/bigwig_merge.sh $tDir $species $i.merged.bw $i $files;
# 		fi
# 	else
# 		if [ ! -s Alignments/$i.merged.1.bw ]; then
# 			files=$(awk '{if ($2=="'${i}'"){print$0}}' Conditions.txt | cut -f1 | tr "\n" " " | sed 's/ /.1.bw /g')
# 			qsub -hold_jid bw.${r} -N bw_merge -e logs -o logs -b y -j y -cwd -V bash $BASE_DIR/bigwig_merge.sh $tDir $species $i.merged.1.bw $i.1 1 $files;		
# 			files=$(awk '{if ($2=="'${i}'"){print$0}}' Conditions.txt | cut -f1 | tr "\n" " " | sed 's/ /.2.bw /g')
# 			qsub -hold_jid bw.${r} -N bw_merge -e logs -o logs -b y -j y -cwd -V bash $BASE_DIR/bigwig_merge.sh $tDir $species $i.merged.2.bw $i.2 2 $files;
# 		fi
# 	fi
# done

#run the HTSeq DE and analysis script as a single event 
#if [ ! -d Results ]; then
	#echo "Launching HTSeq EdgeR DE analysis job..."
	qsub -hold_jid htseq.${r} -e logs -o logs -b y -N jm_pipe_de -j y -cwd -V bash $BASE_DIR/HTSeq_DE_analysis.sh $species $BASE_DIR $design
#fi

#run the RSEM DE and analysis script as a single event 
if [ ! -s $tDir/rsem_data/all.genes.results ]; then
	jobStats="qstat -j \`qstat | grep rc.${r} | cut -d \" \" -f1\` | grep vmem"
	qsub -hold_jid rsem.${r} -e logs -o logs -b y -N rc.${r} -j y -cwd -V "$BASE_DIR/Rsem_combined.sh $BASE_DIR $species $tDir; $jobStats"
fi
