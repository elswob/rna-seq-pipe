#! /bin/bash 

#merge_bw.sh tmp_dir species outfile chrom_sizes bigwig1 bigwig2 ...

tDir=$1
shift

species=$1
shift

out_bw=$1
shift

name=$1
shift

strand=$1
shift

files=$@

if [ $species = Human ]; then
	chromSizes="/share/ClusterShare/biodata/contrib/benels/genome_indices/star/human/hg38_gencode_20_ercc/chrNameLength.txt"
elif [ $species = Mouse ]; then
	chromSizes="/share/ClusterShare/biodata/contrib/benels/genome_indices/star/mouse/GRCm38_M3/chrNameLength.txt"
fi

#n=$(grep $name Conditions.txt | cut -f1 | wc -l)
n=$(echo $files | wc -w)

cd Alignments/
echo `pwd`
echo "files=$files"
echo "n=$n"

# calculate factor for normalization
f=$(echo "1 / $n" | bc -l)
echo "merging $n bigwig files" >&2
echo "norm factor: $f" >&2

# merge bigwig files; this creates a bedgraph file which
# is saved as a temporary file
echo "/share/ClusterShare/software/contrib/benels/ucsc/bigWigMerge $files"
/share/ClusterShare/software/contrib/benels/ucsc/bigWigMerge $files stdout | awk -v f=$f 'BEGIN{OFS="\t"}{$4=f*$4; print}' > $name.tmp || fail

#make values negative for strand 2
if [ $strand == 1 ]; then
	awk '{print $1 "\t" $2 "\t" $3 "\t-" $4}' $name.tmp > $name.tmp2
	cp $name.tmp2 $name.tmp
fi

# create new bigwig file from temporary bedgraph file
NB="/share/ClusterShare/software/contrib/benels/ucsc/bedGraphToBigWig $name.tmp ${chromSizes} ${out_bw} || fail"
echo $NB
eval $NB
#rm $name.tmp
cd ../