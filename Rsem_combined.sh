#!/bin/bash

#set the variables
pipe=$1
species=$2
tDir=$3

#merge rsem files and create files with zero counts lines removed
rsem-generate-data-matrix $tDir/rsem/*.genes.results > rsem_data/all.genes.results
#add gene_id to first column
cd rsem_data
sed -i 's/^[ \t]/\"ens_id\"\t/' all.genes.results
#remove trailing integer
perl -i -plane 's/(ENS.*?)\.\d+/$1/' all.genes.results
#add gene names
R --vanilla --quiet --slave --args all.genes.results $species < $pipe/add_gene_name_to_rsem_gene_matrix.R
cd ../
#merge isoform predictions
rsem-generate-data-matrix $tDir/rsem/*.isoforms.results > rsem_data/all.iso.results
cd rsem_data
#add ens_id to first column
sed -i 's/^[ \t]/\"ens_id\"\t/' all.iso.results
#remove trailing integer
perl -i -plane 's/(ENS.*?)\.\d+/$1/' all.iso.results
#add gene names
R --vanilla --quiet --slave --args all.iso.results $species < $pipe/add_gene_name_to_rsem_gene_matrix.R 
