#R --vanilla --quiet --slave <~/scripts/add_gene_name_to_rsem_gene_matrix.R
Args <- commandArgs();
gFile<-as.character(Args[6])
cat('Counts file =',gFile,'\n')
species<-as.character(Args[7])
cat('Species =',species,'\n')

x <- read.delim(gFile,row.names="ens_id")
#x <- read.delim("all.iso.results",row.names="ens_id")
#remove data with no counts
keep<-rowSums(x>0) >0
x<-x[keep,]
#create ensembl id to gene name map
if (species == "Human"){
	if(gFile == 'all.genes.results'){
		Gencode<-read.table("/share/ClusterShare/biodata/contrib/benels/gencode/human/g20/ensGene_to_geneName.txt",sep='\t',header=T)
	}else{
		Gencode<-read.table("/share/ClusterShare/biodata/contrib/benels/gencode/human/g20/ensTrans_to_geneName.txt",sep='\t',header=T)
	}
}
if (species == "Mouse"){
	if(gFile == 'all.genes.results'){
		Gencode<-read.table("/share/ClusterShare/biodata/contrib/benels/gencode/mouse/m3/ensGene_to_geneName.txt",sep='\t',header=T)
	}else{
		Gencode<-read.table("/share/ClusterShare/biodata/contrib/benels/gencode/mouse/m3/ensTrans_to_geneName.txt",sep='\t',header=T)
	}
}
#create ensembl id column
x$ens_id<-row.names(x)
#use this to add in the gene names
x<-merge(x,Gencode,by.x="ens_id",by.y="ens_id")
#create empty gene name column at start of data frame
x<-cbind(gene_name=0,x)
#add gene name data to new column
x$gene_name<-x$name
#remove old gene name column
x$name<-NULL
#write to file
if(gFile == 'all.genes.results'){
	write.table(x,file="all.genes.results.gt0.names",sep="\t",quote=F,row.names=F)
}else{
	write.table(x,file="all.iso.results.gt0.names",sep="\t",quote=F,row.names=F)
}