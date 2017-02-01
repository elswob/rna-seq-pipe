library(edgeR)
library(RColorBrewer)
library(ggplot2)

args <- commandArgs();
species<-as.character(args[6])
if(species == "Human"){
	cat("Loading human GencodeNames\n")
	Gencode<-read.table("/share/ClusterShare/biodata/contrib/benels/gencode/human/g20/ensGene_to_geneName_default.txt",sep='\t',header=T)
}else if(species == "Mouse"){
	Gencode<-read.table("/share/ClusterShare/biodata/contrib/benels/gencode/mouse/m3/ensGene_to_geneName_default.txt",sep='\t',header=T)
	cat("Loading mouse GencodeNames\n")
}

design_file<-as.character(args[7])

count_table<-read.table("HTSeq/Combined/Allcombined", header=T, sep="\t", row.names=1)
conditions<-read.table("Conditions.txt",sep='\t')
sorterad<-as.data.frame(table(conditions$V2))
sorterad<-sort(sorterad$Freq)
noint<-rownames(count_table) %in% c("__no_feature","__ambiguous","__too_low_aQual","__not_aligned","__alignment_not_unique")

#cpms<-cpm(count_table)
#keep data with at least one read
#keep<-rowSums(x>0) >0
#keep<-rowSums(cpms>1) >=sorterad[1] & !noint
#count_table<-count_table[keep,]
#count_table_temp=count_table
#count_table_temp$Genes<-row.names(count_table)
#merged<-merge(count_table_temp,Gencode,by.x="Genes",by.y="Gene_names_ENSG")

x=count_table
x=x[!noint,]
filter <- apply(x, 1, function(x) length(x[x>5])>=2)
x <- x[filter,]

#print counts with gene names
x$Genes<-row.names(x)
x<-merge(x,Gencode,by.x="Genes",by.y="ens_id")
row.names(x)<-x$Genes
x$Genes<-x$name
x$name<-NULL
x<-cbind(ENGS_ID=rownames(x),x)
rownames(x)<-NULL
write.table(x,file="HTSeq/Combined/Allcombined_gene_names",sep="\t",quote=F,row.names=F)


expt_designEr <- data.frame(row.names = colnames(count_table),conditions=conditions$V2)

group<-expt_designEr[,1]
cds<-DGEList(count_table,group=group)

if (design_file == "n"){
	design <- model.matrix(~0+group,data=cds$samples)
}else{
	dFile<-read.delim(design_file,header=T, sep="\t")
  	dim(dFile)
  	Subject <- factor(dFile$subject)
  	Treat <- factor(dFile$treatment)
  	design <- model.matrix(~Subject+Treat)
}

print("Design for DE")
design

cds <- calcNormFactors(cds, method="upperquartile")
cool<-factor(c(conditions$V2))
uncond<-unique(conditions$V2)
cols<-brewer.pal(length(uncond),"Set1")
namn<-t(levels(conditions$V2))
namn<-namn[1:length(unique(conditions$V2))]
namn<-paste(namn,collapse=" ")
namn<-paste0("MDS plot of ",namn)
dir.create("Results")
dir.create("Results/Graphs")
mds_cords <- function(table) {
     temp <- plotMDS(table)
     temp <- cbind(temp$x, temp$y)
     colnames(temp) <- c("x","y")
     temp <- as.data.frame(temp)
     return(as.data.frame(temp))
     rm(temp)
 }
ata<-mds_cords(cds)


pdf("Results/Graphs/MDS.pdf", height=10, width=15)
ggplot(ata,aes(x, y, label=rownames(ata),colour=conditions$V2)) + geom_text(hjust=1.2)+geom_point(size=3)+ggtitle(namn)
dev.off()


pca<-princomp(cds$counts)
pdf("Results/Graphs/pca_explain.pdf", height=10, width=15)
plot(pca)
dev.off()

cds <- estimateGLMCommonDisp(cds, design, verbose=TRUE)
cds <- estimateGLMTagwiseDisp(cds, design)


pdf("Results/Graphs/BiologicalCoefficientOfVariation.pdf")
plotBCV(cds)
dev.off()

pdf("Results/Graphs/MeanVariance.pdf")
meanVarPlot <- plotMeanVar(cds, show.raw.vars=TRUE, show.tagwise.vars=TRUE, show.binned.common.disp.vars=FALSE, show.ave.raw.vars=FALSE, NBline = TRUE , nbins = 100 , pch = 16 , xlab ="Mean Expression (Log10 Scale)" , ylab = "Variance (Log10 Scale)" , main = "Mean-Variance Plot" )
dev.off()

normCounts<-cpm(cds, normalized.lib.sizes=T)

library(pheatmap)

if (design_file == "n"){
	colnames(design)<-levels(conditions$V2)
	colnames(cds$design)<-levels(conditions$V2)
}

dir.create("Results/DE_Genes")
dir.create("Results/GSEA")
dir.create("Results/Graphs/Heatmaps")
objectVec<-c(unique(levels(conditions$V2)))
#NoCond<-length(unique(conditions$V2)) 
fakeList=list()
for(i in objectVec){fakeList[[i]]=objectVec[i]}
order=combn(objectVec,2)
sampleOrder=combn(objectVec,2)
sig<-as.data.frame(t(sampleOrder))
sig$value<-0
for(j in 1:length(sampleOrder[1,])){
	Prov1<-sampleOrder[[1,j]]
	Prov2<-sampleOrder[[2,j]]
	
	E12<-exactTest(cds,pair=c(Prov1,Prov2))
	T12<- topTags(E12, n=nrow(cds$counts))$table	
	
	T12$Genes<-row.names(T12)
	T12<-merge(T12,Gencode,by.x="Genes",by.y="ens_id")
	T12<-T12[order(T12$FDR),]
	row.names(T12)<-T12$Genes
	T12$Genes<-T12$name
	T12$name<-NULL
	D12<- rownames(T12[T12$FDR<0.05,])
	T12<-T12[D12,]
	T12<-cbind(ENGS_ID=rownames(T12),T12)
	rownames(T12)<-NULL

	EttTva<-paste(Prov2,Prov1)
	EttTva<-gsub(" ","_vs_",EttTva)
	cat("Conducting", EttTva,"comparison","\n")

	sig[(sig$V1==Prov1&sig$V2==Prov2),3]=length(D12)
	if (length(D12)<= 1) {
		cat("No significant genes between", EttTva, '\n')
	} else if (length(D12)>1 ) {
		fout<-paste0("Results/Graphs/Heatmaps/",EttTva,"_Heatmap.pdf")
		pdf(fout)
		pheatmap(log(normCounts[D12,]+1), cluster_cols=T, cluster_rows=T, main=EttTva,show_rownames=F,color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(100))
		dev.off()
	
		fout<-paste0("Results/DE_Genes/EdgeR_",EttTva,"_FDR_0.05.tsv")
		write.table(T12,file=fout,sep="\t",quote=F,row.names=F)
	
	
		T12<- topTags(E12, n=nrow(cds$counts))$table
		T12$Genes<-row.names(T12)
		T12<-merge(T12,Gencode,by.x="Genes",by.y="ens_id")
		T12<-T12[order(T12$FDR),]
		row.names(T12)<-T12$Genes
		T12$Genes<-T12$name
		T12$name<-NULL
		T12<-cbind(ENGS_ID=rownames(T12),T12)
		rownames(T12)<-NULL

		fout<-paste0("Results/GSEA/",EttTva,".tsv")
		write.table(T12,file=fout,sep="\t",quote=F,row.names=F)
		}
	}

library(ggplot2)
pdf("Results/Graphs/SignificantMatrix.pdf",height=8,width=10)
ggplot(sig, aes(V1,V2, fill=value)) + geom_tile()+geom_text( vjust = 0,size=6,aes(label=value))+ggtitle("Significant genes between each condition FDR<0.05")+ylab("Samples")+xlab("Samples")+scale_fill_gradient()
dev.off()

normCounts<-as.data.frame(normCounts)
normCounts$Genes<-row.names(normCounts)
normCounts<-merge(normCounts,Gencode,by.x="Genes",by.y="ens_id")
row.names(normCounts)<-normCounts$Genes
normCounts$Genes<-normCounts$name
normCounts$name<-NULL
normCounts<-cbind(ENGS_ID=rownames(normCounts),normCounts)
rownames(normCounts)<-NULL

write.table(normCounts,file="Results/NormalisedCounts.tsv",sep="\t",quote=F,row.names=F)


