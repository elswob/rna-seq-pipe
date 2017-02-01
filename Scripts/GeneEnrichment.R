library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(GOstats)
library(GSEABase)
library(GO.db)
library(DBI)

args <- commandArgs();
pipeDir<-as.character(args[6])
species<-as.character(args[7])

source(paste(pipeDir,"/Scripts/GOanalysis.R",sep=""))
source(paste(pipeDir,"/Scripts/GeneSetCollection-methods.R",sep=""))
source(paste(pipeDir,"/Scripts/import_superclusters.R",sep=""))


fname=dir(pattern="^[^.].*$")
filetype=gsub("EdgeR_|Downreg_|Upreg_|_FDR_0.05.tsv"," ",fname)
ffname<-split(fname,filetype)
titles <- ffname

dir.create("GO_Analysis")

for(i in 1:length(titles)) {
  titles[[i]] <- sub(names(titles)[i], paste0(names(titles)[i], " - "), titles[[i]])
}
#titles2<-gsub(names("-enrichment.tsv|-|_"titles)[i])
 titles[[i]] <- sub(names(titles)[i], paste0(names(titles)[i], " - "), titles[[i]])

if(species=="Mouse"){
	print("Running mouse")
	for (i in 1:length(ffname)){
		gsm1<-readLines(ffname[[i]][1])
		gsmm1<-sub("[.0-9].$","",gsm1)
		gsmmm1<-sub("\\.", "", gsmm1)
		geneIdType <- ENSEMBLIdentifier("org.Mm.eg.db")
		title<-gsub("_EdgeR_|_FDR0.05.tsv|_|-","_",titles[[i]][1])
		genesetlist<-GeneSet(gsmmm1, geneIdType=geneIdType, setName=title)
		geneIdType(genesetlist) <- EntrezIdentifier("org.Mm.eg.db")
		res <- GeneSetCollection(genesetlist)
		title1<-gsub("Downreg_|EdgeR_|_FDR_0.05.tsv","",titles[[i]][1])
		outdir<-paste0("GO_Analysis/",title1)
		module.enrichment.analysis(res, test.type="BP",outdir=outdir)
		module.enrichment.analysis(res, test.type="CC",outdir=outdir)
		module.enrichment.analysis(res, test.type="MF",outdir=outdir)
	
		gsm1<-readLines(ffname[[i]][2])
		gsmm1<-sub("[.0-9].$","",gsm1)
		gsmmm1<-sub("\\.", "", gsmm1)
		geneIdType <- ENSEMBLIdentifier("org.Mm.eg.db")
		title<-gsub("_EdgeR_|_FDR0.05.tsv|_|-","_",titles[[i]][2])
		genesetlist<-GeneSet(gsmmm1, geneIdType=geneIdType, setName=title)
		geneIdType(genesetlist) <- EntrezIdentifier("org.Mm.eg.db")
		res <- GeneSetCollection(genesetlist)
		title1<-gsub("Downreg_|EdgeR_|_FDR_0.05.tsv","",titles[[i]][1])
		outdir<-paste0("GO_Analysis/",title1)
		module.enrichment.analysis(res, test.type="BP",outdir=outdir)
		module.enrichment.analysis(res, test.type="CC",outdir=outdir)
		module.enrichment.analysis(res, test.type="MF",outdir=outdir)
	}
}else{
	print("Running human")
	for (i in 1:length(ffname)){
		gsm1<-readLines(ffname[[i]][1])
		gsmm1<-sub("[.0-9].$","",gsm1)
		gsmmm1<-sub("\\.", "", gsmm1)
		geneIdType <- ENSEMBLIdentifier("org.Hs.eg.db")
		title<-gsub("_EdgeR_|_FDR0.05.tsv|_|-","_",titles[[i]][1])
		genesetlist<-GeneSet(gsmmm1, geneIdType=geneIdType, setName=title)
		geneIdType(genesetlist) <- EntrezIdentifier("org.Hs.eg.db")
		res <- GeneSetCollection(genesetlist)
		title1<-gsub("Downreg_|EdgeR_|_FDR_0.05.tsv","",titles[[i]][1])
		outdir<-paste0("GO_Analysis/",title1)
		module.enrichment.analysis(res, test.type="BP",outdir=outdir)
		module.enrichment.analysis(res, test.type="CC",outdir=outdir)
		module.enrichment.analysis(res, test.type="MF",outdir=outdir)
	
		gsm1<-readLines(ffname[[i]][2])
		gsmm1<-sub("[.0-9].$","",gsm1)
		gsmmm1<-sub("\\.", "", gsmm1)
		geneIdType <- ENSEMBLIdentifier("org.Hs.eg.db")
		title<-gsub("_EdgeR_|_FDR0.05.tsv|_|-","_",titles[[i]][2])
		genesetlist<-GeneSet(gsmmm1, geneIdType=geneIdType, setName=title)
		geneIdType(genesetlist) <- EntrezIdentifier("org.Hs.eg.db")
		res <- GeneSetCollection(genesetlist)
		title1<-gsub("Downreg_|EdgeR_|_FDR_0.05.tsv","",titles[[i]][1])
		outdir<-paste0("GO_Analysis/",title1)
		module.enrichment.analysis(res, test.type="BP",outdir=outdir)
		module.enrichment.analysis(res, test.type="CC",outdir=outdir)
		module.enrichment.analysis(res, test.type="MF",outdir=outdir)
	}
}	