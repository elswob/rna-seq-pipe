library(reshape2)
library(ggplot2)

args <- commandArgs();
pipeDir<-as.character(args[6])

source(paste(pipeDir,"/Scripts/multiplot.R",sep=""))
source(paste(pipeDir,"/Scripts/SummarizesData.R",sep=""))

dir.create("../../Graphs/TopDEgenes")
fname=dir(pattern="^[^.].*$")
fname=fname[!file.info(fname)$isdir]
#Replaces the end of filenames with "" from fnmaes
filetype=gsub("_FDR0.05.tsv_|EdgeR_|_|miRNA|pseudo|ProteinCoding|lncRNA|_top20_NormalisedMAtrix","",fname)
ffname<-split(fname,filetype)
titles <- ffname
for(i in 1:length(titles)) {
  titles[[i]] <- sub(names(titles)[i], paste0(names(titles)[i], " - "), titles[[i]])
}

for (i in 1:length(ffname)){

	data<-read.table(ffname[[i]][1])
	colnames(data)<-readLines("../Sampelnames.txt")
	conditions<-read.table("../../../Conditions.txt")
	cond<-conditions$V2
	langd<-length(data)
	melteddata<-melt(melt(data, id=1:langd, na.rm=T))
	melteddata[5]<-melteddata$variable
	levels(melteddata$variable)<-cond
	melteddata$value<-log10(melteddata$value+1)
	colnames(melteddata)<-c("Gene_Id", "Gene_names", "Condition","value","Sample")   
	Prot<-summarySE(melteddata, measurevar="value", groupvars=c("Condition","Gene_names"))
	title<-gsub("_FDR_0.05.tsv_|_EdgeR_"," ",titles[[i]][1])
	plnc<-ggplot(data=Prot, aes(x=Gene_names, y=value, fill=Condition)) + geom_bar(stat="identity", position=position_dodge(), colour="black") +theme(axis.text.x = element_text(angle = 75, hjust = 1,size=15))+ geom_errorbar(aes(ymin=value-se, ymax=value+se),size=.4,width=.3,position=position_dodge(.9))+xlab("Gene Names")+ylab("Log10(Cpm +1)")+ggtitle(title)

	if (length(ffname[[i]])>1) {data<-read.table(ffname[[i]][2])
	colnames(data)<-readLines("../Sampelnames.txt")
	conditions<-read.table("../../../Conditions.txt")
	cond<-conditions$V2
	langd<-length(data)
	melteddata<-melt(melt(data, id=1:langd, na.rm=T))
	melteddata[5]<-melteddata$variable
	levels(melteddata$variable)<-cond
	melteddata$value<-log10(melteddata$value+1)
	colnames(melteddata)<-c("Gene_Id", "Gene_names", "Condition","value","Sample")   
	Prot<-summarySE(melteddata, measurevar="value", groupvars=c("Condition","Gene_names"))
	title<-gsub("_FDR_0.05.tsv_|_EdgeR_"," ",titles[[i]][2])
	pprot<-ggplot(data=Prot, aes(x=Gene_names, y=value, fill=Condition)) + geom_bar(stat="identity", position=position_dodge(), colour="black") +theme(axis.text.x = element_text(angle = 75, hjust = 1,size=15))+ geom_errorbar(aes(ymin=value-se, ymax=value+se),size=.4,width=.3,position=position_dodge(.9))+xlab("Gene Names")+ylab("Log10(Cpm +1)")+ggtitle(title)
} 
if (length(ffname[[i]])>2) {data<-read.table(ffname[[i]][3])
	colnames(data)<-readLines("../Sampelnames.txt")
	conditions<-read.table("../../../Conditions.txt")
	cond<-conditions$V2
	langd<-length(data)
	melteddata<-melt(melt(data, id=1:langd, na.rm=T))
	melteddata[5]<-melteddata$variable
	levels(melteddata$variable)<-cond
	melteddata$value<-log10(melteddata$value+1)
	colnames(melteddata)<-c("Gene_Id", "Gene_names", "Condition","value","Sample")   
	Prot<-summarySE(melteddata, measurevar="value", groupvars=c("Condition","Gene_names"))
	title<-gsub("_FDR_0.05.tsv_|_EdgeR_"," ",titles[[i]][3])
	pseudo<-ggplot(data=Prot, aes(x=Gene_names, y=value, fill=Condition)) + geom_bar(stat="identity", position=position_dodge(), colour="black") +theme(axis.text.x = element_text(angle = 75, hjust = 1,size=15))+ geom_errorbar(aes(ymin=value-se, ymax=value+se),size=.4,width=.3,position=position_dodge(.9))+xlab("Gene Names")+ylab("Log10(Cpm +1)")+ggtitle(title)
}
	if (length(ffname[[i]])==1) {Cancer<-gsub("EdgeR_|_lncRNA|_FDR_0.05.tsv","",ffname[[i]][1])
	fout<-paste0("../../Graphs/TopDEgenes/",Cancer,".pdf")
	pdf(file=fout, height=25,width=20)
	multiplot(plnc, cols=1)
	dev.off()
	} else if (length(ffname[[i]])==2) {Cancer<-gsub("EdgeR_|_lncRNA|_FDR_0.05.tsv","",ffname[[i]][1])
	fout<-paste0("../../Graphs/TopDEgenes/",Cancer,".pdf")
	pdf(file=fout, height=25,width=20)
	multiplot(pprot, plnc, cols=1)
	dev.off()
} else if (length(ffname[[i]])==3) {Cancer<-gsub("EdgeR_|_lncRNA|_FDR_0.05.tsv","",ffname[[i]][1])
	fout<-paste0("../../Graphs/TopDEgenes/",Cancer,".pdf")
	pdf(file=fout, height=25,width=20)
	multiplot(pprot, plnc, pseudo, cols=1)
	dev.off()
}
}


