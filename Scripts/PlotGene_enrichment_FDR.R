args <- commandArgs();
pipeDir<-as.character(args[6])

source(paste(pipeDir,"/Scripts/multiplot.R",sep=""))
source(paste(pipeDir,"/Scripts/SummarizesData.R",sep=""))

library(reshape2)
library(ggplot2)

print("Plotting gene enrichment")

fname=dir(pattern="^[^.].*$")
fname=fname[!file.info(fname)$isdir]
#Replaces the end of filenames with "" from fnmaes
filetype=gsub("_down-|_Up|_UP|_Down|-enrichment.tsv|_FDR_0.05.tsv|BP|CC|MF|-|_EdgeR_|_FDR0.05.tsv|Downreg_|Upreg_"," ",fname)
ffname<-split(fname,filetype)
titles <- ffname
#dir.create("Graphs")
dir.create("../../Graphs/Enrichment")
for(i in 1:length(titles)) {
  titles[[i]] <- sub(names(titles)[i], paste0(names(titles)[i], " - "), titles[[i]])
}
#titles2<-gsub(names("-enrichment.tsv|-|_"titles)[i])
 titles[[i]] <- sub(names(titles)[i], paste0(names(titles)[i], " - "), titles[[i]])
for (i in 1:length(ffname)){
   data<-read.table(ffname[[i]][1], sep='\t',header=T)
   melteddata<-melt(data,1:6,na.rm=T)
#   melteddata<-melteddata[melteddata$FDR<0.05,]
   melteddata$FDR<--log10(melteddata$FDR)
   melteddata<-melteddata[1:20,]
   melteddata$Enrichment<-round(melteddata$Enrichment,digits=1)
   title<-gsub("-enrichment.tsv|FDR_0.05.tsv|_|-"," ",titles[[i]][1])
   colnames(melteddata)<-c("Enrichment_Term","Count","Size","Enrichment", "Pvalue","Minus_Log10_FDR")
   pprot1<-ggplot(data=melteddata, aes(x=Enrichment_Term,y=Minus_Log10_FDR))+geom_bar(stat="identity",colour="black", fill="blue")+theme(axis.text.x = element_text(angle = 75, hjust = 1,size=10))+xlab("Enrichment Term")+ylab("-Log10(Adjusted P-value)")+geom_text( vjust = -0.3,size=4,aes(label=Enrichment))+ggtitle(title)


   data<-read.table(ffname[[i]][2], sep='\t',header=T)
   melteddata<-melt(data,1:6,na.rm=T)
#   melteddata<-melteddata[melteddata$FDR<0.05,]
   melteddata$FDR<--log10(melteddata$FDR)
   melteddata<-melteddata[1:20,]
   melteddata$Enrichment<-round(melteddata$Enrichment,digits=1)
   title<-gsub("-enrichment.tsv|FDR_0.05.tsv|_|-"," ",titles[[i]][2])
   colnames(melteddata)<-c("Enrichment_Term","Count","Size","Enrichment", "Pvalue","Minus_Log10_FDR")
   pprot2<-ggplot(data=melteddata, aes(x=Enrichment_Term,y=Minus_Log10_FDR))+geom_bar(stat="identity",colour="black", fill="blue")+theme(axis.text.x = element_text(angle = 75, hjust = 1,size=10))+xlab("Enrichment Term")+ylab("-Log10(Adjusted P-value)")+geom_text( vjust = -0.3,size=4,aes(label=Enrichment))+ggtitle(title)


   data<-read.table(ffname[[i]][3], sep='\t',header=T)
   melteddata<-melt(data,1:6,na.rm=T)
#   melteddata<-melteddata[melteddata$FDR<0.05,]
   melteddata$FDR<--log10(melteddata$FDR)
   melteddata<-melteddata[1:20,]
   melteddata$Enrichment<-round(melteddata$Enrichment,digits=1)
   title<-gsub("-enrichment.tsv|FDR_0.05.tsv|_|-"," ",titles[[i]][3])
   colnames(melteddata)<-c("Enrichment_Term","Count","Size","Enrichment", "Pvalue","Minus_Log10_FDR")
   pprot3<-ggplot(data=melteddata, aes(x=Enrichment_Term,y=Minus_Log10_FDR))+geom_bar(stat="identity",colour="black", fill="blue")+theme(axis.text.x = element_text(angle = 75, hjust = 1,size=10))+xlab("Enrichment Term")+ylab("-Log10(Adjusted P-value)")+geom_text( vjust = -0.3,size=4,aes(label=Enrichment))+ggtitle(title)

  
   if (length(ffname[[i]])==6) {data<-read.table(ffname[[i]][4], sep='\t',header=T)
   melteddata<-melt(data,1:6,na.rm=T)
#   melteddata<-melteddata[melteddata$FDR<0.05,]
   melteddata$FDR<--log10(melteddata$FDR)
   melteddata<-melteddata[1:20,]
   melteddata$Enrichment<-round(melteddata$Enrichment,digits=1)
   title<-gsub("-enrichment.tsv|FDR_0.05.tsv|_|-"," ",titles[[i]][4])
   colnames(melteddata)<-c("Enrichment_Term","Count","Size","Enrichment", "Pvalue","Minus_Log10_FDR")
   pprot4<-ggplot(data=melteddata, aes(x=Enrichment_Term,y=Minus_Log10_FDR))+geom_bar(stat="identity",colour="black", fill="red")+theme(axis.text.x = element_text(angle = 75, hjust = 1,size=10))+xlab("Enrichment Term")+ylab("-Log10(Adjusted P-value)")+geom_text( vjust = -0.3,size=4,aes(label=Enrichment))+ggtitle(title)

   data<-read.table(ffname[[i]][5], sep='\t',header=T)
   melteddata<-melt(data,1:6,na.rm=T)
#   melteddata<-melteddata[melteddata$FDR<0.05,]
   melteddata$FDR<--log10(melteddata$FDR)
   melteddata<-melteddata[1:20,]
   melteddata$Enrichment<-round(melteddata$Enrichment,digits=1)
   title<-gsub("-enrichment.tsv|FDR_0.05.tsv|_|-"," ",titles[[i]][5])
   colnames(melteddata)<-c("Enrichment_Term","Count","Size","Enrichment", "Pvalue","Minus_Log10_FDR")
   pprot5<-ggplot(data=melteddata, aes(x=Enrichment_Term,y=Minus_Log10_FDR))+geom_bar(stat="identity",colour="black", fill="red")+theme(axis.text.x = element_text(angle = 75, hjust = 1,size=10))+xlab("Enrichment Term")+ylab("-Log10(Adjusted P-value)")+geom_text( vjust = -0.3,size=4,aes(label=Enrichment))+ggtitle(title)

   data<-read.table(ffname[[i]][6], sep='\t',header=T)
   melteddata<-melt(data,1:6,na.rm=T)
#   melteddata<-melteddata[melteddata$FDR<0.05,]
   melteddata$FDR<--log10(melteddata$FDR)
   melteddata<-melteddata[1:20,]
   melteddata$Enrichment<-round(melteddata$Enrichment,digits=1)
   title<-gsub("-enrichment.tsv|FDR_0.05.tsv|_|-"," ",titles[[i]][6])
   colnames(melteddata)<-c("Enrichment_Term","Count","Size","Enrichment", "Pvalue","Minus_Log10_FDR")
   pprot6<-ggplot(data=melteddata, aes(x=Enrichment_Term,y=Minus_Log10_FDR))+geom_bar(stat="identity",colour="black", fill="red")+theme(axis.text.x = element_text(angle = 75, hjust = 1,size=10))+xlab("Enrichment Term")+ylab("-Log10(Adjusted P-value)")+geom_text( vjust = -0.3,size=4,aes(label=Enrichment))+ggtitle(title)


  Cancer<-gsub("-BP-enrichment.tsv|_Up-BP-enrichment.tsv|_Down-BP-enrichment.tsv|Downreg_","",ffname[[i]][1])
  fout<-paste0("../../Graphs/Enrichment/",Cancer,"_enrichment.pdf")
  pdf(file=fout, height=25,width=15)
  multiplot(pprot1,pprot2,pprot3,pprot4,pprot5,pprot6, cols=2)
  dev.off()
} else if (length(ffname[[i]])==3) {
	Cancer<-gsub("-BP-enrichment.tsv|_Up-BP-enrichment.tsv|_Down-BP-enrichment.tsv|Downreg_","",ffname[[i]][1])
	  fout<-paste0("../../Graphs/Enrichment/",Cancer,"_enrichment.pdf")
	  pdf(file=fout, height=25,width=15)
	  multiplot(pprot1,pprot2,pprot3,cols=1)
	  dev.off()
} 
}
