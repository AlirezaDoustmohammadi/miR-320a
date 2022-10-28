library(ballgown)
library(devtools)
library(RSkittleBrewer)
library(genefilter)
library(dplyr)
library(limma)
#Import Phenotype Data
pheno_data = read.csv("phenodata.csv")


blg= ballgown(dataDir = "ballgown/", samplePattern = "SRR22015", pData=pheno_data)
#remove low-abundance genes
blg_filt = subset(blg,"rowVars(texpr(blg)) >1",genomesubset=TRUE)

#statistically signifcant differences between groups
results_transcripts = stattest(blg_filt,feature="transcript",covariate="sample.type",getFC=TRUE, meas="FPKM")
results_genes = stattest(blg_filt, feature="gene",covariate="sample.type", getFC=TRUE,meas="FPKM")
#Add gene names and gene IDs
results_transcripts =data.frame(geneNames=ballgown::geneNames(blg_filt),geneIDs=ballgown::geneIDs(blg_filt), results_transcripts)
#sort
results_transcripts = arrange(results_transcripts,pval)
results_genes = arrange(results_genes,pval)

write.csv(results_transcripts, "transcript_results.csv",row.names=FALSE)
write.csv(results_genes, "gene_results.csv",row.names=FALSE)

subset(results_transcripts,results_transcripts$qval<0.05)
subset(results_genes,results_genes$qval<0.05)

tropical= c('darkorange', 'dodgerblue','hotpink', 'limegreen', 'yellow')
palette(tropical)

fpkm = texpr(blg,meas="FPKM")
fpkm = log2(fpkm+1)
boxplot(fpkm,col=as.numeric(pheno_data$sex),las=2,ylab='log2(FPKM+1)')

#transcript 7
plot(fpkm[7,] ~ pheno_data$sample.type, border=c(1,2),main=paste(ballgown::geneNames(blg)[7],' : ',ballgown::transcriptNames(blg)[7]),pch=19, xlab="Sample Type",ylab='log2(FPKM+1)')
points(fpkm[7,] ~ jitter(as.numeric(pheno_data$sample.type)),col=as.numeric(pheno_data$sample.type))

#transcript 9
plot(fpkm[9,] ~ pheno_data$sample.type, border=c(1,2),main=paste(ballgown::geneNames(blg)[9],' : ',ballgown::transcriptNames(blg)[9]),pch=19, xlab="Sample Type",ylab='log2(FPKM+1)')
points(fpkm[9,] ~ jitter(as.numeric(pheno_data$sample.type)),col=as.numeric(pheno_data$sample.type))

#transcript 33
plot(fpkm[33,] ~ pheno_data$sample.type, border=c(1,2),main=paste(ballgown::geneNames(blg)[33],' : ',ballgown::transcriptNames(blg)[33]),pch=19, xlab="Sample Type",ylab='log2(FPKM+1)')
points(fpkm[33,] ~ jitter(as.numeric(pheno_data$sample.type)),col=as.numeric(pheno_data$sample.type))


plotMeans('MSTRG.7', blg_filt,groupvar="sample.type",legend=FALSE)
plotMeans('MSTRG.9', blg_filt,groupvar="sample.type",legend=FALSE)
plotMeans('MSTRG.33', blg_filt,groupvar="sample.type",legend=FALSE)

