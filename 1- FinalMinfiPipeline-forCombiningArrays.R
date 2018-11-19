#!/usr/bin/env Rscript

source("http://bioconductor.org/biocLite.R")
#biocUpgrade()
#biocLite("minfi")
#biocLite("IlluminaHumanMethylation450kanno.ilmn12.hg19")
#biocLite("IlluminaHumanMethylation450kmanifest")
#biocLite("IlluminaHumanMethylationEPICmanifest")
#biocLite("IlluminaHumanMethylationEPICanno.ilm10b2.hg19")
#biocLite("sva")
#biocLite("limma")
library("minfi")
library("IlluminaHumanMethylationEPICmanifest")
library("IlluminaHumanMethylationEPICanno.ilm10b2.hg19")
library("Biobase")
library("RColorBrewer")
library("limma")
library(reshape2)
library(ggplot2)
library(gplots)

###Reading input IDAT SampleSheet###
setwd("C:\\Users\\Ramona Bledea\\Documents\\PROJECT RESULTS\\IDH Tumor Project\\Samples and Raw Data\\Test Samples\\Idats")
baseDir <- getwd()
outDir <- getwd()

##Get IDH mutation status from and the sub-sample sheets for 450k and 850k. You can merge those data as follows with the minfi-package:
targets_450k <- read.metharray.sheet(baseDir, "samplesheet-analysis-450K.csv", recursive=TRUE)
targets_850k <- read.metharray.sheet(baseDir, "samplesheet-analysis-850K.csv", recursive=TRUE)
RGSet_450k <- read.metharray.exp(targets=targets_450k, recursive=TRUE, force=TRUE)
RGSet_850k <- read.metharray.exp(targets=targets_850k, recursive=TRUE, force=TRUE)
# Merge Datasets
RGSet <- combineArrays(RGSet_850k,RGSet_450k,outType="IlluminaHumanMethylation450k")
targets <- read.metharray.sheet(baseDir, "all-samples-order-matched.csv", recursive=TRUE)
rownames(targets) <- targets$Basename
targets <- targets[pData(RGSet)$Basename,]
# The targets object then contains the info for all samples about the mutation status.


###Read IDAT intensity data as RGChannelSet (Red/Green Channel set)
phenoData <- pData(RGSet)
manifest <- getManifest(RGSet)

#####QC and plots####
qcReportpdf = qcReport(RGSet,sampNames = targets$Sample_Name,sampGroups=targets$Sample_Group,pdf="qcReport.pdf")
qcReportpdf

Mset.illumina <- preprocessIllumina(RGSet, bg.correct=TRUE, normalize="controls") ##After Norm
qc <- getQC(Mset.illumina)
png(file=paste(outDir,"/QC.png", sep=""),width=2048,height=2048,pointsize=50)
plotQC(qc)
dev.off()

###Removing poor quality samples###
detP <- detectionP(RGSet)
colnames(detP) <- RGSet$Sample_Group

# examine mean detection p-values across all samples to identify any failed samples
pal <- brewer.pal(8,"Dark2")
par(mfrow=c(1,2))
png(file=paste(outDir,"/pvalue_sample_filter.png", sep=""),width=2048,height=2048,pointsize=70)
barplot(colMeans(detP), col=pal[factor(targets$Sample_Group)], las=2, 
        cex.names=1, cex.axis=1, names.arg=targets$Sample_Name, ylab="Mean Detection p-values")
abline(h=0.05,col="red")
legend("topleft", legend=levels(factor(targets$Sample_Group)), fill=pal,
       bg="white")
dev.off()

# remove poor quality samples
keep <- colMeans(detP) < 0.05
RGSet <- RGSet[,keep]
RGSet

# remove poor quality samples from targets data
targets <- targets[keep,]
detP <- detP[,keep]
dim(detP)
mSetSq <- preprocessQuantile(RGSet)

# ensure probes are in the same order in the mSetSq and detP objects
detP <- detP[match(featureNames(mSetSq),rownames(detP)),] 

# remove any probes that have failed in one or more samples
keep <- rowSums(detP < 0.01) == ncol(mSetSq) 
table(keep)
gset.funnorm <- mSetSq[keep,]
dim(gset.funnorm)

###Normalisation###
gset.funnorm <- addSnpInfo(gset.funnorm) ##add the genomic ranges info to gset
gset.funnorm <- dropLociWithSnps(gset.funnorm,snps=c("SBE", "CpG"), maf=0) ##drop the loci which has snps

###Get annotation###
annot = getAnnotation(gset.funnorm)

###Remove sex probes###
sex_probes = annot$Name[annot$chr %in% c("chrX", "chrY")]
gset.funnorm = gset.funnorm[!(rownames(gset.funnorm) %in% sex_probes),]
dim(gset.funnorm)

## Subset the conserved mouse probes (from paper) from EPIC annot ##
#annot <- read.csv("differential_methylation_withAnnotation.csv")
#mouse_probes <- read.delim("Mouse_Probes_EPIC.csv",header = T, sep =",")
#annot_mouse_probes <- subset(annot, Name %in% mouse_probes$EPICProbeID)
#annot_mouse_probes <-read.csv("probes_from_EPIC.csv")
#dim(annot_mouse_probes)
#gset.funnorm = gset.funnorm[rownames(gset.funnorm) %in% annot_mouse_probes$Name,]
#dim(gset.funnorm)

###Subset Probes by Gene of Interest AND Region of Interest###
#DRD1Bodya <- annot[annot$UCSC_RefGene_Name == "DRD1" & annot$UCSC_RefGene_Group == "Body",]
#DRD1Bodyb <- annot[annot$UCSC_RefGene_Name == "DRD1;DRD1" & annot$UCSC_RefGene_Group == "5'UTR;1stExon",]
#DRD1Body_probes <- rbind(DRD1Bodya, DRD1Bodyb)
#gset.funnorm = gset.funnorm[rownames(gset.funnorm) %in% DRD1Body_probes$Name,]

###Get beta values to obtain the data matrix###
#gset.funnorm.beta <- getBeta(gset.funnorm)
#colnames(gset.funnorm.beta) <- gset.funnorm$Sample_Group

###Get M values to obtain the data matrix###
gset.funnorm.beta <- getM(gset.funnorm)
colnames(gset.funnorm.beta) <- gset.funnorm$Sample_Group

#MDSplots
pal <- brewer.pal(8,"Dark2")
png(file=paste(outDir,"/MDSplotWithSampleNames-bypatients.png", sep=""),width=2048,height=2048,pointsize=60)
plotMDS(gset.funnorm.beta, top=1000, labels = targets$Sample_Name, gene.selection="common",col=pal[factor(targets$Sample_Group)])
legend("topleft", legend=levels(factor(targets$Sample_Group)), text.col=pal,bg="white", cex=1.2)
dev.off()

pal <- brewer.pal(8,"Dark2")
png(file=paste(outDir,"/MDSplotTop1000.png", sep=""),width=2048,height=2048,pointsize=50)
plotMDS(gset.funnorm.beta, top=1000, gene.selection="common",col=pal[factor(targets$Sample_Group)],pch=16)
legend("topleft", legend=levels(factor(targets$Sample_Group)), text.col=pal,bg="white", cex=1.5)
dev.off()

write.csv(gset.funnorm.beta,file="beta.csv") ##write to a csv the beta value datamatrix

###Differential Methylation###
condition <- pData(gset.funnorm)$Sample_Group #provide the sample groups for diff.methylation (if all groups)
dmp <- dmpFinder(gset.funnorm.beta, pheno=condition, type="categorical") #run dmp on the data
dmp <- cbind(dmp, ID=rownames(dmp))
#write.csv(dmp,file="differential_methylation.csv")

#annotation <- getAnnotation(RGSet) #get the probe annotation for reference
dmp_annot_combined <- cbind(annot[row.names(dmp),],dmp) #combine the dmp pvalue data with the annotation.
annot<-write.csv(dmp_annot_combined,file="differential_methylation_withAnnotation.csv")

#annotation <- getAnnotation(RGSet) #get the probe annotation for reference
dmp_annot_combined <- cbind(annot[row.names(dmp),],dmp) #combine the dmp pvalue data with the annotation.

###Plot the heatmap###
cell_colors = colorRampPalette( c("#010F57", "#010F57", "#FAFAFA", "#B21212", "#B21212") )(300)
f <- factor(targets$Sample_Group)
#Can include colnames as well (remove labCol) and cexCol = 0.5
#how do we get these to reflect beta values instead of z score?
png(file=paste(outDir,"/Heatmap(Top 50).png", sep=""),width=2050,height=2048,pointsize=50)
heatmap.2(gset.funnorm.beta[row.names(dmp[1:50,]),],trace = 'none',key.title="Methylation", labRow=FALSE, labCol = targets$Sample_Name, cexCol = 1.2, scale = 'column',col = cell_colors,cexRow = 0.49,key.xlab = "BetaValue",main = "Heatmap (Top 50)",ColSideColors = pal[factor(targets$Sample_Group)])
legend("topright", legend=levels(f), col=pal[factor(levels(f))], pch=15,cex = 0.7)
dev.off()
png(file=paste(outDir,"/Heatmap(Top 100).png", sep=""),width=2050,height=2048,pointsize=50)
heatmap.2(gset.funnorm.beta[row.names(dmp[1:100,]),],trace = 'none',key.title="Methylation",labRow=FALSE, labCol = targets$Sample_Name, cexCol = 1.2, scale = 'column',col = cell_colors,cexRow = 0.35,key.xlab = "BetaValue",main = "Heatmap (Top 100)",ColSideColors = pal[factor(targets$Sample_Group)])
legend("topright", legend=levels(f), col=pal[factor(levels(f))], pch=15,cex = 0.7)
dev.off()
png(file=paste(outDir,"/Heatmap(Top 1000).png", sep=""),width=2050,height=2048,pointsize=50)
heatmap.2(gset.funnorm.beta[row.names(dmp[1:1000,]),],trace = 'none',key.title="Methylation", labRow=FALSE, labCol = targets$Sample_Name, cexCol = 1.2, scale = 'column',col = cell_colors,cexRow = 0.52,key.xlab = "BetaValue",main = "Heatmap (Top 1000)",ColSideColors = pal[factor(targets$Sample_Group)])
legend("topright", legend=levels(f), col=pal[factor(levels(f))], pch=15,cex = 0.7)
dev.off()
png(file=paste(outDir,"/Heatmap(Top 5000) - IDH Mut, legend bottomright.png", sep=""),width=2050,height=2048,pointsize=50)
heatmap.2(gset.funnorm.beta[row.names(dmp[1:5000,]),],trace = 'none',key.title="Methylation", labRow=FALSE, labCol=FALSE, cexCol = 1.2, scale = 'column',col = cell_colors,cexRow = 0.52,key.xlab = "Z-Score",main = "Heatmap (Top 5000)",ColSideColors = pal[factor(targets$Sample_Group)])
legend("bottomright", legend=levels(f), col=pal[factor(levels(f))], pch=15,cex = 0.5)
dev.off()
png(file=paste(outDir,"/Heatmap(Top 10000) - IDH Mut, legend topleft.png", sep=""),width=2050,height=2048,pointsize=50)
heatmap.2(gset.funnorm.beta[row.names(dmp[1:10000,]),],trace = 'none',key.title="Methylation", labRow = FALSE, labCol=FALSE, cexCol = 1.2, scale = 'column',col = cell_colors,cexRow = 0.52,key.xlab = "Z-Score",main = "Heatmap (Top 10000)",ColSideColors = pal[factor(targets$Sample_Group)])
legend("topleft", legend=levels(f), col=pal[factor(levels(f))], pch=15,cex = 0.5)
dev.off()





























###Plot Boxplot and Violin plots of global methylation###
df <- melt(gset.funnorm.beta,value.name = "Methylation")
colnames(df) <- c("Probes","Samples","Methylation")
png(file=paste(outDir,"/GlobalMethylation_Boxplot.png", sep=""),width=800,height=800,pointsize=50)
ggplot(aes(y=Methylation,x=Samples,fill=Samples),data = df)+geom_boxplot()
dev.off()

png(file=paste(outDir,"/GlobalMethylation_Violinplot.png", sep=""),width=800,height=800,pointsize=50)
ggplot(aes(y=Methylation,x=Samples,fill=Samples),data = df)+geom_violin()
dev.off()

###For Pathway Analysis###

#In beta file, change the first column name to cpgisland. The Illumina reference files should already have it changed.
#In beta file, rename the colums with the sample names
#****for Chiang project, subset for only the classic/pure samples in each carcinoma subtype - AC removed, and CCC 8,9,10,11 - "beta2"
#Load new beta sheet
beta <- read.csv("C:/Users/Ramona Bledea/Documents/PROJECT RESULTS/ChiangMSKCC-forMinfi/all samples results - Functional Norm/beta2.csv")

#Load Illumina Reference Annotation File
annotation850K <- read.csv("C:/Users/Ramona Bledea/Documents/PROJECT RESULTS/Illumina EPIC 850K Methylation Annotation.csv")

#Examine beta sheet to see up to what columns your beta values extend
names(beta)
names(annotation850K)

##Merge the generated beta sheet with the Illumina annotation files
mergeddata<-merge(annotation850K,beta,by='cpgisland')
names(mergeddata)

#Find the variance of the beta values for each probe across all samples while ignoring any NAs. 
mergeddata$var <- apply(mergeddata[,48:78],1,var, na.rm=TRUE)

#Sort the beta values by the ones with the most variance
#beta_variance_sorted <- mergeddata[order(sort(mergeddata$var)),]
#^sorting by variance by column (sample), and not by rows (probes)
#this does not seem right. In the variance sorted merged data file, the cpgisland column is in numerical order and var column is not. Ended up sorting manually, by var column largest to smallest

#Use your MsigDB access to find the pathway(s) most involved with the top 100 genes sorted by maximum variance found in CpG shores from our beta data set
setwd("C:/Users/Ramona Bledea/Documents/PROJECT RESULTS/ChiangMSKCC-forMinfi/all samples results - Functional Norm")
write.csv(beta_variance_sorted, "mergeddata_variance_sorted.csv")
#manually sort by the var column
#http://software.broadinstitute.org/gsea/msigdb/annotate.jsp
#copy and paste the top 100 UCSC Gene Names into MSigDb
#check off H, C1, and C2 (bone sarcoma paper seems to only have used C2 KEGG, REACTOME, and BIOCARTA?) Any others?
#how many top genesets to use, and is the FDR q-value being below 0.05 ok?
#click computer overlaps? and then how to interpret the results
