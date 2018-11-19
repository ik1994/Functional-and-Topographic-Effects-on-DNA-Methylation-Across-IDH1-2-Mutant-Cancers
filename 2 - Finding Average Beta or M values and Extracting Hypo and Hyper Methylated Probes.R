############FIND AVERAGE BETA or M VALUES FOR EACH PROBE, FOR EACH SAMPLE GROUP############
##subset each list for mean betas < 0.2 and mean betas > 0.8, and for M <-2 and >2
##merge the PROBES from each list
test_M_values <- read.csv("C:/Users.csv")
test_probes <- read.csv("C:/.csv")
control_M_values <- read.csv("C:/.csv")
control_probes <- read.csv("C:/.csv")
test_samples <- merge(test_M_values, test_probes, by="ID")
control_samples <- merge(control_M_values, control_probes, by="ID")
names(test_samples)
names(control_samples)

################EXTRACTING HYPOMETHYLATED PROBES####################
####TEST SAMPLES####
##AML
MeanAML = NULL
MeanAML$AMLprobes <- test_samples[ ,1]
#MeanAML$MeanBeta <- rowMeans(test_samples[ ,26:46], na.rm=TRUE)
MeanAML$MeanM <- rowMeans(test_samples[ ,26:46], na.rm=TRUE)
MeanAML <- as.data.frame(MeanAML)
#AMLSubsettedB <- subset(MeanAML, MeanBeta<=0.2)
AMLSubsettedM <- subset(MeanAML, MeanM<=-2)
##Astro
MeanAstro = NULL
MeanAstro$Astroprobes <- test_samples[ ,1]
#MeanAstro$MeanBeta <- rowMeans(test_samples[ ,c(2:7, 47:71)], na.rm=TRUE)
MeanAstro$MeanM <- rowMeans(test_samples[ ,c(2:7, 47:71)], na.rm=TRUE)
MeanAstro <- as.data.frame(MeanAstro)
#AstroSubsettedB <- subset(MeanAstro, MeanBeta<=0.2)
AstroSubsettedM <- subset(MeanAstro, MeanM<=-2)
##Oligo
MeanOligo = NULL
MeanOligo$Oligoprobes <- test_samples[ ,1]
#MeanOligo$MeanBeta <- rowMeans(test_samples[ ,87:111], na.rm=TRUE)
MeanOligo$MeanM <- rowMeans(test_samples[ ,c(13:17, 81:95)], na.rm=TRUE)
MeanOligo <- as.data.frame(MeanOligo)
#OligoSubsettedB <- subset(MeanOligo, MeanBeta<=0.2)
OligoSubsettedM <- subset(MeanOligo, MeanM<=-2)
##Breast
MeanBreast = NULL
MeanBreast$Breastprobes <- test_samples[ ,1]
#MeanBreast$MeanBeta <- rowMeans(test_samples[ ,77:81], na.rm=TRUE)
MeanBreast$MeanM <- rowMeans(test_samples[ ,8:12], na.rm=TRUE)
MeanBreast <- as.data.frame(MeanBreast)
#BreastSubsettedB <- subset(MeanBreast, MeanBeta<=0.2)
BreastSubsettedM <- subset(MeanBreast, MeanM<=-2)
##SNUC
MeanSNUC = NULL
MeanSNUC$SNUCprobes <- test_samples[ ,1]
#MeanSNUC$MeanBeta <- rowMeans(test_samples[ ,112:119], na.rm=TRUE)
MeanSNUC$MeanM <- rowMeans(test_samples[ ,18:25], na.rm=TRUE)
MeanSNUC <- as.data.frame(MeanSNUC)
#SNUCSubsettedB <- subset(MeanSNUC, MeanBeta<=0.2)
SNUCSubsettedM <- subset(MeanSNUC, MeanM<=-2)
##Cholangio
MeanChol = NULL
MeanChol$Cholprobes <- test_samples[ ,1]
#MeanChol$MeanBeta <- rowMeans(test_samples[ ,2:13], na.rm=TRUE)
MeanChol$MeanM <- rowMeans(test_samples[ ,72:80], na.rm=TRUE)
MeanChol <- as.data.frame(MeanChol)
#CholSubsettedB <- subset(MeanChol, MeanBeta<=0.2)
CholSubsettedM <- subset(MeanChol, MeanM<=-2)

##Subset the means to create a separate list
#install.packages("rowr")
library(rowr)
subsetted_probes_list<-cbind.fill(AMLSubsettedM[ ,1], AstroSubsettedM[ ,1], BreastSubsettedM[ ,1], CholSubsettedM[ ,1], OligoSubsettedM[ ,1], SNUCSubsettedM[ ,1], fill="")
MeanM_probes_list<-cbind.fill(AMLSubsettedM[ ,2], AstroSubsettedM[ ,2], BreastSubsettedM[ ,2], CholSubsettedM[ ,2], OligoSubsettedM[ ,2], SNUCSubsettedM[ ,2], fill="")
setwd("C:/")
write.csv(subsetted_probes_list, file=".csv")
write.csv(MeanM_probes_list, file="t.csv")
##manually rename columns





################EXTRACTING HYPERMETHYLATED PROBES####################
####TEST SAMPLES####
##AML
MeanAML = NULL
MeanAML$AMLprobes <- test_samples[ ,1]
#MeanAML$MeanBeta <- rowMeans(test_samples[ ,26:46], na.rm=TRUE)
MeanAML$MeanM <- rowMeans(test_samples[ ,26:46], na.rm=TRUE)
MeanAML <- as.data.frame(MeanAML)
#AMLSubsettedB <- subset(MeanAML, MeanBeta>=0.8)
AMLSubsettedM <- subset(MeanAML, MeanM>=2)
##Astro
MeanAstro = NULL
MeanAstro$Astroprobes <- test_samples[ ,1]
#MeanAstro$MeanBeta <- rowMeans(test_samples[ ,c(2:7, 47:71)], na.rm=TRUE)
MeanAstro$MeanM <- rowMeans(test_samples[ ,c(2:7, 47:71)], na.rm=TRUE)
MeanAstro <- as.data.frame(MeanAstro)
#AstroSubsettedB <- subset(MeanAstro, MeanBeta>=0.8)
AstroSubsettedM <- subset(MeanAstro, MeanM>=2)
##Oligo
MeanOligo = NULL
MeanOligo$Oligoprobes <- test_samples[ ,1]
#MeanOligo$MeanBeta <- rowMeans(test_samples[ ,87:111], na.rm=TRUE)
MeanOligo$MeanM <- rowMeans(test_samples[ ,c(13:17, 81:95)], na.rm=TRUE)
MeanOligo <- as.data.frame(MeanOligo)
#OligoSubsettedB <- subset(MeanOligo, MeanBeta>=0.8)
OligoSubsettedM <- subset(MeanOligo, MeanM>=2)
##Breast
MeanBreast = NULL
MeanBreast$Breastprobes <- test_samples[ ,1]
#MeanBreast$MeanBeta <- rowMeans(test_samples[ ,77:81], na.rm=TRUE)
MeanBreast$MeanM <- rowMeans(test_samples[ ,8:12], na.rm=TRUE)
MeanBreast <- as.data.frame(MeanBreast)
#BreastSubsettedB <- subset(MeanBreast, MeanBeta>=0.8)
BreastSubsettedM <- subset(MeanBreast, MeanM>=2)
##SNUC
MeanSNUC = NULL
MeanSNUC$SNUCprobes <- test_samples[ ,1]
#MeanSNUC$MeanBeta <- rowMeans(test_samples[ ,112:119], na.rm=TRUE)
MeanSNUC$MeanM <- rowMeans(test_samples[ ,18:25], na.rm=TRUE)
MeanSNUC <- as.data.frame(MeanSNUC)
#SNUCSubsettedB <- subset(MeanSNUC, MeanBeta>=0.8)
SNUCSubsettedM <- subset(MeanSNUC, MeanM>=2)
##Cholangio
MeanChol = NULL
MeanChol$Cholprobes <- test_samples[ ,1]
#MeanChol$MeanBeta <- rowMeans(test_samples[ ,2:13], na.rm=TRUE)
MeanChol$MeanM <- rowMeans(test_samples[ ,72:80], na.rm=TRUE)
MeanChol <- as.data.frame(MeanChol)
#CholSubsettedB <- subset(MeanChol, MeanBeta>=0.8)
CholSubsettedM <- subset(MeanChol, MeanM>=2)

##Subset the means to create a separate list
#install.packages("rowr")
library(rowr)
subsetted_probes_list<-cbind.fill(AMLSubsettedM[ ,1], AstroSubsettedM[ ,1], BreastSubsettedM[ ,1], CholSubsettedM[ ,1], OligoSubsettedM[ ,1], SNUCSubsettedM[ ,1], fill="")
MeanM_probes_list<-cbind.fill(AMLSubsettedM[ ,2], AstroSubsettedM[ ,2], BreastSubsettedM[ ,2], CholSubsettedM[ ,2], OligoSubsettedM[ ,2], SNUCSubsettedM[ ,2], fill="")
setwd("C:/")
write.csv(subsetted_probes_list, file=".csv")
write.csv(MeanM_probes_list, file=".csv")
##manually rename columns
