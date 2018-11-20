source("https://bioconductor.org/biocLite.R")
#biocLite("karyoploteR")
#biocLite("regioneR")
#biocLite("stringi")
library(stringi)
library(karyoploteR)
library(regioneR)


#####LOAD IN THE UNIQUE AND COMMON PROBE ANNOTATED FILES####
common_probes_test_hypo_annot <- read.csv("C:/.csv")
common_probes_test_hyper_annot <- read.csv("C:/.csv")
AML_hypo_unique <- read.csv("C:/.csv")
Astro_hypo_unique <- read.csv("C:/.csv")
Breast_hypo_unique <- read.csv("C:/.csv")
Chol_hypo_unique <- read.csv("C:/.csv")
SNUC_hypo_unique <- read.csv("C:/.csv")
Oligo_hypo_unique <- read.csv("C:/.csv")
AML_hyper_unique <- read.csv("C:/.csv")
Astro_hyper_unique <- read.csv("C:/.csv")
Breast_hyper_unique <- read.csv("C:/.csv")
Chol_hyper_unique <- read.csv("C:/.csv")
SNUC_hyper_unique <- read.csv("C:/.csv")
Oligo_hyper_unique <- read.csv("C:/.csv")




#####EXTRACT ENHANCERS, PROMOTERS, AND BODIES FROM EACH DATAFRAME#####
AML_hypo_unique_enhancers <- subset(AML_hypo_unique, AML_hypo_unique$Enhancer == "TRUE")
AML_hypo_unique_bodies <- AML_hypo_unique[grep("Body", AML_hypo_unique$UCSC_RefGene_Group), ]
AML_hypo_unique_promoters <- AML_hypo_unique[grep("TSS200", AML_hypo_unique$UCSC_RefGene_Group), ]

Astro_hypo_unique_enhancers <- subset(Astro_hypo_unique, Astro_hypo_unique$Enhancer == "TRUE")
Astro_hypo_unique_bodies <- Astro_hypo_unique[grep("Body", Astro_hypo_unique$UCSC_RefGene_Group), ]
Astro_hypo_unique_promoters <- Astro_hypo_unique[grep("TSS200", Astro_hypo_unique$UCSC_RefGene_Group), ]

Breast_hypo_unique_enhancers <- subset(Breast_hypo_unique, Breast_hypo_unique$Enhancer == "TRUE")
Breast_hypo_unique_bodies <- Breast_hypo_unique[grep("Body", Breast_hypo_unique$UCSC_RefGene_Group), ]
Breast_hypo_unique_promoters <- Breast_hypo_unique[grep("TSS200", Breast_hypo_unique$UCSC_RefGene_Group), ]

Chol_hypo_unique_enhancers <- subset(Chol_hypo_unique, Chol_hypo_unique$Enhancer == "TRUE")
Chol_hypo_unique_bodies <- Chol_hypo_unique[grep("Body", Chol_hypo_unique$UCSC_RefGene_Group), ]
Chol_hypo_unique_promoters <- Chol_hypo_unique[grep("TSS200", Chol_hypo_unique$UCSC_RefGene_Group), ]

Oligo_hypo_unique_enhancers <- subset(Oligo_hypo_unique, Oligo_hypo_unique$Enhancer == "TRUE")
Oligo_hypo_unique_bodies <- Oligo_hypo_unique[grep("Body", Oligo_hypo_unique$UCSC_RefGene_Group), ]
Oligo_hypo_unique_promoters <- Oligo_hypo_unique[grep("TSS200", Oligo_hypo_unique$UCSC_RefGene_Group), ]

SNUC_hypo_unique_enhancers <- subset(SNUC_hypo_unique, SNUC_hypo_unique$Enhancer == "TRUE")
SNUC_hypo_unique_bodies <- SNUC_hypo_unique[grep("Body", SNUC_hypo_unique$UCSC_RefGene_Group), ]
SNUC_hypo_unique_promoters <- SNUC_hypo_unique[grep("TSS200", SNUC_hypo_unique$UCSC_RefGene_Group), ]

common_probes_test_hypo_annot_enhancers <- subset(common_probes_test_hypo_annot, common_probes_test_hypo_annot$Enhancer == "TRUE")
common_probes_test_hypo_annot_bodies <- common_probes_test_hypo_annot[grep("Body", common_probes_test_hypo_annot$UCSC_RefGene_Group), ]
common_probes_test_hypo_annot_promoters <- common_probes_test_hypo_annot[grep("TSS200", common_probes_test_hypo_annot$UCSC_RefGene_Group), ]





AML_hyper_unique_enhancers <- subset(AML_hyper_unique, AML_hyper_unique$Enhancer == "TRUE")
AML_hyper_unique_bodies <- AML_hyper_unique[grep("Body", AML_hyper_unique$UCSC_RefGene_Group), ]
AML_hyper_unique_promoters <- AML_hyper_unique[grep("TSS200", AML_hyper_unique$UCSC_RefGene_Group), ]

Astro_hyper_unique_enhancers <- subset(Astro_hyper_unique, Astro_hyper_unique$Enhancer == "TRUE")
Astro_hyper_unique_bodies <- Astro_hyper_unique[grep("Body", Astro_hyper_unique$UCSC_RefGene_Group), ]
Astro_hyper_unique_promoters <- Astro_hyper_unique[grep("TSS200", Astro_hyper_unique$UCSC_RefGene_Group), ]

Breast_hyper_unique_enhancers <- subset(Breast_hyper_unique, Breast_hyper_unique$Enhancer == "TRUE")
Breast_hyper_unique_bodies <- Breast_hyper_unique[grep("Body", Breast_hyper_unique$UCSC_RefGene_Group), ]
Breast_hyper_unique_promoters <- Breast_hyper_unique[grep("TSS200", Breast_hyper_unique$UCSC_RefGene_Group), ]

Chol_hyper_unique_enhancers <- subset(Chol_hyper_unique, Chol_hyper_unique$Enhancer == "TRUE")
Chol_hyper_unique_bodies <- Chol_hyper_unique[grep("Body", Chol_hyper_unique$UCSC_RefGene_Group), ]
Chol_hyper_unique_promoters <- Chol_hyper_unique[grep("TSS200", Chol_hyper_unique$UCSC_RefGene_Group), ]

Oligo_hyper_unique_enhancers <- subset(Oligo_hyper_unique, Oligo_hyper_unique$Enhancer == "TRUE")
Oligo_hyper_unique_bodies <- Oligo_hyper_unique[grep("Body", Oligo_hyper_unique$UCSC_RefGene_Group), ]
Oligo_hyper_unique_promoters <- Oligo_hyper_unique[grep("TSS200", Oligo_hyper_unique$UCSC_RefGene_Group), ]

SNUC_hyper_unique_enhancers <- subset(SNUC_hyper_unique, SNUC_hyper_unique$Enhancer == "TRUE")
SNUC_hyper_unique_bodies <- SNUC_hyper_unique[grep("Body", SNUC_hyper_unique$UCSC_RefGene_Group), ]
SNUC_hyper_unique_promoters <- SNUC_hyper_unique[grep("TSS200", SNUC_hyper_unique$UCSC_RefGene_Group), ]

common_probes_test_hyper_annot_enhancers <- subset(common_probes_test_hyper_annot, common_probes_test_hyper_annot$Enhancer == "TRUE")
common_probes_test_hyper_annot_bodies <- common_probes_test_hyper_annot[grep("Body", common_probes_test_hyper_annot$UCSC_RefGene_Group), ]
common_probes_test_hyper_annot_promoters <- common_probes_test_hyper_annot[grep("TSS200", common_probes_test_hyper_annot$UCSC_RefGene_Group), ]

##Find the rowmeans for each##
##Enhancers
AML_E_hypo=NULL
AML_E_hypo$Chromosome <- AML_hypo_unique_enhancers[ ,99]
AML_E_hypo$ID <- AML_hypo_unique_enhancers[ ,2]
AML_E_hypo$POS <- AML_hypo_unique_enhancers[ ,100]
AML_E_hypo$UCSC_Gene <- AML_hypo_unique_enhancers[ ,122]
AML_E_hypo$Mean_Methylation <- rowMeans(AML_hypo_unique_enhancers[ ,28:48], na.rm=TRUE)
AML_E_hypo$Tumor_Type <- "AML"
AML_E_hypo$Gene_Region <- "Enhancer"
AML_E_hypo <- as.data.frame(AML_E_hypo)
Astro_E_hypo=NULL
Astro_E_hypo$Chromosome <- Astro_hypo_unique_enhancers[ ,99]
Astro_E_hypo$ID <- Astro_hypo_unique_enhancers[ ,2]
Astro_E_hypo$POS <- Astro_hypo_unique_enhancers[ ,100]
Astro_E_hypo$UCSC_Gene <- Astro_hypo_unique_enhancers[ ,122]
Astro_E_hypo$Mean_Methylation <- rowMeans(Astro_hypo_unique_enhancers[ ,c(4:9, 49:73)], na.rm=TRUE)
Astro_E_hypo$Tumor_Type <- "Astrocytoma"
Astro_E_hypo$Gene_Region <- "Enhancer"
Astro_E_hypo <- as.data.frame(Astro_E_hypo)
Chol_E_hypo=NULL
Chol_E_hypo$Chromosome <- Chol_hypo_unique_enhancers[ ,99]
Chol_E_hypo$ID <- Chol_hypo_unique_enhancers[ ,2]
Chol_E_hypo$POS <- Chol_hypo_unique_enhancers[ ,100]
Chol_E_hypo$UCSC_Gene <- Chol_hypo_unique_enhancers[ ,122]
Chol_E_hypo$Mean_Methylation <- rowMeans(Chol_hypo_unique_enhancers[ ,74:82], na.rm=TRUE)
Chol_E_hypo$Tumor_Type <- "Cholangiocarcinoma"
Chol_E_hypo$Gene_Region <- "Enhancer"
Chol_E_hypo <- as.data.frame(Chol_E_hypo)
Breast_E_hypo=NULL
Breast_E_hypo$Chromosome <- Breast_hypo_unique_enhancers[ ,99]
Breast_E_hypo$ID <- Breast_hypo_unique_enhancers[ ,2]
Breast_E_hypo$POS <- Breast_hypo_unique_enhancers[ ,100]
Breast_E_hypo$UCSC_Gene <- Breast_hypo_unique_enhancers[ ,122]
Breast_E_hypo$Mean_Methylation <- rowMeans(Breast_hypo_unique_enhancers[ ,10:14], na.rm=TRUE)
Breast_E_hypo$Tumor_Type <- "Breast Cancer"
Breast_E_hypo$Gene_Region <- "Enhancer"
Breast_E_hypo <- as.data.frame(Breast_E_hypo)
Oligo_E_hypo=NULL
Oligo_E_hypo$Chromosome <- Oligo_hypo_unique_enhancers[ ,99]
Oligo_E_hypo$ID <- Oligo_hypo_unique_enhancers[ ,2]
Oligo_E_hypo$POS <- Oligo_hypo_unique_enhancers[ ,100]
Oligo_E_hypo$UCSC_Gene <- Oligo_hypo_unique_enhancers[ ,122]
Oligo_E_hypo$Mean_Methylation <- rowMeans(Oligo_hypo_unique_enhancers[ ,c(15:19, 83:97)], na.rm=TRUE)
Oligo_E_hypo$Tumor_Type <- "Oligodendroglioma"
Oligo_E_hypo$Gene_Region <- "Enhancer"
Oligo_E_hypo <- as.data.frame(Oligo_E_hypo)
SNUC_E_hypo=NULL
SNUC_E_hypo$Chromosome <- SNUC_hypo_unique_enhancers[ ,99]
SNUC_E_hypo$ID <- SNUC_hypo_unique_enhancers[ ,2]
SNUC_E_hypo$POS <- SNUC_hypo_unique_enhancers[ ,100]
SNUC_E_hypo$UCSC_Gene <- SNUC_hypo_unique_enhancers[ ,122]
SNUC_E_hypo$Mean_Methylation <- rowMeans(SNUC_hypo_unique_enhancers[ ,20:27], na.rm=TRUE)
SNUC_E_hypo$Tumor_Type <- "SNUC"
SNUC_E_hypo$Gene_Region <- "Enhancer"
SNUC_E_hypo <- as.data.frame(SNUC_E_hypo)
common_test_E_hypo=NULL
common_test_E_hypo$Chromosome <- common_probes_test_hypo_annot_enhancers[ ,98]
common_test_E_hypo$ID <- common_probes_test_hypo_annot_enhancers[ ,2]
common_test_E_hypo$POS <- common_probes_test_hypo_annot_enhancers[ ,99]
common_test_E_hypo$UCSC_Gene <- common_probes_test_hypo_annot_enhancers[ ,121]
common_test_E_hypo$Mean_Methylation <- rowMeans(common_probes_test_hypo_annot_enhancers[ ,3:96], na.rm=TRUE)
common_test_E_hypo$Tumor_Type <- "Common Probes, IDH Mutant"
common_test_E_hypo$Gene_Region <- "Enhancer"
common_test_E_hypo <- as.data.frame(common_test_E_hypo)


AML_E_hyper=NULL
AML_E_hyper$Chromosome <- AML_hyper_unique_enhancers[ ,99]
AML_E_hyper$ID <- AML_hyper_unique_enhancers[ ,2]
AML_E_hyper$POS <- AML_hyper_unique_enhancers[ ,100]
AML_E_hyper$UCSC_Gene <- AML_hyper_unique_enhancers[ ,122]
AML_E_hyper$Mean_Methylation <- rowMeans(AML_hyper_unique_enhancers[ ,28:48], na.rm=TRUE)
AML_E_hyper$Tumor_Type <- "AML"
AML_E_hyper$Gene_Region <- "Enhancer"
AML_E_hyper <- as.data.frame(AML_E_hyper)
Astro_E_hyper=NULL
Astro_E_hyper$Chromosome <- Astro_hyper_unique_enhancers[ ,99]
Astro_E_hyper$ID <- Astro_hyper_unique_enhancers[ ,2]
Astro_E_hyper$POS <- Astro_hyper_unique_enhancers[ ,100]
Astro_E_hyper$UCSC_Gene <- Astro_hyper_unique_enhancers[ ,122]
Astro_E_hyper$Mean_Methylation <- rowMeans(Astro_hyper_unique_enhancers[ ,c(4:9, 49:73)], na.rm=TRUE)
Astro_E_hyper$Tumor_Type <- "Astrocytoma"
Astro_E_hyper$Gene_Region <- "Enhancer"
Astro_E_hyper <- as.data.frame(Astro_E_hyper)
Chol_E_hyper=NULL
Chol_E_hyper$Chromosome <- Chol_hyper_unique_enhancers[ ,99]
Chol_E_hyper$ID <- Chol_hyper_unique_enhancers[ ,2]
Chol_E_hyper$POS <- Chol_hyper_unique_enhancers[ ,100]
Chol_E_hyper$UCSC_Gene <- Chol_hyper_unique_enhancers[ ,122]
Chol_E_hyper$Mean_Methylation <- rowMeans(Chol_hyper_unique_enhancers[ ,74:82], na.rm=TRUE)
Chol_E_hyper$Tumor_Type <- "Cholangiocarcinoma"
Chol_E_hyper$Gene_Region <- "Enhancer"
Chol_E_hyper <- as.data.frame(Chol_E_hyper)
Breast_E_hyper=NULL
Breast_E_hyper$Chromosome <- Breast_hyper_unique_enhancers[ ,99]
Breast_E_hyper$ID <- Breast_hyper_unique_enhancers[ ,2]
Breast_E_hyper$POS <- Breast_hyper_unique_enhancers[ ,100]
Breast_E_hyper$UCSC_Gene <- Breast_hyper_unique_enhancers[ ,122]
Breast_E_hyper$Mean_Methylation <- rowMeans(Breast_hyper_unique_enhancers[ ,10:14], na.rm=TRUE)
Breast_E_hyper$Tumor_Type <- "Breast Cancer"
Breast_E_hyper$Gene_Region <- "Enhancer"
Breast_E_hyper <- as.data.frame(Breast_E_hyper)
Oligo_E_hyper=NULL
Oligo_E_hyper$Chromosome <- Oligo_hyper_unique_enhancers[ ,99]
Oligo_E_hyper$ID <- Oligo_hyper_unique_enhancers[ ,2]
Oligo_E_hyper$POS <- Oligo_hyper_unique_enhancers[ ,100]
Oligo_E_hyper$UCSC_Gene <- Oligo_hyper_unique_enhancers[ ,122]
Oligo_E_hyper$Mean_Methylation <- rowMeans(Oligo_hyper_unique_enhancers[ ,c(15:19, 83:97)], na.rm=TRUE)
Oligo_E_hyper$Tumor_Type <- "Oligodendroglioma"
Oligo_E_hyper$Gene_Region <- "Enhancer"
Oligo_E_hyper <- as.data.frame(Oligo_E_hyper)
SNUC_E_hyper=NULL
SNUC_E_hyper$Chromosome <- SNUC_hyper_unique_enhancers[ ,99]
SNUC_E_hyper$ID <- SNUC_hyper_unique_enhancers[ ,2]
SNUC_E_hyper$POS <- SNUC_hyper_unique_enhancers[ ,100]
SNUC_E_hyper$UCSC_Gene <- SNUC_hyper_unique_enhancers[ ,122]
SNUC_E_hyper$Mean_Methylation <- rowMeans(SNUC_hyper_unique_enhancers[ ,20:27], na.rm=TRUE)
SNUC_E_hyper$Tumor_Type <- "SNUC"
SNUC_E_hyper$Gene_Region <- "Enhancer"
SNUC_E_hyper <- as.data.frame(SNUC_E_hyper)
common_test_E_hyper=NULL
common_test_E_hyper$Chromosome <- common_probes_test_hyper_annot_enhancers[ ,98]
common_test_E_hyper$ID <- common_probes_test_hyper_annot_enhancers[ ,2]
common_test_E_hyper$POS <- common_probes_test_hyper_annot_enhancers[ ,99]
common_test_E_hyper$UCSC_Gene <- common_probes_test_hyper_annot_enhancers[ ,121]
common_test_E_hyper$Mean_Methylation <- rowMeans(common_probes_test_hyper_annot_enhancers[ ,3:96], na.rm=TRUE)
common_test_E_hyper$Tumor_Type <- "Common Probes, IDH Mutant"
common_test_E_hyper$Gene_Region <- "Enhancer"
common_test_E_hyper <- as.data.frame(common_test_E_hyper)

##Promoters
AML_P_hypo=NULL
AML_P_hypo$Chromosome <- AML_hypo_unique_promoters[ ,99]
AML_P_hypo$ID <- AML_hypo_unique_promoters[ ,2]
AML_P_hypo$POS <- AML_hypo_unique_promoters[ ,100]
AML_P_hypo$UCSC_Gene <- AML_hypo_unique_promoters[ ,122]
AML_P_hypo$Mean_Methylation <- rowMeans(AML_hypo_unique_promoters[ ,28:48], na.rm=TRUE)
AML_P_hypo$Tumor_Type <- "AML"
AML_P_hypo$Gene_Region <- "TSS200"
AML_P_hypo <- as.data.frame(AML_P_hypo)
Astro_P_hypo=NULL
Astro_P_hypo$Chromosome <- Astro_hypo_unique_promoters[ ,99]
Astro_P_hypo$ID <- Astro_hypo_unique_promoters[ ,2]
Astro_P_hypo$POS <- Astro_hypo_unique_promoters[ ,100]
Astro_P_hypo$UCSC_Gene <- Astro_hypo_unique_promoters[ ,122]
Astro_P_hypo$Mean_Methylation <- rowMeans(Astro_hypo_unique_promoters[ ,c(4:9, 49:73)], na.rm=TRUE)
Astro_P_hypo$Tumor_Type <- "Astrocytoma"
Astro_P_hypo$Gene_Region <- "TSS200"
Astro_P_hypo <- as.data.frame(Astro_P_hypo)
Chol_P_hypo=NULL
Chol_P_hypo$Chromosome <- Chol_hypo_unique_promoters[ ,99]
Chol_P_hypo$ID <- Chol_hypo_unique_promoters[ ,2]
Chol_P_hypo$POS <- Chol_hypo_unique_promoters[ ,100]
Chol_P_hypo$UCSC_Gene <- Chol_hypo_unique_promoters[ ,122]
Chol_P_hypo$Mean_Methylation <- rowMeans(Chol_hypo_unique_promoters[ ,74:82], na.rm=TRUE)
Chol_P_hypo$Tumor_Type <- "Cholangiocarcinoma"
Chol_P_hypo$Gene_Region <- "TSS200"
Chol_P_hypo <- as.data.frame(Chol_P_hypo)
Breast_P_hypo=NULL
Breast_P_hypo$Chromosome <- Breast_hypo_unique_promoters[ ,99]
Breast_P_hypo$ID <- Breast_hypo_unique_promoters[ ,2]
Breast_P_hypo$POS <- Breast_hypo_unique_promoters[ ,100]
Breast_P_hypo$UCSC_Gene <- Breast_hypo_unique_promoters[ ,122]
Breast_P_hypo$Mean_Methylation <- rowMeans(Breast_hypo_unique_promoters[ ,10:14], na.rm=TRUE)
Breast_P_hypo$Tumor_Type <- "Breast Cancer"
Breast_P_hypo$Gene_Region <- "TSS200"
Breast_P_hypo <- as.data.frame(Breast_P_hypo)
Oligo_P_hypo=NULL
Oligo_P_hypo$Chromosome <- Oligo_hypo_unique_promoters[ ,99]
Oligo_P_hypo$ID <- Oligo_hypo_unique_promoters[ ,2]
Oligo_P_hypo$POS <- Oligo_hypo_unique_promoters[ ,100]
Oligo_P_hypo$UCSC_Gene <- Oligo_hypo_unique_promoters[ ,122]
Oligo_P_hypo$Mean_Methylation <- rowMeans(Oligo_hypo_unique_promoters[ ,c(15:19, 83:97)], na.rm=TRUE)
Oligo_P_hypo$Tumor_Type <- "Oligodendroglioma"
Oligo_P_hypo$Gene_Region <- "TSS200"
Oligo_P_hypo <- as.data.frame(Oligo_P_hypo)
SNUC_P_hypo=NULL
SNUC_P_hypo$Chromosome <- SNUC_hypo_unique_promoters[ ,99]
SNUC_P_hypo$ID <- SNUC_hypo_unique_promoters[ ,2]
SNUC_P_hypo$POS <- SNUC_hypo_unique_promoters[ ,100]
SNUC_P_hypo$UCSC_Gene <- SNUC_hypo_unique_promoters[ ,122]
SNUC_P_hypo$Mean_Methylation <- rowMeans(SNUC_hypo_unique_promoters[ ,20:27], na.rm=TRUE)
SNUC_P_hypo$Tumor_Type <- "SNUC"
SNUC_P_hypo$Gene_Region <- "TSS200"
SNUC_P_hypo <- as.data.frame(SNUC_P_hypo)
common_test_P_hypo=NULL
common_test_P_hypo$Chromosome <- common_probes_test_hypo_annot_promoters[ ,98]
common_test_P_hypo$ID <- common_probes_test_hypo_annot_promoters[ ,2]
common_test_P_hypo$POS <- common_probes_test_hypo_annot_promoters[ ,99]
common_test_P_hypo$UCSC_Gene <- common_probes_test_hypo_annot_promoters[ ,121]
common_test_P_hypo$Mean_Methylation <- rowMeans(common_probes_test_hypo_annot_promoters[ ,3:96], na.rm=TRUE)
common_test_P_hypo$Tumor_Type <- "Common Probes, IDH Mutant"
common_test_P_hypo$Gene_Region <- "TSS200"
common_test_P_hypo <- as.data.frame(common_test_P_hypo)

AML_P_hyper=NULL
AML_P_hyper$Chromosome <- AML_hyper_unique_promoters[ ,99]
AML_P_hyper$ID <- AML_hyper_unique_promoters[ ,2]
AML_P_hyper$POS <- AML_hyper_unique_promoters[ ,100]
AML_P_hyper$UCSC_Gene <- AML_hyper_unique_promoters[ ,122]
AML_P_hyper$Mean_Methylation <- rowMeans(AML_hyper_unique_promoters[ ,28:48], na.rm=TRUE)
AML_P_hyper$Tumor_Type <- "AML"
AML_P_hyper$Gene_Region <- "TSS200"
AML_P_hyper <- as.data.frame(AML_P_hyper)
Astro_P_hyper=NULL
Astro_P_hyper$Chromosome <- Astro_hyper_unique_promoters[ ,99]
Astro_P_hyper$ID <- Astro_hyper_unique_promoters[ ,2]
Astro_P_hyper$POS <- Astro_hyper_unique_promoters[ ,100]
Astro_P_hyper$UCSC_Gene <- Astro_hyper_unique_promoters[ ,122]
Astro_P_hyper$Mean_Methylation <- rowMeans(Astro_hyper_unique_promoters[ ,c(4:9, 49:73)], na.rm=TRUE)
Astro_P_hyper$Tumor_Type <- "Astrocytoma"
Astro_P_hyper$Gene_Region <- "TSS200"
Astro_P_hyper <- as.data.frame(Astro_P_hyper)
Chol_P_hyper=NULL
Chol_P_hyper$Chromosome <- Chol_hyper_unique_promoters[ ,99]
Chol_P_hyper$ID <- Chol_hyper_unique_promoters[ ,2]
Chol_P_hyper$POS <- Chol_hyper_unique_promoters[ ,100]
Chol_P_hyper$UCSC_Gene <- Chol_hyper_unique_promoters[ ,122]
Chol_P_hyper$Mean_Methylation <- rowMeans(Chol_hyper_unique_promoters[ ,74:82], na.rm=TRUE)
Chol_P_hyper$Tumor_Type <- "Cholangiocarcinoma"
Chol_P_hyper$Gene_Region <- "TSS200"
Chol_P_hyper <- as.data.frame(Chol_P_hyper)
Breast_P_hyper=NULL
Breast_P_hyper$Chromosome <- Breast_hyper_unique_promoters[ ,99]
Breast_P_hyper$ID <- Breast_hyper_unique_promoters[ ,2]
Breast_P_hyper$POS <- Breast_hyper_unique_promoters[ ,100]
Breast_P_hyper$UCSC_Gene <- Breast_hyper_unique_promoters[ ,122]
Breast_P_hyper$Mean_Methylation <- rowMeans(Breast_hyper_unique_promoters[ ,10:14], na.rm=TRUE)
Breast_P_hyper$Tumor_Type <- "Breast Cancer"
Breast_P_hyper$Gene_Region <- "TSS200"
Breast_P_hyper <- as.data.frame(Breast_P_hyper)
Oligo_P_hyper=NULL
Oligo_P_hyper$Chromosome <- Oligo_hyper_unique_promoters[ ,99]
Oligo_P_hyper$ID <- Oligo_hyper_unique_promoters[ ,2]
Oligo_P_hyper$POS <- Oligo_hyper_unique_promoters[ ,100]
Oligo_P_hyper$UCSC_Gene <- Oligo_hyper_unique_promoters[ ,122]
Oligo_P_hyper$Mean_Methylation <- rowMeans(Oligo_hyper_unique_promoters[ ,c(15:19, 83:97)], na.rm=TRUE)
Oligo_P_hyper$Tumor_Type <- "Oligodendroglioma"
Oligo_P_hyper$Gene_Region <- "TSS200"
Oligo_P_hyper <- as.data.frame(Oligo_P_hyper)
SNUC_P_hyper=NULL
SNUC_P_hyper$Chromosome <- SNUC_hyper_unique_promoters[ ,99]
SNUC_P_hyper$ID <- SNUC_hyper_unique_promoters[ ,2]
SNUC_P_hyper$POS <- SNUC_hyper_unique_promoters[ ,100]
SNUC_P_hyper$UCSC_Gene <- SNUC_hyper_unique_promoters[ ,122]
SNUC_P_hyper$Mean_Methylation <- rowMeans(SNUC_hyper_unique_promoters[ ,20:27], na.rm=TRUE)
SNUC_P_hyper$Tumor_Type <- "SNUC"
SNUC_P_hyper$Gene_Region <- "TSS200"
SNUC_P_hyper <- as.data.frame(SNUC_P_hyper)
common_test_P_hyper=NULL
common_test_P_hyper$Chromosome <- common_probes_test_hyper_annot_promoters[ ,98]
common_test_P_hyper$ID <- common_probes_test_hyper_annot_promoters[ ,2]
common_test_P_hyper$POS <- common_probes_test_hyper_annot_promoters[ ,99]
common_test_P_hyper$UCSC_Gene <- common_probes_test_hyper_annot_promoters[ ,121]
common_test_P_hyper$Mean_Methylation <- rowMeans(common_probes_test_hyper_annot_promoters[ ,3:96], na.rm=TRUE)
common_test_P_hyper$Tumor_Type <- "Common Probes, IDH Mutant"
common_test_P_hyper$Gene_Region <- "TSS200"
common_test_P_hyper <- as.data.frame(common_test_P_hyper)


##Bodies
AML_B_hypo=NULL
AML_B_hypo$Chromosome <- AML_hypo_unique_bodies[ ,99]
AML_B_hypo$ID <- AML_hypo_unique_bodies[ ,2]
AML_B_hypo$POS <- AML_hypo_unique_bodies[ ,100]
AML_B_hypo$UCSC_Gene <- AML_hypo_unique_bodies[ ,122]
AML_B_hypo$Mean_Methylation <- rowMeans(AML_hypo_unique_bodies[ ,28:48], na.rm=TRUE)
AML_B_hypo$Tumor_Type <- "AML"
AML_B_hypo$Gene_Region <- "Body"
AML_B_hypo <- as.data.frame(AML_B_hypo)
Astro_B_hypo=NULL
Astro_B_hypo$Chromosome <- Astro_hypo_unique_bodies[ ,99]
Astro_B_hypo$ID <- Astro_hypo_unique_bodies[ ,2]
Astro_B_hypo$POS <- Astro_hypo_unique_bodies[ ,100]
Astro_B_hypo$UCSC_Gene <- Astro_hypo_unique_bodies[ ,122]
Astro_B_hypo$Mean_Methylation <- rowMeans(Astro_hypo_unique_bodies[ ,c(4:9, 49:73)], na.rm=TRUE)
Astro_B_hypo$Tumor_Type <- "Astrocytoma"
Astro_B_hypo$Gene_Region <- "Body"
Astro_B_hypo <- as.data.frame(Astro_B_hypo)
Chol_B_hypo=NULL
Chol_B_hypo$Chromosome <- Chol_hypo_unique_bodies[ ,99]
Chol_B_hypo$ID <- Chol_hypo_unique_bodies[ ,2]
Chol_B_hypo$POS <- Chol_hypo_unique_bodies[ ,100]
Chol_B_hypo$UCSC_Gene <- Chol_hypo_unique_bodies[ ,122]
Chol_B_hypo$Mean_Methylation <- rowMeans(Chol_hypo_unique_bodies[ ,74:82], na.rm=TRUE)
Chol_B_hypo$Tumor_Type <- "Cholangiocarcinoma"
Chol_B_hypo$Gene_Region <- "Body"
Chol_B_hypo <- as.data.frame(Chol_B_hypo)
Breast_B_hypo=NULL
Breast_B_hypo$Chromosome <- Breast_hypo_unique_bodies[ ,99]
Breast_B_hypo$ID <- Breast_hypo_unique_bodies[ ,2]
Breast_B_hypo$POS <- Breast_hypo_unique_bodies[ ,100]
Breast_B_hypo$UCSC_Gene <- Breast_hypo_unique_bodies[ ,122]
Breast_B_hypo$Mean_Methylation <- rowMeans(Breast_hypo_unique_bodies[ ,10:14], na.rm=TRUE)
Breast_B_hypo$Tumor_Type <- "Breast Cancer"
Breast_B_hypo$Gene_Region <- "Body"
Breast_B_hypo <- as.data.frame(Breast_B_hypo)
Oligo_B_hypo=NULL
Oligo_B_hypo$Chromosome <- Oligo_hypo_unique_bodies[ ,99]
Oligo_B_hypo$ID <- Oligo_hypo_unique_bodies[ ,2]
Oligo_B_hypo$POS <- Oligo_hypo_unique_bodies[ ,100]
Oligo_B_hypo$UCSC_Gene <- Oligo_hypo_unique_bodies[ ,122]
Oligo_B_hypo$Mean_Methylation <- rowMeans(Oligo_hypo_unique_bodies[ ,c(15:19, 83:97)], na.rm=TRUE)
Oligo_B_hypo$Tumor_Type <- "Oligodendroglioma"
Oligo_B_hypo$Gene_Region <- "Body"
Oligo_B_hypo <- as.data.frame(Oligo_B_hypo)
SNUC_B_hypo=NULL
SNUC_B_hypo$Chromosome <- SNUC_hypo_unique_bodies[ ,99]
SNUC_B_hypo$ID <- SNUC_hypo_unique_bodies[ ,2]
SNUC_B_hypo$POS <- SNUC_hypo_unique_bodies[ ,100]
SNUC_B_hypo$UCSC_Gene <- SNUC_hypo_unique_bodies[ ,122]
SNUC_B_hypo$Mean_Methylation <- rowMeans(SNUC_hypo_unique_bodies[ ,20:27], na.rm=TRUE)
SNUC_B_hypo$Tumor_Type <- "SNUC"
SNUC_B_hypo$Gene_Region <- "Body"
SNUC_B_hypo <- as.data.frame(SNUC_B_hypo)
common_test_B_hypo=NULL
common_test_B_hypo$Chromosome <- common_probes_test_hypo_annot_bodies[ ,98]
common_test_B_hypo$ID <- common_probes_test_hypo_annot_bodies[ ,2]
common_test_B_hypo$POS <- common_probes_test_hypo_annot_bodies[ ,99]
common_test_B_hypo$UCSC_Gene <- common_probes_test_hypo_annot_bodies[ ,121]
common_test_B_hypo$Mean_Methylation <- rowMeans(common_probes_test_hypo_annot_bodies[ ,3:96], na.rm=TRUE)
common_test_B_hypo$Tumor_Type <- "Common Probes, IDH Mutant"
common_test_B_hypo$Gene_Region <- "Body"
common_test_B_hypo <- as.data.frame(common_test_B_hypo)




AML_B_hyper=NULL
AML_B_hyper$Chromosome <- AML_hyper_unique_bodies[ ,99]
AML_B_hyper$ID <- AML_hyper_unique_bodies[ ,2]
AML_B_hyper$POS <- AML_hyper_unique_bodies[ ,100]
AML_B_hyper$UCSC_Gene <- AML_hyper_unique_bodies[ ,122]
AML_B_hyper$Mean_Methylation <- rowMeans(AML_hyper_unique_bodies[ ,28:48], na.rm=TRUE)
AML_B_hyper$Tumor_Type <- "AML"
AML_B_hyper$Gene_Region <- "Body"
AML_B_hyper <- as.data.frame(AML_B_hyper)
Astro_B_hyper=NULL
Astro_B_hyper$Chromosome <- Astro_hyper_unique_bodies[ ,99]
Astro_B_hyper$ID <- Astro_hyper_unique_bodies[ ,2]
Astro_B_hyper$POS <- Astro_hyper_unique_bodies[ ,100]
Astro_B_hyper$UCSC_Gene <- Astro_hyper_unique_bodies[ ,122]
Astro_B_hyper$Mean_Methylation <- rowMeans(Astro_hyper_unique_bodies[ ,c(4:9, 49:73)], na.rm=TRUE)
Astro_B_hyper$Tumor_Type <- "Astrocytoma"
Astro_B_hyper$Gene_Region <- "Body"
Astro_B_hyper <- as.data.frame(Astro_B_hyper)
Chol_B_hyper=NULL
Chol_B_hyper$Chromosome <- Chol_hyper_unique_bodies[ ,99]
Chol_B_hyper$ID <- Chol_hyper_unique_bodies[ ,2]
Chol_B_hyper$POS <- Chol_hyper_unique_bodies[ ,100]
Chol_B_hyper$UCSC_Gene <- Chol_hyper_unique_bodies[ ,122]
Chol_B_hyper$Mean_Methylation <- rowMeans(Chol_hyper_unique_bodies[ ,74:82], na.rm=TRUE)
Chol_B_hyper$Tumor_Type <- "Cholangiocarcinoma"
Chol_B_hyper$Gene_Region <- "Body"
Chol_B_hyper <- as.data.frame(Chol_B_hyper)
Breast_B_hyper=NULL
Breast_B_hyper$Chromosome <- Breast_hyper_unique_bodies[ ,99]
Breast_B_hyper$ID <- Breast_hyper_unique_bodies[ ,2]
Breast_B_hyper$POS <- Breast_hyper_unique_bodies[ ,100]
Breast_B_hyper$UCSC_Gene <- Breast_hyper_unique_bodies[ ,122]
Breast_B_hyper$Mean_Methylation <- rowMeans(Breast_hyper_unique_bodies[ ,10:14], na.rm=TRUE)
Breast_B_hyper$Tumor_Type <- "Breast Cancer"
Breast_B_hyper$Gene_Region <- "Body"
Breast_B_hyper <- as.data.frame(Breast_B_hyper)
Oligo_B_hyper=NULL
Oligo_B_hyper$Chromosome <- Oligo_hyper_unique_bodies[ ,99]
Oligo_B_hyper$ID <- Oligo_hyper_unique_bodies[ ,2]
Oligo_B_hyper$POS <- Oligo_hyper_unique_bodies[ ,100]
Oligo_B_hyper$UCSC_Gene <- Oligo_hyper_unique_bodies[ ,122]
Oligo_B_hyper$Mean_Methylation <- rowMeans(Oligo_hyper_unique_bodies[ ,c(15:19, 83:97)], na.rm=TRUE)
Oligo_B_hyper$Tumor_Type <- "Oligodendroglioma"
Oligo_B_hyper$Gene_Region <- "Body"
Oligo_B_hyper <- as.data.frame(Oligo_B_hyper)
SNUC_B_hyper=NULL
SNUC_B_hyper$Chromosome <- SNUC_hyper_unique_bodies[ ,99]
SNUC_B_hyper$ID <- SNUC_hyper_unique_bodies[ ,2]
SNUC_B_hyper$POS <- SNUC_hyper_unique_bodies[ ,100]
SNUC_B_hyper$UCSC_Gene <- SNUC_hyper_unique_bodies[ ,122]
SNUC_B_hyper$Mean_Methylation <- rowMeans(SNUC_hyper_unique_bodies[ ,20:27], na.rm=TRUE)
SNUC_B_hyper$Tumor_Type <- "SNUC"
SNUC_B_hyper$Gene_Region <- "Body"
SNUC_B_hyper <- as.data.frame(SNUC_B_hyper)
common_test_B_hyper=NULL
common_test_B_hyper$Chromosome <- common_probes_test_hyper_annot_bodies[ ,98]
common_test_B_hyper$ID <- common_probes_test_hyper_annot_bodies[ ,2]
common_test_B_hyper$POS <- common_probes_test_hyper_annot_bodies[ ,99]
common_test_B_hyper$UCSC_Gene <- common_probes_test_hyper_annot_bodies[ ,121]
common_test_B_hyper$Mean_Methylation <- rowMeans(common_probes_test_hyper_annot_bodies[ ,3:96], na.rm=TRUE)
common_test_B_hyper$Tumor_Type <- "Common Probes, IDH Mutant"
common_test_B_hyper$Gene_Region <- "Body"
common_test_B_hyper <- as.data.frame(common_test_B_hyper)


################################MUTANT UNIQUE PROBES#######################################

setwd("C:\\")
png(file=paste("KP_MutantEnhancers_Unique.png", sep=""),width=20000,height=3000)
pp <- getDefaultPlotParams(plot.type=3)
pp$data1max <- 4.5
pp$data1min <- 1.9
pp$data2min <- -1.9
pp$data2max <- -4.5
pp$ideogramlateralmargin <- 0.005
pp$leftmargin <- 0.05
pp$rightmargin <- 0.01
kp <- plotKaryotype(plot.type=3, main="IDH Mutant Samples - Unique Gene Enhancer Regions", genome = "hg19",chromosomes="autosomal", cex=7.5, plot.params = pp)

kpDataBackground(kp, data.panel=1, color = "#FFB7B7", r0=1.9, r1=2.29)
kpAxis(kp, labels=c(2, 3.2, 4.5), cex = 2, r0=1.9, r1=2.29)
kpPoints(kp, chr=AML_E_hyper$Chromosome, x=AML_E_hyper$POS, y=AML_E_hyper$Mean_Methylation, cex=1.5, data.panel=1, col="black", r0=1.9, r1=2.29)
#kpArea(kp, chr=AML_E_hyper$Chromosome, x=AML_E_hyper$POS, y=AML_E_hyper$Mean_Methylation, cex=1.5, data.panel=1, col="black", r0=1.9, r1=2.29)
kpAddLabels(kp, labels="AML", cex = 5.5, r0=1.9, r1=2.29, data.panel=1)

kpDataBackground(kp, data.panel=1, color = "#FFC9F3", r0=2.33, r1=2.72)
kpAxis(kp, labels=c(2, 3.2, 4.5), cex = 2, r0=2.33, r1=2.72)
kpPoints(kp, chr=Astro_E_hyper$Chromosome, x=Astro_E_hyper$POS, y=Astro_E_hyper$Mean_Methylation, cex=1.5, data.panel=1, col="black", r0=2.33, r1=2.72)
kpAddLabels(kp, labels="Astrocytoma", cex = 5.5, r0=2.33, r1=2.72, data.panel=1)

kpDataBackground(kp, data.panel=1, color = "#D9C6FF", r0=2.76, r1=3.15)
kpAxis(kp, labels=c(2, 3.2, 4.5), cex = 2, r0=2.76, r1=3.15)
kpPoints(kp, chr=Breast_E_hyper$Chromosome, x=Breast_E_hyper$POS, y=Breast_E_hyper$Mean_Methylation, cex=1.5, data.panel=1, col="black", r0=2.76, r1=3.15)
kpAddLabels(kp, labels="Breast Cancer", cex = 5.5, r0=2.76, r1=3.15, data.panel=1)

kpDataBackground(kp, data.panel=1, color = "#B9CEFF", r0=3.19, r1=3.58)
kpAxis(kp, labels=c(2, 3.2, 4.5), r0=3.19, cex = 2, r1=3.58)
kpPoints(kp, chr=Chol_E_hyper$Chromosome, x=Chol_E_hyper$POS, y=Chol_E_hyper$Mean_Methylation, cex=1.5, data.panel=1, col="black", r0=3.19, r1=3.58)
kpAddLabels(kp, labels="Cholangiocarcinoma", cex = 5.5, r0=3.19, r1=3.58, data.panel=1)

kpDataBackground(kp, data.panel=1, color = "#BAFFBF", r0=3.62, r1=4.01)
kpAxis(kp, labels=c(2, 3.2, 4.5), r0=3.62, r1=4.01, cex = 2)
kpPoints(kp, chr=SNUC_E_hyper$Chromosome, x=SNUC_E_hyper$POS, y=SNUC_E_hyper$Mean_Methylation, cex=1.5, data.panel=1, col="black", r0=3.62, r1=4.01)
kpAddLabels(kp, labels="SNUC", cex = 5.5, r0=3.62, r1=4.01, data.panel=1)

kpDataBackground(kp, data.panel=1, color = "#F7FFBA", r0=4.05, r1=4.44)
kpAxis(kp, labels=c(2, 3.2, 4.5), r0=4.05, r1=4.44, cex = 2)
kpPoints(kp, chr=Oligo_E_hyper$Chromosome, x=Oligo_E_hyper$POS, y=Oligo_E_hyper$Mean_Methylation, cex=1.5, data.panel=1, col="black", r0=4.05, r1=4.44)
kpAddLabels(kp, labels="Oligodendroglioma", cex = 5.5, r0=4.05, r1=4.44, data.panel=1)



kpDataBackground(kp, data.panel=2, color = "#FFB7B7", r0=-1.9, r1=-2.29)
kpAxis(kp, labels=c(-2, -3.2, -4.5), cex = 2, r0=-1.9, r1=-2.29)
kpPoints(kp, chr=AML_E_hypo$Chromosome, x=AML_E_hypo$POS, y=AML_E_hypo$Mean_Methylation, cex=1.5, data.panel=2, col="black", r0=-1.9, r1=-2.29)
#kpArea(kp, chr=AML_E_hypo$Chromosome, x=AML_E_hypo$POS, y=AML_E_hypo$Mean_Methylation, cex=1.5, data.panel=2, col="black", r0=-1.9, r1=-2.29)
kpAddLabels(kp, labels="AML", cex = 5.5, r0=-1.9, r1=-2.29, data.panel=2)

kpDataBackground(kp, data.panel=2, color = "#FFC9F3", r0=-2.33, r1=-2.72)
kpAxis(kp, labels=c(-2, -3.2, -4.5), cex = 2, r0=-2.33, r1=-2.72)
kpPoints(kp, chr=Astro_E_hypo$Chromosome, x=Astro_E_hypo$POS, y=Astro_E_hypo$Mean_Methylation, cex=1.5, data.panel=2, col="black", r0=-2.33, r1=-2.72)
kpAddLabels(kp, labels="Astrocytoma", cex = 5.5, r0=-2.33, r1=-2.72, data.panel=2)

kpDataBackground(kp, data.panel=2, color = "#D9C6FF", r0=-2.76, r1=-3.15)
kpAxis(kp, labels=c(-2, -3.2, -4.5), cex = 2, r0=-2.76, r1=-3.15)
kpPoints(kp, chr=Breast_E_hypo$Chromosome, x=Breast_E_hypo$POS, y=Breast_E_hypo$Mean_Methylation, cex=1.5, data.panel=2, col="black", r0=-2.76, r1=-3.15)
kpAddLabels(kp, labels="Breast Cancer", cex = 5.5, r0=-2.76, r1=-3.15, data.panel=2)

kpDataBackground(kp, data.panel=2, color = "#B9CEFF", r0=-3.19, r1=-3.58)
kpAxis(kp, labels=c(-2, -3.2, -4.5), r0=-3.19, cex = 2, r1=-3.58)
kpPoints(kp, chr=Chol_E_hypo$Chromosome, x=Chol_E_hypo$POS, y=Chol_E_hypo$Mean_Methylation, cex=1.5, data.panel=2, col="black", r0=-3.19, r1=-3.58)
kpAddLabels(kp, labels="Cholangiocarcinoma", cex = 5.5, r0=-3.19, r1=-3.58, data.panel=2)

kpDataBackground(kp, data.panel=2, color = "#BAFFBF", r0=-3.62, r1=-4.01)
kpAxis(kp, labels=c(-2, -3.2, -4.5), r0=-3.62, r1=-4.01, cex = 2)
kpPoints(kp, chr=SNUC_E_hypo$Chromosome, x=SNUC_E_hypo$POS, y=SNUC_E_hypo$Mean_Methylation, cex=1.5, data.panel=2, col="black", r0=-3.62, r1=-4.01)
kpAddLabels(kp, labels="SNUC", cex = 5.5, r0=-3.62, r1=-4.01, data.panel=2)

kpDataBackground(kp, data.panel=2, color = "#F7FFBA", r0=-4.05, r1=-4.44)
kpAxis(kp, labels=c(-2, -3.2, -4.5), r0=-4.05, r1=-4.44, cex = 2)
kpPoints(kp, chr=Oligo_E_hypo$Chromosome, x=Oligo_E_hypo$POS, y=Oligo_E_hypo$Mean_Methylation, cex=1.5, data.panel=2, col="black", r0=-4.05, r1=-4.44)
kpAddLabels(kp, labels="Oligodendroglioma", cex = 5.5, r0=-4.05, r1=-4.44, data.panel=2)

kpAddCytobandLabels(kp, cex=2.5)
dev.off()






################################COMMON MUTANT PROBES#######################################

setwd("C:\\")
png(file=paste("KP_BodiesvsPromotersvsEnhancers_Mutant_Common.png", sep=""),width=20000,height=3000)
pp <- getDefaultPlotParams(plot.type=3)
pp$data1max <- 7
pp$data1min <- 1.9
pp$data2min <- -1.9
pp$data2max <- -7
pp$ideogramlateralmargin <- 0.005
pp$leftmargin <- 0.05
pp$rightmargin <- 0.01
kp <- plotKaryotype(plot.type=3, main="IDH Mutant Samples - Common Probes in Body vs. Promoter vs. Enhancer Gene Regions", genome = "hg19",chromosomes="autosomal", cex=7.5, plot.params = pp)

kpDataBackground(kp, data.panel=1, color = "#FFB7B7", r0=1.9, r1=3.2)
kpAxis(kp, labels=c(2, 4.5, 7), cex = 2, r0=1.9, r1=3.2)
kpPoints(kp, chr=common_test_B_hyper$Chromosome, x=common_test_B_hyper$POS, y=common_test_B_hyper$Mean_Methylation, cex=1.5, data.panel=1, col="black", r0=1.9, r1=3.2)
kpAddLabels(kp, labels="Body", cex = 5.5, r0=1.9, r1=3.2, data.panel=1)

kpDataBackground(kp, data.panel=1, color = "#FFC9F3", r0=3.6, r1=4.9)
kpAxis(kp, labels=c(2, 4.5, 7), cex = 2, r0=3.6, r1=4.9)
kpPoints(kp, chr=common_test_P_hyper$Chromosome, x=common_test_P_hyper$POS, y=common_test_P_hyper$Mean_Methylation, cex=1.5, data.panel=1, col="black", r0=3.6, r1=4.9)
kpAddLabels(kp, labels="Promoter", cex = 5.5, r0=3.6, r1=4.9, data.panel=1)

kpDataBackground(kp, data.panel=1, color = "#D9C6FF", r0=5.3, r1=7)
kpAxis(kp, labels=c(2, 4.5, 7), cex = 2, r0=5.3, r1=7)
kpPoints(kp, chr=common_test_E_hyper$Chromosome, x=common_test_E_hyper$POS, y=common_test_E_hyper$Mean_Methylation, cex=1.5, data.panel=1, col="black", r0=5.3, r1=7)
kpAddLabels(kp, labels="Enhancer", cex = 5.5, r0=5.3, r1=7, data.panel=1)



kpDataBackground(kp, data.panel=2, color = "#FFB7B7", r0=-1.9, r1=-3.2)
kpAxis(kp, labels=c(-2, -4.5, -7), cex = 2, r0=-1.9, r1=-3.2)
kpPoints(kp, chr=common_test_B_hypo$Chromosome, x=common_test_B_hypo$POS, y=common_test_B_hypo$Mean_Methylation, cex=1.5, data.panel=2, col="black", r0=-1.9, r1=-3.2)
kpAddLabels(kp, labels="Body", cex = 5.5, r0=-1.9, r1=-3.2, data.panel=2)

kpDataBackground(kp, data.panel=2, color = "#FFC9F3", r0=-3.6, r1=-4.9)
kpAxis(kp, labels=c(-2, -4.5, -7), cex = 2, r0=-3.6, r1=-4.9)
kpPoints(kp, chr=common_test_P_hypo$Chromosome, x=common_test_P_hypo$POS, y=common_test_P_hypo$Mean_Methylation, cex=1.5, data.panel=2, col="black", r0=-3.6, r1=-4.9)
kpAddLabels(kp, labels="Promoter", cex = 5.5, r0=-3.6, r1=-4.9, data.panel=2)

kpDataBackground(kp, data.panel=2, color = "#D9C6FF", r0=-5.3, r1=-7)
kpAxis(kp, labels=c(-2, -4.5, -7), cex = 2, r0=-5.3, r1=-7)
kpPoints(kp, chr=common_test_E_hypo$Chromosome, x=common_test_E_hypo$POS, y=common_test_E_hypo$Mean_Methylation, cex=1.5, data.panel=2, col="black", r0=-5.3, r1=-7)
kpAddLabels(kp, labels="Enhancer", cex = 5.5, r0=-5.3, r1=-7, data.panel=2)

kpAddCytobandLabels(kp, cex=2.5)
dev.off()
