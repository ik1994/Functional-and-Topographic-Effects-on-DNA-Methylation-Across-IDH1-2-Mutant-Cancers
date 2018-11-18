source("https://bioconductor.org/biocLite.R")
#biocLite("karyoploteR")
#biocLite("regioneR")
#biocLite("stringi")
library(stringi)
library(karyoploteR)
library(regioneR)


#####LOAD IN THE UNIQUE AND COMMON PROBE ANNOTATED FILES####
common_probes_control_hypo_annot <- read.csv("C:/Users/Ramona Bledea/Documents/PROJECT RESULTS/IDH Tumor Project/Results/Control Samples - Differential Methylation/Using M Values/Hypomethylated Probes/control_common_probes_hypo_annot.csv")
common_probes_control_hyper_annot <- read.csv("C:/Users/Ramona Bledea/Documents/PROJECT RESULTS/IDH Tumor Project/Results/Control Samples - Differential Methylation/Using M Values/Hypermethylated Probes/control_common_probes_hyper_annot.csv")
common_probes_test_hypo_annot <- read.csv("C:/Users/Ramona Bledea/Documents/PROJECT RESULTS/IDH Tumor Project/Results/Test Samples - Differential Methylation/Using M Values/Hypomethylated Probes/test_common_probes_hypo_annot.csv")
common_probes_test_hyper_annot <- read.csv("C:/Users/Ramona Bledea/Documents/PROJECT RESULTS/IDH Tumor Project/Results/Test Samples - Differential Methylation/Using M Values/Hypermethylated Probes/test_common_probes_hyper_annot.csv")
AML_hypo_unique <- read.csv("C:/Users/Ramona Bledea/Documents/PROJECT RESULTS/IDH Tumor Project/Results/Test Samples - Differential Methylation/Using M Values/Hypomethylated Probes/AML_unique_hypo_probes.csv")
Astro_hypo_unique <- read.csv("C:/Users/Ramona Bledea/Documents/PROJECT RESULTS/IDH Tumor Project/Results/Test Samples - Differential Methylation/Using M Values/Hypomethylated Probes/Astro_unique_hypo_probes.csv")
Breast_hypo_unique <- read.csv("C:/Users/Ramona Bledea/Documents/PROJECT RESULTS/IDH Tumor Project/Results/Test Samples - Differential Methylation/Using M Values/Hypomethylated Probes/Breast_unique_hypo_probes.csv")
Chol_hypo_unique <- read.csv("C:/Users/Ramona Bledea/Documents/PROJECT RESULTS/IDH Tumor Project/Results/Test Samples - Differential Methylation/Using M Values/Hypomethylated Probes/Chol_unique_hypo_probes.csv")
SNUC_hypo_unique <- read.csv("C:/Users/Ramona Bledea/Documents/PROJECT RESULTS/IDH Tumor Project/Results/Test Samples - Differential Methylation/Using M Values/Hypomethylated Probes/SNUC_unique_hypo_probes.csv")
Oligo_hypo_unique <- read.csv("C:/Users/Ramona Bledea/Documents/PROJECT RESULTS/IDH Tumor Project/Results/Test Samples - Differential Methylation/Using M Values/Hypomethylated Probes/Oligo_unique_hypo_probes.csv")
AML_hyper_unique <- read.csv("C:/Users/Ramona Bledea/Documents/PROJECT RESULTS/IDH Tumor Project/Results/Test Samples - Differential Methylation/Using M Values/hypermethylated Probes/AML_unique_hyper_probes.csv")
Astro_hyper_unique <- read.csv("C:/Users/Ramona Bledea/Documents/PROJECT RESULTS/IDH Tumor Project/Results/Test Samples - Differential Methylation/Using M Values/hypermethylated Probes/Astro_unique_hyper_probes.csv")
Breast_hyper_unique <- read.csv("C:/Users/Ramona Bledea/Documents/PROJECT RESULTS/IDH Tumor Project/Results/Test Samples - Differential Methylation/Using M Values/hypermethylated Probes/Breast_unique_hyper_probes.csv")
Chol_hyper_unique <- read.csv("C:/Users/Ramona Bledea/Documents/PROJECT RESULTS/IDH Tumor Project/Results/Test Samples - Differential Methylation/Using M Values/hypermethylated Probes/Chol_unique_hyper_probes.csv")
SNUC_hyper_unique <- read.csv("C:/Users/Ramona Bledea/Documents/PROJECT RESULTS/IDH Tumor Project/Results/Test Samples - Differential Methylation/Using M Values/hypermethylated Probes/SNUC_unique_hyper_probes.csv")
Oligo_hyper_unique <- read.csv("C:/Users/Ramona Bledea/Documents/PROJECT RESULTS/IDH Tumor Project/Results/Test Samples - Differential Methylation/Using M Values/hypermethylated Probes/Oligo_unique_hyper_probes.csv")
GBM_hypo_unique <- read.csv("C:/Users/Ramona Bledea/Documents/PROJECT RESULTS/IDH Tumor Project/Results/Control Samples - Differential Methylation/Using M Values/hypomethylated Probes/GBM_unique_hypo_probes.csv")
Blood_hypo_unique <- read.csv("C:/Users/Ramona Bledea/Documents/PROJECT RESULTS/IDH Tumor Project/Results/Control Samples - Differential Methylation/Using M Values/hypomethylated Probes/Blood_unique_hypo_probes.csv")
Neuro_hypo_unique <- read.csv("C:/Users/Ramona Bledea/Documents/PROJECT RESULTS/IDH Tumor Project/Results/Control Samples - Differential Methylation/Using M Values/hypomethylated Probes/Neuro_unique_hypo_probes.csv")
SUDEP_hypo_unique <- read.csv("C:/Users/Ramona Bledea/Documents/PROJECT RESULTS/IDH Tumor Project/Results/Control Samples - Differential Methylation/Using M Values/hypomethylated Probes/SUDEP_unique_hypo_probes.csv")
PAstro_hypo_unique <- read.csv("C:/Users/Ramona Bledea/Documents/PROJECT RESULTS/IDH Tumor Project/Results/Control Samples - Differential Methylation/Using M Values/hypomethylated Probes/PAstro_unique_hypo_probes.csv")
GBM_hyper_unique <- read.csv("C:/Users/Ramona Bledea/Documents/PROJECT RESULTS/IDH Tumor Project/Results/Control Samples - Differential Methylation/Using M Values/hypermethylated Probes/GBM_unique_hyper_probes.csv")
Blood_hyper_unique <- read.csv("C:/Users/Ramona Bledea/Documents/PROJECT RESULTS/IDH Tumor Project/Results/Control Samples - Differential Methylation/Using M Values/hypermethylated Probes/Blood_unique_hyper_probes.csv")
Neuro_hyper_unique <- read.csv("C:/Users/Ramona Bledea/Documents/PROJECT RESULTS/IDH Tumor Project/Results/Control Samples - Differential Methylation/Using M Values/hypermethylated Probes/Neuro_unique_hyper_probes.csv")
SUDEP_hyper_unique <- read.csv("C:/Users/Ramona Bledea/Documents/PROJECT RESULTS/IDH Tumor Project/Results/Control Samples - Differential Methylation/Using M Values/hypermethylated Probes/SUDEP_unique_hyper_probes.csv")
PAstro_hyper_unique <- read.csv("C:/Users/Ramona Bledea/Documents/PROJECT RESULTS/IDH Tumor Project/Results/Control Samples - Differential Methylation/Using M Values/hypermethylated Probes/PAstro_unique_hyper_probes.csv")





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



Blood_hypo_unique_enhancers <- subset(Blood_hypo_unique, Blood_hypo_unique$Enhancer == "TRUE")
Blood_hypo_unique_bodies <- Blood_hypo_unique[grep("Body", Blood_hypo_unique$UCSC_RefGene_Group), ]
Blood_hypo_unique_promoters <- Blood_hypo_unique[grep("TSS200", Blood_hypo_unique$UCSC_RefGene_Group), ]

GBM_hypo_unique_enhancers <- subset(GBM_hypo_unique, GBM_hypo_unique$Enhancer == "TRUE")
GBM_hypo_unique_bodies <- GBM_hypo_unique[grep("Body", GBM_hypo_unique$UCSC_RefGene_Group), ]
GBM_hypo_unique_promoters <- GBM_hypo_unique[grep("TSS200", GBM_hypo_unique$UCSC_RefGene_Group), ]

Neuro_hypo_unique_enhancers <- subset(Neuro_hypo_unique, Neuro_hypo_unique$Enhancer == "TRUE")
Neuro_hypo_unique_bodies <- Neuro_hypo_unique[grep("Body", Neuro_hypo_unique$UCSC_RefGene_Group), ]
Neuro_hypo_unique_promoters <- Neuro_hypo_unique[grep("TSS200", Neuro_hypo_unique$UCSC_RefGene_Group), ]

PAstro_hypo_unique_enhancers <- subset(PAstro_hypo_unique, PAstro_hypo_unique$Enhancer == "TRUE")
PAstro_hypo_unique_bodies <- PAstro_hypo_unique[grep("Body", PAstro_hypo_unique$UCSC_RefGene_Group), ]
PAstro_hypo_unique_promoters <- PAstro_hypo_unique[grep("TSS200", PAstro_hypo_unique$UCSC_RefGene_Group), ]

SUDEP_hypo_unique_enhancers <- subset(SUDEP_hypo_unique, SUDEP_hypo_unique$Enhancer == "TRUE")
SUDEP_hypo_unique_bodies <- SUDEP_hypo_unique[grep("Body", SUDEP_hypo_unique$UCSC_RefGene_Group), ]
SUDEP_hypo_unique_promoters <- SUDEP_hypo_unique[grep("TSS200", SUDEP_hypo_unique$UCSC_RefGene_Group), ]

common_probes_control_hypo_annot_enhancers <- subset(common_probes_control_hypo_annot, common_probes_control_hypo_annot$Enhancer == "TRUE")
common_probes_control_hypo_annot_bodies <- common_probes_control_hypo_annot[grep("Body", common_probes_control_hypo_annot$UCSC_RefGene_Group), ]
common_probes_control_hypo_annot_promoters <- common_probes_control_hypo_annot[grep("TSS200", common_probes_control_hypo_annot$UCSC_RefGene_Group), ]

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

Blood_hyper_unique_enhancers <- subset(Blood_hyper_unique, Blood_hyper_unique$Enhancer == "TRUE")
Blood_hyper_unique_bodies <- Blood_hyper_unique[grep("Body", Blood_hyper_unique$UCSC_RefGene_Group), ]
Blood_hyper_unique_promoters <- Blood_hyper_unique[grep("TSS200", Blood_hyper_unique$UCSC_RefGene_Group), ]

GBM_hyper_unique_enhancers <- subset(GBM_hyper_unique, GBM_hyper_unique$Enhancer == "TRUE")
GBM_hyper_unique_bodies <- GBM_hyper_unique[grep("Body", GBM_hyper_unique$UCSC_RefGene_Group), ]
GBM_hyper_unique_promoters <- GBM_hyper_unique[grep("TSS200", GBM_hyper_unique$UCSC_RefGene_Group), ]

Neuro_hyper_unique_enhancers <- subset(Neuro_hyper_unique, Neuro_hyper_unique$Enhancer == "TRUE")
Neuro_hyper_unique_bodies <- Neuro_hyper_unique[grep("Body", Neuro_hyper_unique$UCSC_RefGene_Group), ]
Neuro_hyper_unique_promoters <- Neuro_hyper_unique[grep("TSS200", Neuro_hyper_unique$UCSC_RefGene_Group), ]

PAstro_hyper_unique_enhancers <- subset(PAstro_hyper_unique, PAstro_hyper_unique$Enhancer == "TRUE")
PAstro_hyper_unique_bodies <- PAstro_hyper_unique[grep("Body", PAstro_hyper_unique$UCSC_RefGene_Group), ]
PAstro_hyper_unique_promoters <- PAstro_hyper_unique[grep("TSS200", PAstro_hyper_unique$UCSC_RefGene_Group), ]

SUDEP_hyper_unique_enhancers <- subset(SUDEP_hyper_unique, SUDEP_hyper_unique$Enhancer == "TRUE")
SUDEP_hyper_unique_bodies <- SUDEP_hyper_unique[grep("Body", SUDEP_hyper_unique$UCSC_RefGene_Group), ]
SUDEP_hyper_unique_promoters <- SUDEP_hyper_unique[grep("TSS200", SUDEP_hyper_unique$UCSC_RefGene_Group), ]

common_probes_control_hyper_annot_enhancers <- subset(common_probes_control_hyper_annot, common_probes_control_hyper_annot$Enhancer == "TRUE")
common_probes_control_hyper_annot_bodies <- common_probes_control_hyper_annot[grep("Body", common_probes_control_hyper_annot$UCSC_RefGene_Group), ]
common_probes_control_hyper_annot_promoters <- common_probes_control_hyper_annot[grep("TSS200", common_probes_control_hyper_annot$UCSC_RefGene_Group), ]

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
Blood_E_hypo=NULL
Blood_E_hypo$Chromosome <- Blood_hypo_unique_enhancers[ ,56]
Blood_E_hypo$ID <- Blood_hypo_unique_enhancers[ ,2]
Blood_E_hypo$POS <- Blood_hypo_unique_enhancers[ ,57]
Blood_E_hypo$UCSC_Gene <- Blood_hypo_unique_enhancers[ ,79]
Blood_E_hypo$Mean_Methylation <- rowMeans(Blood_hypo_unique_enhancers[ ,4:35], na.rm=TRUE)
Blood_E_hypo$Tumor_Type <- "Normal Blood"
Blood_E_hypo$Gene_Region <- "Enhancer"
Blood_E_hypo <- as.data.frame(Blood_E_hypo)
PAstro_E_hypo=NULL
PAstro_E_hypo$Chromosome <- PAstro_hypo_unique_enhancers[ ,56]
PAstro_E_hypo$ID <- PAstro_hypo_unique_enhancers[ ,2]
PAstro_E_hypo$POS <- PAstro_hypo_unique_enhancers[ ,57]
PAstro_E_hypo$UCSC_Gene <- PAstro_hypo_unique_enhancers[ ,79]
PAstro_E_hypo$Mean_Methylation <- rowMeans(PAstro_hypo_unique_enhancers[ ,36:41], na.rm=TRUE)
PAstro_E_hypo$Tumor_Type <- "Pilocytic Astrocytoma"
PAstro_E_hypo$Gene_Region <- "Enhancer"
PAstro_E_hypo <- as.data.frame(PAstro_E_hypo)
Neuro_E_hypo=NULL
Neuro_E_hypo$Chromosome <- Neuro_hypo_unique_enhancers[ ,56]
Neuro_E_hypo$ID <- Neuro_hypo_unique_enhancers[ ,2]
Neuro_E_hypo$POS <- Neuro_hypo_unique_enhancers[ ,57]
Neuro_E_hypo$UCSC_Gene <- Neuro_hypo_unique_enhancers[ ,79]
Neuro_E_hypo$Mean_Methylation <- rowMeans(Neuro_hypo_unique_enhancers[ ,42:44], na.rm=TRUE)
Neuro_E_hypo$Tumor_Type <- "Neurocytoma"
Neuro_E_hypo$Gene_Region <- "Enhancer"
Neuro_E_hypo <- as.data.frame(Neuro_E_hypo)
GBM_E_hypo=NULL
GBM_E_hypo$Chromosome <- GBM_hypo_unique_enhancers[ ,56]
GBM_E_hypo$ID <- GBM_hypo_unique_enhancers[ ,2]
GBM_E_hypo$POS <- GBM_hypo_unique_enhancers[ ,57]
GBM_E_hypo$UCSC_Gene <- GBM_hypo_unique_enhancers[ ,79]
GBM_E_hypo$Mean_Methylation <- rowMeans(GBM_hypo_unique_enhancers[ ,45:50], na.rm=TRUE)
GBM_E_hypo$Tumor_Type <- "GBM"
GBM_E_hypo$Gene_Region <- "Enhancer"
GBM_E_hypo <- as.data.frame(GBM_E_hypo)
SUDEP_E_hypo=NULL
SUDEP_E_hypo$Chromosome <- SUDEP_hypo_unique_enhancers[ ,56]
SUDEP_E_hypo$ID <- SUDEP_hypo_unique_enhancers[ ,2]
SUDEP_E_hypo$POS <- SUDEP_hypo_unique_enhancers[ ,57]
SUDEP_E_hypo$UCSC_Gene <- SUDEP_hypo_unique_enhancers[ ,79]
SUDEP_E_hypo$Mean_Methylation <- rowMeans(SUDEP_hypo_unique_enhancers[ ,51:54], na.rm=TRUE)
SUDEP_E_hypo$Tumor_Type <- "SUDEP"
SUDEP_E_hypo$Gene_Region <- "Enhancer"
SUDEP_E_hypo <- as.data.frame(SUDEP_E_hypo)
common_control_E_hypo=NULL
common_control_E_hypo$Chromosome <- common_probes_control_hypo_annot_enhancers[ ,56]
common_control_E_hypo$ID <- common_probes_control_hypo_annot_enhancers[ ,2]
common_control_E_hypo$POS <- common_probes_control_hypo_annot_enhancers[ ,57]
common_control_E_hypo$UCSC_Gene <- common_probes_control_hypo_annot_enhancers[ ,79]
common_control_E_hypo$Mean_Methylation <- rowMeans(common_probes_control_hypo_annot_enhancers[ ,4:54], na.rm=TRUE)
common_control_E_hypo$Tumor_Type <- "Common Probes, IDH WT"
common_control_E_hypo$Gene_Region <- "Enhancer"
common_control_E_hypo <- as.data.frame(common_control_E_hypo)
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
Blood_E_hyper=NULL
Blood_E_hyper$Chromosome <- Blood_hyper_unique_enhancers[ ,56]
Blood_E_hyper$ID <- Blood_hyper_unique_enhancers[ ,2]
Blood_E_hyper$POS <- Blood_hyper_unique_enhancers[ ,57]
Blood_E_hyper$UCSC_Gene <- Blood_hyper_unique_enhancers[ ,79]
Blood_E_hyper$Mean_Methylation <- rowMeans(Blood_hyper_unique_enhancers[ ,4:35], na.rm=TRUE)
Blood_E_hyper$Tumor_Type <- "Normal Blood"
Blood_E_hyper$Gene_Region <- "Enhancer"
Blood_E_hyper <- as.data.frame(Blood_E_hyper)
PAstro_E_hyper=NULL
PAstro_E_hyper$Chromosome <- PAstro_hyper_unique_enhancers[ ,56]
PAstro_E_hyper$ID <- PAstro_hyper_unique_enhancers[ ,2]
PAstro_E_hyper$POS <- PAstro_hyper_unique_enhancers[ ,57]
PAstro_E_hyper$UCSC_Gene <- PAstro_hyper_unique_enhancers[ ,79]
PAstro_E_hyper$Mean_Methylation <- rowMeans(PAstro_hyper_unique_enhancers[ ,36:41], na.rm=TRUE)
PAstro_E_hyper$Tumor_Type <- "Pilocytic Astrocytoma"
PAstro_E_hyper$Gene_Region <- "Enhancer"
PAstro_E_hyper <- as.data.frame(PAstro_E_hyper)
Neuro_E_hyper=NULL
Neuro_E_hyper$Chromosome <- Neuro_hyper_unique_enhancers[ ,56]
Neuro_E_hyper$ID <- Neuro_hyper_unique_enhancers[ ,2]
Neuro_E_hyper$POS <- Neuro_hyper_unique_enhancers[ ,57]
Neuro_E_hyper$UCSC_Gene <- Neuro_hyper_unique_enhancers[ ,79]
Neuro_E_hyper$Mean_Methylation <- rowMeans(Neuro_hyper_unique_enhancers[ ,42:44], na.rm=TRUE)
Neuro_E_hyper$Tumor_Type <- "Neurocytoma"
Neuro_E_hyper$Gene_Region <- "Enhancer"
Neuro_E_hyper <- as.data.frame(Neuro_E_hyper)
GBM_E_hyper=NULL
GBM_E_hyper$Chromosome <- GBM_hyper_unique_enhancers[ ,56]
GBM_E_hyper$ID <- GBM_hyper_unique_enhancers[ ,2]
GBM_E_hyper$POS <- GBM_hyper_unique_enhancers[ ,57]
GBM_E_hyper$UCSC_Gene <- GBM_hyper_unique_enhancers[ ,79]
GBM_E_hyper$Mean_Methylation <- rowMeans(GBM_hyper_unique_enhancers[ ,45:50], na.rm=TRUE)
GBM_E_hyper$Tumor_Type <- "GBM"
GBM_E_hyper$Gene_Region <- "Enhancer"
GBM_E_hyper <- as.data.frame(GBM_E_hyper)
SUDEP_E_hyper=NULL
SUDEP_E_hyper$Chromosome <- SUDEP_hyper_unique_enhancers[ ,56]
SUDEP_E_hyper$ID <- SUDEP_hyper_unique_enhancers[ ,2]
SUDEP_E_hyper$POS <- SUDEP_hyper_unique_enhancers[ ,57]
SUDEP_E_hyper$UCSC_Gene <- SUDEP_hyper_unique_enhancers[ ,79]
SUDEP_E_hyper$Mean_Methylation <- rowMeans(SUDEP_hyper_unique_enhancers[ ,51:54], na.rm=TRUE)
SUDEP_E_hyper$Tumor_Type <- "SUDEP"
SUDEP_E_hyper$Gene_Region <- "Enhancer"
SUDEP_E_hyper <- as.data.frame(SUDEP_E_hyper)
common_control_E_hyper=NULL
common_control_E_hyper$Chromosome <- common_probes_control_hyper_annot_enhancers[ ,56]
common_control_E_hyper$ID <- common_probes_control_hyper_annot_enhancers[ ,2]
common_control_E_hyper$POS <- common_probes_control_hyper_annot_enhancers[ ,57]
common_control_E_hyper$UCSC_Gene <- common_probes_control_hyper_annot_enhancers[ ,79]
common_control_E_hyper$Mean_Methylation <- rowMeans(common_probes_control_hyper_annot_enhancers[ ,4:54], na.rm=TRUE)
common_control_E_hyper$Tumor_Type <- "Common Probes, IDH WT"
common_control_E_hyper$Gene_Region <- "Enhancer"
common_control_E_hyper <- as.data.frame(common_control_E_hyper)
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
Blood_P_hypo=NULL
Blood_P_hypo$Chromosome <- Blood_hypo_unique_promoters[ ,56]
Blood_P_hypo$ID <- Blood_hypo_unique_promoters[ ,2]
Blood_P_hypo$POS <- Blood_hypo_unique_promoters[ ,57]
Blood_P_hypo$UCSC_Gene <- Blood_hypo_unique_promoters[ ,79]
Blood_P_hypo$Mean_Methylation <- rowMeans(Blood_hypo_unique_promoters[ ,4:35], na.rm=TRUE)
Blood_P_hypo$Tumor_Type <- "Normal Blood"
Blood_P_hypo$Gene_Region <- "TSS200"
Blood_P_hypo <- as.data.frame(Blood_P_hypo)
PAstro_P_hypo=NULL
PAstro_P_hypo$Chromosome <- PAstro_hypo_unique_promoters[ ,56]
PAstro_P_hypo$ID <- PAstro_hypo_unique_promoters[ ,2]
PAstro_P_hypo$POS <- PAstro_hypo_unique_promoters[ ,57]
PAstro_P_hypo$UCSC_Gene <- PAstro_hypo_unique_promoters[ ,79]
PAstro_P_hypo$Mean_Methylation <- rowMeans(PAstro_hypo_unique_promoters[ ,36:41], na.rm=TRUE)
PAstro_P_hypo$Tumor_Type <- "Pilocytic Astrocytoma"
PAstro_P_hypo$Gene_Region <- "TSS200"
PAstro_P_hypo <- as.data.frame(PAstro_P_hypo)
Neuro_P_hypo=NULL
Neuro_P_hypo$Chromosome <- Neuro_hypo_unique_promoters[ ,56]
Neuro_P_hypo$ID <- Neuro_hypo_unique_promoters[ ,2]
Neuro_P_hypo$POS <- Neuro_hypo_unique_promoters[ ,57]
Neuro_P_hypo$UCSC_Gene <- Neuro_hypo_unique_promoters[ ,79]
Neuro_P_hypo$Mean_Methylation <- rowMeans(Neuro_hypo_unique_promoters[ ,42:44], na.rm=TRUE)
Neuro_P_hypo$Tumor_Type <- "Neurocytoma"
Neuro_P_hypo$Gene_Region <- "TSS200"
Neuro_P_hypo <- as.data.frame(Neuro_P_hypo)
GBM_P_hypo=NULL
GBM_P_hypo$Chromosome <- GBM_hypo_unique_promoters[ ,56]
GBM_P_hypo$ID <- GBM_hypo_unique_promoters[ ,2]
GBM_P_hypo$POS <- GBM_hypo_unique_promoters[ ,57]
GBM_P_hypo$UCSC_Gene <- GBM_hypo_unique_promoters[ ,79]
GBM_P_hypo$Mean_Methylation <- rowMeans(GBM_hypo_unique_promoters[ ,45:50], na.rm=TRUE)
GBM_P_hypo$Tumor_Type <- "GBM"
GBM_P_hypo$Gene_Region <- "TSS200"
GBM_P_hypo <- as.data.frame(GBM_P_hypo)
SUDEP_P_hypo=NULL
SUDEP_P_hypo$Chromosome <- SUDEP_hypo_unique_promoters[ ,56]
SUDEP_P_hypo$ID <- SUDEP_hypo_unique_promoters[ ,2]
SUDEP_P_hypo$POS <- SUDEP_hypo_unique_promoters[ ,57]
SUDEP_P_hypo$UCSC_Gene <- SUDEP_hypo_unique_promoters[ ,79]
SUDEP_P_hypo$Mean_Methylation <- rowMeans(SUDEP_hypo_unique_promoters[ ,51:54], na.rm=TRUE)
SUDEP_P_hypo$Tumor_Type <- "SUDEP"
SUDEP_P_hypo$Gene_Region <- "TSS200"
SUDEP_P_hypo <- as.data.frame(SUDEP_P_hypo)
common_control_P_hypo=NULL
common_control_P_hypo$Chromosome <- common_probes_control_hypo_annot_promoters[ ,56]
common_control_P_hypo$ID <- common_probes_control_hypo_annot_promoters[ ,2]
common_control_P_hypo$POS <- common_probes_control_hypo_annot_promoters[ ,57]
common_control_P_hypo$UCSC_Gene <- common_probes_control_hypo_annot_promoters[ ,79]
common_control_P_hypo$Mean_Methylation <- rowMeans(common_probes_control_hypo_annot_promoters[ ,4:54], na.rm=TRUE)
common_control_P_hypo$Tumor_Type <- "Common Probes, IDH WT"
common_control_P_hypo$Gene_Region <- "TSS200"
common_control_P_hypo <- as.data.frame(common_control_P_hypo)
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
Blood_P_hyper=NULL
Blood_P_hyper$Chromosome <- Blood_hyper_unique_promoters[ ,56]
Blood_P_hyper$ID <- Blood_hyper_unique_promoters[ ,2]
Blood_P_hyper$POS <- Blood_hyper_unique_promoters[ ,57]
Blood_P_hyper$UCSC_Gene <- Blood_hyper_unique_promoters[ ,79]
Blood_P_hyper$Mean_Methylation <- rowMeans(Blood_hyper_unique_promoters[ ,4:35], na.rm=TRUE)
Blood_P_hyper$Tumor_Type <- "Normal Blood"
Blood_P_hyper$Gene_Region <- "TSS200"
Blood_P_hyper <- as.data.frame(Blood_P_hyper)
PAstro_P_hyper=NULL
PAstro_P_hyper$Chromosome <- PAstro_hyper_unique_promoters[ ,56]
PAstro_P_hyper$ID <- PAstro_hyper_unique_promoters[ ,2]
PAstro_P_hyper$POS <- PAstro_hyper_unique_promoters[ ,57]
PAstro_P_hyper$UCSC_Gene <- PAstro_hyper_unique_promoters[ ,79]
PAstro_P_hyper$Mean_Methylation <- rowMeans(PAstro_hyper_unique_promoters[ ,36:41], na.rm=TRUE)
PAstro_P_hyper$Tumor_Type <- "Pilocytic Astrocytoma"
PAstro_P_hyper$Gene_Region <- "TSS200"
PAstro_P_hyper <- as.data.frame(PAstro_P_hyper)
Neuro_P_hyper=NULL
Neuro_P_hyper$Chromosome <- Neuro_hyper_unique_promoters[ ,56]
Neuro_P_hyper$ID <- Neuro_hyper_unique_promoters[ ,2]
Neuro_P_hyper$POS <- Neuro_hyper_unique_promoters[ ,57]
Neuro_P_hyper$UCSC_Gene <- Neuro_hyper_unique_promoters[ ,79]
Neuro_P_hyper$Mean_Methylation <- rowMeans(Neuro_hyper_unique_promoters[ ,42:44], na.rm=TRUE)
Neuro_P_hyper$Tumor_Type <- "Neurocytoma"
Neuro_P_hyper$Gene_Region <- "TSS200"
Neuro_P_hyper <- as.data.frame(Neuro_P_hyper)
GBM_P_hyper=NULL
GBM_P_hyper$Chromosome <- GBM_hyper_unique_promoters[ ,56]
GBM_P_hyper$ID <- GBM_hyper_unique_promoters[ ,2]
GBM_P_hyper$POS <- GBM_hyper_unique_promoters[ ,57]
GBM_P_hyper$UCSC_Gene <- GBM_hyper_unique_promoters[ ,79]
GBM_P_hyper$Mean_Methylation <- rowMeans(GBM_hyper_unique_promoters[ ,45:50], na.rm=TRUE)
GBM_P_hyper$Tumor_Type <- "GBM"
GBM_P_hyper$Gene_Region <- "TSS200"
GBM_P_hyper <- as.data.frame(GBM_P_hyper)
SUDEP_P_hyper=NULL
SUDEP_P_hyper$Chromosome <- SUDEP_hyper_unique_promoters[ ,56]
SUDEP_P_hyper$ID <- SUDEP_hyper_unique_promoters[ ,2]
SUDEP_P_hyper$POS <- SUDEP_hyper_unique_promoters[ ,57]
SUDEP_P_hyper$UCSC_Gene <- SUDEP_hyper_unique_promoters[ ,79]
SUDEP_P_hyper$Mean_Methylation <- rowMeans(SUDEP_hyper_unique_promoters[ ,51:54], na.rm=TRUE)
SUDEP_P_hyper$Tumor_Type <- "SUDEP"
SUDEP_P_hyper$Gene_Region <- "TSS200"
SUDEP_P_hyper <- as.data.frame(SUDEP_P_hyper)
common_control_P_hyper=NULL
common_control_P_hyper$Chromosome <- common_probes_control_hyper_annot_promoters[ ,56]
common_control_P_hyper$ID <- common_probes_control_hyper_annot_promoters[ ,2]
common_control_P_hyper$POS <- common_probes_control_hyper_annot_promoters[ ,57]
common_control_P_hyper$UCSC_Gene <- common_probes_control_hyper_annot_promoters[ ,79]
common_control_P_hyper$Mean_Methylation <- rowMeans(common_probes_control_hyper_annot_promoters[ ,4:54], na.rm=TRUE)
common_control_P_hyper$Tumor_Type <- "Common Probes, IDH WT"
common_control_P_hyper$Gene_Region <- "TSS200"
common_control_P_hyper <- as.data.frame(common_control_P_hyper)
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
Blood_B_hypo=NULL
Blood_B_hypo$Chromosome <- Blood_hypo_unique_bodies[ ,56]
Blood_B_hypo$ID <- Blood_hypo_unique_bodies[ ,2]
Blood_B_hypo$POS <- Blood_hypo_unique_bodies[ ,57]
Blood_B_hypo$UCSC_Gene <- Blood_hypo_unique_bodies[ ,79]
Blood_B_hypo$Mean_Methylation <- rowMeans(Blood_hypo_unique_bodies[ ,4:35], na.rm=TRUE)
Blood_B_hypo$Tumor_Type <- "Normal Blood"
Blood_B_hypo$Gene_Region <- "Body"
Blood_B_hypo <- as.data.frame(Blood_B_hypo)
PAstro_B_hypo=NULL
PAstro_B_hypo$Chromosome <- PAstro_hypo_unique_bodies[ ,56]
PAstro_B_hypo$ID <- PAstro_hypo_unique_bodies[ ,2]
PAstro_B_hypo$POS <- PAstro_hypo_unique_bodies[ ,57]
PAstro_B_hypo$UCSC_Gene <- PAstro_hypo_unique_bodies[ ,79]
PAstro_B_hypo$Mean_Methylation <- rowMeans(PAstro_hypo_unique_bodies[ ,36:41], na.rm=TRUE)
PAstro_B_hypo$Tumor_Type <- "Pilocytic Astrocytoma"
PAstro_B_hypo$Gene_Region <- "Body"
PAstro_B_hypo <- as.data.frame(PAstro_B_hypo)
Neuro_B_hypo=NULL
Neuro_B_hypo$Chromosome <- Neuro_hypo_unique_bodies[ ,56]
Neuro_B_hypo$ID <- Neuro_hypo_unique_bodies[ ,2]
Neuro_B_hypo$POS <- Neuro_hypo_unique_bodies[ ,57]
Neuro_B_hypo$UCSC_Gene <- Neuro_hypo_unique_bodies[ ,79]
Neuro_B_hypo$Mean_Methylation <- rowMeans(Neuro_hypo_unique_bodies[ ,42:44], na.rm=TRUE)
Neuro_B_hypo$Tumor_Type <- "Neurocytoma"
Neuro_B_hypo$Gene_Region <- "Body"
Neuro_B_hypo <- as.data.frame(Neuro_B_hypo)
GBM_B_hypo=NULL
GBM_B_hypo$Chromosome <- GBM_hypo_unique_bodies[ ,56]
GBM_B_hypo$ID <- GBM_hypo_unique_bodies[ ,2]
GBM_B_hypo$POS <- GBM_hypo_unique_bodies[ ,57]
GBM_B_hypo$UCSC_Gene <- GBM_hypo_unique_bodies[ ,79]
GBM_B_hypo$Mean_Methylation <- rowMeans(GBM_hypo_unique_bodies[ ,45:50], na.rm=TRUE)
GBM_B_hypo$Tumor_Type <- "GBM"
GBM_B_hypo$Gene_Region <- "Body"
GBM_B_hypo <- as.data.frame(GBM_B_hypo)
SUDEP_B_hypo=NULL
SUDEP_B_hypo$Chromosome <- SUDEP_hypo_unique_bodies[ ,56]
SUDEP_B_hypo$ID <- SUDEP_hypo_unique_bodies[ ,2]
SUDEP_B_hypo$POS <- SUDEP_hypo_unique_bodies[ ,57]
SUDEP_B_hypo$UCSC_Gene <- SUDEP_hypo_unique_bodies[ ,79]
SUDEP_B_hypo$Mean_Methylation <- rowMeans(SUDEP_hypo_unique_bodies[ ,51:54], na.rm=TRUE)
SUDEP_B_hypo$Tumor_Type <- "SUDEP"
SUDEP_B_hypo$Gene_Region <- "Body"
SUDEP_B_hypo <- as.data.frame(SUDEP_B_hypo)
common_control_B_hypo=NULL
common_control_B_hypo$Chromosome <- common_probes_control_hypo_annot_bodies[ ,56]
common_control_B_hypo$ID <- common_probes_control_hypo_annot_bodies[ ,2]
common_control_B_hypo$POS <- common_probes_control_hypo_annot_bodies[ ,57]
common_control_B_hypo$UCSC_Gene <- common_probes_control_hypo_annot_bodies[ ,79]
common_control_B_hypo$Mean_Methylation <- rowMeans(common_probes_control_hypo_annot_bodies[ ,4:54], na.rm=TRUE)
common_control_B_hypo$Tumor_Type <- "Common Probes, IDH WT"
common_control_B_hypo$Gene_Region <- "Body"
common_control_B_hypo <- as.data.frame(common_control_B_hypo)
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
Blood_B_hyper=NULL
Blood_B_hyper$Chromosome <- Blood_hyper_unique_bodies[ ,56]
Blood_B_hyper$ID <- Blood_hyper_unique_bodies[ ,2]
Blood_B_hyper$POS <- Blood_hyper_unique_bodies[ ,57]
Blood_B_hyper$UCSC_Gene <- Blood_hyper_unique_bodies[ ,79]
Blood_B_hyper$Mean_Methylation <- rowMeans(Blood_hyper_unique_bodies[ ,4:35], na.rm=TRUE)
Blood_B_hyper$Tumor_Type <- "Normal Blood"
Blood_B_hyper$Gene_Region <- "Body"
Blood_B_hyper <- as.data.frame(Blood_B_hyper)
PAstro_B_hyper=NULL
PAstro_B_hyper$Chromosome <- PAstro_hyper_unique_bodies[ ,56]
PAstro_B_hyper$ID <- PAstro_hyper_unique_bodies[ ,2]
PAstro_B_hyper$POS <- PAstro_hyper_unique_bodies[ ,57]
PAstro_B_hyper$UCSC_Gene <- PAstro_hyper_unique_bodies[ ,79]
PAstro_B_hyper$Mean_Methylation <- rowMeans(PAstro_hyper_unique_bodies[ ,36:41], na.rm=TRUE)
PAstro_B_hyper$Tumor_Type <- "Pilocytic Astrocytoma"
PAstro_B_hyper$Gene_Region <- "Body"
PAstro_B_hyper <- as.data.frame(PAstro_B_hyper)
Neuro_B_hyper=NULL
Neuro_B_hyper$Chromosome <- Neuro_hyper_unique_bodies[ ,56]
Neuro_B_hyper$ID <- Neuro_hyper_unique_bodies[ ,2]
Neuro_B_hyper$POS <- Neuro_hyper_unique_bodies[ ,57]
Neuro_B_hyper$UCSC_Gene <- Neuro_hyper_unique_bodies[ ,79]
Neuro_B_hyper$Mean_Methylation <- rowMeans(Neuro_hyper_unique_bodies[ ,42:44], na.rm=TRUE)
Neuro_B_hyper$Tumor_Type <- "Neurocytoma"
Neuro_B_hyper$Gene_Region <- "Body"
Neuro_B_hyper <- as.data.frame(Neuro_B_hyper)
GBM_B_hyper=NULL
GBM_B_hyper$Chromosome <- GBM_hyper_unique_bodies[ ,56]
GBM_B_hyper$ID <- GBM_hyper_unique_bodies[ ,2]
GBM_B_hyper$POS <- GBM_hyper_unique_bodies[ ,57]
GBM_B_hyper$UCSC_Gene <- GBM_hyper_unique_bodies[ ,79]
GBM_B_hyper$Mean_Methylation <- rowMeans(GBM_hyper_unique_bodies[ ,45:50], na.rm=TRUE)
GBM_B_hyper$Tumor_Type <- "GBM"
GBM_B_hyper$Gene_Region <- "Body"
GBM_B_hyper <- as.data.frame(GBM_B_hyper)
SUDEP_B_hyper=NULL
SUDEP_B_hyper$Chromosome <- SUDEP_hyper_unique_bodies[ ,56]
SUDEP_B_hyper$ID <- SUDEP_hyper_unique_bodies[ ,2]
SUDEP_B_hyper$POS <- SUDEP_hyper_unique_bodies[ ,57]
SUDEP_B_hyper$UCSC_Gene <- SUDEP_hyper_unique_bodies[ ,79]
SUDEP_B_hyper$Mean_Methylation <- rowMeans(SUDEP_hyper_unique_bodies[ ,51:54], na.rm=TRUE)
SUDEP_B_hyper$Tumor_Type <- "SUDEP"
SUDEP_B_hyper$Gene_Region <- "Body"
SUDEP_B_hyper <- as.data.frame(SUDEP_B_hyper)
common_control_B_hyper=NULL
common_control_B_hyper$Chromosome <- common_probes_control_hyper_annot_bodies[ ,56]
common_control_B_hyper$ID <- common_probes_control_hyper_annot_bodies[ ,2]
common_control_B_hyper$POS <- common_probes_control_hyper_annot_bodies[ ,57]
common_control_B_hyper$UCSC_Gene <- common_probes_control_hyper_annot_bodies[ ,79]
common_control_B_hyper$Mean_Methylation <- rowMeans(common_probes_control_hyper_annot_bodies[ ,4:54], na.rm=TRUE)
common_control_B_hyper$Tumor_Type <- "Common Probes, IDH WT"
common_control_B_hyper$Gene_Region <- "Body"
common_control_B_hyper <- as.data.frame(common_control_B_hyper)
common_test_B_hyper=NULL
common_test_B_hyper$Chromosome <- common_probes_test_hyper_annot_bodies[ ,98]
common_test_B_hyper$ID <- common_probes_test_hyper_annot_bodies[ ,2]
common_test_B_hyper$POS <- common_probes_test_hyper_annot_bodies[ ,99]
common_test_B_hyper$UCSC_Gene <- common_probes_test_hyper_annot_bodies[ ,121]
common_test_B_hyper$Mean_Methylation <- rowMeans(common_probes_test_hyper_annot_bodies[ ,3:96], na.rm=TRUE)
common_test_B_hyper$Tumor_Type <- "Common Probes, IDH Mutant"
common_test_B_hyper$Gene_Region <- "Body"
common_test_B_hyper <- as.data.frame(common_test_B_hyper)

##Whole Genome - dont include??
AML_WG_hypo=NULL
AML_WG_hypo$Chromosome <- AML_hypo_unique[ ,99]
AML_WG_hypo$Mean_Methylation <- rowMeans(AML_hypo_unique[ ,28:48], na.rm=TRUE)
AML_WG_hypo$Tumor_Type <- "AML"
AML_WG_hypo$Gene_Region <- "All Hypomethylated"
AML_WG_hypo <- as.data.frame(AML_WG_hypo)
Astro_WG_hypo=NULL
Astro_WG_hypo$Chromosome <- Astro_hypo_unique[ ,99]
Astro_WG_hypo$Mean_Methylation <- rowMeans(Astro_hypo_unique[ ,c(4:9, 49:73)], na.rm=TRUE)
Astro_WG_hypo$Tumor_Type <- "Astrocytoma"
Astro_WG_hypo$Gene_Region <- "All Hypomethylated"
Astro_WG_hypo <- as.data.frame(Astro_WG_hypo)
Chol_WG_hypo=NULL
Chol_WG_hypo$Chromosome <- Chol_hypo_unique[ ,99]
Chol_WG_hypo$Mean_Methylation <- rowMeans(Chol_hypo_unique[ ,74:82], na.rm=TRUE)
Chol_WG_hypo$Tumor_Type <- "Cholangiocarcinoma"
Chol_WG_hypo$Gene_Region <- "All Hypomethylated"
Chol_WG_hypo <- as.data.frame(Chol_WG_hypo)
Breast_WG_hypo=NULL
Breast_WG_hypo$Chromosome <- Breast_hypo_unique[ ,99]
Breast_WG_hypo$Mean_Methylation <- rowMeans(Breast_hypo_unique[ ,10:14], na.rm=TRUE)
Breast_WG_hypo$Tumor_Type <- "Breast Cancer"
Breast_WG_hypo$Gene_Region <- "All Hypomethylated"
Breast_WG_hypo <- as.data.frame(Breast_WG_hypo)
Oligo_WG_hypo=NULL
Oligo_WG_hypo$Chromosome <- Oligo_hypo_unique[ ,99]
Oligo_WG_hypo$Mean_Methylation <- rowMeans(Oligo_hypo_unique[ ,c(15:19, 83:97)], na.rm=TRUE)
Oligo_WG_hypo$Tumor_Type <- "Oligodendroglioma"
Oligo_WG_hypo$Gene_Region <- "All Hypomethylated"
Oligo_WG_hypo <- as.data.frame(Oligo_WG_hypo)
SNUC_WG_hypo=NULL
SNUC_WG_hypo$Chromosome <- SNUC_hypo_unique[ ,99]
SNUC_WG_hypo$Mean_Methylation <- rowMeans(SNUC_hypo_unique[ ,20:27], na.rm=TRUE)
SNUC_WG_hypo$Tumor_Type <- "SNUC"
SNUC_WG_hypo$Gene_Region <- "All Hypomethylated"
SNUC_WG_hypo <- as.data.frame(SNUC_WG_hypo)
Blood_WG_hypo=NULL
Blood_WG_hypo$Chromosome <- Blood_hypo_unique[ ,56]
Blood_WG_hypo$Mean_Methylation <- rowMeans(Blood_hypo_unique[ ,4:35], na.rm=TRUE)
Blood_WG_hypo$Tumor_Type <- "Normal Blood"
Blood_WG_hypo$Gene_Region <- "All Hypomethylated"
Blood_WG_hypo <- as.data.frame(Blood_WG_hypo)
PAstro_WG_hypo=NULL
PAstro_WG_hypo$Chromosome <- PAstro_hypo_unique[ ,56]
PAstro_WG_hypo$Mean_Methylation <- rowMeans(PAstro_hypo_unique[ ,36:41], na.rm=TRUE)
PAstro_WG_hypo$Tumor_Type <- "Pilocytic Astrocytoma"
PAstro_WG_hypo$Gene_Region <- "All Hypomethylated"
PAstro_WG_hypo <- as.data.frame(PAstro_WG_hypo)
Neuro_WG_hypo=NULL
Neuro_WG_hypo$Chromosome <- Neuro_hypo_unique[ ,56]
Neuro_WG_hypo$Mean_Methylation <- rowMeans(Neuro_hypo_unique[ ,42:44], na.rm=TRUE)
Neuro_WG_hypo$Tumor_Type <- "Neurocytoma"
Neuro_WG_hypo$Gene_Region <- "All Hypomethylated"
Neuro_WG_hypo <- as.data.frame(Neuro_WG_hypo)
GBM_WG_hypo=NULL
GBM_WG_hypo$Chromosome <- GBM_hypo_unique[ ,56]
GBM_WG_hypo$Mean_Methylation <- rowMeans(GBM_hypo_unique[ ,45:50], na.rm=TRUE)
GBM_WG_hypo$Tumor_Type <- "GBM"
GBM_WG_hypo$Gene_Region <- "All Hypomethylated"
GBM_WG_hypo <- as.data.frame(GBM_WG_hypo)
SUDEP_WG_hypo=NULL
SUDEP_WG_hypo$Chromosome <- SUDEP_hypo_unique[ ,56]
SUDEP_WG_hypo$Mean_Methylation <- rowMeans(SUDEP_hypo_unique[ ,51:54], na.rm=TRUE)
SUDEP_WG_hypo$Tumor_Type <- "SUDEP"
SUDEP_WG_hypo$Gene_Region <- "All Hypomethylated"
SUDEP_WG_hypo <- as.data.frame(SUDEP_WG_hypo)
common_control_WG_hypo=NULL
common_control_WG_hypo$Chromosome <- common_probes_control_hypo_annot[ ,56]
common_control_WG_hypo$Mean_Methylation <- rowMeans(common_probes_control_hypo_annot[ ,4:54], na.rm=TRUE)
common_control_WG_hypo$Tumor_Type <- "Common Probes, IDH WT"
common_control_WG_hypo$Gene_Region <- "All Hypomethylated"
common_control_WG_hypo <- as.data.frame(common_control_WG_hypo)
common_test_WG_hypo=NULL
common_test_WG_hypo$Chromosome <- common_probes_test_hypo_annot[ ,98]
common_test_WG_hypo$Mean_Methylation <- rowMeans(common_probes_test_hypo_annot[ ,3:96], na.rm=TRUE)
common_test_WG_hypo$Tumor_Type <- "Common Probes, IDH Mutant"
common_test_WG_hypo$Gene_Region <- "All Hypomethylated"
common_test_WG_hypo <- as.data.frame(common_test_WG_hypo)

AML_WG_hyper=NULL
AML_WG_hyper$Chromosome <- AML_hyper_unique[ ,99]
AML_WG_hyper$Mean_Methylation <- rowMeans(AML_hyper_unique[ ,28:48], na.rm=TRUE)
AML_WG_hyper$Tumor_Type <- "AML"
AML_WG_hyper$Gene_Region <- "All Hypermethylated"
AML_WG_hyper <- as.data.frame(AML_WG_hyper)
Astro_WG_hyper=NULL
Astro_WG_hyper$Chromosome <- Astro_hyper_unique[ ,99]
Astro_WG_hyper$Mean_Methylation <- rowMeans(Astro_hyper_unique[ ,c(4:9), (49:73)], na.rm=TRUE)
Astro_WG_hyper$Tumor_Type <- "Astrocytoma"
Astro_WG_hyper$Gene_Region <- "All Hypermethylated"
Astro_WG_hyper <- as.data.frame(Astro_WG_hyper)
Chol_WG_hyper=NULL
Chol_WG_hyper$Chromosome <- Chol_hyper_unique[ ,99]
Chol_WG_hyper$Mean_Methylation <- rowMeans(Chol_hyper_unique[ ,74:82], na.rm=TRUE)
Chol_WG_hyper$Tumor_Type <- "Cholangiocarcinoma"
Chol_WG_hyper$Gene_Region <- "All Hypermethylated"
Chol_WG_hyper <- as.data.frame(Chol_WG_hyper)
Breast_WG_hyper=NULL
Breast_WG_hyper$Chromosome <- Breast_hyper_unique[ ,99]
Breast_WG_hyper$Mean_Methylation <- rowMeans(Breast_hyper_unique[ ,10:14], na.rm=TRUE)
Breast_WG_hyper$Tumor_Type <- "Breast Cancer"
Breast_WG_hyper$Gene_Region <- "All Hypermethylated"
Breast_WG_hyper <- as.data.frame(Breast_WG_hyper)
Oligo_WG_hyper=NULL
Oligo_WG_hyper$Chromosome <- Oligo_hyper_unique[ ,99]
Oligo_WG_hyper$Mean_Methylation <- rowMeans(Oligo_hyper_unique[ ,c(15:19, 83:97)], na.rm=TRUE)
Oligo_WG_hyper$Tumor_Type <- "Oligodendroglioma"
Oligo_WG_hyper$Gene_Region <- "All Hypermethylated"
Oligo_WG_hyper <- as.data.frame(Oligo_WG_hyper)
SNUC_WG_hyper=NULL
SNUC_WG_hyper$Chromosome <- SNUC_hyper_unique[ ,99]
SNUC_WG_hyper$Mean_Methylation <- rowMeans(SNUC_hyper_unique[ ,20:27], na.rm=TRUE)
SNUC_WG_hyper$Tumor_Type <- "SNUC"
SNUC_WG_hyper$Gene_Region <- "All Hypermethylated"
SNUC_WG_hyper <- as.data.frame(SNUC_WG_hyper)
Blood_WG_hyper=NULL
Blood_WG_hyper$Chromosome <- Blood_hyper_unique[ ,56]
Blood_WG_hyper$Mean_Methylation <- rowMeans(Blood_hyper_unique[ ,4:35], na.rm=TRUE)
Blood_WG_hyper$Tumor_Type <- "Normal Blood"
Blood_WG_hyper$Gene_Region <- "All Hypermethylated"
Blood_WG_hyper <- as.data.frame(Blood_WG_hyper)
PAstro_WG_hyper=NULL
PAstro_WG_hyper$Chromosome <- PAstro_hyper_unique[ ,56]
PAstro_WG_hyper$Mean_Methylation <- rowMeans(PAstro_hyper_unique[ ,36:41], na.rm=TRUE)
PAstro_WG_hyper$Tumor_Type <- "Pilocytic Astrocytoma"
PAstro_WG_hyper$Gene_Region <- "All Hypermethylated"
PAstro_WG_hyper <- as.data.frame(PAstro_WG_hyper)
Neuro_WG_hyper=NULL
Neuro_WG_hyper$Chromosome <- Neuro_hyper_unique[ ,56]
Neuro_WG_hyper$Mean_Methylation <- rowMeans(Neuro_hyper_unique[ ,42:44], na.rm=TRUE)
Neuro_WG_hyper$Tumor_Type <- "Neurocytoma"
Neuro_WG_hyper$Gene_Region <- "All Hypermethylated"
Neuro_WG_hyper <- as.data.frame(Neuro_WG_hyper)
GBM_WG_hyper=NULL
GBM_WG_hyper$Chromosome <- GBM_hyper_unique[ ,56]
GBM_WG_hyper$Mean_Methylation <- rowMeans(GBM_hyper_unique[ ,45:50], na.rm=TRUE)
GBM_WG_hyper$Tumor_Type <- "GBM"
GBM_WG_hyper$Gene_Region <- "All Hypermethylated"
GBM_WG_hyper <- as.data.frame(GBM_WG_hyper)
SUDEP_WG_hyper=NULL
SUDEP_WG_hyper$Chromosome <- SUDEP_hyper_unique[ ,56]
SUDEP_WG_hyper$Mean_Methylation <- rowMeans(SUDEP_hyper_unique[ ,51:54], na.rm=TRUE)
SUDEP_WG_hyper$Tumor_Type <- "SUDEP"
SUDEP_WG_hyper$Gene_Region <- "All Hypermethylated"
SUDEP_WG_hyper <- as.data.frame(SUDEP_WG_hyper)
common_control_WG_hyper=NULL
common_control_WG_hyper$Chromosome <- common_probes_control_hyper_annot[ ,56]
common_control_WG_hyper$Mean_Methylation <- rowMeans(common_probes_control_hyper_annot[ ,4:54], na.rm=TRUE)
common_control_WG_hyper$Tumor_Type <- "Common Probes, IDH WT"
common_control_WG_hyper$Gene_Region <- "All Hypermethylated"
common_control_WG_hyper <- as.data.frame(common_control_WG_hyper)
common_test_WG_hyper=NULL
common_test_WG_hyper$Chromosome <- common_probes_test_hyper_annot[ ,98]
common_test_WG_hyper$Mean_Methylation <- rowMeans(common_probes_test_hyper_annot[ ,3:96], na.rm=TRUE)
common_test_WG_hyper$Tumor_Type <- "Common Probes, IDH Mutant"
common_test_WG_hyper$Gene_Region <- "All Hypermethylated"
common_test_WG_hyper <- as.data.frame(common_test_WG_hyper)




################################MUTANT UNIQUE PROBES#######################################

setwd("C:\\Users\\Ramona Bledea\\Documents\\PROJECT RESULTS\\IDH Tumor Project\\Results\\All Figures\\Mean Methylation vs Position")
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






################################WT UNIQUE PROBES#######################################

setwd("C:\\Users\\Ramona Bledea\\Documents\\PROJECT RESULTS\\IDH Tumor Project\\Results\\All Figures\\Mean Methylation vs Position")
png(file=paste("KP_WTPromoters_Unique.png", sep=""),width=20000,height=3000)
pp <- getDefaultPlotParams(plot.type=3)
pp$data1max <- 4.5
pp$data1min <- 1.9
pp$data2min <- -1.9
pp$data2max <- -4.5
pp$ideogramlateralmargin <- 0.005
pp$leftmargin <- 0.05
pp$rightmargin <- 0.01
kp <- plotKaryotype(plot.type=3, main="IDH WT Samples - Unique Gene Promoter Regions", genome = "hg19",chromosomes="autosomal", cex=7.5, plot.params = pp)

kpDataBackground(kp, data.panel=1, color = "#FFB7B7", r0=1.9, r1=2.29)
kpAxis(kp, labels=c(2, 3.2, 4.5), cex = 2, r0=1.9, r1=2.29)
kpPoints(kp, chr=Blood_P_hyper$Chromosome, x=Blood_P_hyper$POS, y=Blood_P_hyper$Mean_Methylation, cex=1.5, data.panel=1, col="black", r0=1.9, r1=2.29)
kpAddLabels(kp, labels="Normal Blood", cex = 5.5, r0=1.9, r1=2.29, data.panel=1)

kpDataBackground(kp, data.panel=1, color = "#FFC9F3", r0=2.33, r1=2.72)
kpAxis(kp, labels=c(2, 3.2, 4.5), cex = 2, r0=2.33, r1=2.72)
kpPoints(kp, chr=PAstro_P_hyper$Chromosome, x=PAstro_P_hyper$POS, y=PAstro_P_hyper$Mean_Methylation, cex=1.5, data.panel=1, col="black", r0=2.33, r1=2.72)
kpAddLabels(kp, labels="Pilocytic Astrocytoma", cex = 5.5, r0=2.33, r1=2.72, data.panel=1)

kpDataBackground(kp, data.panel=1, color = "#D9C6FF", r0=2.76, r1=3.15)
kpAxis(kp, labels=c(2, 3.2, 4.5), cex = 2, r0=2.76, r1=3.15)
kpPoints(kp, chr=Neuro_P_hyper$Chromosome, x=Neuro_P_hyper$POS, y=Neuro_P_hyper$Mean_Methylation, cex=1.5, data.panel=1, col="black", r0=2.76, r1=3.15)
kpAddLabels(kp, labels="Neurocytoma", cex = 5.5, r0=2.76, r1=3.15, data.panel=1)

kpDataBackground(kp, data.panel=1, color = "#B9CEFF", r0=3.19, r1=3.58)
kpAxis(kp, labels=c(2, 3.2, 4.5), r0=3.19, cex = 2, r1=3.58)
kpPoints(kp, chr=GBM_P_hyper$Chromosome, x=GBM_P_hyper$POS, y=GBM_P_hyper$Mean_Methylation, cex=1.5, data.panel=1, col="black", r0=3.19, r1=3.58)
kpAddLabels(kp, labels="GBM", cex = 5.5, r0=3.19, r1=3.58, data.panel=1)

kpDataBackground(kp, data.panel=1, color = "#BAFFBF", r0=3.62, r1=4.01)
kpAxis(kp, labels=c(2, 3.2, 4.5), r0=3.62, r1=4.01, cex = 2)
kpPoints(kp, chr=SUDEP_P_hyper$Chromosome, x=SUDEP_P_hyper$POS, y=SUDEP_P_hyper$Mean_Methylation, cex=1.5, data.panel=1, col="black", r0=3.62, r1=4.01)
kpAddLabels(kp, labels="SUDEP", cex = 5.5, r0=3.62, r1=4.01, data.panel=1)



kpDataBackground(kp, data.panel=2, color = "#FFB7B7", r0=-1.9, r1=-2.29)
kpAxis(kp, labels=c(-2, -3.2, -4.5), cex = 2, r0=-1.9, r1=-2.29)
kpPoints(kp, chr=Blood_P_hypo$Chromosome, x=Blood_P_hypo$POS, y=Blood_P_hypo$Mean_Methylation, cex=1.5, data.panel=2, col="black", r0=-1.9, r1=-2.29)
kpAddLabels(kp, labels="Normal Blood", cex = 5.5, r0=-1.9, r1=-2.29, data.panel=2)

kpDataBackground(kp, data.panel=2, color = "#FFC9F3", r0=-2.33, r1=-2.72)
kpAxis(kp, labels=c(-2, -3.2, -4.5), cex = 2, r0=-2.33, r1=-2.72)
kpPoints(kp, chr=PAstro_P_hypo$Chromosome, x=PAstro_P_hypo$POS, y=PAstro_P_hypo$Mean_Methylation, cex=1.5, data.panel=2, col="black", r0=-2.33, r1=-2.72)
kpAddLabels(kp, labels="Pilocytic Astrocytoma", cex = 5.5, r0=-2.33, r1=-2.72, data.panel=2)

kpDataBackground(kp, data.panel=2, color = "#D9C6FF", r0=-2.76, r1=-3.15)
kpAxis(kp, labels=c(-2, -3.2, -4.5), cex = 2, r0=-2.76, r1=-3.15)
kpPoints(kp, chr=Neuro_P_hypo$Chromosome, x=Neuro_P_hypo$POS, y=Neuro_P_hypo$Mean_Methylation, cex=1.5, data.panel=2, col="black", r0=-2.76, r1=-3.15)
kpAddLabels(kp, labels="Neurocytoma", cex = 5.5, r0=-2.76, r1=-3.15, data.panel=2)

kpDataBackground(kp, data.panel=2, color = "#B9CEFF", r0=-3.19, r1=-3.58)
kpAxis(kp, labels=c(-2, -3.2, -4.5), r0=-3.19, cex = 2, r1=-3.58)
kpPoints(kp, chr=GBM_P_hypo$Chromosome, x=GBM_P_hypo$POS, y=GBM_P_hypo$Mean_Methylation, cex=1.5, data.panel=2, col="black", r0=-3.19, r1=-3.58)
kpAddLabels(kp, labels="GBM", cex = 5.5, r0=-3.19, r1=-3.58, data.panel=2)

kpDataBackground(kp, data.panel=2, color = "#BAFFBF", r0=-3.62, r1=-4.01)
kpAxis(kp, labels=c(-2, -3.2, -4.5), r0=-3.62, r1=-4.01, cex = 2)
kpPoints(kp, chr=SUDEP_P_hypo$Chromosome, x=SUDEP_P_hypo$POS, y=SUDEP_P_hypo$Mean_Methylation, cex=1.5, data.panel=2, col="black", r0=-3.62, r1=-4.01)
kpAddLabels(kp, labels="SUDEP", cex = 5.5, r0=-3.62, r1=-4.01, data.panel=2)

kpAddCytobandLabels(kp, cex=2.5)
dev.off()






################################COMMON MUTANT vs WT PROBES#######################################

setwd("C:\\Users\\Ramona Bledea\\Documents\\PROJECT RESULTS\\IDH Tumor Project\\Results\\All Figures\\Mean Methylation vs Position")
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