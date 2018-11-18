###########MEAN GLOBAL METHYLATION BOXPLOTS AND HYPO/HYPER METHYLATED PROBE DISTRIBUTION BARPLOTS###########

##Load in the M-value files and merge them with the annotated differential meth files
test_M_values <- read.csv("C:/Users/Ramona Bledea/Documents/PROJECT RESULTS/IDH Tumor Project/Results/Test Samples - Differential Methylation/Using M Values/M-values.csv")
test_probes <- read.csv("C:/Users/Ramona Bledea/Documents/PROJECT RESULTS/IDH Tumor Project/Results/Test Samples - Differential Methylation/Using M Values/differential_methylation_withAnnotation-M.csv")
control_M_values <- read.csv("C:/Users/Ramona Bledea/Documents/PROJECT RESULTS/IDH Tumor Project/Results/Control Samples - Differential Methylation/Using M Values/M-values.csv")
control_probes <- read.csv("C:/Users/Ramona Bledea/Documents/PROJECT RESULTS/IDH Tumor Project/Results/Control Samples - Differential Methylation/Using M Values/differential_methylation_withAnnotation-M.csv")
test_samples <- merge(test_M_values, test_probes, by="ID")
control_samples <- merge(control_M_values, control_probes, by="ID")
##Write to csvs. May need to manually resort each by q value.
setwd("C:/Users/Ramona Bledea/Documents/PROJECT RESULTS/IDH Tumor Project/Results/Test Samples - Differential Methylation/Using M Values")
write.csv(test_samples, file="test_samples_diffmeth_annot_mvalues_masterfile.csv")
setwd("C:/Users/Ramona Bledea/Documents/PROJECT RESULTS/IDH Tumor Project/Results/Control Samples - Differential Methylation/Using M Values")
write.csv(control_samples, file="control_samples_diffmeth_annot_mvalues_masterfile.csv")

##Open each file and filter in excel for Enhancers, Bodies, and TSS200. Save each separately and reload (total of 6 files - 3 for test data, 3 for control)
test_WG <- read.csv("C:/Users/Ramona Bledea/Documents/PROJECT RESULTS/IDH Tumor Project/Results/Test Samples - Differential Methylation/Using M Values/test_samples_diffmeth_annot_mvalues_masterfile.csv")
control_WG <-read.csv("C:/Users/Ramona Bledea/Documents/PROJECT RESULTS/IDH Tumor Project/Results/Control Samples - Differential Methylation/Using M Values/control_samples_diffmeth_annot_mvalues_masterfile.csv") 
test_enhancers <- read.csv("C:/Users/Ramona Bledea/Documents/PROJECT RESULTS/IDH Tumor Project/Results/Test Samples - Differential Methylation/Using M Values/test_samples_diffmeth_annot_mvalues_masterfile-ENHANCERS.csv")
test_promoters <- read.csv("C:/Users/Ramona Bledea/Documents/PROJECT RESULTS/IDH Tumor Project/Results/Test Samples - Differential Methylation/Using M Values/test_samples_diffmeth_annot_mvalues_masterfile-PROMOTERS.csv")
test_bodies <- read.csv("C:/Users/Ramona Bledea/Documents/PROJECT RESULTS/IDH Tumor Project/Results/Test Samples - Differential Methylation/Using M Values/test_samples_diffmeth_annot_mvalues_masterfile-BODIES.csv")
control_enhancers <- read.csv("C:/Users/Ramona Bledea/Documents/PROJECT RESULTS/IDH Tumor Project/Results/Control Samples - Differential Methylation/Using M Values/control_samples_diffmeth_annot_mvalues_masterfile-ENHANCERS.csv")
control_promoters <- read.csv("C:/Users/Ramona Bledea/Documents/PROJECT RESULTS/IDH Tumor Project/Results/Control Samples - Differential Methylation/Using M Values/control_samples_diffmeth_annot_mvalues_masterfile-PROMOTERS.csv")
control_bodies <- read.csv("C:/Users/Ramona Bledea/Documents/PROJECT RESULTS/IDH Tumor Project/Results/Control Samples - Differential Methylation/Using M Values/control_samples_diffmeth_annot_mvalues_masterfile-BODIES.csv")


##Enhancers
names(test_enhancers)
names(control_enhancers)

AML_E=NULL
AML_E$Chromosome <- test_enhancers[ ,98]
AML_E$Mean_Methylation <- rowMeans(test_enhancers[ ,27:47], na.rm=TRUE)
AML_E$Tumor_Type <- "AML"
AML_E$Gene_Region <- "Enhancer"
AML_E <- as.data.frame(AML_E)
AML_E_hypo <- subset(AML_E, Mean_Methylation<=-2)
AML_E_hyper <- subset(AML_E, Mean_Methylation>=2)
Astro_E=NULL
Astro_E$Chromosome <- test_enhancers[ ,98]
Astro_E$Mean_Methylation <- rowMeans(test_enhancers[ ,c(3:8, 48:72)], na.rm=TRUE)
Astro_E$Tumor_Type <- "Astrocytoma"
Astro_E$Gene_Region <- "Enhancer"
Astro_E <- as.data.frame(Astro_E)
Astro_E_hypo <- subset(Astro_E, Mean_Methylation<=-2)
Astro_E_hyper <- subset(Astro_E, Mean_Methylation>=2)
Chol_E=NULL
Chol_E$Chromosome <- test_enhancers[ ,98]
Chol_E$Mean_Methylation <- rowMeans(test_enhancers[ ,9:13], na.rm=TRUE)
Chol_E$Tumor_Type <- "Breast Cancer"
Chol_E$Gene_Region <- "Enhancer"
Chol_E <- as.data.frame(Chol_E)
Chol_E_hypo <- subset(Chol_E, Mean_Methylation<=-2)
Chol_E_hyper <- subset(Chol_E, Mean_Methylation>=2)
Chol_E=NULL
Chol_E$Chromosome <- test_enhancers[ ,98]
Chol_E$Mean_Methylation <- rowMeans(test_enhancers[ ,73:81], na.rm=TRUE)
Chol_E$Tumor_Type <- "Cholangiocarcinoma"
Chol_E$Gene_Region <- "Enhancer"
Chol_E <- as.data.frame(Chol_E)
Chol_E_hypo <- subset(Chol_E, Mean_Methylation<=-2)
Chol_E_hyper <- subset(Chol_E, Mean_Methylation>=2)
Oligo_E=NULL
Oligo_E$Chromosome <- test_enhancers[ ,98]
Oligo_E$Mean_Methylation <- rowMeans(test_enhancers[ ,c(14:18, 82:96)], na.rm=TRUE)
Oligo_E$Tumor_Type <- "Oligodendroglioma"
Oligo_E$Gene_Region <- "Enhancer"
Oligo_E <- as.data.frame(Oligo_E)
Oligo_E_hypo <- subset(Oligo_E, Mean_Methylation<=-2)
Oligo_E_hyper <- subset(Oligo_E, Mean_Methylation>=2)
SNUC_E=NULL
SNUC_E$Chromosome <- test_enhancers[ ,98]
SNUC_E$Mean_Methylation <- rowMeans(test_enhancers[ ,19:26], na.rm=TRUE)
SNUC_E$Tumor_Type <- "SNUC"
SNUC_E$Gene_Region <- "Enhancer"
SNUC_E <- as.data.frame(SNUC_E)
SNUC_E_hypo <- subset(SNUC_E, Mean_Methylation<=-2)
SNUC_E_hyper <- subset(SNUC_E, Mean_Methylation>=2)
Blood_E=NULL
Blood_E$Chromosome <- control_enhancers[ ,55]
Blood_E$Mean_Methylation <- rowMeans(control_enhancers[ ,3:34], na.rm=TRUE)
Blood_E$Tumor_Type <- "Normal Blood"
Blood_E$Gene_Region <- "Enhancer"
Blood_E <- as.data.frame(Blood_E)
Blood_E_hypo <- subset(Blood_E, Mean_Methylation<=-2)
Blood_E_hyper <- subset(Blood_E, Mean_Methylation>=2)
GBM_E=NULL
GBM_E$Chromosome <- control_enhancers[ ,55]
GBM_E$Mean_Methylation <- rowMeans(control_enhancers[ ,35:40], na.rm=TRUE)
GBM_E$Tumor_Type <- "Pilocytic Astrocytoma"
GBM_E$Gene_Region <- "Enhancer"
GBM_E <- as.data.frame(GBM_E)
GBM_E_hypo <- subset(GBM_E, Mean_Methylation<=-2)
GBM_E_hyper <- subset(GBM_E, Mean_Methylation>=2)
Neuro_E=NULL
Neuro_E$Chromosome <- control_enhancers[ ,55]
Neuro_E$Mean_Methylation <- rowMeans(control_enhancers[ ,41:43], na.rm=TRUE)
Neuro_E$Tumor_Type <- "Neurocytoma"
Neuro_E$Gene_Region <- "Enhancer"
Neuro_E <- as.data.frame(Neuro_E)
Neuro_E_hypo <- subset(Neuro_E, Mean_Methylation<=-2)
Neuro_E_hyper <- subset(Neuro_E, Mean_Methylation>=2)
GBM_E=NULL
GBM_E$Chromosome <- control_enhancers[ ,55]
GBM_E$Mean_Methylation <- rowMeans(control_enhancers[ ,44:49], na.rm=TRUE)
GBM_E$Tumor_Type <- "GBM"
GBM_E$Gene_Region <- "Enhancer"
GBM_E <- as.data.frame(GBM_E)
GBM_E_hypo <- subset(GBM_E, Mean_Methylation<=-2)
GBM_E_hyper <- subset(GBM_E, Mean_Methylation>=2)
SUDEP_E=NULL
SUDEP_E$Chromosome <- control_enhancers[ ,55]
SUDEP_E$Mean_Methylation <- rowMeans(control_enhancers[ ,50:53], na.rm=TRUE)
SUDEP_E$Tumor_Type <- "SUDEP"
SUDEP_E$Gene_Region <- "Enhancer"
SUDEP_E <- as.data.frame(SUDEP_E)
SUDEP_E_hypo <- subset(SUDEP_E, Mean_Methylation<=-2)
SUDEP_E_hyper <- subset(SUDEP_E, Mean_Methylation>=2)

##Promoters
names(test_promoters)
names(control_promoters)

AML_P=NULL
AML_P$Chromosome <- test_promoters[ ,98]
AML_P$Mean_Methylation <- rowMeans(test_promoters[ ,27:47], na.rm=TRUE)
AML_P$Tumor_Type <- "AML"
AML_P$Gene_Region <- "TSS200"
AML_P <- as.data.frame(AML_P)
AML_P_hypo <- subset(AML_P, Mean_Methylation<=-2)
AML_P_hyper <- subset(AML_P, Mean_Methylation>=2)
Astro_P=NULL
Astro_P$Chromosome <- test_promoters[ ,98]
Astro_P$Mean_Methylation <- rowMeans(test_promoters[ ,c(3:8, 48:72)], na.rm=TRUE)
Astro_P$Tumor_Type <- "Astrocytoma"
Astro_P$Gene_Region <- "TSS200"
Astro_P <- as.data.frame(Astro_P)
Astro_P_hypo <- subset(Astro_P, Mean_Methylation<=-2)
Astro_P_hyper <- subset(Astro_P, Mean_Methylation>=2)
Breast_P=NULL
Breast_P$Chromosome <- test_promoters[ ,98]
Breast_P$Mean_Methylation <- rowMeans(test_promoters[ ,9:13], na.rm=TRUE)
Breast_P$Tumor_Type <- "Breast Cancer"
Breast_P$Gene_Region <- "TSS200"
Breast_P <- as.data.frame(Breast_P)
Breast_P_hypo <- subset(Breast_P, Mean_Methylation<=-2)
Breast_P_hyper <- subset(Breast_P, Mean_Methylation>=2)
Chol_P=NULL
Chol_P$Chromosome <- test_promoters[ ,98]
Chol_P$Mean_Methylation <- rowMeans(test_promoters[ ,73:81], na.rm=TRUE)
Chol_P$Tumor_Type <- "Cholangiocarcinoma"
Chol_P$Gene_Region <- "TSS200"
Chol_P <- as.data.frame(Chol_P)
Chol_P_hypo <- subset(Chol_P, Mean_Methylation<=-2)
Chol_P_hyper <- subset(Chol_P, Mean_Methylation>=2)
Oligo_P=NULL
Oligo_P$Chromosome <- test_promoters[ ,98]
Oligo_P$Mean_Methylation <- rowMeans(test_promoters[ ,c(14:18, 82:96)], na.rm=TRUE)
Oligo_P$Tumor_Type <- "Oligodendroglioma"
Oligo_P$Gene_Region <- "TSS200"
Oligo_P <- as.data.frame(Oligo_P)
Oligo_P_hypo <- subset(Oligo_P, Mean_Methylation<=-2)
Oligo_P_hyper <- subset(Oligo_P, Mean_Methylation>=2)
SNUC_P=NULL
SNUC_P$Chromosome <- test_promoters[ ,98]
SNUC_P$Mean_Methylation <- rowMeans(test_promoters[ ,19:26], na.rm=TRUE)
SNUC_P$Tumor_Type <- "SNUC"
SNUC_P$Gene_Region <- "TSS200"
SNUC_P <- as.data.frame(SNUC_P)
SNUC_P_hypo <- subset(SNUC_P, Mean_Methylation<=-2)
SNUC_P_hyper <- subset(SNUC_P, Mean_Methylation>=2)
Blood_P=NULL
Blood_P$Chromosome <- control_promoters[ ,55]
Blood_P$Mean_Methylation <- rowMeans(control_promoters[ ,3:34], na.rm=TRUE)
Blood_P$Tumor_Type <- "Normal Blood"
Blood_P$Gene_Region <- "TSS200"
Blood_P <- as.data.frame(Blood_P)
Blood_P_hypo <- subset(Blood_P, Mean_Methylation<=-2)
Blood_P_hyper <- subset(Blood_P, Mean_Methylation>=2)
PAstro_P=NULL
PAstro_P$Chromosome <- control_promoters[ ,55]
PAstro_P$Mean_Methylation <- rowMeans(control_promoters[ ,35:40], na.rm=TRUE)
PAstro_P$Tumor_Type <- "Pilocytic Astrocytoma"
PAstro_P$Gene_Region <- "TSS200"
PAstro_P <- as.data.frame(PAstro_P)
PAstro_P_hypo <- subset(PAstro_P, Mean_Methylation<=-2)
PAstro_P_hyper <- subset(PAstro_P, Mean_Methylation>=2)
Neuro_P=NULL
Neuro_P$Chromosome <- control_promoters[ ,55]
Neuro_P$Mean_Methylation <- rowMeans(control_promoters[ ,41:43], na.rm=TRUE)
Neuro_P$Tumor_Type <- "Neurocytoma"
Neuro_P$Gene_Region <- "TSS200"
Neuro_P <- as.data.frame(Neuro_P)
Neuro_P_hypo <- subset(Neuro_P, Mean_Methylation<=-2)
Neuro_P_hyper <- subset(Neuro_P, Mean_Methylation>=2)
GBM_P=NULL
GBM_P$Chromosome <- control_promoters[ ,55]
GBM_P$Mean_Methylation <- rowMeans(control_promoters[ ,44:49], na.rm=TRUE)
GBM_P$Tumor_Type <- "GBM"
GBM_P$Gene_Region <- "TSS200"
GBM_P <- as.data.frame(GBM_P)
GBM_P_hypo <- subset(GBM_P, Mean_Methylation<=-2)
GBM_P_hyper <- subset(GBM_P, Mean_Methylation>=2)
SUDEP_P=NULL
SUDEP_P$Chromosome <- control_promoters[ ,55]
SUDEP_P$Mean_Methylation <- rowMeans(control_promoters[ ,50:53], na.rm=TRUE)
SUDEP_P$Tumor_Type <- "SUDEP"
SUDEP_P$Gene_Region <- "TSS200"
SUDEP_P <- as.data.frame(SUDEP_P)
SUDEP_P_hypo <- subset(SUDEP_P, Mean_Methylation<=-2)
SUDEP_P_hyper <- subset(SUDEP_P, Mean_Methylation>=2)

##Bodies
names(test_bodies)
names(control_bodies)

AML_B=NULL
AML_B$Chromosome <- test_bodies[ ,98]
AML_B$Mean_Methylation <- rowMeans(test_bodies[ ,27:47], na.rm=TRUE)
AML_B$Tumor_Type <- "AML"
AML_B$Gene_Region <- "Body"
AML_B <- as.data.frame(AML_B)
AML_B_hypo <- subset(AML_B, Mean_Methylation<=-2)
AML_B_hyper <- subset(AML_B, Mean_Methylation>=2)
Astro_B=NULL
Astro_B$Chromosome <- test_bodies[ ,98]
Astro_B$Mean_Methylation <- rowMeans(test_bodies[ ,c(3:8, 48:72)], na.rm=TRUE)
Astro_B$Tumor_Type <- "Astrocytoma"
Astro_B$Gene_Region <- "Body"
Astro_B <- as.data.frame(Astro_B)
Astro_B_hypo <- subset(Astro_B, Mean_Methylation<=-2)
Astro_B_hyper <- subset(Astro_B, Mean_Methylation>=2)
Breast_B=NULL
Breast_B$Chromosome <- test_bodies[ ,98]
Breast_B$Mean_Methylation <- rowMeans(test_bodies[ ,9:13], na.rm=TRUE)
Breast_B$Tumor_Type <- "Breast Cancer"
Breast_B$Gene_Region <- "Body"
Breast_B <- as.data.frame(Breast_B)
Breast_B_hypo <- subset(Breast_B, Mean_Methylation<=-2)
Breast_B_hyper <- subset(Breast_B, Mean_Methylation>=2)
Chol_B=NULL
Chol_B$Chromosome <- test_bodies[ ,98]
Chol_B$Mean_Methylation <- rowMeans(test_bodies[ ,73:81], na.rm=TRUE)
Chol_B$Tumor_Type <- "Cholangiocarcinoma"
Chol_B$Gene_Region <- "Body"
Chol_B <- as.data.frame(Chol_B)
Chol_B_hypo <- subset(Chol_B, Mean_Methylation<=-2)
Chol_B_hyper <- subset(Chol_B, Mean_Methylation>=2)
Oligo_B=NULL
Oligo_B$Chromosome <- test_bodies[ ,98]
Oligo_B$Mean_Methylation <- rowMeans(test_bodies[ ,c(14:18, 82:96)], na.rm=TRUE)
Oligo_B$Tumor_Type <- "Oligodendroglioma"
Oligo_B$Gene_Region <- "Body"
Oligo_B <- as.data.frame(Oligo_B)
Oligo_B_hypo <- subset(Oligo_B, Mean_Methylation<=-2)
Oligo_B_hyper <- subset(Oligo_B, Mean_Methylation>=2)
SNUC_B=NULL
SNUC_B$Chromosome <- test_bodies[ ,98]
SNUC_B$Mean_Methylation <- rowMeans(test_bodies[ ,19:26], na.rm=TRUE)
SNUC_B$Tumor_Type <- "SNUC"
SNUC_B$Gene_Region <- "Body"
SNUC_B <- as.data.frame(SNUC_B)
SNUC_B_hypo <- subset(SNUC_B, Mean_Methylation<=-2)
SNUC_B_hyper <- subset(SNUC_B, Mean_Methylation>=2)
Blood_B=NULL
Blood_B$Chromosome <- control_bodies[ ,55]
Blood_B$Mean_Methylation <- rowMeans(control_bodies[ ,3:34], na.rm=TRUE)
Blood_B$Tumor_Type <- "Normal Blood"
Blood_B$Gene_Region <- "Body"
Blood_B <- as.data.frame(Blood_B)
Blood_B_hypo <- subset(Blood_B, Mean_Methylation<=-2)
Blood_B_hyper <- subset(Blood_B, Mean_Methylation>=2)
PAstro_B=NULL
PAstro_B$Chromosome <- control_bodies[ ,55]
PAstro_B$Mean_Methylation <- rowMeans(control_bodies[ ,35:40], na.rm=TRUE)
PAstro_B$Tumor_Type <- "Pilocytic Astrocytoma"
PAstro_B$Gene_Region <- "Body"
PAstro_B <- as.data.frame(PAstro_B)
PAstro_B_hypo <- subset(PAstro_B, Mean_Methylation<=-2)
PAstro_B_hyper <- subset(PAstro_B, Mean_Methylation>=2)
Neuro_B=NULL
Neuro_B$Chromosome <- control_bodies[ ,55]
Neuro_B$Mean_Methylation <- rowMeans(control_bodies[ ,41:43], na.rm=TRUE)
Neuro_B$Tumor_Type <- "Neurocytoma"
Neuro_B$Gene_Region <- "Body"
Neuro_B <- as.data.frame(Neuro_B)
Neuro_B_hypo <- subset(Neuro_B, Mean_Methylation<=-2)
Neuro_B_hyper <- subset(Neuro_B, Mean_Methylation>=2)
GBM_B=NULL
GBM_B$Chromosome <- control_bodies[ ,55]
GBM_B$Mean_Methylation <- rowMeans(control_bodies[ ,44:49], na.rm=TRUE)
GBM_B$Tumor_Type <- "GBM"
GBM_B$Gene_Region <- "Body"
GBM_B <- as.data.frame(GBM_B)
GBM_B_hypo <- subset(GBM_B, Mean_Methylation<=-2)
GBM_B_hyper <- subset(GBM_B, Mean_Methylation>=2)
SUDEP_B=NULL
SUDEP_B$Chromosome <- control_bodies[ ,55]
SUDEP_B$Mean_Methylation <- rowMeans(control_bodies[ ,50:53], na.rm=TRUE)
SUDEP_B$Tumor_Type <- "SUDEP"
SUDEP_B$Gene_Region <- "Body"
SUDEP_B <- as.data.frame(SUDEP_B)
SUDEP_B_hypo <- subset(SUDEP_B, Mean_Methylation<=-2)
SUDEP_B_hyper <- subset(SUDEP_B, Mean_Methylation>=2)


##Whole Genome
names(test_samples)
names(control_samples)

AML_WG=NULL
AML_WG$Chromosome <- test_samples[ ,97]
AML_WG$Mean_Methylation <- rowMeans(test_samples[ ,26:46], na.rm=TRUE)
AML_WG$Tumor_Type <- "AML"
AML_WG$Gene_Region <- "Whole Genome"
AML_WG <- as.data.frame(AML_WG)
AML_WG_hypo <- subset(AML_WG, Mean_Methylation<=-2)
AML_WG_hyper <- subset(AML_WG, Mean_Methylation>=2)
Astro_WG=NULL
Astro_WG$Chromosome <- test_samples[ ,97]
Astro_WG$Mean_Methylation <- rowMeans(test_samples[ ,c(2:7, 47:71)], na.rm=TRUE)
Astro_WG$Tumor_Type <- "Astrocytoma"
Astro_WG$Gene_Region <- "Whole Genome"
Astro_WG <- as.data.frame(Astro_WG)
Astro_WG_hypo <- subset(Astro_WG, Mean_Methylation<=-2)
Astro_WG_hyper <- subset(Astro_WG, Mean_Methylation>=2)
Chol_WG=NULL
Chol_WG$Chromosome <- test_samples[ ,97]
Chol_WG$Mean_Methylation <- rowMeans(test_samples[ ,8:12], na.rm=TRUE)
Chol_WG$Tumor_Type <- "Breast Cancer"
Chol_WG$Gene_Region <- "Whole Genome"
Chol_WG <- as.data.frame(Chol_WG)
Chol_WG_hypo <- subset(Chol_WG, Mean_Methylation<=-2)
Chol_WG_hyper <- subset(Chol_WG, Mean_Methylation>=2)
Chol_WG=NULL
Chol_WG$Chromosome <- test_samples[ ,97]
Chol_WG$Mean_Methylation <- rowMeans(test_samples[ ,72:80], na.rm=TRUE)
Chol_WG$Tumor_Type <- "Cholangiocarcinoma"
Chol_WG$Gene_Region <- "Whole Genome"
Chol_WG <- as.data.frame(Chol_WG)
Chol_WG_hypo <- subset(Chol_WG, Mean_Methylation<=-2)
Chol_WG_hyper <- subset(Chol_WG, Mean_Methylation>=2)
Oligo_WG=NULL
Oligo_WG$Chromosome <- test_samples[ ,97]
Oligo_WG$Mean_Methylation <- rowMeans(test_samples[ ,c(13:17, 81:95)], na.rm=TRUE)
Oligo_WG$Tumor_Type <- "Oligodendroglioma"
Oligo_WG$Gene_Region <- "Whole Genome"
Oligo_WG <- as.data.frame(Oligo_WG)
Oligo_WG_hypo <- subset(Oligo_WG, Mean_Methylation<=-2)
Oligo_WG_hyper <- subset(Oligo_WG, Mean_Methylation>=2)
SNUC_WG=NULL
SNUC_WG$Chromosome <- test_samples[ ,97]
SNUC_WG$Mean_Methylation <- rowMeans(test_samples[ ,18:25], na.rm=TRUE)
SNUC_WG$Tumor_Type <- "SNUC"
SNUC_WG$Gene_Region <- "Whole Genome"
SNUC_WG <- as.data.frame(SNUC_WG)
SNUC_WG_hypo <- subset(SNUC_WG, Mean_Methylation<=-2)
SNUC_WG_hyper <- subset(SNUC_WG, Mean_Methylation>=2)
Blood_WG=NULL
Blood_WG$Chromosome <- control_samples[ ,54]
Blood_WG$Mean_Methylation <- rowMeans(control_samples[ ,2:33], na.rm=TRUE)
Blood_WG$Tumor_Type <- "Normal Blood"
Blood_WG$Gene_Region <- "Whole Genome"
Blood_WG <- as.data.frame(Blood_WG)
Blood_WG_hypo <- subset(Blood_WG, Mean_Methylation<=-2)
Blood_WG_hyper <- subset(Blood_WG, Mean_Methylation>=2)
PAstro_WG=NULL
PAstro_WG$Chromosome <- control_samples[ ,54]
PAstro_WG$Mean_Methylation <- rowMeans(control_samples[ ,34:39], na.rm=TRUE)
PAstro_WG$Tumor_Type <- "Pilocytic Astrocytoma"
PAstro_WG$Gene_Region <- "Whole Genome"
PAstro_WG <- as.data.frame(PAstro_WG)
PAstro_WG_hypo <- subset(PAstro_WG, Mean_Methylation<=-2)
PAstro_WG_hyper <- subset(PAstro_WG, Mean_Methylation>=2)
Neuro_WG=NULL
Neuro_WG$Chromosome <- control_samples[ ,54]
Neuro_WG$Mean_Methylation <- rowMeans(control_samples[ ,40:42], na.rm=TRUE)
Neuro_WG$Tumor_Type <- "Neurocytoma"
Neuro_WG$Gene_Region <- "Whole Genome"
Neuro_WG <- as.data.frame(Neuro_WG)
Neuro_WG_hypo <- subset(Neuro_WG, Mean_Methylation<=-2)
Neuro_WG_hyper <- subset(Neuro_WG, Mean_Methylation>=2)
GBM_WG=NULL
GBM_WG$Chromosome <- control_samples[ ,54]
GBM_WG$Mean_Methylation <- rowMeans(control_samples[ ,43:48], na.rm=TRUE)
GBM_WG$Tumor_Type <- "GBM"
GBM_WG$Gene_Region <- "Whole Genome"
GBM_WG <- as.data.frame(GBM_WG)
GBM_WG_hypo <- subset(GBM_WG, Mean_Methylation<=-2)
GBM_WG_hyper <- subset(GBM_WG, Mean_Methylation>=2)
SUDEP_WG=NULL
SUDEP_WG$Chromosome <- control_samples[ ,54]
SUDEP_WG$Mean_Methylation <- rowMeans(control_samples[ ,49:52], na.rm=TRUE)
SUDEP_WG$Tumor_Type <- "SUDEP"
SUDEP_WG$Gene_Region <- "Whole Genome"
SUDEP_WG <- as.data.frame(SUDEP_WG)
SUDEP_WG_hypo <- subset(SUDEP_WG, Mean_Methylation<=-2)
SUDEP_WG_hyper <- subset(SUDEP_WG, Mean_Methylation>=2)


######Create the Mean Methylation Boxplots#####
library("ggplot2")
library("RColorBrewer")
library("limma")
library("reshape2")
library("gplots")
library("extrafont")

##Bind the rows of all the data frames
test_final_df=NULL
test_final_df <- rbind(AML_E, Astro_E, Chol_E, Chol_E, Oligo_E, SNUC_E)
test_final_df <- rbind(test_final_df,AML_P, Astro_P, Breast_P, Chol_P, Oligo_P, SNUC_P)
test_final_df <- rbind(test_final_df, AML_B, Astro_B, Breast_B, Chol_B, Oligo_B, SNUC_B)
test_final_df <- rbind(test_final_df, AML_WG, Astro_WG, Chol_WG, Chol_WG, Oligo_WG, SNUC_WG)
setwd("C:\\Users\\Ramona Bledea\\Documents\\PROJECT RESULTS\\IDH Tumor Project\\Results\\Test Samples - Differential Methylation\\Using M Values\\Global Methylation")
write.csv(final_df,file="Test_Datasheet_for_MeanMeth_Boxplots.csv")
control_final_df=NULL
control_final_df <- rbind(Blood_E, GBM_E, Neuro_E, GBM_E, SUDEP_E)
control_final_df <- rbind(control_final_df,Blood_P, PAstro_P, Neuro_P, GBM_P, SUDEP_P)
control_final_df <- rbind(control_final_df, Blood_B, PAstro_B, Neuro_B, GBM_B, SUDEP_B)
control_final_df <- rbind(control_final_df, Blood_WG, PAstro_WG, Neuro_WG, GBM_WG, SUDEP_WG)
setwd("C:\\Users\\Ramona Bledea\\Documents\\PROJECT RESULTS\\IDH Tumor Project\\Results\\Control Samples - Differential Methylation\\Using M Values\\Global Methylation")
write.csv(final_df,file="Control_Datasheet_for_MeanMeth_Boxplots.csv")

##Generate plots
setwd("C:\\Users\\Ramona Bledea\\Documents\\PROJECT RESULTS\\IDH Tumor Project\\Results\\Test Samples - Differential Methylation\\Using M Values\\Global Methylation")
png(file=paste("Test_MeanMeth_Boxplot.png", sep=""),width=1400,height=600,pointsize=80)
ggplot(test_final_df, aes(x=Tumor_Type, y=Mean_Methylation, fill=Gene_Region)) +
  geom_boxplot() +
  ggtitle("IDH Mutant Samples Mean Global Methylation") +
  xlab("Tumor Type") + ylab("Mean Methylation (M)") +
  theme(text=element_text(size=20, family="Arial"))
dev.off()
setwd("C:\\Users\\Ramona Bledea\\Documents\\PROJECT RESULTS\\IDH Tumor Project\\Results\\Control Samples - Differential Methylation\\Using M Values\\Global Methylation")
png(file=paste("Control_MeanMeth_Boxplot.png", sep=""),width=1400,height=600,pointsize=80)
ggplot(control_final_df, aes(x=Tumor_Type, y=Mean_Methylation, fill=Gene_Region)) +
  geom_boxplot() +
  ggtitle("IDH WT Control Samples Mean Global Methylation") +
  xlab("Tumor Type") + ylab("Mean Methylation (M)") +
  theme(text=element_text(size=20, family="Arial"))
dev.off()


######Generate the Probe Gene Region Distribution Barplots#####
##Take the dim of the all the subsetted dataframes and create an excel sheet with the info. Create 4 files total.
#HYPER#
dim(AML_E_hyper)
dim(Astro_E_hyper)
dim(Chol_E_hyper)
dim(Chol_E_hyper)
dim(Oligo_E_hyper)
dim(SNUC_E_hyper)
dim(Blood_E_hyper)
dim(GBM_E_hyper)
dim(Neuro_E_hyper)
dim(GBM_E_hyper)
dim(SUDEP_E_hyper)
dim(AML_P_hyper)
dim(Astro_P_hyper)
dim(Breast_P_hyper)
dim(Chol_P_hyper)
dim(Oligo_P_hyper)
dim(SNUC_P_hyper)
dim(Blood_P_hyper)
dim(PAstro_P_hyper)
dim(Neuro_P_hyper)
dim(GBM_P_hyper)
dim(SUDEP_P_hyper)
dim(AML_B_hyper)
dim(Astro_B_hyper)
dim(Breast_B_hyper)
dim(Chol_B_hyper)
dim(Oligo_B_hyper)
dim(SNUC_B_hyper)
dim(Blood_B_hyper)
dim(PAstro_B_hyper)
dim(Neuro_B_hyper)
dim(GBM_B_hyper)
dim(SUDEP_B_hyper)
dim(AML_WG_hyper)
dim(Astro_WG_hyper)
dim(Chol_WG_hyper)
dim(Chol_WG_hyper)
dim(Oligo_WG_hyper)
dim(SNUC_WG_hyper)
dim(Blood_WG_hyper)
dim(PAstro_WG_hyper)
dim(Neuro_WG_hyper)
dim(GBM_WG_hyper)
dim(SUDEP_WG_hyper)

#HYPO#
dim(AML_E_hypo)
dim(Astro_E_hypo)
dim(Chol_E_hypo)
dim(Chol_E_hypo)
dim(Oligo_E_hypo)
dim(SNUC_E_hypo)
dim(Blood_E_hypo)
dim(GBM_E_hypo)
dim(Neuro_E_hypo)
dim(GBM_E_hypo)
dim(SUDEP_E_hypo)
dim(AML_P_hypo)
dim(Astro_P_hypo)
dim(Breast_P_hypo)
dim(Chol_P_hypo)
dim(Oligo_P_hypo)
dim(SNUC_P_hypo)
dim(Blood_P_hypo)
dim(PAstro_P_hypo)
dim(Neuro_P_hypo)
dim(GBM_P_hypo)
dim(SUDEP_P_hypo)
dim(AML_B_hypo)
dim(Astro_B_hypo)
dim(Breast_B_hypo)
dim(Chol_B_hypo)
dim(Oligo_B_hypo)
dim(SNUC_B_hypo)
dim(Blood_B_hypo)
dim(PAstro_B_hypo)
dim(Neuro_B_hypo)
dim(GBM_B_hypo)
dim(SUDEP_B_hypo)
dim(AML_WG_hypo)
dim(Astro_WG_hypo)
dim(Chol_WG_hypo)
dim(Chol_WG_hypo)
dim(Oligo_WG_hypo)
dim(SNUC_WG_hypo)
dim(Blood_WG_hypo)
dim(PAstro_WG_hypo)
dim(Neuro_WG_hypo)
dim(GBM_WG_hypo)
dim(SUDEP_WG_hypo)

##Read in data files (or create dataframes using data.frame)
test_hypo_probes <-read.csv("C:\\Users\\Ramona Bledea\\Documents\\PROJECT RESULTS\\IDH Tumor Project\\Results\\Test Samples - Differential Methylation\\Using M Values\\Hypomethylated Probes\\Test_Hypo_Probes_Distribution.csv")
test_hyper_probes <-read.csv("C:\\Users\\Ramona Bledea\\Documents\\PROJECT RESULTS\\IDH Tumor Project\\Results\\Test Samples - Differential Methylation\\Using M Values\\Hypermethylated Probes\\Test_Hyper_Probes_Distribution.csv")
control_hypo_probes <-read.csv("C:\\Users\\Ramona Bledea\\Documents\\PROJECT RESULTS\\IDH Tumor Project\\Results\\Control Samples - Differential Methylation\\Using M Values\\Hypomethylated Probes\\Control_Hypo_Probes_Distribution.csv")
control_hyper_probes <-read.csv("C:\\Users\\Ramona Bledea\\Documents\\PROJECT RESULTS\\IDH Tumor Project\\Results\\Control Samples - Differential Methylation\\Using M Values\\Hypermethylated Probes\\Control_Hyper_Probes_Distribution.csv")

##Generate Plots
setwd("C:\\Users\\Ramona Bledea\\Documents\\PROJECT RESULTS\\IDH Tumor Project\\Results\\Test Samples - Differential Methylation\\Using M Values\\Hypomethylated Probes")
png(file=paste("Test_Hypo_Probe_Distribution_Barplot.png", sep=""),width=1000,height=700,pointsize=30)
ggplot(data=test_hypo_probes, aes(x=Tumor_Type, y=Number_of_Probes, fill=Gene_Region)) +
  geom_bar(stat="identity", position=position_dodge()) +
  ggtitle("IDH Mutant Samples Hypomethylated Probe Distribution by Region") +
  xlab("Tumor Type") + ylab("Number of Probes") +
  theme(text=element_text(size=15, family="Arial"))
dev.off()
setwd("C:\\Users\\Ramona Bledea\\Documents\\PROJECT RESULTS\\IDH Tumor Project\\Results\\Test Samples - Differential Methylation\\Using M Values\\Hypermethylated Probes")
png(file=paste("Test_Hyper_Probe_Distribution_Barplot.png", sep=""),width=1000,height=700,pointsize=30)
ggplot(data=test_hyper_probes, aes(x=Tumor_Type, y=Number_of_Probes, fill=Gene_Region)) +
  geom_bar(stat="identity", position=position_dodge())+
  ggtitle("IDH Mutant Samples Hypermethylated Probe Distribution by Region") +
  xlab("Tumor Type") + ylab("Number of Probes") +
  theme(text=element_text(size=15, family="Arial"))
dev.off()
setwd("C:\\Users\\Ramona Bledea\\Documents\\PROJECT RESULTS\\IDH Tumor Project\\Results\\Control Samples - Differential Methylation\\Using M Values\\Hypomethylated Probes")
png(file=paste("Control_Hypo_Probe_Distribution_Barplot.png", sep=""),width=1000,height=700,pointsize=30)
ggplot(data=control_hypo_probes, aes(x=Tumor_Type, y=Number_of_Probes, fill=Gene_Region)) +
  geom_bar(stat="identity", position=position_dodge())+
  ggtitle("IDH WT Samples Hypomethylated Probe Distribution by Region") +
  xlab("Tumor Type") + ylab("Number of Probes") +
  theme(text=element_text(size=15, family="Arial"))
dev.off()
setwd("C:\\Users\\Ramona Bledea\\Documents\\PROJECT RESULTS\\IDH Tumor Project\\Results\\Control Samples - Differential Methylation\\Using M Values\\Hypermethylated Probes")
png(file=paste("Control_Hyper_Probe_Distribution_Barplot.png", sep=""),width=1000,height=700,pointsize=30)
ggplot(data=control_hyper_probes, aes(x=Tumor_Type, y=Number_of_Probes, fill=Gene_Region)) +
  geom_bar(stat="identity", position=position_dodge())+
  ggtitle("IDH WT Samples Hypermethylated Probe Distribution by Region") +
  xlab("Tumor Type") + ylab("Number of Probes") +
  theme(text=element_text(size=15, family="Arial"))
dev.off()







