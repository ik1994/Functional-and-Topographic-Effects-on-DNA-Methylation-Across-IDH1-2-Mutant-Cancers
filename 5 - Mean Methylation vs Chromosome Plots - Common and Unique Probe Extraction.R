#########MEAN METHYLATION VS CHROMOSOME PLOTS (PER REGION AND WHOLE GENOME)
######EXTRACTING COMMON PROBES FOR HYPO AND HYPER#####

##Read in masterfiles##
test_annot <- read.csv("C:/.csv")
control_annot <- read.csv("C:/.csv")

######HYPOMETHYLATED PROBES ONLY#######
####TEST PROBES#####
##Read in all hypomethylated test probes lists##
subsetted_probes_list <- read.csv("C:/.csv")
subsetted_probes_list_Astro <- subsetted_probes_list$Astrocytoma
subsetted_probes_list_AML <- subsetted_probes_list$AML
subsetted_probes_list_Chol <- subsetted_probes_list$Cholangiocarcinoma
subsetted_probes_list_SNUC <- subsetted_probes_list$SNUC
subsetted_probes_list_Breast <- subsetted_probes_list$Breast
subsetted_probes_list_Oligo <- subsetted_probes_list$Oligodendroglioma
##Use intersect function to extract probes that are common to all 6 lists
test_probes_common_hypo=NULL
test_probes_common_hypo$ID <- intersect(subsetted_probes_list_Astro, subsetted_probes_list_AML)
test_probes_common_hypo$ID <- intersect(test_probes_common_hypo$ID, subsetted_probes_list_Chol)
test_probes_common_hypo$ID <- intersect(test_probes_common_hypo$ID, subsetted_probes_list_SNUC) 
test_probes_common_hypo$ID <- intersect(test_probes_common_hypo$ID, subsetted_probes_list_Breast)
test_probes_common_hypo$ID <- intersect(test_probes_common_hypo$ID, subsetted_probes_list_Oligo)
#DONT need to use brain for intersections above
##Merge with annotated M-value files and write to csv
test_probes_common_hypo_annot <- merge(test_probes_common_hypo, test_WG, by="ID")
dim(test_probes_common_hypo_annot)
setwd("C:/")
write.csv(test_probes_common_hypo_annot, file=".csv")



######HYPERMETHYLATED PROBES ONLY#######

####TEST PROBES#####
##Read in all hypermethylated test probes lists##
subsetted_probes_list <- read.csv("C:/.csv")
subsetted_probes_list_Astro <- subsetted_probes_list$Astrocytoma
subsetted_probes_list_AML <- subsetted_probes_list$AML
subsetted_probes_list_Chol <- subsetted_probes_list$Cholangiocarcinoma
subsetted_probes_list_SNUC <- subsetted_probes_list$SNUC
subsetted_probes_list_Breast <- subsetted_probes_list$Breast
subsetted_probes_list_Oligo <- subsetted_probes_list$Oligodendroglioma
##Use intersect function to extract probes that are common to all 6 lists
test_probes_common_hyper=NULL
test_probes_common_hyper$ID <- intersect(subsetted_probes_list_Astro, subsetted_probes_list_AML)
test_probes_common_hyper$ID <- intersect(test_probes_common_hyper$ID, subsetted_probes_list_Chol)
test_probes_common_hyper$ID <- intersect(test_probes_common_hyper$ID, subsetted_probes_list_SNUC) 
test_probes_common_hyper$ID <- intersect(test_probes_common_hyper$ID, subsetted_probes_list_Breast)
test_probes_common_hyper$ID <- intersect(test_probes_common_hyper$ID, subsetted_probes_list_Oligo) 
#DONT need to use brain for intersections above
##Merge with annotated M-value files and write to csv
test_probes_common_hyper_annot <- merge(test_probes_common_hyper, test_WG, by="ID")
setwd("C:/")
write.csv(test_probes_common_hyper_annot, file=".csv")




######IDENTIFY UNIQUE PROBES FOR EACH DISEASE FOR HYPER AND HYPO#######

##Read in probe lists and annotated m-value masterfiles
test_probes_list_hypo <- read.csv("C:/Users/Ramona Bledea/Documents/PROJECT RESULTS/IDH Tumor Project/Results/Test Samples - Differential Methylation/Using M Values/Hypomethylated Probes/test_samples_hypomethylated_probes_list.csv")
test_probes_list_hyper <- read.csv("C:/Users/Ramona Bledea/Documents/PROJECT RESULTS/IDH Tumor Project/Results/Test Samples - Differential Methylation/Using M Values/Hypermethylated Probes/test_samples_hypermethylated_probes_list.csv")
names(test_probes_list_hypo)
test_annot <- read.csv("C:/Users/Ramona Bledea/Documents/PROJECT RESULTS/IDH Tumor Project/Results/Test Samples - Differential Methylation/Using M Values/test_samples_diffmeth_annot_mvalues_masterfile.csv")

#Use setdiff(x,y) to determined probes that are in x but not in y
GBM_hypo_unique=NULL
GBM_hypo_unique$ID <- setdiff(control_probes_list_hypo$GBM, control_probes_list_hypo$SUDEP)
GBM_hypo_unique$ID <- setdiff(GBM_hypo_unique$ID, control_probes_list_hypo$Neurocytoma)
GBM_hypo_unique$ID <- setdiff(GBM_hypo_unique$ID, control_probes_list_hypo$Blood)
GBM_hypo_unique$ID <- setdiff(GBM_hypo_unique$ID, control_probes_list_hypo$Pilocytic_Astrocytoma)
GBM_hypo_unique <- merge(GBM_hypo_unique, control_annot, by="ID")
setwd("C:/Users/Ramona Bledea/Documents/PROJECT RESULTS/IDH Tumor Project/Results/Control Samples - Differential Methylation/Using M Values/hypomethylated Probes")
write.csv(GBM_hypo_unique, file="GBM_unique_hypo_probes.csv")

#####LOAD IN THE UNIQUE AND COMMON PROBE FILES CREATED ABOVE####
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
SNUC_E_hypo$ID <- SNUC_hypo_unique_enhancers[ ,2]
SNUC_E_hypo$POS <- SNUC_hypo_unique_enhancers[ ,57]
SNUC_E_hypo$UCSC_Gene <- SNUC_hypo_unique_enhancers[ ,79]
Blood_E_hypo$Mean_Methylation <- rowMeans(Blood_hypo_unique_enhancers[ ,4:35], na.rm=TRUE)
Blood_E_hypo$Tumor_Type <- "Normal Blood"
Blood_E_hypo$Gene_Region <- "Enhancer"
Blood_E_hypo <- as.data.frame(Blood_E_hypo)
PAstro_E_hypo=NULL
PAstro_E_hypo$Chromosome <- PAstro_hypo_unique_enhancers[ ,56]
SNUC_E_hypo$ID <- SNUC_hypo_unique_enhancers[ ,2]
SNUC_E_hypo$POS <- SNUC_hypo_unique_enhancers[ ,57]
SNUC_E_hypo$UCSC_Gene <- SNUC_hypo_unique_enhancers[ ,79]
PAstro_E_hypo$Mean_Methylation <- rowMeans(PAstro_hypo_unique_enhancers[ ,36:41], na.rm=TRUE)
PAstro_E_hypo$Tumor_Type <- "Pilocytic Astrocytoma"
PAstro_E_hypo$Gene_Region <- "Enhancer"
PAstro_E_hypo <- as.data.frame(PAstro_E_hypo)
Neuro_E_hypo=NULL
Neuro_E_hypo$Chromosome <- Neuro_hypo_unique_enhancers[ ,56]
SNUC_E_hypo$ID <- SNUC_hypo_unique_enhancers[ ,2]
SNUC_E_hypo$POS <- SNUC_hypo_unique_enhancers[ ,57]
SNUC_E_hypo$UCSC_Gene <- SNUC_hypo_unique_enhancers[ ,79]
Neuro_E_hypo$Mean_Methylation <- rowMeans(Neuro_hypo_unique_enhancers[ ,42:44], na.rm=TRUE)
Neuro_E_hypo$Tumor_Type <- "Neurocytoma"
Neuro_E_hypo$Gene_Region <- "Enhancer"
Neuro_E_hypo <- as.data.frame(Neuro_E_hypo)
GBM_E_hypo=NULL
GBM_E_hypo$Chromosome <- GBM_hypo_unique_enhancers[ ,56]
SNUC_E_hypo$ID <- SNUC_hypo_unique_enhancers[ ,2]
SNUC_E_hypo$POS <- SNUC_hypo_unique_enhancers[ ,57]
SNUC_E_hypo$UCSC_Gene <- SNUC_hypo_unique_enhancers[ ,79]
GBM_E_hypo$Mean_Methylation <- rowMeans(GBM_hypo_unique_enhancers[ ,45:50], na.rm=TRUE)
GBM_E_hypo$Tumor_Type <- "GBM"
GBM_E_hypo$Gene_Region <- "Enhancer"
GBM_E_hypo <- as.data.frame(GBM_E_hypo)
SUDEP_E_hypo=NULL
SUDEP_E_hypo$Chromosome <- SUDEP_hypo_unique_enhancers[ ,56]
SNUC_E_hypo$ID <- SNUC_hypo_unique_enhancers[ ,2]
SNUC_E_hypo$POS <- SNUC_hypo_unique_enhancers[ ,57]
SNUC_E_hypo$UCSC_Gene <- SNUC_hypo_unique_enhancers[ ,79]
SUDEP_E_hypo$Mean_Methylation <- rowMeans(SUDEP_hypo_unique_enhancers[ ,51:54], na.rm=TRUE)
SUDEP_E_hypo$Tumor_Type <- "SUDEP"
SUDEP_E_hypo$Gene_Region <- "Enhancer"
SUDEP_E_hypo <- as.data.frame(SUDEP_E_hypo)
common_control_E_hypo=NULL
common_control_E_hypo$Chromosome <- common_probes_control_hypo_annot_enhancers[ ,56]
SNUC_E_hypo$ID <- SNUC_hypo_unique_enhancers[ ,2]
SNUC_E_hypo$POS <- SNUC_hypo_unique_enhancers[ ,57]
SNUC_E_hypo$UCSC_Gene <- SNUC_hypo_unique_enhancers[ ,79]
common_control_E_hypo$Mean_Methylation <- rowMeans(common_probes_control_hypo_annot_enhancers[ ,4:54], na.rm=TRUE)
common_control_E_hypo$Tumor_Type <- "Common Probes, IDH WT"
common_control_E_hypo$Gene_Region <- "Enhancer"
common_control_E_hypo <- as.data.frame(common_control_E_hypo)
common_test_E_hypo=NULL
common_test_E_hypo$Chromosome <- common_probes_test_hypo_annot_enhancers[ ,98]
common_test_E_hypo$Mean_Methylation <- rowMeans(common_probes_test_hypo_annot_enhancers[ ,3:96], na.rm=TRUE)
common_test_E_hypo$Tumor_Type <- "Common Probes, IDH Mutant"
common_test_E_hypo$Gene_Region <- "Enhancer"
common_test_E_hypo <- as.data.frame(common_test_E_hypo)

AML_E_hyper=NULL
AML_E_hyper$Chromosome <- AML_hyper_unique_enhancers[ ,99]
AML_E_hyper$Mean_Methylation <- rowMeans(AML_hyper_unique_enhancers[ ,28:48], na.rm=TRUE)
AML_E_hyper$Tumor_Type <- "AML"
AML_E_hyper$Gene_Region <- "Enhancer"
AML_E_hyper <- as.data.frame(AML_E_hyper)
Astro_E_hyper=NULL
Astro_E_hyper$Chromosome <- Astro_hyper_unique_enhancers[ ,99]
Astro_E_hyper$Mean_Methylation <- rowMeans(Astro_hyper_unique_enhancers[ ,c(4:9, 49:73)], na.rm=TRUE)
Astro_E_hyper$Tumor_Type <- "Astrocytoma"
Astro_E_hyper$Gene_Region <- "Enhancer"
Astro_E_hyper <- as.data.frame(Astro_E_hyper)
Chol_E_hyper=NULL
Chol_E_hyper$Chromosome <- Chol_hyper_unique_enhancers[ ,99]
Chol_E_hyper$Mean_Methylation <- rowMeans(Chol_hyper_unique_enhancers[ ,74:82], na.rm=TRUE)
Chol_E_hyper$Tumor_Type <- "Cholangiocarcinoma"
Chol_E_hyper$Gene_Region <- "Enhancer"
Chol_E_hyper <- as.data.frame(Chol_E_hyper)
Breast_E_hyper=NULL
Breast_E_hyper$Chromosome <- Breast_hyper_unique_enhancers[ ,99]
Breast_E_hyper$Mean_Methylation <- rowMeans(Breast_hyper_unique_enhancers[ ,10:14], na.rm=TRUE)
Breast_E_hyper$Tumor_Type <- "Breast Cancer"
Breast_E_hyper$Gene_Region <- "Enhancer"
Breast_E_hyper <- as.data.frame(Breast_E_hyper)
Oligo_E_hyper=NULL
Oligo_E_hyper$Chromosome <- Oligo_hyper_unique_enhancers[ ,99]
Oligo_E_hyper$Mean_Methylation <- rowMeans(Oligo_hyper_unique_enhancers[ ,c(15:19, 83:97)], na.rm=TRUE)
Oligo_E_hyper$Tumor_Type <- "Oligodendroglioma"
Oligo_E_hyper$Gene_Region <- "Enhancer"
Oligo_E_hyper <- as.data.frame(Oligo_E_hyper)
SNUC_E_hyper=NULL
SNUC_E_hyper$Chromosome <- SNUC_hyper_unique_enhancers[ ,99]
SNUC_E_hyper$Mean_Methylation <- rowMeans(SNUC_hyper_unique_enhancers[ ,20:27], na.rm=TRUE)
SNUC_E_hyper$Tumor_Type <- "SNUC"
SNUC_E_hyper$Gene_Region <- "Enhancer"
SNUC_E_hyper <- as.data.frame(SNUC_E_hyper)
Blood_E_hyper=NULL
Blood_E_hyper$Chromosome <- Blood_hyper_unique_enhancers[ ,56]
Blood_E_hyper$Mean_Methylation <- rowMeans(Blood_hyper_unique_enhancers[ ,4:35], na.rm=TRUE)
Blood_E_hyper$Tumor_Type <- "Normal Blood"
Blood_E_hyper$Gene_Region <- "Enhancer"
Blood_E_hyper <- as.data.frame(Blood_E_hyper)
PAstro_E_hyper=NULL
PAstro_E_hyper$Chromosome <- PAstro_hyper_unique_enhancers[ ,56]
PAstro_E_hyper$Mean_Methylation <- rowMeans(PAstro_hyper_unique_enhancers[ ,36:41], na.rm=TRUE)
PAstro_E_hyper$Tumor_Type <- "Pilocytic Astrocytoma"
PAstro_E_hyper$Gene_Region <- "Enhancer"
PAstro_E_hyper <- as.data.frame(PAstro_E_hyper)
Neuro_E_hyper=NULL
Neuro_E_hyper$Chromosome <- Neuro_hyper_unique_enhancers[ ,56]
Neuro_E_hyper$Mean_Methylation <- rowMeans(Neuro_hyper_unique_enhancers[ ,42:44], na.rm=TRUE)
Neuro_E_hyper$Tumor_Type <- "Neurocytoma"
Neuro_E_hyper$Gene_Region <- "Enhancer"
Neuro_E_hyper <- as.data.frame(Neuro_E_hyper)
GBM_E_hyper=NULL
GBM_E_hyper$Chromosome <- GBM_hyper_unique_enhancers[ ,56]
GBM_E_hyper$Mean_Methylation <- rowMeans(GBM_hyper_unique_enhancers[ ,45:50], na.rm=TRUE)
GBM_E_hyper$Tumor_Type <- "GBM"
GBM_E_hyper$Gene_Region <- "Enhancer"
GBM_E_hyper <- as.data.frame(GBM_E_hyper)
SUDEP_E_hyper=NULL
SUDEP_E_hyper$Chromosome <- SUDEP_hyper_unique_enhancers[ ,56]
SUDEP_E_hyper$Mean_Methylation <- rowMeans(SUDEP_hyper_unique_enhancers[ ,51:54], na.rm=TRUE)
SUDEP_E_hyper$Tumor_Type <- "SUDEP"
SUDEP_E_hyper$Gene_Region <- "Enhancer"
SUDEP_E_hyper <- as.data.frame(SUDEP_E_hyper)
common_control_E_hyper=NULL
common_control_E_hyper$Chromosome <- common_probes_control_hyper_annot_enhancers[ ,56]
common_control_E_hyper$Mean_Methylation <- rowMeans(common_probes_control_hyper_annot_enhancers[ ,4:54], na.rm=TRUE)
common_control_E_hyper$Tumor_Type <- "Common Probes, IDH WT"
common_control_E_hyper$Gene_Region <- "Enhancer"
common_control_E_hyper <- as.data.frame(common_control_E_hyper)
common_test_E_hyper=NULL
common_test_E_hyper$Chromosome <- common_probes_test_hyper_annot_enhancers[ ,98]
common_test_E_hyper$Mean_Methylation <- rowMeans(common_probes_test_hyper_annot_enhancers[ ,3:96], na.rm=TRUE)
common_test_E_hyper$Tumor_Type <- "Common Probes, IDH Mutant"
common_test_E_hyper$Gene_Region <- "Enhancer"
common_test_E_hyper <- as.data.frame(common_test_E_hyper)

##Promoters
AML_P_hypo=NULL
AML_P_hypo$Chromosome <- AML_hypo_unique_promoters[ ,99]
AML_P_hypo$Mean_Methylation <- rowMeans(AML_hypo_unique_promoters[ ,28:48], na.rm=TRUE)
AML_P_hypo$Tumor_Type <- "AML"
AML_P_hypo$Gene_Region <- "TSS200"
AML_P_hypo <- as.data.frame(AML_P_hypo)
Astro_P_hypo=NULL
Astro_P_hypo$Chromosome <- Astro_hypo_unique_promoters[ ,99]
Astro_P_hypo$Mean_Methylation <- rowMeans(Astro_hypo_unique_promoters[ ,c(4:9, 49:73)], na.rm=TRUE)
Astro_P_hypo$Tumor_Type <- "Astrocytoma"
Astro_P_hypo$Gene_Region <- "TSS200"
Astro_P_hypo <- as.data.frame(Astro_P_hypo)
Chol_P_hypo=NULL
Chol_P_hypo$Chromosome <- Chol_hypo_unique_promoters[ ,99]
Chol_P_hypo$Mean_Methylation <- rowMeans(Chol_hypo_unique_promoters[ ,74:82], na.rm=TRUE)
Chol_P_hypo$Tumor_Type <- "Cholangiocarcinoma"
Chol_P_hypo$Gene_Region <- "TSS200"
Chol_P_hypo <- as.data.frame(Chol_P_hypo)
Breast_P_hypo=NULL
Breast_P_hypo$Chromosome <- Breast_hypo_unique_promoters[ ,99]
Breast_P_hypo$Mean_Methylation <- rowMeans(Breast_hypo_unique_promoters[ ,10:14], na.rm=TRUE)
Breast_P_hypo$Tumor_Type <- "Breast Cancer"
Breast_P_hypo$Gene_Region <- "TSS200"
Breast_P_hypo <- as.data.frame(Breast_P_hypo)
Oligo_P_hypo=NULL
Oligo_P_hypo$Chromosome <- Oligo_hypo_unique_promoters[ ,99]
Oligo_P_hypo$Mean_Methylation <- rowMeans(Oligo_hypo_unique_promoters[ ,c(15:19, 83:97)], na.rm=TRUE)
Oligo_P_hypo$Tumor_Type <- "Oligodendroglioma"
Oligo_P_hypo$Gene_Region <- "TSS200"
Oligo_P_hypo <- as.data.frame(Oligo_P_hypo)
SNUC_P_hypo=NULL
SNUC_P_hypo$Chromosome <- SNUC_hypo_unique_promoters[ ,99]
SNUC_P_hypo$Mean_Methylation <- rowMeans(SNUC_hypo_unique_promoters[ ,20:27], na.rm=TRUE)
SNUC_P_hypo$Tumor_Type <- "SNUC"
SNUC_P_hypo$Gene_Region <- "TSS200"
SNUC_P_hypo <- as.data.frame(SNUC_P_hypo)
Blood_P_hypo=NULL
Blood_P_hypo$Chromosome <- Blood_hypo_unique_promoters[ ,56]
Blood_P_hypo$Mean_Methylation <- rowMeans(Blood_hypo_unique_promoters[ ,4:35], na.rm=TRUE)
Blood_P_hypo$Tumor_Type <- "Normal Blood"
Blood_P_hypo$Gene_Region <- "TSS200"
Blood_P_hypo <- as.data.frame(Blood_P_hypo)
PAstro_P_hypo=NULL
PAstro_P_hypo$Chromosome <- PAstro_hypo_unique_promoters[ ,56]
PAstro_P_hypo$Mean_Methylation <- rowMeans(PAstro_hypo_unique_promoters[ ,36:41], na.rm=TRUE)
PAstro_P_hypo$Tumor_Type <- "Pilocytic Astrocytoma"
PAstro_P_hypo$Gene_Region <- "TSS200"
PAstro_P_hypo <- as.data.frame(PAstro_P_hypo)
Neuro_P_hypo=NULL
Neuro_P_hypo$Chromosome <- Neuro_hypo_unique_promoters[ ,56]
Neuro_P_hypo$Mean_Methylation <- rowMeans(Neuro_hypo_unique_promoters[ ,42:44], na.rm=TRUE)
Neuro_P_hypo$Tumor_Type <- "Neurocytoma"
Neuro_P_hypo$Gene_Region <- "TSS200"
Neuro_P_hypo <- as.data.frame(Neuro_P_hypo)
GBM_P_hypo=NULL
GBM_P_hypo$Chromosome <- GBM_hypo_unique_promoters[ ,56]
GBM_P_hypo$Mean_Methylation <- rowMeans(GBM_hypo_unique_promoters[ ,45:50], na.rm=TRUE)
GBM_P_hypo$Tumor_Type <- "GBM"
GBM_P_hypo$Gene_Region <- "TSS200"
GBM_P_hypo <- as.data.frame(GBM_P_hypo)
SUDEP_P_hypo=NULL
SUDEP_P_hypo$Chromosome <- SUDEP_hypo_unique_promoters[ ,56]
SUDEP_P_hypo$Mean_Methylation <- rowMeans(SUDEP_hypo_unique_promoters[ ,51:54], na.rm=TRUE)
SUDEP_P_hypo$Tumor_Type <- "SUDEP"
SUDEP_P_hypo$Gene_Region <- "TSS200"
SUDEP_P_hypo <- as.data.frame(SUDEP_P_hypo)
common_control_P_hypo=NULL
common_control_P_hypo$Chromosome <- common_probes_control_hypo_annot_promoters[ ,56]
common_control_P_hypo$Mean_Methylation <- rowMeans(common_probes_control_hypo_annot_promoters[ ,4:54], na.rm=TRUE)
common_control_P_hypo$Tumor_Type <- "Common Probes, IDH WT"
common_control_P_hypo$Gene_Region <- "TSS200"
common_control_P_hypo <- as.data.frame(common_control_P_hypo)
common_test_P_hypo=NULL
common_test_P_hypo$Chromosome <- common_probes_test_hypo_annot_promoters[ ,98]
common_test_P_hypo$Mean_Methylation <- rowMeans(common_probes_test_hypo_annot_promoters[ ,3:96], na.rm=TRUE)
common_test_P_hypo$Tumor_Type <- "Common Probes, IDH Mutant"
common_test_P_hypo$Gene_Region <- "TSS200"
common_test_P_hypo <- as.data.frame(common_test_P_hypo)

AML_P_hyper=NULL
AML_P_hyper$Chromosome <- AML_hyper_unique_promoters[ ,99]
AML_P_hyper$Mean_Methylation <- rowMeans(AML_hyper_unique_promoters[ ,28:48], na.rm=TRUE)
AML_P_hyper$Tumor_Type <- "AML"
AML_P_hyper$Gene_Region <- "TSS200"
AML_P_hyper <- as.data.frame(AML_P_hyper)
Astro_P_hyper=NULL
Astro_P_hyper$Chromosome <- Astro_hyper_unique_promoters[ ,99]
Astro_P_hyper$Mean_Methylation <- rowMeans(Astro_hyper_unique_promoters[ ,c(4:9, 49:73)], na.rm=TRUE)
Astro_P_hyper$Tumor_Type <- "Astrocytoma"
Astro_P_hyper$Gene_Region <- "TSS200"
Astro_P_hyper <- as.data.frame(Astro_P_hyper)
Chol_P_hyper=NULL
Chol_P_hyper$Chromosome <- Chol_hyper_unique_promoters[ ,99]
Chol_P_hyper$Mean_Methylation <- rowMeans(Chol_hyper_unique_promoters[ ,74:82], na.rm=TRUE)
Chol_P_hyper$Tumor_Type <- "Cholangiocarcinoma"
Chol_P_hyper$Gene_Region <- "TSS200"
Chol_P_hyper <- as.data.frame(Chol_P_hyper)
Breast_P_hyper=NULL
Breast_P_hyper$Chromosome <- Breast_hyper_unique_promoters[ ,99]
Breast_P_hyper$Mean_Methylation <- rowMeans(Breast_hyper_unique_promoters[ ,10:14], na.rm=TRUE)
Breast_P_hyper$Tumor_Type <- "Breast Cancer"
Breast_P_hyper$Gene_Region <- "TSS200"
Breast_P_hyper <- as.data.frame(Breast_P_hyper)
Oligo_P_hyper=NULL
Oligo_P_hyper$Chromosome <- Oligo_hyper_unique_promoters[ ,99]
Oligo_P_hyper$Mean_Methylation <- rowMeans(Oligo_hyper_unique_promoters[ ,c(15:19, 83:97)], na.rm=TRUE)
Oligo_P_hyper$Tumor_Type <- "Oligodendroglioma"
Oligo_P_hyper$Gene_Region <- "TSS200"
Oligo_P_hyper <- as.data.frame(Oligo_P_hyper)
SNUC_P_hyper=NULL
SNUC_P_hyper$Chromosome <- SNUC_hyper_unique_promoters[ ,99]
SNUC_P_hyper$Mean_Methylation <- rowMeans(SNUC_hyper_unique_promoters[ ,20:27], na.rm=TRUE)
SNUC_P_hyper$Tumor_Type <- "SNUC"
SNUC_P_hyper$Gene_Region <- "TSS200"
SNUC_P_hyper <- as.data.frame(SNUC_P_hyper)
Blood_P_hyper=NULL
Blood_P_hyper$Chromosome <- Blood_hyper_unique_promoters[ ,56]
Blood_P_hyper$Mean_Methylation <- rowMeans(Blood_hyper_unique_promoters[ ,4:35], na.rm=TRUE)
Blood_P_hyper$Tumor_Type <- "Normal Blood"
Blood_P_hyper$Gene_Region <- "TSS200"
Blood_P_hyper <- as.data.frame(Blood_P_hyper)
PAstro_P_hyper=NULL
PAstro_P_hyper$Chromosome <- PAstro_hyper_unique_promoters[ ,56]
PAstro_P_hyper$Mean_Methylation <- rowMeans(PAstro_hyper_unique_promoters[ ,36:41], na.rm=TRUE)
PAstro_P_hyper$Tumor_Type <- "Pilocytic Astrocytoma"
PAstro_P_hyper$Gene_Region <- "TSS200"
PAstro_P_hyper <- as.data.frame(PAstro_P_hyper)
Neuro_P_hyper=NULL
Neuro_P_hyper$Chromosome <- Neuro_hyper_unique_promoters[ ,56]
Neuro_P_hyper$Mean_Methylation <- rowMeans(Neuro_hyper_unique_promoters[ ,42:44], na.rm=TRUE)
Neuro_P_hyper$Tumor_Type <- "Neurocytoma"
Neuro_P_hyper$Gene_Region <- "TSS200"
Neuro_P_hyper <- as.data.frame(Neuro_P_hyper)
GBM_P_hyper=NULL
GBM_P_hyper$Chromosome <- GBM_hyper_unique_promoters[ ,56]
GBM_P_hyper$Mean_Methylation <- rowMeans(GBM_hyper_unique_promoters[ ,45:50], na.rm=TRUE)
GBM_P_hyper$Tumor_Type <- "GBM"
GBM_P_hyper$Gene_Region <- "TSS200"
GBM_P_hyper <- as.data.frame(GBM_P_hyper)
SUDEP_P_hyper=NULL
SUDEP_P_hyper$Chromosome <- SUDEP_hyper_unique_promoters[ ,56]
SUDEP_P_hyper$Mean_Methylation <- rowMeans(SUDEP_hyper_unique_promoters[ ,51:54], na.rm=TRUE)
SUDEP_P_hyper$Tumor_Type <- "SUDEP"
SUDEP_P_hyper$Gene_Region <- "TSS200"
SUDEP_P_hyper <- as.data.frame(SUDEP_P_hyper)
common_control_P_hyper=NULL
common_control_P_hyper$Chromosome <- common_probes_control_hyper_annot_promoters[ ,56]
common_control_P_hyper$Mean_Methylation <- rowMeans(common_probes_control_hyper_annot_promoters[ ,4:54], na.rm=TRUE)
common_control_P_hyper$Tumor_Type <- "Common Probes, IDH WT"
common_control_P_hyper$Gene_Region <- "TSS200"
common_control_P_hyper <- as.data.frame(common_control_P_hyper)
common_test_P_hyper=NULL
common_test_P_hyper$Chromosome <- common_probes_test_hyper_annot_promoters[ ,98]
common_test_P_hyper$Mean_Methylation <- rowMeans(common_probes_test_hyper_annot_promoters[ ,3:96], na.rm=TRUE)
common_test_P_hyper$Tumor_Type <- "Common Probes, IDH Mutant"
common_test_P_hyper$Gene_Region <- "TSS200"
common_test_P_hyper <- as.data.frame(common_test_P_hyper)


##Bodies
AML_B_hypo=NULL
AML_B_hypo$Chromosome <- AML_hypo_unique_bodies[ ,99]
AML_B_hypo$Mean_Methylation <- rowMeans(AML_hypo_unique_bodies[ ,28:48], na.rm=TRUE)
AML_B_hypo$Tumor_Type <- "AML"
AML_B_hypo$Gene_Region <- "Body"
AML_B_hypo <- as.data.frame(AML_B_hypo)
Astro_B_hypo=NULL
Astro_B_hypo$Chromosome <- Astro_hypo_unique_bodies[ ,99]
Astro_B_hypo$Mean_Methylation <- rowMeans(Astro_hypo_unique_bodies[ ,c(4:9, 49:73)], na.rm=TRUE)
Astro_B_hypo$Tumor_Type <- "Astrocytoma"
Astro_B_hypo$Gene_Region <- "Body"
Astro_B_hypo <- as.data.frame(Astro_B_hypo)
Chol_B_hypo=NULL
Chol_B_hypo$Chromosome <- Chol_hypo_unique_bodies[ ,99]
Chol_B_hypo$Mean_Methylation <- rowMeans(Chol_hypo_unique_bodies[ ,74:82], na.rm=TRUE)
Chol_B_hypo$Tumor_Type <- "Cholangiocarcinoma"
Chol_B_hypo$Gene_Region <- "Body"
Chol_B_hypo <- as.data.frame(Chol_B_hypo)
Breast_B_hypo=NULL
Breast_B_hypo$Chromosome <- Breast_hypo_unique_bodies[ ,99]
Breast_B_hypo$Mean_Methylation <- rowMeans(Breast_hypo_unique_bodies[ ,10:14], na.rm=TRUE)
Breast_B_hypo$Tumor_Type <- "Breast Cancer"
Breast_B_hypo$Gene_Region <- "Body"
Breast_B_hypo <- as.data.frame(Breast_B_hypo)
Oligo_B_hypo=NULL
Oligo_B_hypo$Chromosome <- Oligo_hypo_unique_bodies[ ,99]
Oligo_B_hypo$Mean_Methylation <- rowMeans(Oligo_hypo_unique_bodies[ ,c(15:19, 83:97)], na.rm=TRUE)
Oligo_B_hypo$Tumor_Type <- "Oligodendroglioma"
Oligo_B_hypo$Gene_Region <- "Body"
Oligo_B_hypo <- as.data.frame(Oligo_B_hypo)
SNUC_B_hypo=NULL
SNUC_B_hypo$Chromosome <- SNUC_hypo_unique_bodies[ ,99]
SNUC_B_hypo$Mean_Methylation <- rowMeans(SNUC_hypo_unique_bodies[ ,20:27], na.rm=TRUE)
SNUC_B_hypo$Tumor_Type <- "SNUC"
SNUC_B_hypo$Gene_Region <- "Body"
SNUC_B_hypo <- as.data.frame(SNUC_B_hypo)
Blood_B_hypo=NULL
Blood_B_hypo$Chromosome <- Blood_hypo_unique_bodies[ ,56]
Blood_B_hypo$Mean_Methylation <- rowMeans(Blood_hypo_unique_bodies[ ,4:35], na.rm=TRUE)
Blood_B_hypo$Tumor_Type <- "Normal Blood"
Blood_B_hypo$Gene_Region <- "Body"
Blood_B_hypo <- as.data.frame(Blood_B_hypo)
PAstro_B_hypo=NULL
PAstro_B_hypo$Chromosome <- PAstro_hypo_unique_bodies[ ,56]
PAstro_B_hypo$Mean_Methylation <- rowMeans(PAstro_hypo_unique_bodies[ ,36:41], na.rm=TRUE)
PAstro_B_hypo$Tumor_Type <- "Pilocytic Astrocytoma"
PAstro_B_hypo$Gene_Region <- "Body"
PAstro_B_hypo <- as.data.frame(PAstro_B_hypo)
Neuro_B_hypo=NULL
Neuro_B_hypo$Chromosome <- Neuro_hypo_unique_bodies[ ,56]
Neuro_B_hypo$Mean_Methylation <- rowMeans(Neuro_hypo_unique_bodies[ ,42:44], na.rm=TRUE)
Neuro_B_hypo$Tumor_Type <- "Neurocytoma"
Neuro_B_hypo$Gene_Region <- "Body"
Neuro_B_hypo <- as.data.frame(Neuro_B_hypo)
GBM_B_hypo=NULL
GBM_B_hypo$Chromosome <- GBM_hypo_unique_bodies[ ,56]
GBM_B_hypo$Mean_Methylation <- rowMeans(GBM_hypo_unique_bodies[ ,45:50], na.rm=TRUE)
GBM_B_hypo$Tumor_Type <- "GBM"
GBM_B_hypo$Gene_Region <- "Body"
GBM_B_hypo <- as.data.frame(GBM_B_hypo)
SUDEP_B_hypo=NULL
SUDEP_B_hypo$Chromosome <- SUDEP_hypo_unique_bodies[ ,56]
SUDEP_B_hypo$Mean_Methylation <- rowMeans(SUDEP_hypo_unique_bodies[ ,51:54], na.rm=TRUE)
SUDEP_B_hypo$Tumor_Type <- "SUDEP"
SUDEP_B_hypo$Gene_Region <- "Body"
SUDEP_B_hypo <- as.data.frame(SUDEP_B_hypo)
common_control_B_hypo=NULL
common_control_B_hypo$Chromosome <- common_probes_control_hypo_annot_bodies[ ,56]
common_control_B_hypo$Mean_Methylation <- rowMeans(common_probes_control_hypo_annot_bodies[ ,4:54], na.rm=TRUE)
common_control_B_hypo$Tumor_Type <- "Common Probes, IDH WT"
common_control_B_hypo$Gene_Region <- "Body"
common_control_B_hypo <- as.data.frame(common_control_B_hypo)
common_test_B_hypo=NULL
common_test_B_hypo$Chromosome <- common_probes_test_hypo_annot_bodies[ ,98]
common_test_B_hypo$Mean_Methylation <- rowMeans(common_probes_test_hypo_annot_bodies[ ,3:96], na.rm=TRUE)
common_test_B_hypo$Tumor_Type <- "Common Probes, IDH Mutant"
common_test_B_hypo$Gene_Region <- "Body"
common_test_B_hypo <- as.data.frame(common_test_B_hypo)

AML_B_hyper=NULL
AML_B_hyper$Chromosome <- AML_hyper_unique_bodies[ ,99]
AML_B_hyper$Mean_Methylation <- rowMeans(AML_hyper_unique_bodies[ ,28:48], na.rm=TRUE)
AML_B_hyper$Tumor_Type <- "AML"
AML_B_hyper$Gene_Region <- "Body"
AML_B_hyper <- as.data.frame(AML_B_hyper)
Astro_B_hyper=NULL
Astro_B_hyper$Chromosome <- Astro_hyper_unique_bodies[ ,99]
Astro_B_hyper$Mean_Methylation <- rowMeans(Astro_hyper_unique_bodies[ ,c(4:9, 49:73)], na.rm=TRUE)
Astro_B_hyper$Tumor_Type <- "Astrocytoma"
Astro_B_hyper$Gene_Region <- "Body"
Astro_B_hyper <- as.data.frame(Astro_B_hyper)
Chol_B_hyper=NULL
Chol_B_hyper$Chromosome <- Chol_hyper_unique_bodies[ ,99]
Chol_B_hyper$Mean_Methylation <- rowMeans(Chol_hyper_unique_bodies[ ,74:82], na.rm=TRUE)
Chol_B_hyper$Tumor_Type <- "Cholangiocarcinoma"
Chol_B_hyper$Gene_Region <- "Body"
Chol_B_hyper <- as.data.frame(Chol_B_hyper)
Breast_B_hyper=NULL
Breast_B_hyper$Chromosome <- Breast_hyper_unique_bodies[ ,99]
Breast_B_hyper$Mean_Methylation <- rowMeans(Breast_hyper_unique_bodies[ ,10:14], na.rm=TRUE)
Breast_B_hyper$Tumor_Type <- "Breast Cancer"
Breast_B_hyper$Gene_Region <- "Body"
Breast_B_hyper <- as.data.frame(Breast_B_hyper)
Oligo_B_hyper=NULL
Oligo_B_hyper$Chromosome <- Oligo_hyper_unique_bodies[ ,99]
Oligo_B_hyper$Mean_Methylation <- rowMeans(Oligo_hyper_unique_bodies[ ,c(15:19, 83:97)], na.rm=TRUE)
Oligo_B_hyper$Tumor_Type <- "Oligodendroglioma"
Oligo_B_hyper$Gene_Region <- "Body"
Oligo_B_hyper <- as.data.frame(Oligo_B_hyper)
SNUC_B_hyper=NULL
SNUC_B_hyper$Chromosome <- SNUC_hyper_unique_bodies[ ,99]
SNUC_B_hyper$Mean_Methylation <- rowMeans(SNUC_hyper_unique_bodies[ ,20:27], na.rm=TRUE)
SNUC_B_hyper$Tumor_Type <- "SNUC"
SNUC_B_hyper$Gene_Region <- "Body"
SNUC_B_hyper <- as.data.frame(SNUC_B_hyper)
Blood_B_hyper=NULL
Blood_B_hyper$Chromosome <- Blood_hyper_unique_bodies[ ,56]
Blood_B_hyper$Mean_Methylation <- rowMeans(Blood_hyper_unique_bodies[ ,4:35], na.rm=TRUE)
Blood_B_hyper$Tumor_Type <- "Normal Blood"
Blood_B_hyper$Gene_Region <- "Body"
Blood_B_hyper <- as.data.frame(Blood_B_hyper)
PAstro_B_hyper=NULL
PAstro_B_hyper$Chromosome <- PAstro_hyper_unique_bodies[ ,56]
PAstro_B_hyper$Mean_Methylation <- rowMeans(PAstro_hyper_unique_bodies[ ,36:41], na.rm=TRUE)
PAstro_B_hyper$Tumor_Type <- "Pilocytic Astrocytoma"
PAstro_B_hyper$Gene_Region <- "Body"
PAstro_B_hyper <- as.data.frame(PAstro_B_hyper)
Neuro_B_hyper=NULL
Neuro_B_hyper$Chromosome <- Neuro_hyper_unique_bodies[ ,56]
Neuro_B_hyper$Mean_Methylation <- rowMeans(Neuro_hyper_unique_bodies[ ,42:44], na.rm=TRUE)
Neuro_B_hyper$Tumor_Type <- "Neurocytoma"
Neuro_B_hyper$Gene_Region <- "Body"
Neuro_B_hyper <- as.data.frame(Neuro_B_hyper)
GBM_B_hyper=NULL
GBM_B_hyper$Chromosome <- GBM_hyper_unique_bodies[ ,56]
GBM_B_hyper$Mean_Methylation <- rowMeans(GBM_hyper_unique_bodies[ ,45:50], na.rm=TRUE)
GBM_B_hyper$Tumor_Type <- "GBM"
GBM_B_hyper$Gene_Region <- "Body"
GBM_B_hyper <- as.data.frame(GBM_B_hyper)
SUDEP_B_hyper=NULL
SUDEP_B_hyper$Chromosome <- SUDEP_hyper_unique_bodies[ ,56]
SUDEP_B_hyper$Mean_Methylation <- rowMeans(SUDEP_hyper_unique_bodies[ ,51:54], na.rm=TRUE)
SUDEP_B_hyper$Tumor_Type <- "SUDEP"
SUDEP_B_hyper$Gene_Region <- "Body"
SUDEP_B_hyper <- as.data.frame(SUDEP_B_hyper)
common_control_B_hyper=NULL
common_control_B_hyper$Chromosome <- common_probes_control_hyper_annot_bodies[ ,56]
common_control_B_hyper$Mean_Methylation <- rowMeans(common_probes_control_hyper_annot_bodies[ ,4:54], na.rm=TRUE)
common_control_B_hyper$Tumor_Type <- "Common Probes, IDH WT"
common_control_B_hyper$Gene_Region <- "Body"
common_control_B_hyper <- as.data.frame(common_control_B_hyper)
common_test_B_hyper=NULL
common_test_B_hyper$Chromosome <- common_probes_test_hyper_annot_bodies[ ,98]
common_test_B_hyper$Mean_Methylation <- rowMeans(common_probes_test_hyper_annot_bodies[ ,3:96], na.rm=TRUE)
common_test_B_hyper$Tumor_Type <- "Common Probes, IDH Mutant"
common_test_B_hyper$Gene_Region <- "Body"
common_test_B_hyper <- as.data.frame(common_test_B_hyper)


##Find the mean of the Mean_Methylation column values in each data frame, for each chromosome
#subset for each chromosome

#####################################TEST#############################################
E_1 <- AML_B_hyper[AML_B_hyper$Chromosome == "chr1",]
E_2 <- AML_B_hyper[AML_B_hyper$Chromosome == "chr2",]
E_3 <- AML_B_hyper[AML_B_hyper$Chromosome == "chr3",]
E_4 <- AML_B_hyper[AML_B_hyper$Chromosome == "chr4",]
E_5 <- AML_B_hyper[AML_B_hyper$Chromosome == "chr5",]
E_6 <- AML_B_hyper[AML_B_hyper$Chromosome == "chr6",]
E_7 <- AML_B_hyper[AML_B_hyper$Chromosome == "chr7",]
E_8 <- AML_B_hyper[AML_B_hyper$Chromosome == "chr8",]
E_9 <- AML_B_hyper[AML_B_hyper$Chromosome == "chr9",]
E_10 <- AML_B_hyper[AML_B_hyper$Chromosome == "chr10",]
E_11 <- AML_B_hyper[AML_B_hyper$Chromosome == "chr11",]
E_12 <- AML_B_hyper[AML_B_hyper$Chromosome == "chr12",]
E_13 <- AML_B_hyper[AML_B_hyper$Chromosome == "chr13",]
E_14 <- AML_B_hyper[AML_B_hyper$Chromosome == "chr14",]
E_15 <- AML_B_hyper[AML_B_hyper$Chromosome == "chr15",]
E_16 <- AML_B_hyper[AML_B_hyper$Chromosome == "chr16",]
E_17 <- AML_B_hyper[AML_B_hyper$Chromosome == "chr17",]
E_18 <- AML_B_hyper[AML_B_hyper$Chromosome == "chr18",]
E_19 <- AML_B_hyper[AML_B_hyper$Chromosome == "chr19",]
E_20 <- AML_B_hyper[AML_B_hyper$Chromosome == "chr20",]
E_21 <- AML_B_hyper[AML_B_hyper$Chromosome == "chr21",]
E_22 <- AML_B_hyper[AML_B_hyper$Chromosome == "chr22",]

E_23 <- Astro_B_hyper[Astro_B_hyper$Chromosome == "chr1",]
E_24 <- Astro_B_hyper[Astro_B_hyper$Chromosome == "chr2",]
E_25 <- Astro_B_hyper[Astro_B_hyper$Chromosome == "chr3",]
E_26 <- Astro_B_hyper[Astro_B_hyper$Chromosome == "chr4",]
E_27 <- Astro_B_hyper[Astro_B_hyper$Chromosome == "chr5",]
E_28 <- Astro_B_hyper[Astro_B_hyper$Chromosome == "chr6",]
E_29 <- Astro_B_hyper[Astro_B_hyper$Chromosome == "chr7",]
E_30 <- Astro_B_hyper[Astro_B_hyper$Chromosome == "chr8",]
E_31 <- Astro_B_hyper[Astro_B_hyper$Chromosome == "chr9",]
E_32 <- Astro_B_hyper[Astro_B_hyper$Chromosome == "chr10",]
E_33 <- Astro_B_hyper[Astro_B_hyper$Chromosome == "chr11",]
E_34 <- Astro_B_hyper[Astro_B_hyper$Chromosome == "chr12",]
E_35 <- Astro_B_hyper[Astro_B_hyper$Chromosome == "chr13",]
E_36 <- Astro_B_hyper[Astro_B_hyper$Chromosome == "chr14",]
E_37 <- Astro_B_hyper[Astro_B_hyper$Chromosome == "chr15",]
E_38 <- Astro_B_hyper[Astro_B_hyper$Chromosome == "chr16",]
E_39 <- Astro_B_hyper[Astro_B_hyper$Chromosome == "chr17",]
E_40 <- Astro_B_hyper[Astro_B_hyper$Chromosome == "chr18",]
E_41 <- Astro_B_hyper[Astro_B_hyper$Chromosome == "chr19",]
E_42 <- Astro_B_hyper[Astro_B_hyper$Chromosome == "chr20",]
E_43 <- Astro_B_hyper[Astro_B_hyper$Chromosome == "chr21",]
E_44 <- Astro_B_hyper[Astro_B_hyper$Chromosome == "chr22",]

E_45 <- Breast_B_hyper[Breast_B_hyper$Chromosome == "chr1",]
E_46 <- Breast_B_hyper[Breast_B_hyper$Chromosome == "chr2",]
E_47 <- Breast_B_hyper[Breast_B_hyper$Chromosome == "chr3",]
E_48 <- Breast_B_hyper[Breast_B_hyper$Chromosome == "chr4",]
E_49 <- Breast_B_hyper[Breast_B_hyper$Chromosome == "chr5",]
E_50 <- Breast_B_hyper[Breast_B_hyper$Chromosome == "chr6",]
E_51 <- Breast_B_hyper[Breast_B_hyper$Chromosome == "chr7",]
E_52 <- Breast_B_hyper[Breast_B_hyper$Chromosome == "chr8",]
E_53 <- Breast_B_hyper[Breast_B_hyper$Chromosome == "chr9",]
E_54 <- Breast_B_hyper[Breast_B_hyper$Chromosome == "chr10",]
E_55 <- Breast_B_hyper[Breast_B_hyper$Chromosome == "chr11",]
E_56 <- Breast_B_hyper[Breast_B_hyper$Chromosome == "chr12",]
E_57 <- Breast_B_hyper[Breast_B_hyper$Chromosome == "chr13",]
E_58 <- Breast_B_hyper[Breast_B_hyper$Chromosome == "chr14",]
E_59 <- Breast_B_hyper[Breast_B_hyper$Chromosome == "chr15",]
E_60 <- Breast_B_hyper[Breast_B_hyper$Chromosome == "chr16",]
E_61 <- Breast_B_hyper[Breast_B_hyper$Chromosome == "chr17",]
E_62 <- Breast_B_hyper[Breast_B_hyper$Chromosome == "chr18",]
E_63 <- Breast_B_hyper[Breast_B_hyper$Chromosome == "chr19",]
E_64 <- Breast_B_hyper[Breast_B_hyper$Chromosome == "chr20",]
E_65 <- Breast_B_hyper[Breast_B_hyper$Chromosome == "chr21",]
E_66 <- Breast_B_hyper[Breast_B_hyper$Chromosome == "chr22",]

E_67 <- Chol_B_hyper[Chol_B_hyper$Chromosome == "chr1",]
E_68 <- Chol_B_hyper[Chol_B_hyper$Chromosome == "chr2",]
E_69 <- Chol_B_hyper[Chol_B_hyper$Chromosome == "chr3",]
E_70 <- Chol_B_hyper[Chol_B_hyper$Chromosome == "chr4",]
E_71 <- Chol_B_hyper[Chol_B_hyper$Chromosome == "chr5",]
E_72 <- Chol_B_hyper[Chol_B_hyper$Chromosome == "chr6",]
E_73 <- Chol_B_hyper[Chol_B_hyper$Chromosome == "chr7",]
E_74 <- Chol_B_hyper[Chol_B_hyper$Chromosome == "chr8",]
E_75 <- Chol_B_hyper[Chol_B_hyper$Chromosome == "chr9",]
E_76 <- Chol_B_hyper[Chol_B_hyper$Chromosome == "chr10",]
E_77 <- Chol_B_hyper[Chol_B_hyper$Chromosome == "chr11",]
E_78 <- Chol_B_hyper[Chol_B_hyper$Chromosome == "chr12",]
E_79 <- Chol_B_hyper[Chol_B_hyper$Chromosome == "chr13",]
E_80 <- Chol_B_hyper[Chol_B_hyper$Chromosome == "chr14",]
E_81 <- Chol_B_hyper[Chol_B_hyper$Chromosome == "chr15",]
E_82 <- Chol_B_hyper[Chol_B_hyper$Chromosome == "chr16",]
E_83 <- Chol_B_hyper[Chol_B_hyper$Chromosome == "chr17",]
E_84 <- Chol_B_hyper[Chol_B_hyper$Chromosome == "chr18",]
E_85 <- Chol_B_hyper[Chol_B_hyper$Chromosome == "chr19",]
E_86 <- Chol_B_hyper[Chol_B_hyper$Chromosome == "chr20",]
E_87 <- Chol_B_hyper[Chol_B_hyper$Chromosome == "chr21",]
E_88 <- Chol_B_hyper[Chol_B_hyper$Chromosome == "chr22",]

E_89 <- Oligo_B_hyper[Oligo_B_hyper$Chromosome == "chr1",]
E_90 <- Oligo_B_hyper[Oligo_B_hyper$Chromosome == "chr2",]
E_91 <- Oligo_B_hyper[Oligo_B_hyper$Chromosome == "chr3",]
E_92 <- Oligo_B_hyper[Oligo_B_hyper$Chromosome == "chr4",]
E_93 <- Oligo_B_hyper[Oligo_B_hyper$Chromosome == "chr5",]
E_94 <- Oligo_B_hyper[Oligo_B_hyper$Chromosome == "chr6",]
E_95 <- Oligo_B_hyper[Oligo_B_hyper$Chromosome == "chr7",]
E_96 <- Oligo_B_hyper[Oligo_B_hyper$Chromosome == "chr8",]
E_97 <- Oligo_B_hyper[Oligo_B_hyper$Chromosome == "chr9",]
E_98 <- Oligo_B_hyper[Oligo_B_hyper$Chromosome == "chr10",]
E_99 <- Oligo_B_hyper[Oligo_B_hyper$Chromosome == "chr11",]
E_100 <- Oligo_B_hyper[Oligo_B_hyper$Chromosome == "chr12",]
E_101 <- Oligo_B_hyper[Oligo_B_hyper$Chromosome == "chr13",]
E_102 <- Oligo_B_hyper[Oligo_B_hyper$Chromosome == "chr14",]
E_103 <- Oligo_B_hyper[Oligo_B_hyper$Chromosome == "chr15",]
E_104 <- Oligo_B_hyper[Oligo_B_hyper$Chromosome == "chr16",]
E_105 <- Oligo_B_hyper[Oligo_B_hyper$Chromosome == "chr17",]
E_106 <- Oligo_B_hyper[Oligo_B_hyper$Chromosome == "chr18",]
E_107 <- Oligo_B_hyper[Oligo_B_hyper$Chromosome == "chr19",]
E_108 <- Oligo_B_hyper[Oligo_B_hyper$Chromosome == "chr20",]
E_109 <- Oligo_B_hyper[Oligo_B_hyper$Chromosome == "chr21",]
E_110 <- Oligo_B_hyper[Oligo_B_hyper$Chromosome == "chr22",]

E_111 <- SNUC_B_hyper[SNUC_B_hyper$Chromosome == "chr1",]
E_112 <- SNUC_B_hyper[SNUC_B_hyper$Chromosome == "chr2",]
E_113 <- SNUC_B_hyper[SNUC_B_hyper$Chromosome == "chr3",]
E_114 <- SNUC_B_hyper[SNUC_B_hyper$Chromosome == "chr4",]
E_115 <- SNUC_B_hyper[SNUC_B_hyper$Chromosome == "chr5",]
E_116 <- SNUC_B_hyper[SNUC_B_hyper$Chromosome == "chr6",]
E_117 <- SNUC_B_hyper[SNUC_B_hyper$Chromosome == "chr7",]
E_118 <- SNUC_B_hyper[SNUC_B_hyper$Chromosome == "chr8",]
E_119 <- SNUC_B_hyper[SNUC_B_hyper$Chromosome == "chr9",]
E_120 <- SNUC_B_hyper[SNUC_B_hyper$Chromosome == "chr10",]
E_121 <- SNUC_B_hyper[SNUC_B_hyper$Chromosome == "chr11",]
E_122 <- SNUC_B_hyper[SNUC_B_hyper$Chromosome == "chr12",]
E_123 <- SNUC_B_hyper[SNUC_B_hyper$Chromosome == "chr13",]
E_124 <- SNUC_B_hyper[SNUC_B_hyper$Chromosome == "chr14",]
E_125 <- SNUC_B_hyper[SNUC_B_hyper$Chromosome == "chr15",]
E_126 <- SNUC_B_hyper[SNUC_B_hyper$Chromosome == "chr16",]
E_127 <- SNUC_B_hyper[SNUC_B_hyper$Chromosome == "chr17",]
E_128 <- SNUC_B_hyper[SNUC_B_hyper$Chromosome == "chr18",]
E_129 <- SNUC_B_hyper[SNUC_B_hyper$Chromosome == "chr19",]
E_130 <- SNUC_B_hyper[SNUC_B_hyper$Chromosome == "chr20",]
E_131 <- SNUC_B_hyper[SNUC_B_hyper$Chromosome == "chr21",]
E_132 <- SNUC_B_hyper[SNUC_B_hyper$Chromosome == "chr22",]

E_1 <- common_test_E_hyper[common_test_E_hyper$Chromosome == "chr1",]
E_2 <- common_test_E_hyper[common_test_E_hyper$Chromosome == "chr2",]
E_3 <- common_test_E_hyper[common_test_E_hyper$Chromosome == "chr3",]
E_4 <- common_test_E_hyper[common_test_E_hyper$Chromosome == "chr4",]
E_5 <- common_test_E_hyper[common_test_E_hyper$Chromosome == "chr5",]
E_6 <- common_test_E_hyper[common_test_E_hyper$Chromosome == "chr6",]
E_7 <- common_test_E_hyper[common_test_E_hyper$Chromosome == "chr7",]
E_8 <- common_test_E_hyper[common_test_E_hyper$Chromosome == "chr8",]
E_9 <- common_test_E_hyper[common_test_E_hyper$Chromosome == "chr9",]
E_10 <- common_test_E_hyper[common_test_E_hyper$Chromosome == "chr10",]
E_11 <- common_test_E_hyper[common_test_E_hyper$Chromosome == "chr11",]
E_12 <- common_test_E_hyper[common_test_E_hyper$Chromosome == "chr12",]
E_13 <- common_test_E_hyper[common_test_E_hyper$Chromosome == "chr13",]
E_14 <- common_test_E_hyper[common_test_E_hyper$Chromosome == "chr14",]
E_15 <- common_test_E_hyper[common_test_E_hyper$Chromosome == "chr15",]
E_16 <- common_test_E_hyper[common_test_E_hyper$Chromosome == "chr16",]
E_17 <- common_test_E_hyper[common_test_E_hyper$Chromosome == "chr17",]
E_18 <- common_test_E_hyper[common_test_E_hyper$Chromosome == "chr18",]
E_19 <- common_test_E_hyper[common_test_E_hyper$Chromosome == "chr19",]
E_20 <- common_test_E_hyper[common_test_E_hyper$Chromosome == "chr20",]
E_21 <- common_test_E_hyper[common_test_E_hyper$Chromosome == "chr21",]
E_22 <- common_test_E_hyper[common_test_E_hyper$Chromosome == "chr22",]

E_23 <- common_test_E_hypo[common_test_E_hypo$Chromosome == "chr1",]
E_24 <- common_test_E_hypo[common_test_E_hypo$Chromosome == "chr2",]
E_25 <- common_test_E_hypo[common_test_E_hypo$Chromosome == "chr3",]
E_26 <- common_test_E_hypo[common_test_E_hypo$Chromosome == "chr4",]
E_27 <- common_test_E_hypo[common_test_E_hypo$Chromosome == "chr5",]
E_28 <- common_test_E_hypo[common_test_E_hypo$Chromosome == "chr6",]
E_29 <- common_test_E_hypo[common_test_E_hypo$Chromosome == "chr7",]
E_30 <- common_test_E_hypo[common_test_E_hypo$Chromosome == "chr8",]
E_31 <- common_test_E_hypo[common_test_E_hypo$Chromosome == "chr9",]
E_32 <- common_test_E_hypo[common_test_E_hypo$Chromosome == "chr10",]
E_33 <- common_test_E_hypo[common_test_E_hypo$Chromosome == "chr11",]
E_34 <- common_test_E_hypo[common_test_E_hypo$Chromosome == "chr12",]
E_35 <- common_test_E_hypo[common_test_E_hypo$Chromosome == "chr13",]
E_36 <- common_test_E_hypo[common_test_E_hypo$Chromosome == "chr14",]
E_37 <- common_test_E_hypo[common_test_E_hypo$Chromosome == "chr15",]
E_38 <- common_test_E_hypo[common_test_E_hypo$Chromosome == "chr16",]
E_39 <- common_test_E_hypo[common_test_E_hypo$Chromosome == "chr17",]
E_40 <- common_test_E_hypo[common_test_E_hypo$Chromosome == "chr18",]
E_41 <- common_test_E_hypo[common_test_E_hypo$Chromosome == "chr19",]
E_42 <- common_test_E_hypo[common_test_E_hypo$Chromosome == "chr20",]
E_43 <- common_test_E_hypo[common_test_E_hypo$Chromosome == "chr21",]
E_44 <- common_test_E_hypo[common_test_E_hypo$Chromosome == "chr22",]

################################################CONTROL###############################################3
E_1 <- Blood_WG_hyper[Blood_WG_hyper$Chromosome == "chr1",]
E_2 <- Blood_WG_hyper[Blood_WG_hyper$Chromosome == "chr2",]
E_3 <- Blood_WG_hyper[Blood_WG_hyper$Chromosome == "chr3",]
E_4 <- Blood_WG_hyper[Blood_WG_hyper$Chromosome == "chr4",]
E_5 <- Blood_WG_hyper[Blood_WG_hyper$Chromosome == "chr5",]
E_6 <- Blood_WG_hyper[Blood_WG_hyper$Chromosome == "chr6",]
E_7 <- Blood_WG_hyper[Blood_WG_hyper$Chromosome == "chr7",]
E_8 <- Blood_WG_hyper[Blood_WG_hyper$Chromosome == "chr8",]
E_9 <- Blood_WG_hyper[Blood_WG_hyper$Chromosome == "chr9",]
E_10 <- Blood_WG_hyper[Blood_WG_hyper$Chromosome == "chr10",]
E_11 <- Blood_WG_hyper[Blood_WG_hyper$Chromosome == "chr11",]
E_12 <- Blood_WG_hyper[Blood_WG_hyper$Chromosome == "chr12",]
E_13 <- Blood_WG_hyper[Blood_WG_hyper$Chromosome == "chr13",]
E_14 <- Blood_WG_hyper[Blood_WG_hyper$Chromosome == "chr14",]
E_15 <- Blood_WG_hyper[Blood_WG_hyper$Chromosome == "chr15",]
E_16 <- Blood_WG_hyper[Blood_WG_hyper$Chromosome == "chr16",]
E_17 <- Blood_WG_hyper[Blood_WG_hyper$Chromosome == "chr17",]
E_18 <- Blood_WG_hyper[Blood_WG_hyper$Chromosome == "chr18",]
E_19 <- Blood_WG_hyper[Blood_WG_hyper$Chromosome == "chr19",]
E_20 <- Blood_WG_hyper[Blood_WG_hyper$Chromosome == "chr20",]
E_21 <- Blood_WG_hyper[Blood_WG_hyper$Chromosome == "chr21",]
E_22 <- Blood_WG_hyper[Blood_WG_hyper$Chromosome == "chr22",]

E_23 <- Neuro_WG_hyper[Neuro_WG_hyper$Chromosome == "chr1",]
E_24 <- Neuro_WG_hyper[Neuro_WG_hyper$Chromosome == "chr2",]
E_25 <- Neuro_WG_hyper[Neuro_WG_hyper$Chromosome == "chr3",]
E_26 <- Neuro_WG_hyper[Neuro_WG_hyper$Chromosome == "chr4",]
E_27 <- Neuro_WG_hyper[Neuro_WG_hyper$Chromosome == "chr5",]
E_28 <- Neuro_WG_hyper[Neuro_WG_hyper$Chromosome == "chr6",]
E_29 <- Neuro_WG_hyper[Neuro_WG_hyper$Chromosome == "chr7",]
E_30 <- Neuro_WG_hyper[Neuro_WG_hyper$Chromosome == "chr8",]
E_31 <- Neuro_WG_hyper[Neuro_WG_hyper$Chromosome == "chr9",]
E_32 <- Neuro_WG_hyper[Neuro_WG_hyper$Chromosome == "chr10",]
E_33 <- Neuro_WG_hyper[Neuro_WG_hyper$Chromosome == "chr11",]
E_34 <- Neuro_WG_hyper[Neuro_WG_hyper$Chromosome == "chr12",]
E_35 <- Neuro_WG_hyper[Neuro_WG_hyper$Chromosome == "chr13",]
E_36 <- Neuro_WG_hyper[Neuro_WG_hyper$Chromosome == "chr14",]
E_37 <- Neuro_WG_hyper[Neuro_WG_hyper$Chromosome == "chr15",]
E_38 <- Neuro_WG_hyper[Neuro_WG_hyper$Chromosome == "chr16",]
E_39 <- Neuro_WG_hyper[Neuro_WG_hyper$Chromosome == "chr17",]
E_40 <- Neuro_WG_hyper[Neuro_WG_hyper$Chromosome == "chr18",]
E_41 <- Neuro_WG_hyper[Neuro_WG_hyper$Chromosome == "chr19",]
E_42 <- Neuro_WG_hyper[Neuro_WG_hyper$Chromosome == "chr20",]
E_43 <- Neuro_WG_hyper[Neuro_WG_hyper$Chromosome == "chr21",]
E_44 <- Neuro_WG_hyper[Neuro_WG_hyper$Chromosome == "chr22",]

E_45 <- GBM_WG_hyper[GBM_WG_hyper$Chromosome == "chr1",]
E_46 <- GBM_WG_hyper[GBM_WG_hyper$Chromosome == "chr2",]
E_47 <- GBM_WG_hyper[GBM_WG_hyper$Chromosome == "chr3",]
E_48 <- GBM_WG_hyper[GBM_WG_hyper$Chromosome == "chr4",]
E_49 <- GBM_WG_hyper[GBM_WG_hyper$Chromosome == "chr5",]
E_50 <- GBM_WG_hyper[GBM_WG_hyper$Chromosome == "chr6",]
E_51 <- GBM_WG_hyper[GBM_WG_hyper$Chromosome == "chr7",]
E_52 <- GBM_WG_hyper[GBM_WG_hyper$Chromosome == "chr8",]
E_53 <- GBM_WG_hyper[GBM_WG_hyper$Chromosome == "chr9",]
E_54 <- GBM_WG_hyper[GBM_WG_hyper$Chromosome == "chr10",]
E_55 <- GBM_WG_hyper[GBM_WG_hyper$Chromosome == "chr11",]
E_56 <- GBM_WG_hyper[GBM_WG_hyper$Chromosome == "chr12",]
E_57 <- GBM_WG_hyper[GBM_WG_hyper$Chromosome == "chr13",]
E_58 <- GBM_WG_hyper[GBM_WG_hyper$Chromosome == "chr14",]
E_59 <- GBM_WG_hyper[GBM_WG_hyper$Chromosome == "chr15",]
E_60 <- GBM_WG_hyper[GBM_WG_hyper$Chromosome == "chr16",]
E_61 <- GBM_WG_hyper[GBM_WG_hyper$Chromosome == "chr17",]
E_62 <- GBM_WG_hyper[GBM_WG_hyper$Chromosome == "chr18",]
E_63 <- GBM_WG_hyper[GBM_WG_hyper$Chromosome == "chr19",]
E_64 <- GBM_WG_hyper[GBM_WG_hyper$Chromosome == "chr20",]
E_65 <- GBM_WG_hyper[GBM_WG_hyper$Chromosome == "chr21",]
E_66 <- GBM_WG_hyper[GBM_WG_hyper$Chromosome == "chr22",]

E_67 <- PAstro_WG_hyper[PAstro_WG_hyper$Chromosome == "chr1",]
E_68 <- PAstro_WG_hyper[PAstro_WG_hyper$Chromosome == "chr2",]
E_69 <- PAstro_WG_hyper[PAstro_WG_hyper$Chromosome == "chr3",]
E_70 <- PAstro_WG_hyper[PAstro_WG_hyper$Chromosome == "chr4",]
E_71 <- PAstro_WG_hyper[PAstro_WG_hyper$Chromosome == "chr5",]
E_72 <- PAstro_WG_hyper[PAstro_WG_hyper$Chromosome == "chr6",]
E_73 <- PAstro_WG_hyper[PAstro_WG_hyper$Chromosome == "chr7",]
E_74 <- PAstro_WG_hyper[PAstro_WG_hyper$Chromosome == "chr8",]
E_75 <- PAstro_WG_hyper[PAstro_WG_hyper$Chromosome == "chr9",]
E_76 <- PAstro_WG_hyper[PAstro_WG_hyper$Chromosome == "chr10",]
E_77 <- PAstro_WG_hyper[PAstro_WG_hyper$Chromosome == "chr11",]
E_78 <- PAstro_WG_hyper[PAstro_WG_hyper$Chromosome == "chr12",]
E_79 <- PAstro_WG_hyper[PAstro_WG_hyper$Chromosome == "chr13",]
E_80 <- PAstro_WG_hyper[PAstro_WG_hyper$Chromosome == "chr14",]
E_81 <- PAstro_WG_hyper[PAstro_WG_hyper$Chromosome == "chr15",]
E_82 <- PAstro_WG_hyper[PAstro_WG_hyper$Chromosome == "chr16",]
E_83 <- PAstro_WG_hyper[PAstro_WG_hyper$Chromosome == "chr17",]
E_84 <- PAstro_WG_hyper[PAstro_WG_hyper$Chromosome == "chr18",]
E_85 <- PAstro_WG_hyper[PAstro_WG_hyper$Chromosome == "chr19",]
E_86 <- PAstro_WG_hyper[PAstro_WG_hyper$Chromosome == "chr20",]
E_87 <- PAstro_WG_hyper[PAstro_WG_hyper$Chromosome == "chr21",]
E_88 <- PAstro_WG_hyper[PAstro_WG_hyper$Chromosome == "chr22",]

E_89 <- SUDEP_WG_hyper[SUDEP_WG_hyper$Chromosome == "chr1",]
E_90 <- SUDEP_WG_hyper[SUDEP_WG_hyper$Chromosome == "chr2",]
E_91 <- SUDEP_WG_hyper[SUDEP_WG_hyper$Chromosome == "chr3",]
E_92 <- SUDEP_WG_hyper[SUDEP_WG_hyper$Chromosome == "chr4",]
E_93 <- SUDEP_WG_hyper[SUDEP_WG_hyper$Chromosome == "chr5",]
E_94 <- SUDEP_WG_hyper[SUDEP_WG_hyper$Chromosome == "chr6",]
E_95 <- SUDEP_WG_hyper[SUDEP_WG_hyper$Chromosome == "chr7",]
E_96 <- SUDEP_WG_hyper[SUDEP_WG_hyper$Chromosome == "chr8",]
E_97 <- SUDEP_WG_hyper[SUDEP_WG_hyper$Chromosome == "chr9",]
E_98 <- SUDEP_WG_hyper[SUDEP_WG_hyper$Chromosome == "chr10",]
E_99 <- SUDEP_WG_hyper[SUDEP_WG_hyper$Chromosome == "chr11",]
E_100 <- SUDEP_WG_hyper[SUDEP_WG_hyper$Chromosome == "chr12",]
E_101 <- SUDEP_WG_hyper[SUDEP_WG_hyper$Chromosome == "chr13",]
E_102 <- SUDEP_WG_hyper[SUDEP_WG_hyper$Chromosome == "chr14",]
E_103 <- SUDEP_WG_hyper[SUDEP_WG_hyper$Chromosome == "chr15",]
E_104 <- SUDEP_WG_hyper[SUDEP_WG_hyper$Chromosome == "chr16",]
E_105 <- SUDEP_WG_hyper[SUDEP_WG_hyper$Chromosome == "chr17",]
E_106 <- SUDEP_WG_hyper[SUDEP_WG_hyper$Chromosome == "chr18",]
E_107 <- SUDEP_WG_hyper[SUDEP_WG_hyper$Chromosome == "chr19",]
E_108 <- SUDEP_WG_hyper[SUDEP_WG_hyper$Chromosome == "chr20",]
E_109 <- SUDEP_WG_hyper[SUDEP_WG_hyper$Chromosome == "chr21",]
E_110 <- SUDEP_WG_hyper[SUDEP_WG_hyper$Chromosome == "chr22",]

E_1 <- common_control_P_hyper[common_control_P_hyper$Chromosome == "chr1",]
E_2 <- common_control_P_hyper[common_control_P_hyper$Chromosome == "chr2",]
E_3 <- common_control_P_hyper[common_control_P_hyper$Chromosome == "chr3",]
E_4 <- common_control_P_hyper[common_control_P_hyper$Chromosome == "chr4",]
E_5 <- common_control_P_hyper[common_control_P_hyper$Chromosome == "chr5",]
E_6 <- common_control_P_hyper[common_control_P_hyper$Chromosome == "chr6",]
E_7 <- common_control_P_hyper[common_control_P_hyper$Chromosome == "chr7",]
E_8 <- common_control_P_hyper[common_control_P_hyper$Chromosome == "chr8",]
E_9 <- common_control_P_hyper[common_control_P_hyper$Chromosome == "chr9",]
E_10 <- common_control_P_hyper[common_control_P_hyper$Chromosome == "chr10",]
E_11 <- common_control_P_hyper[common_control_P_hyper$Chromosome == "chr11",]
E_12 <- common_control_P_hyper[common_control_P_hyper$Chromosome == "chr12",]
E_13 <- common_control_P_hyper[common_control_P_hyper$Chromosome == "chr13",]
E_14 <- common_control_P_hyper[common_control_P_hyper$Chromosome == "chr14",]
E_15 <- common_control_P_hyper[common_control_P_hyper$Chromosome == "chr15",]
E_16 <- common_control_P_hyper[common_control_P_hyper$Chromosome == "chr16",]
E_17 <- common_control_P_hyper[common_control_P_hyper$Chromosome == "chr17",]
E_18 <- common_control_P_hyper[common_control_P_hyper$Chromosome == "chr18",]
E_19 <- common_control_P_hyper[common_control_P_hyper$Chromosome == "chr19",]
E_20 <- common_control_P_hyper[common_control_P_hyper$Chromosome == "chr20",]
E_21 <- common_control_P_hyper[common_control_P_hyper$Chromosome == "chr21",]
E_22 <- common_control_P_hyper[common_control_P_hyper$Chromosome == "chr22",]

E_23 <- common_control_P_hypo[common_control_P_hypo$Chromosome == "chr1",]
E_24 <- common_control_P_hypo[common_control_P_hypo$Chromosome == "chr2",]
E_25 <- common_control_P_hypo[common_control_P_hypo$Chromosome == "chr3",]
E_26 <- common_control_P_hypo[common_control_P_hypo$Chromosome == "chr4",]
E_27 <- common_control_P_hypo[common_control_P_hypo$Chromosome == "chr5",]
E_28 <- common_control_P_hypo[common_control_P_hypo$Chromosome == "chr6",]
E_29 <- common_control_P_hypo[common_control_P_hypo$Chromosome == "chr7",]
E_30 <- common_control_P_hypo[common_control_P_hypo$Chromosome == "chr8",]
E_31 <- common_control_P_hypo[common_control_P_hypo$Chromosome == "chr9",]
E_32 <- common_control_P_hypo[common_control_P_hypo$Chromosome == "chr10",]
E_33 <- common_control_P_hypo[common_control_P_hypo$Chromosome == "chr11",]
E_34 <- common_control_P_hypo[common_control_P_hypo$Chromosome == "chr12",]
E_35 <- common_control_P_hypo[common_control_P_hypo$Chromosome == "chr13",]
E_36 <- common_control_P_hypo[common_control_P_hypo$Chromosome == "chr14",]
E_37 <- common_control_P_hypo[common_control_P_hypo$Chromosome == "chr15",]
E_38 <- common_control_P_hypo[common_control_P_hypo$Chromosome == "chr16",]
E_39 <- common_control_P_hypo[common_control_P_hypo$Chromosome == "chr17",]
E_40 <- common_control_P_hypo[common_control_P_hypo$Chromosome == "chr18",]
E_41 <- common_control_P_hypo[common_control_P_hypo$Chromosome == "chr19",]
E_42 <- common_control_P_hypo[common_control_P_hypo$Chromosome == "chr20",]
E_43 <- common_control_P_hypo[common_control_P_hypo$Chromosome == "chr21",]
E_44 <- common_control_P_hypo[common_control_P_hypo$Chromosome == "chr22",]



#find the means of the mean meth values in each dataframe
chr_means = NULL
chr_means$means <- rbind(
  mean(E_1$Mean_Methylation),
  mean(E_2$Mean_Methylation),
  mean(E_3$Mean_Methylation),
  mean(E_4$Mean_Methylation),
  mean(E_5$Mean_Methylation),
  mean(E_6$Mean_Methylation),
  mean(E_7$Mean_Methylation),
  mean(E_8$Mean_Methylation),
  mean(E_9$Mean_Methylation),
  mean(E_10$Mean_Methylation),
  mean(E_11$Mean_Methylation),
  mean(E_12$Mean_Methylation),
  mean(E_13$Mean_Methylation),
  mean(E_14$Mean_Methylation),
  mean(E_15$Mean_Methylation),
  mean(E_16$Mean_Methylation),
  mean(E_17$Mean_Methylation),
  mean(E_18$Mean_Methylation),
  mean(E_19$Mean_Methylation),
  mean(E_20$Mean_Methylation),
  mean(E_21$Mean_Methylation),
  mean(E_22$Mean_Methylation),
  mean(E_23$Mean_Methylation),
  mean(E_24$Mean_Methylation),
  mean(E_25$Mean_Methylation),
  mean(E_26$Mean_Methylation),
  mean(E_27$Mean_Methylation),
  mean(E_28$Mean_Methylation),
  mean(E_29$Mean_Methylation),
  mean(E_30$Mean_Methylation),
  mean(E_31$Mean_Methylation),
  mean(E_32$Mean_Methylation),
  mean(E_33$Mean_Methylation),
  mean(E_34$Mean_Methylation),
  mean(E_35$Mean_Methylation),
  mean(E_36$Mean_Methylation),
  mean(E_37$Mean_Methylation),
  mean(E_38$Mean_Methylation),
  mean(E_39$Mean_Methylation),
  mean(E_40$Mean_Methylation),
  mean(E_41$Mean_Methylation),
  mean(E_42$Mean_Methylation),
  mean(E_43$Mean_Methylation),
  mean(E_44$Mean_Methylation),
  mean(E_45$Mean_Methylation),
  mean(E_46$Mean_Methylation),
  mean(E_47$Mean_Methylation),
  mean(E_48$Mean_Methylation),
  mean(E_49$Mean_Methylation),
  mean(E_50$Mean_Methylation),
  mean(E_51$Mean_Methylation),
  mean(E_52$Mean_Methylation),
  mean(E_53$Mean_Methylation),
  mean(E_54$Mean_Methylation),
  mean(E_55$Mean_Methylation),
  mean(E_56$Mean_Methylation),
  mean(E_57$Mean_Methylation),
  mean(E_58$Mean_Methylation),
  mean(E_59$Mean_Methylation),
  mean(E_60$Mean_Methylation),
  mean(E_61$Mean_Methylation),
  mean(E_62$Mean_Methylation),
  mean(E_63$Mean_Methylation),
  mean(E_64$Mean_Methylation),
  mean(E_65$Mean_Methylation),
  mean(E_66$Mean_Methylation),
  mean(E_67$Mean_Methylation),
  mean(E_68$Mean_Methylation),
  mean(E_69$Mean_Methylation),
  mean(E_70$Mean_Methylation),
  mean(E_71$Mean_Methylation),
  mean(E_72$Mean_Methylation),
  mean(E_73$Mean_Methylation),  
  mean(E_74$Mean_Methylation),
  mean(E_75$Mean_Methylation),  
  mean(E_76$Mean_Methylation),
  mean(E_77$Mean_Methylation),
  mean(E_78$Mean_Methylation),
  mean(E_79$Mean_Methylation),
  mean(E_80$Mean_Methylation),
  mean(E_81$Mean_Methylation),
  mean(E_82$Mean_Methylation),
  mean(E_83$Mean_Methylation),
  mean(E_84$Mean_Methylation),
  mean(E_85$Mean_Methylation),
  mean(E_86$Mean_Methylation),
  mean(E_87$Mean_Methylation),
  mean(E_88$Mean_Methylation),
  mean(E_89$Mean_Methylation),
  mean(E_90$Mean_Methylation),
  mean(E_91$Mean_Methylation),
  mean(E_92$Mean_Methylation),
  mean(E_93$Mean_Methylation),
  mean(E_94$Mean_Methylation),
  mean(E_95$Mean_Methylation),
  mean(E_96$Mean_Methylation),
  mean(E_97$Mean_Methylation),
  mean(E_98$Mean_Methylation),
  mean(E_99$Mean_Methylation),
  mean(E_100$Mean_Methylation),
  mean(E_101$Mean_Methylation),
  mean(E_102$Mean_Methylation),
  mean(E_103$Mean_Methylation),
  mean(E_104$Mean_Methylation),
  mean(E_105$Mean_Methylation),
  mean(E_106$Mean_Methylation),
  mean(E_107$Mean_Methylation),
  mean(E_108$Mean_Methylation),
  mean(E_109$Mean_Methylation),
  mean(E_110$Mean_Methylation),
  mean(E_111$Mean_Methylation),
  mean(E_112$Mean_Methylation),
  mean(E_113$Mean_Methylation),
  mean(E_114$Mean_Methylation),
  mean(E_115$Mean_Methylation),
  mean(E_116$Mean_Methylation),
  mean(E_117$Mean_Methylation),
  mean(E_118$Mean_Methylation),
  mean(E_119$Mean_Methylation),
  mean(E_120$Mean_Methylation),
  mean(E_121$Mean_Methylation),
  mean(E_122$Mean_Methylation),
  mean(E_123$Mean_Methylation),
  mean(E_124$Mean_Methylation),
  mean(E_125$Mean_Methylation),
  mean(E_126$Mean_Methylation),
  mean(E_127$Mean_Methylation),
  mean(E_128$Mean_Methylation),
  mean(E_129$Mean_Methylation),
  mean(E_130$Mean_Methylation),
  mean(E_131$Mean_Methylation),
  mean(E_132$Mean_Methylation)
)
chr_means
#copy and paste above values into excel
#to remove n characters from end of value in a cell,  formula== RIGHT(G1, 10)

##Read in dataframes
WG_test_hypo_final_df <-read.csv("C:\\Users\\Ramona Bledea\\Documents\\PROJECT RESULTS\\IDH Tumor Project\\Results\\Test Samples - Differential Methylation\\Using M Values\\Hypomethylated Probes\\Test_ChromvsMethPlot_Hypo_Datasheet_WHOLEGENOME.csv")
B_test_hypo_final_df <-read.csv("C:\\Users\\Ramona Bledea\\Documents\\PROJECT RESULTS\\IDH Tumor Project\\Results\\Test Samples - Differential Methylation\\Using M Values\\Hypomethylated Probes\\Test_ChromvsMethPlot_Hypo_Datasheet_BODIES.csv")
P_test_hypo_final_df <-read.csv("C:\\Users\\Ramona Bledea\\Documents\\PROJECT RESULTS\\IDH Tumor Project\\Results\\Test Samples - Differential Methylation\\Using M Values\\Hypomethylated Probes\\Test_ChromvsMethPlot_Hypo_Datasheet_PROMOTERS.csv")
E_test_hypo_final_df <-read.csv("C:\\Users\\Ramona Bledea\\Documents\\PROJECT RESULTS\\IDH Tumor Project\\Results\\Test Samples - Differential Methylation\\Using M Values\\Hypomethylated Probes\\Test_ChromvsMethPlot_Hypo_Datasheet_ENHANCERS.csv")
A_test_hypo_final_df <-read.csv("C:\\Users\\Ramona Bledea\\Documents\\PROJECT RESULTS\\IDH Tumor Project\\Results\\Test Samples - Differential Methylation\\Using M Values\\Hypomethylated Probes\\Test_ChromvsMethPlot_Hypo_Datasheet_ALL.csv")
WG_control_hypo_final_df <-read.csv("C:\\Users\\Ramona Bledea\\Documents\\PROJECT RESULTS\\IDH Tumor Project\\Results\\Control Samples - Differential Methylation\\Using M Values\\Hypomethylated Probes\\control_ChromvsMethPlot_Hypo_Datasheet_WHOLEGENOME.csv")
B_control_hypo_final_df <-read.csv("C:\\Users\\Ramona Bledea\\Documents\\PROJECT RESULTS\\IDH Tumor Project\\Results\\Control Samples - Differential Methylation\\Using M Values\\Hypomethylated Probes\\control_ChromvsMethPlot_Hypo_Datasheet_BODIES.csv")
P_control_hypo_final_df <-read.csv("C:\\Users\\Ramona Bledea\\Documents\\PROJECT RESULTS\\IDH Tumor Project\\Results\\Control Samples - Differential Methylation\\Using M Values\\Hypomethylated Probes\\control_ChromvsMethPlot_Hypo_Datasheet_PROMOTERS.csv")
E_control_hypo_final_df <-read.csv("C:\\Users\\Ramona Bledea\\Documents\\PROJECT RESULTS\\IDH Tumor Project\\Results\\Control Samples - Differential Methylation\\Using M Values\\Hypomethylated Probes\\control_ChromvsMethPlot_Hypo_Datasheet_ENHANCERS.csv")
A_control_hypo_final_df <-read.csv("C:\\Users\\Ramona Bledea\\Documents\\PROJECT RESULTS\\IDH Tumor Project\\Results\\Control Samples - Differential Methylation\\Using M Values\\Hypomethylated Probes\\control_ChromvsMethPlot_Hypo_Datasheet_ALL.csv")

WG_test_hyper_final_df <-read.csv("C:\\Users\\Ramona Bledea\\Documents\\PROJECT RESULTS\\IDH Tumor Project\\Results\\Test Samples - Differential Methylation\\Using M Values\\hypermethylated Probes\\Test_ChromvsMethPlot_hyper_Datasheet_WHOLEGENOME.csv")
B_test_hyper_final_df <-read.csv("C:\\Users\\Ramona Bledea\\Documents\\PROJECT RESULTS\\IDH Tumor Project\\Results\\Test Samples - Differential Methylation\\Using M Values\\hypermethylated Probes\\Test_ChromvsMethPlot_hyper_Datasheet_BODIES.csv")
P_test_hyper_final_df <-read.csv("C:\\Users\\Ramona Bledea\\Documents\\PROJECT RESULTS\\IDH Tumor Project\\Results\\Test Samples - Differential Methylation\\Using M Values\\hypermethylated Probes\\Test_ChromvsMethPlot_hyper_Datasheet_PROMOTERS.csv")
E_test_hyper_final_df <-read.csv("C:\\Users\\Ramona Bledea\\Documents\\PROJECT RESULTS\\IDH Tumor Project\\Results\\Test Samples - Differential Methylation\\Using M Values\\hypermethylated Probes\\Test_ChromvsMethPlot_hyper_Datasheet_ENHANCERS.csv")
A_test_hyper_final_df <-read.csv("C:\\Users\\Ramona Bledea\\Documents\\PROJECT RESULTS\\IDH Tumor Project\\Results\\Test Samples - Differential Methylation\\Using M Values\\Hypermethylated Probes\\Test_ChromvsMethPlot_Unique_Hyper_Datasheet_ALL.csv")
WG_control_hyper_final_df <-read.csv("C:\\Users\\Ramona Bledea\\Documents\\PROJECT RESULTS\\IDH Tumor Project\\Results\\Control Samples - Differential Methylation\\Using M Values\\hypermethylated Probes\\control_ChromvsMethPlot_hyper_Datasheet_WHOLEGENOME.csv")
B_control_hyper_final_df <-read.csv("C:\\Users\\Ramona Bledea\\Documents\\PROJECT RESULTS\\IDH Tumor Project\\Results\\Control Samples - Differential Methylation\\Using M Values\\hypermethylated Probes\\control_ChromvsMethPlot_hyper_Datasheet_BODIES.csv")
P_control_hyper_final_df <-read.csv("C:\\Users\\Ramona Bledea\\Documents\\PROJECT RESULTS\\IDH Tumor Project\\Results\\Control Samples - Differential Methylation\\Using M Values\\hypermethylated Probes\\control_ChromvsMethPlot_hyper_Datasheet_PROMOTERS.csv")
E_control_hyper_final_df <-read.csv("C:\\Users\\Ramona Bledea\\Documents\\PROJECT RESULTS\\IDH Tumor Project\\Results\\Control Samples - Differential Methylation\\Using M Values\\hypermethylated Probes\\control_ChromvsMethPlot_hyper_Datasheet_ENHANCERS.csv")
A_control_hyper_final_df <-read.csv("C:\\Users\\Ramona Bledea\\Documents\\PROJECT RESULTS\\IDH Tumor Project\\Results\\Control Samples - Differential Methylation\\Using M Values\\hypermethylated Probes\\control_ChromvsMethPlot_hyper_Datasheet_ALL.csv")

hypo_common_df <- read.csv("C:\\Users\\Ramona Bledea\\Documents\\PROJECT RESULTS\\IDH Tumor Project\\Results\\ChromvsMethPlot_CommonProbes_Datasheet_ALL_hypo.csv")
hyper_common_df <- read.csv("C:\\Users\\Ramona Bledea\\Documents\\PROJECT RESULTS\\IDH Tumor Project\\Results\\ChromvsMethPlot_CommonProbes_Datasheet_ALL_hyper.csv")
common_df <- read.csv("C:\\Users\\Ramona Bledea\\Documents\\PROJECT RESULTS\\IDH Tumor Project\\Results\\ChromvsMethPlot_CommonProbes_Datasheet_ALL.csv")

##Generate the plots
library(ggplot2)

#UNIQUE TEST SAMPLES
setwd("C:\\Users\\Ramona Bledea\\Documents\\PROJECT RESULTS\\IDH Tumor Project\\Results\\Test Samples - Differential Methylation\\Using M Values\\hypermethylated Probes")
png(file=paste("Test_All_hyper_ChromvsMeanMeth2.png", sep=""),width=3000,height=3000,pointsize=60)
a<-ggplot(data=A_test_hyper_final_df, aes(x=Chromosome, y=Mean_Methylation, fill=Tumor_Type)) +
  geom_bar(stat="identity", position=position_dodge(), width=0.7) +
  ggtitle("IDH Mutant Samples Unique Hypermethylated Probes") +
  xlab("Chromosome") + ylab("Mean Methylation (M)") +
  theme(text=element_text(size=35, family="Arial")) +
  coord_cartesian(ylim=c(2,2.6))
a + xlim("1", "2", "3","4", "5", "6","7", "8", "9","10", "11", "12","13", "14", "15","16", "17", "18","19", "20", "21","22")
a + facet_wrap(~Gene_Region, ncol = 2)
dev.off()
setwd("C:\\Users\\Ramona Bledea\\Documents\\PROJECT RESULTS\\IDH Tumor Project\\Results\\Test Samples - Differential Methylation\\Using M Values\\hypomethylated Probes")
png(file=paste("Test_All_hypo_ChromvsMeanMeth2.png", sep=""),width=3000,height=3000,pointsize=60)
a<-ggplot(data=A_test_hypo_final_df, aes(x=Chromosome, y=Mean_Methylation, fill=Tumor_Type)) +
  geom_bar(stat="identity", position=position_dodge(), width=0.7) +
  ggtitle("IDH Mutant Samples Unique Hypomethylated Probes") +
  xlab("Chromosome") + ylab("Mean Methylation (M)") +
  theme(text=element_text(size=35, family="Arial")) +
  coord_cartesian(ylim=c(-3.6,-2))
a + xlim("1", "2", "3","4", "5", "6","7", "8", "9","10", "11", "12","13", "14", "15","16", "17", "18","19", "20", "21","22")
a + facet_wrap(~Gene_Region, ncol = 2)
dev.off()

#UNIQUE CONTROL SAMPLES
setwd("C:\\Users\\Ramona Bledea\\Documents\\PROJECT RESULTS\\IDH Tumor Project\\Results\\Control Samples - Differential Methylation\\Using M Values\\hypermethylated Probes")
png(file=paste("Control_All_hyper_ChromvsMeanMeth2.png", sep=""),width=3000,height=3000,pointsize=60)
a<-ggplot(data=A_control_hyper_final_df, aes(x=Chromosome, y=Mean_Methylation, fill=Tumor_Type)) +
  geom_bar(stat="identity", position=position_dodge(), width=0.7) +
  ggtitle("IDH WT Samples Unique Hypermethylated Probes") +
  xlab("Chromosome") + ylab("Mean Methylation (M)") +
  theme(text=element_text(size=35, family="Arial")) +
  coord_cartesian(ylim=c(2,2.6))
a + xlim("1", "2", "3","4", "5", "6","7", "8", "9","10", "11", "12","13", "14", "15","16", "17", "18","19", "20", "21","22")
a + facet_wrap(~Gene_Region, ncol = 2)
dev.off()
setwd("C:\\Users\\Ramona Bledea\\Documents\\PROJECT RESULTS\\IDH Tumor Project\\Results\\Control Samples - Differential Methylation\\Using M Values\\hypomethylated Probes")
png(file=paste("Control_All_hypo_ChromvsMeanMeth2.png", sep=""),width=3000,height=3000,pointsize=60)
a<-ggplot(data=A_control_hypo_final_df, aes(x=Chromosome, y=Mean_Methylation, fill=Tumor_Type)) +
  geom_bar(stat="identity", position=position_dodge(), width=0.7) +
  ggtitle("IDH WT Samples Unique Hypomethylated Probes") +
  xlab("Chromosome") + ylab("Mean Methylation (M)") +
  theme(text=element_text(size=35, family="Arial")) +
  coord_cartesian(ylim=c(-3.6,-2))
a + xlim("1", "2", "3","4", "5", "6","7", "8", "9","10", "11", "12","13", "14", "15","16", "17", "18","19", "20", "21","22")
a + facet_wrap(~Gene_Region, ncol = 2)
dev.off()


####COMMON TEST VS CONTROL####
setwd("C:\\Users\\Ramona Bledea\\Documents\\PROJECT RESULTS\\IDH Tumor Project\\Results")
png(file=paste("TestvsControl_All__hyper_Common_ChromvsMeanMeth-yaxiszoom.png", sep=""),width=4000,height=2000,pointsize=30)
a<-ggplot(data=hyper_common_df, aes(x=Chromosome, y=Mean_Methylation, fill=Gene_Region)) +
  geom_bar(stat="identity", position=position_dodge(), width=0.5) +
  ggtitle("All Common Probes Mean Methylation vs. Chromosome") +
  xlab("Chromosome") + ylab("Mean Methylation (M)") +
  theme(text=element_text(size=35, family="Arial")) +
  coord_cartesian(ylim=c(2.5,3))
a + xlim("1", "2", "3","4", "5", "6","7", "8", "9","10", "11", "12","13", "14", "15","16", "17", "18","19", "20", "21","22")
a + facet_wrap(~Tumor_Type, ncol = 1)
dev.off()
setwd("C:\\Users\\Ramona Bledea\\Documents\\PROJECT RESULTS\\IDH Tumor Project\\Results")
png(file=paste("TestvsControl_All__hypo_Common_ChromvsMeanMeth-yaxiszoom.png", sep=""),width=4000,height=2000,pointsize=30)
a<-ggplot(data=hypo_common_df, aes(x=Chromosome, y=Mean_Methylation, fill=Gene_Region)) +
  geom_bar(stat="identity", position=position_dodge(), width=0.5) +
  ggtitle("All Common Probes Mean Methylation vs. Chromosome") +
  xlab("Chromosome") + ylab("Mean Methylation (M)") +
  theme(text=element_text(size=35, family="Arial")) +
  coord_cartesian(ylim=c(-3.75,-3))
a + xlim("1", "2", "3","4", "5", "6","7", "8", "9","10", "11", "12","13", "14", "15","16", "17", "18","19", "20", "21","22")
a + facet_wrap(~Tumor_Type, ncol = 1)
dev.off()
