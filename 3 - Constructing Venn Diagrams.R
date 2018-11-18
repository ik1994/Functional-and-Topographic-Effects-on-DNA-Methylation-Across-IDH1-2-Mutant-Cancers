############VENN DIAGRAM ANALYSIS################
##########TEST SAMPLES###########
########HYPOMETHYLATED PROBES#######
project_dir="C:/Users/Ramona Bledea/Documents/PROJECT RESULTS/IDH Tumor Project/Results/Test Samples - Differential Methylation/Using M Values"
venn_output_dir <- file.path(project_dir, "Hypomethylated Probes")
#dir.create(path = venn_output_dir, recursive = TRUE)
subsetted_probes_list <- read.csv("C:/Users/Ramona Bledea/Documents/PROJECT RESULTS/IDH Tumor Project/Results/Test Samples - Differential Methylation/Using M Values/Hypomethylated Probes/test_samples_hypomethylated_probes_list.csv")
subsetted_probes_list_Astro <- subsetted_probes_list$Astrocytoma
subsetted_probes_list_AML <- subsetted_probes_list$AML
subsetted_probes_list_Chol <- subsetted_probes_list$Cholangiocarcinoma
subsetted_probes_list_SNUC <- subsetted_probes_list$SNUC
subsetted_probes_list_Breast <- subsetted_probes_list$Breast
subsetted_probes_list_Oligo <- subsetted_probes_list$Oligodendroglioma
subsetted_probes_list_Brain <- subsetted_probes_list$Brain

library(VennDiagram)
library(rowr)

#####FOR NON-BRAIN TUMORS ONLY######
probe_list_combined <-cbind.fill(subsetted_probes_list_AML, subsetted_probes_list_Chol, subsetted_probes_list_SNUC, subsetted_probes_list_Breast, subsetted_probes_list_Brain, fill="")
probe_list_combined <- as.list(probe_list_combined)
make_venn_from_annotation_files <- function(output_file = FALSE){
  if(output_file != FALSE) pdf(file = output_file, width = 10, height = 10, onefile = TRUE)
  #alpha_values <- rep(0.2, dim(probe_list_combined))
  #venn_colors <- get_plot_colors(probe_list_combined)
  venn_plot <- venn.diagram(x = probe_list_combined, category.names = c("AML", "Cholangiocarcinoma", "SNUC", "Breast Cancer", "Brain"),
                            filename = NULL,
                            fill = c("blue","red","black", "purple", "dark green"),
                            col = c("blue","red","black", "purple", "dark green"),
                            cat.col = c("blue","red","black", "purple", "dark green"),
                            cat.cex = 0.8,
                            cat.dist = 0.07,
                            #alpha = alpha_values, 
                            scaled = TRUE)
  grid.draw(venn_plot)
  
  if(output_file != FALSE) dev.off()
}
make_venn_from_annotation_files(output_file = file.path(venn_output_dir, "test_hypo_venn_all_tumors - AML, Brain, Chol, SNUC, Breast.pdf"))

#########FOR BRAIN TUMORS ONLY###########
probe_list_combined <-cbind.fill(subsetted_probes_list_Astro, subsetted_probes_list_Oligo, fill="")
probe_list_combined <- as.list(probe_list_combined)
make_venn_from_annotation_files <- function(output_file = FALSE){
  if(output_file != FALSE) pdf(file = output_file, width = 10, height = 10, onefile = TRUE)
  #alpha_values <- rep(0.2, dim(probe_list_combined))
  #venn_colors <- get_plot_colors(probe_list_combined)
  venn_plot <- venn.diagram(x = probe_list_combined, category.names = c("Astrocytoma", "Oligodendroglioma"),
                            filename = NULL,
                            fill = c("blue","red"),
                            col = c("blue","red"),
                            cat.col = c("blue","red"),
                            cat.cex = 0.8,
                            cat.dist = 0.07,
                            #alpha = alpha_values, 
                            scaled = FALSE)
  grid.draw(venn_plot)
  
  if(output_file != FALSE) dev.off()
}
make_venn_from_annotation_files(output_file = file.path(venn_output_dir, "test_hypo_venn_brain_tumors - Astro, Oligo.pdf"))

##########CONTROL SAMPLES###########
########HYPOMETHYLATED PROBES#######
project_dir="C:/Users/Ramona Bledea/Documents/PROJECT RESULTS/IDH Tumor Project/Results/Control Samples - Differential Methylation/Using M Values"
venn_output_dir <- file.path(project_dir, "Hypomethylated Probes")
#dir.create(path = venn_output_dir, recursive = TRUE)
subsetted_probes_list_2 <- read.csv("C:/Users/Ramona Bledea/Documents/PROJECT RESULTS/IDH Tumor Project/Results/Control Samples - Differential Methylation/Using M Values/Hypomethylated Probes/control_samples_hypomethylated_probes_list.csv")
subsetted_probes_list_Blood <- subsetted_probes_list_2$Blood
subsetted_probes_list_GBM <- subsetted_probes_list_2$GBM
subsetted_probes_list_Neuro <- subsetted_probes_list_2$Neurocytoma
subsetted_probes_list_PAstro <- subsetted_probes_list_2$Pilocytic_Astrocytoma
subsetted_probes_list_SUDEP <- subsetted_probes_list_2$SUDEP

#####ALL TUMORS######
probe_list_combined <-cbind.fill(subsetted_probes_list_Blood, subsetted_probes_list_GBM, subsetted_probes_list_Neuro, subsetted_probes_list_PAstro, subsetted_probes_list_SUDEP, fill="")
probe_list_combined <- as.list(probe_list_combined)
make_venn_from_annotation_files <- function(output_file = FALSE){
  if(output_file != FALSE) pdf(file = output_file, width = 10, height = 10, onefile = TRUE)
  #alpha_values <- rep(0.2, dim(probe_list_combined))
  #venn_colors <- get_plot_colors(probe_list_combined)
  venn_plot <- venn.diagram(x = probe_list_combined, category.names = c("Blood", "GBM", "Neurocytoma", "Pilocytic Astrocytoma", "SUDEP"),
                            filename = NULL,
                            fill = c("blue","red","black", "purple", "dark green"),
                            col = c("blue","red","black", "purple", "dark green"),
                            cat.col = c("blue","red","black", "purple", "dark green"),
                            cat.cex = 0.8,
                            cat.dist = 0.07,
                            #alpha = alpha_values, 
                            scaled = TRUE)
  grid.draw(venn_plot)
  
  if(output_file != FALSE) dev.off()
}
make_venn_from_annotation_files(output_file = file.path(venn_output_dir, "control_hypo_venn_all_tumors_2 - Blood, GBM, Neuro, Astro, SUDEP.pdf"))

#####CONTROL VS TEST#####
venn_output_dir="C:/Users/Ramona Bledea/Documents/PROJECT RESULTS/IDH Tumor Project/Results/All Figures/Control vs Test Samples Venn Diagrams/Hypomethylated"
probe_list_combined <-cbind.fill(subsetted_probes_list_Astro, subsetted_probes_list_Oligo, subsetted_probes_list_PAstro, subsetted_probes_list_Neuro, subsetted_probes_list_GBM, fill="")
probe_list_combined <- as.list(probe_list_combined)
make_venn_from_annotation_files <- function(output_file = FALSE){
  if(output_file != FALSE) pdf(file = output_file, width = 10, height = 10, onefile = TRUE)
  #alpha_values <- rep(0.2, dim(probe_list_combined))
  #venn_colors <- get_plot_colors(probe_list_combined)
  venn_plot <- venn.diagram(x = probe_list_combined, category.names = c("Astrocytoma, IDH Mutant", "Oligodendroglioma, IDH Mutant", "Pilocytic Astrocytoma, IDH WT", "Neurocytoma, IDH WT", "GBM, IDH WT"),
                            filename = NULL,
                            fill = c("blue","red","black", "purple", "dark green"),
                            col = c("blue","red","black", "purple", "dark green"),
                            cat.col = c("blue","red","black", "purple", "dark green"),
                            cat.cex = 0.8,
                            cat.dist = 0.07,
                            #alpha = alpha_values, 
                            scaled = FALSE)
  grid.draw(venn_plot)
  
  if(output_file != FALSE) dev.off()
}
make_venn_from_annotation_files(output_file = file.path(venn_output_dir, "brain_control_test_hypo_venn_2.pdf"))


##########TEST SAMPLES###########
########HYPERMETHYLATED PROBES#######
project_dir="C:/Users/Ramona Bledea/Documents/PROJECT RESULTS/IDH Tumor Project/Results/Test Samples - Differential Methylation/Using M Values"
venn_output_dir <- file.path(project_dir, "Hypermethylated Probes")
#dir.create(path = venn_output_dir, recursive = TRUE)
subsetted_probes_list <- read.csv("C:/Users/Ramona Bledea/Documents/PROJECT RESULTS/IDH Tumor Project/Results/Test Samples - Differential Methylation/Using M Values/Hypermethylated Probes/test_samples_hypermethylated_probes_list.csv")
subsetted_probes_list_Astro <- subsetted_probes_list$Astrocytoma
subsetted_probes_list_AML <- subsetted_probes_list$AML
subsetted_probes_list_Chol <- subsetted_probes_list$Cholangiocarcinoma
subsetted_probes_list_SNUC <- subsetted_probes_list$SNUC
subsetted_probes_list_Breast <- subsetted_probes_list$Breast
subsetted_probes_list_Oligo <- subsetted_probes_list$Oligodendroglioma
subsetted_probes_list_Brain <- subsetted_probes_list$Brain

library(VennDiagram)
library(rowr)

#####FOR NON-BRAIN TUMORS ONLY######
probe_list_combined <-cbind.fill(subsetted_probes_list_AML, subsetted_probes_list_Chol, subsetted_probes_list_SNUC, subsetted_probes_list_Breast, subsetted_probes_list_Brain, fill="")
probe_list_combined <- as.list(probe_list_combined)
make_venn_from_annotation_files <- function(output_file = FALSE){
  if(output_file != FALSE) pdf(file = output_file, width = 10, height = 10, onefile = TRUE)
  #alpha_values <- rep(0.2, dim(probe_list_combined))
  #venn_colors <- get_plot_colors(probe_list_combined)
  venn_plot <- venn.diagram(x = probe_list_combined, category.names = c("AML", "Cholangiocarcinoma", "SNUC", "Breast Cancer", "Brain"),
                            filename = NULL,
                            fill = c("blue","red","black", "purple", "dark green"),
                            col = c("blue","red","black", "purple", "dark green"),
                            cat.col = c("blue","red","black", "purple", "dark green"),
                            cat.cex = 0.8,
                            cat.dist = 0.07,
                            #alpha = alpha_values, 
                            scaled = TRUE)
  grid.draw(venn_plot)
  
  if(output_file != FALSE) dev.off()
}
make_venn_from_annotation_files(output_file = file.path(venn_output_dir, "test_hyper_venn_all_tumors - AML, Chol, SNUC, Breast, Brain.pdf"))

#########FOR BRAIN TUMORS ONLY###########
probe_list_combined <-cbind.fill(subsetted_probes_list_Astro, subsetted_probes_list_Oligo, fill="")
probe_list_combined <- as.list(probe_list_combined)
make_venn_from_annotation_files <- function(output_file = FALSE){
  if(output_file != FALSE) pdf(file = output_file, width = 10, height = 10, onefile = TRUE)
  #alpha_values <- rep(0.2, dim(probe_list_combined))
  #venn_colors <- get_plot_colors(probe_list_combined)
  venn_plot <- venn.diagram(x = probe_list_combined, category.names = c("Astrocytoma", "Oligodendroglioma"),
                            filename = NULL,
                            fill = c("blue","red"),
                            col = c("blue","red"),
                            cat.col = c("blue","red"),
                            cat.cex = 0.8,
                            cat.dist = 0.07,
                            #alpha = alpha_values, 
                            scaled = FALSE)
  grid.draw(venn_plot)
  
  if(output_file != FALSE) dev.off()
}
make_venn_from_annotation_files(output_file = file.path(venn_output_dir, "test_hyper_venn_brain_tumors - Astro, Oligo.pdf"))

##########CONTROL SAMPLES###########
########HYPERMETHYLATED PROBES#######
project_dir="C:/Users/Ramona Bledea/Documents/PROJECT RESULTS/IDH Tumor Project/Results/Control Samples - Differential Methylation/Using M Values"
venn_output_dir <- file.path(project_dir, "Hypermethylated Probes")
#dir.create(path = venn_output_dir, recursive = TRUE)
subsetted_probes_list_2 <- read.csv("C:/Users/Ramona Bledea/Documents/PROJECT RESULTS/IDH Tumor Project/Results/Control Samples - Differential Methylation/Using M Values/Hypermethylated Probes/control_samples_hypermethylated_probes_list.csv")
subsetted_probes_list_Blood <- subsetted_probes_list_2$Blood
subsetted_probes_list_GBM <- subsetted_probes_list_2$GBM
subsetted_probes_list_Neuro <- subsetted_probes_list_2$Neurocytoma
subsetted_probes_list_PAstro <- subsetted_probes_list_2$Pilocytic_Astrocytoma
subsetted_probes_list_SUDEP <- subsetted_probes_list_2$SUDEP

#####ALL TUMORS######
probe_list_combined <-cbind.fill(subsetted_probes_list_Blood, subsetted_probes_list_GBM, subsetted_probes_list_Neuro, subsetted_probes_list_PAstro, subsetted_probes_list_SUDEP, fill="")
probe_list_combined <- as.list(probe_list_combined)
make_venn_from_annotation_files <- function(output_file = FALSE){
  if(output_file != FALSE) pdf(file = output_file, width = 10, height = 10, onefile = TRUE)
  #alpha_values <- rep(0.2, dim(probe_list_combined))
  #venn_colors <- get_plot_colors(probe_list_combined)
  venn_plot <- venn.diagram(x = probe_list_combined, category.names = c("Blood", "GBM", "Neurocytoma", "Pilocytic Astrocytoma", "SUDEP"),
                            filename = NULL,
                            fill = c("blue","red","black", "purple", "dark green"),
                            col = c("blue","red","black", "purple", "dark green"),
                            cat.col = c("blue","red","black", "purple", "dark green"),
                            cat.cex = 0.8,
                            cat.dist = 0.07,
                            #alpha = alpha_values, 
                            scaled = TRUE)
  grid.draw(venn_plot)
  
  if(output_file != FALSE) dev.off()
}
make_venn_from_annotation_files(output_file = file.path(venn_output_dir, "control_hyper_venn_all_tumors_2 - Blood, GBM, Neuro, Astro, SUDEP.pdf"))

#####CONTROL VS TEST#####
project_dir="C:/Users/Ramona Bledea/Documents/PROJECT RESULTS/IDH Tumor Project/Results"
venn_output_dir <- file.path(project_dir, "Control vs Test Samples Venn Diagrams")
probe_list_combined <-cbind.fill(subsetted_probes_list_Astro, subsetted_probes_list_PAstro, fill="")
probe_list_combined <- as.list(probe_list_combined)
make_venn_from_annotation_files <- function(output_file = FALSE){
  if(output_file != FALSE) pdf(file = output_file, width = 10, height = 10, onefile = TRUE)
  alpha_values <- rep(0.2, dim(probe_list_combined))
  #venn_colors <- get_plot_colors(probe_list_combined)
  venn_plot <- venn.diagram(x = probe_list_combined, category.names = c("Astrocytoma, IDH Mutant", "Astrocytoma, IDH WT"),
                            filename = NULL,
                            fill = c("blue","red"),
                            col = c("blue","red"),
                            cat.col = c("blue","red"),
                            cat.cex = 0.8,
                            cat.dist = 0.07,
                            alpha = alpha_values, 
                            scaled = FALSE)
  grid.draw(venn_plot)
  
  if(output_file != FALSE) dev.off()
}
make_venn_from_annotation_files(output_file = file.path(venn_output_dir, "hyper_venn_Astro_WTvsMut.pdf"))

##double check the numbers on the plots with the dim of each column in list










