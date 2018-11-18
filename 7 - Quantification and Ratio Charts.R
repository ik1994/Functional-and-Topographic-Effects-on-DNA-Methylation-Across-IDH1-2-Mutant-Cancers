####Ratio Plots####

####STEPS####
#1. run the #5 rscript (without creating methylation plots)
#2. use subsetting code below to extract lists with just the chromosome of interest
#3. Find list dimensions

subset_1 <- subset(SNUC_B_hypo, SNUC_B_hypo$Chromosome == "chr1")
subset_2 <- subset(SNUC_B_hypo, SNUC_B_hypo$Chromosome == "chr2")
subset_3 <- subset(SNUC_B_hypo, SNUC_B_hypo$Chromosome == "chr3")
subset_4 <- subset(SNUC_B_hypo, SNUC_B_hypo$Chromosome == "chr4")
subset_5 <- subset(SNUC_B_hypo, SNUC_B_hypo$Chromosome == "chr5")
subset_6 <- subset(SNUC_B_hypo, SNUC_B_hypo$Chromosome == "chr6")
subset_7 <- subset(SNUC_B_hypo, SNUC_B_hypo$Chromosome == "chr7")
subset_8 <- subset(SNUC_B_hypo, SNUC_B_hypo$Chromosome == "chr8")
subset_9 <- subset(SNUC_B_hypo, SNUC_B_hypo$Chromosome == "chr9")
subset_10 <- subset(SNUC_B_hypo, SNUC_B_hypo$Chromosome == "chr10")
subset_11 <- subset(SNUC_B_hypo, SNUC_B_hypo$Chromosome == "chr11")
subset_12 <- subset(SNUC_B_hypo, SNUC_B_hypo$Chromosome == "chr12")
subset_13 <- subset(SNUC_B_hypo, SNUC_B_hypo$Chromosome == "chr13")
subset_14 <- subset(SNUC_B_hypo, SNUC_B_hypo$Chromosome == "chr14")
subset_15 <- subset(SNUC_B_hypo, SNUC_B_hypo$Chromosome == "chr15")
subset_16 <- subset(SNUC_B_hypo, SNUC_B_hypo$Chromosome == "chr16")
subset_17 <- subset(SNUC_B_hypo, SNUC_B_hypo$Chromosome == "chr17")
subset_18 <- subset(SNUC_B_hypo, SNUC_B_hypo$Chromosome == "chr18")
subset_19 <- subset(SNUC_B_hypo, SNUC_B_hypo$Chromosome == "chr19")
subset_20 <- subset(SNUC_B_hypo, SNUC_B_hypo$Chromosome == "chr20")
subset_21 <- subset(SNUC_B_hypo, SNUC_B_hypo$Chromosome == "chr21")
subset_22 <- subset(SNUC_B_hypo, SNUC_B_hypo$Chromosome == "chr22")


subset_1 <- as.data.frame(subset_1)
subset_2 <- as.data.frame(subset_2)
subset_3 <- as.data.frame(subset_3)
subset_4 <- as.data.frame(subset_4)
subset_5 <- as.data.frame(subset_5)
subset_6 <- as.data.frame(subset_6)
subset_7 <- as.data.frame(subset_7)
subset_8 <- as.data.frame(subset_8)
subset_9 <- as.data.frame(subset_9)
subset_10 <- as.data.frame(subset_10)
subset_11 <- as.data.frame(subset_11)
subset_12 <- as.data.frame(subset_12)
subset_13 <- as.data.frame(subset_13)
subset_14 <- as.data.frame(subset_14)
subset_15 <- as.data.frame(subset_15)
subset_16 <- as.data.frame(subset_16)
subset_17 <- as.data.frame(subset_17)
subset_18 <- as.data.frame(subset_18)
subset_19 <- as.data.frame(subset_19)
subset_20 <- as.data.frame(subset_20)
subset_21 <- as.data.frame(subset_21)
subset_22 <- as.data.frame(subset_22)


dim(subset_1)
dim(subset_2)
dim(subset_3)
dim(subset_4)
dim(subset_5)
dim(subset_6)
dim(subset_7)
dim(subset_8)
dim(subset_9)
dim(subset_10)
dim(subset_11)
dim(subset_12)
dim(subset_13)
dim(subset_14)
dim(subset_15)
dim(subset_16)
dim(subset_17)
dim(subset_18)
dim(subset_19)
dim(subset_20)
dim(subset_21)
dim(subset_22)





####to find total number of probes in each arm for each tumor type:

####STEPS####
#1. run the #5 rscript (without creating methylation plots)
#2. use subsetting code below to extract lists with just the chromosome of interest
#3. subset the subsets for arm of interest
#4. Find list dimensions

subset_1 <- subset(SNUC_E_hypo, SNUC_E_hypo$Chromosome == "chr1")
subset_1p <- subset(subset_1, subset_1$POS < 121535434)
subset_1q <- subset(subset_1, subset_1$POS > 124535434)

subset_2 <- subset(SNUC_E_hypo, SNUC_E_hypo$Chromosome == "chr2")
subset_2p <- subset(subset_2, subset_2$POS<92326171)
subset_2q <- subset(subset_2, subset_2$POS>95326171)

subset_3 <- subset(SNUC_E_hypo, SNUC_E_hypo$Chromosome == "chr3")
subset_3p <- subset(subset_3, subset_3$POS<90504854)
subset_3q <- subset(subset_3, subset_3$POS>93504854)

subset_4 <- subset(SNUC_E_hypo, SNUC_E_hypo$Chromosome == "chr4")
subset_4p <- subset(subset_4, subset_4$POS<49660117)
subset_4q <- subset(subset_4, subset_4$POS>52660117)

subset_5 <- subset(SNUC_E_hypo, SNUC_E_hypo$Chromosome == "chr5")
subset_5p <- subset(subset_5, subset_5$POS<46405641)
subset_5q <- subset(subset_5, subset_5$POS>49405641)

subset_6 <- subset(SNUC_E_hypo, SNUC_E_hypo$Chromosome == "chr6")
subset_6p <- subset(subset_6, subset_6$POS<58830166)
subset_6q <- subset(subset_6, subset_6$POS>61830166)

subset_7 <- subset(SNUC_E_hypo, SNUC_E_hypo$Chromosome == "chr7")
subset_7p <- subset(subset_7, subset_7$POS<58054331)
subset_7q <- subset(subset_7, subset_7$POS>61054331)

subset_8 <- subset(SNUC_E_hypo, SNUC_E_hypo$Chromosome == "chr8")
subset_8p <- subset(subset_8, subset_8$POS<43838887)
subset_8q <- subset(subset_8, subset_8$POS>46838887)

subset_9 <- subset(SNUC_E_hypo, SNUC_E_hypo$Chromosome == "chr9")
subset_9p <- subset(subset_9, subset_9$POS<47367679)
subset_9q <- subset(subset_9, subset_9$POS>50367679)

subset_10 <- subset(SNUC_E_hypo, SNUC_E_hypo$Chromosome == "chr10")
subset_10p <- subset(subset_10, subset_10$POS<39254935)
subset_10q <- subset(subset_10, subset_10$POS>42254935)

subset_11 <- subset(SNUC_E_hypo, SNUC_E_hypo$Chromosome == "chr11")
subset_11p <- subset(subset_11, subset_11$POS<51644205)
subset_11q <- subset(subset_11, subset_11$POS>54644205)

subset_12 <- subset(SNUC_E_hypo, SNUC_E_hypo$Chromosome == "chr12")
subset_12p <- subset(subset_12, subset_12$POS<34856694)
subset_12q <- subset(subset_12, subset_12$POS>37856694)

subset_13 <- subset(SNUC_E_hypo, SNUC_E_hypo$Chromosome == "chr13")
subset_13p <- subset(subset_13, subset_13$POS<16000000)
subset_13q <- subset(subset_13, subset_13$POS>19000000)

subset_14 <- subset(SNUC_E_hypo, SNUC_E_hypo$Chromosome == "chr14")
subset_14p <- subset(subset_14, subset_14$POS<16000000)
subset_14q <- subset(subset_14, subset_14$POS>19000000)

subset_15 <- subset(SNUC_E_hypo, SNUC_E_hypo$Chromosome == "chr15")
subset_15p <- subset(subset_15, subset_15$POS<17000000)
subset_15q <- subset(subset_15, subset_15$POS>20000000)

subset_16 <- subset(SNUC_E_hypo, SNUC_E_hypo$Chromosome == "chr16")
subset_16p <- subset(subset_16, subset_16$POS<35335801)
subset_16q <- subset(subset_16, subset_16$POS>38335801)

subset_17 <- subset(SNUC_E_hypo, SNUC_E_hypo$Chromosome == "chr17")
subset_17p <- subset(subset_17, subset_17$POS<22263006)
subset_17q <- subset(subset_17, subset_17$POS>25263006)

subset_18 <- subset(SNUC_E_hypo, SNUC_E_hypo$Chromosome == "chr18")
subset_18p <- subset(subset_18, subset_18$POS<15460898)
subset_18q <- subset(subset_18, subset_18$POS>18460898)

subset_19 <- subset(SNUC_E_hypo, SNUC_E_hypo$Chromosome == "chr19")
subset_19p <- subset(subset_19, subset_19$POS<24681782)
subset_19q <- subset(subset_19, subset_19$POS>27681782)

subset_20 <- subset(SNUC_E_hypo, SNUC_E_hypo$Chromosome == "chr20")
subset_20p <- subset(subset_20, subset_20$POS<26369569)
subset_20q <- subset(subset_20, subset_20$POS>29369569)

subset_21 <- subset(SNUC_E_hypo, SNUC_E_hypo$Chromosome == "chr21")
subset_21p <- subset(subset_21, subset_21$POS<11288129)
subset_21q <- subset(subset_21, subset_21$POS>14288129)

subset_22 <- subset(SNUC_E_hypo, SNUC_E_hypo$Chromosome == "chr22")
subset_22p <- subset(subset_22, subset_22$POS<13000000)
subset_22q <- subset(subset_22, subset_22$POS>16000000)



subset_1p <- as.data.frame(subset_1p)
subset_1q <- as.data.frame(subset_1q)

subset_2p <- as.data.frame(subset_2p)
subset_2q <- as.data.frame(subset_2q)

subset_3p <- as.data.frame(subset_3p)
subset_3q <- as.data.frame(subset_3q)

subset_4p <- as.data.frame(subset_4p)
subset_4q <- as.data.frame(subset_4q)

subset_5p <- as.data.frame(subset_5p)
subset_5q <- as.data.frame(subset_5q)

subset_6p <- as.data.frame(subset_6p)
subset_6q <- as.data.frame(subset_6q)

subset_7p <- as.data.frame(subset_7p)
subset_7q <- as.data.frame(subset_7q)

subset_8p <- as.data.frame(subset_8p)
subset_8q <- as.data.frame(subset_8q)

subset_9p <- as.data.frame(subset_9p)
subset_9q <- as.data.frame(subset_9q)

subset_10p <- as.data.frame(subset_10p)
subset_10q <- as.data.frame(subset_10q)

subset_11p <- as.data.frame(subset_11p)
subset_11q <- as.data.frame(subset_11q)

subset_12p <- as.data.frame(subset_12p)
subset_12q <- as.data.frame(subset_12q)

subset_13p <- as.data.frame(subset_13p)
subset_13q <- as.data.frame(subset_13q)

subset_14p <- as.data.frame(subset_14p)
subset_14q <- as.data.frame(subset_14q)

subset_15p <- as.data.frame(subset_15p)
subset_15q <- as.data.frame(subset_15q)

subset_16p <- as.data.frame(subset_16p)
subset_16q <- as.data.frame(subset_16q)

subset_17p <- as.data.frame(subset_17p)
subset_17q <- as.data.frame(subset_17q)

subset_18p <- as.data.frame(subset_18p)
subset_18q <- as.data.frame(subset_18q)

subset_19p <- as.data.frame(subset_19p)
subset_19q <- as.data.frame(subset_19q)

subset_20p <- as.data.frame(subset_20p)
subset_20q <- as.data.frame(subset_20q)

subset_21p <- as.data.frame(subset_21p)
subset_21q <- as.data.frame(subset_21q)

subset_22p <- as.data.frame(subset_22p)
subset_22q <- as.data.frame(subset_22q)


df <- data.frame(
  dim(subset_1p),
  dim(subset_1q),
  dim(subset_2p),
  dim(subset_2q),
  dim(subset_3p),
  dim(subset_3q),
  dim(subset_4p),
  dim(subset_4q),
  dim(subset_5p),
  dim(subset_5q),
  dim(subset_6p),
  dim(subset_6q),
  dim(subset_7p),
  dim(subset_7q),
  dim(subset_8p),
  dim(subset_8q),
  dim(subset_9p),
  dim(subset_9q),
  dim(subset_10p),
  dim(subset_10q),
  dim(subset_11p),
  dim(subset_11q),
  dim(subset_12p),
  dim(subset_12q),
  dim(subset_13p),
  dim(subset_13q),
  dim(subset_14p),
  dim(subset_14q),
  dim(subset_15p),
  dim(subset_15q),
  dim(subset_16p),
  dim(subset_16q),
  dim(subset_17p),
  dim(subset_17q),
  dim(subset_18p),
  dim(subset_18q),
  dim(subset_19p),
  dim(subset_19q),
  dim(subset_20p),
  dim(subset_20q),
  dim(subset_21p),
  dim(subset_21q),
  dim(subset_22p),
  dim(subset_22q))

write.excel <- function(x,row.names=FALSE,col.names=TRUE,...) {
  write.table(x,"clipboard",sep="\t",row.names=row.names,col.names=col.names,...)
}

write.excel(df)






####################################################################################################################################################









#### to find total number of probes in each chromosome:

####steps:
#1. load in the reference masterfile
#2. subset for a specific chrom
#3. Find dimensions

annot <- read.csv("C:/Users/rbled/Documents/PROJECT RESULTS/IDH Tumor Project/Results/Test Samples - Differential Methylation/Using M Values/test_samples_diffmeth_annot_mvalues_masterfile.csv")
subset_1 <- subset(annot, annot$chr == "chr1")
subset_2 <- subset(annot, annot$chr == "chr2")
subset_3 <- subset(annot, annot$chr == "chr3")
subset_4 <- subset(annot, annot$chr == "chr4")
subset_5 <- subset(annot, annot$chr == "chr5")
subset_6 <- subset(annot, annot$chr == "chr6")
subset_7 <- subset(annot, annot$chr == "chr7")
subset_8 <- subset(annot, annot$chr == "chr8")
subset_9 <- subset(annot, annot$chr == "chr9")
subset_10 <- subset(annot, annot$chr == "chr10")
subset_11 <- subset(annot, annot$chr == "chr11")
subset_12 <- subset(annot, annot$chr == "chr12")
subset_13 <- subset(annot, annot$chr == "chr13")
subset_14 <- subset(annot, annot$chr == "chr14")
subset_15 <- subset(annot, annot$chr == "chr15")
subset_16 <- subset(annot, annot$chr == "chr16")
subset_17 <- subset(annot, annot$chr == "chr17")
subset_18 <- subset(annot, annot$chr == "chr18")
subset_19 <- subset(annot, annot$chr == "chr19")
subset_20 <- subset(annot, annot$chr == "chr20")
subset_21 <- subset(annot, annot$chr == "chr21")
subset_22 <- subset(annot, annot$chr == "chr22")


subset_1 <- as.data.frame(subset_1)
subset_2 <- as.data.frame(subset_2)
subset_3 <- as.data.frame(subset_3)
subset_4 <- as.data.frame(subset_4)
subset_5 <- as.data.frame(subset_5)
subset_6 <- as.data.frame(subset_6)
subset_7 <- as.data.frame(subset_7)
subset_8 <- as.data.frame(subset_8)
subset_9 <- as.data.frame(subset_9)
subset_10 <- as.data.frame(subset_10)
subset_11 <- as.data.frame(subset_11)
subset_12 <- as.data.frame(subset_12)
subset_13 <- as.data.frame(subset_13)
subset_14 <- as.data.frame(subset_14)
subset_15 <- as.data.frame(subset_15)
subset_16 <- as.data.frame(subset_16)
subset_17 <- as.data.frame(subset_17)
subset_18 <- as.data.frame(subset_18)
subset_19 <- as.data.frame(subset_19)
subset_20 <- as.data.frame(subset_20)
subset_21 <- as.data.frame(subset_21)
subset_22 <- as.data.frame(subset_22)


dim(subset_1)
dim(subset_2)
dim(subset_3)
dim(subset_4)
dim(subset_5)
dim(subset_6)
dim(subset_7)
dim(subset_8)
dim(subset_9)
dim(subset_10)
dim(subset_11)
dim(subset_12)
dim(subset_13)
dim(subset_14)
dim(subset_15)
dim(subset_16)
dim(subset_17)
dim(subset_18)
dim(subset_19)
dim(subset_20)
dim(subset_21)
dim(subset_22)




####to find number of probes in each chromosome arm:

####steps:
#1. load in the reference masterfile
#2. subset for a specific chrom
#3. subset the subset for positions corresponding to each arm
#4. find dimensions

annot <- read.csv("C:/Users/rbled/Documents/PROJECT RESULTS/IDH Tumor Project/Results/Test Samples - Differential Methylation/Using M Values/test_samples_diffmeth_annot_mvalues_masterfile.csv")

subset_1 <- subset(annot, annot$chr == "chr1")
#subset_1 <- subset_1[grep("TSS200", subset_1$UCSC_RefGene_Group), ]
subset_1 <- subset(subset_1, subset_1$Enhancer == "TRUE")
subset_1p <- subset(subset_1, subset_1$pos<121535434)
subset_1q <- subset(subset_1, subset_1$pos>124535434)

subset_2 <- subset(annot, annot$chr == "chr2")
#subset_2 <- subset_2[grep("TSS200", subset_2$UCSC_RefGene_Group), ]
subset_2 <- subset(subset_2, subset_2$Enhancer == "TRUE")
subset_2p <- subset(subset_2, subset_2$pos<92326171)
subset_2q <- subset(subset_2, subset_2$pos>95326171)

subset_3 <- subset(annot, annot$chr == "chr3")
#subset_3 <- subset_3[grep("TSS200", subset_3$UCSC_RefGene_Group), ]
subset_3 <- subset(subset_3, subset_3$Enhancer == "TRUE")
subset_3p <- subset(subset_3, subset_3$pos<90504854)
subset_3q <- subset(subset_3, subset_3$pos>93504854)

subset_4 <- subset(annot, annot$chr == "chr4")
#subset_4 <- subset_4[grep("TSS200", subset_4$UCSC_RefGene_Group), ]
subset_4 <- subset(subset_4, subset_4$Enhancer == "TRUE")
subset_4p <- subset(subset_4, subset_4$pos<49660117)
subset_4q <- subset(subset_4, subset_4$pos>52660117)

subset_5 <- subset(annot, annot$chr == "chr5")
#subset_5 <- subset_5[grep("TSS200", subset_5$UCSC_RefGene_Group), ]
subset_5 <- subset(subset_5, subset_5$Enhancer == "TRUE")
subset_5p <- subset(subset_5, subset_5$pos<46405641)
subset_5q <- subset(subset_5, subset_5$pos>49405641)

subset_6 <- subset(annot, annot$chr == "chr6")
#subset_6 <- subset_6[grep("TSS200", subset_6$UCSC_RefGene_Group), ]
subset_6 <- subset(subset_6, subset_6$Enhancer == "TRUE")
subset_6p <- subset(subset_6, subset_6$pos<58830166)
subset_6q <- subset(subset_6, subset_6$pos>61830166)

subset_7 <- subset(annot, annot$chr == "chr7")
#subset_7 <- subset_7[grep("TSS200", subset_7$UCSC_RefGene_Group), ]
subset_7 <- subset(subset_7, subset_7$Enhancer == "TRUE")
subset_7p <- subset(subset_7, subset_7$pos<58054331)
subset_7q <- subset(subset_7, subset_7$pos>61054331)

subset_8 <- subset(annot, annot$chr == "chr8")
#subset_8 <- subset_8[grep("TSS200", subset_8$UCSC_RefGene_Group), ]
subset_8 <- subset(subset_8, subset_8$Enhancer == "TRUE")
subset_8p <- subset(subset_8, subset_8$pos<43838887)
subset_8q <- subset(subset_8, subset_8$pos>46838887)

subset_9 <- subset(annot, annot$chr == "chr9")
#subset_9 <- subset_9[grep("TSS200", subset_9$UCSC_RefGene_Group), ]
subset_9 <- subset(subset_9, subset_9$Enhancer == "TRUE")
subset_9p <- subset(subset_9, subset_9$pos<47367679)
subset_9q <- subset(subset_9, subset_9$pos>50367679)

subset_10 <- subset(annot, annot$chr == "chr10")
#subset_10 <- subset_10[grep("TSS200", subset_10$UCSC_RefGene_Group), ]
subset_10 <- subset(subset_10, subset_10$Enhancer == "TRUE")
subset_10p <- subset(subset_10, subset_10$pos<39254935)
subset_10q <- subset(subset_10, subset_10$pos>42254935)

subset_11 <- subset(annot, annot$chr == "chr11")
#subset_11 <- subset_11[grep("TSS200", subset_11$UCSC_RefGene_Group), ]
subset_11 <- subset(subset_11, subset_11$Enhancer == "TRUE")
subset_11p <- subset(subset_11, subset_11$pos<51644205)
subset_11q <- subset(subset_11, subset_11$pos>54644205)

subset_12 <- subset(annot, annot$chr == "chr12")
#subset_12 <- subset_12[grep("TSS200", subset_12$UCSC_RefGene_Group), ]
subset_12 <- subset(subset_12, subset_12$Enhancer == "TRUE")
subset_12p <- subset(subset_12, subset_12$pos<34856694)
subset_12q <- subset(subset_12, subset_12$pos>37856694)

subset_13 <- subset(annot, annot$chr == "chr13")
#subset_13 <- subset_13[grep("TSS200", subset_13$UCSC_RefGene_Group), ]
subset_13 <- subset(subset_13, subset_13$Enhancer == "TRUE")
subset_13p <- subset(subset_13, subset_13$pos<16000000)
subset_13q <- subset(subset_13, subset_13$pos>19000000)

subset_14 <- subset(annot, annot$chr == "chr14")
#subset_14 <- subset_14[grep("TSS200", subset_14$UCSC_RefGene_Group), ]
subset_14 <- subset(subset_14, subset_14$Enhancer == "TRUE")
subset_14p <- subset(subset_14, subset_14$pos<16000000)
subset_14q <- subset(subset_14, subset_14$pos>19000000)

subset_15 <- subset(annot, annot$chr == "chr15")
#subset_15 <- subset_15[grep("TSS200", subset_15$UCSC_RefGene_Group), ]
subset_15 <- subset(subset_15, subset_15$Enhancer == "TRUE")
subset_15p <- subset(subset_15, subset_15$pos<17000000)
subset_15q <- subset(subset_15, subset_15$pos>20000000)

subset_16 <- subset(annot, annot$chr == "chr16")
#subset_16 <- subset_16[grep("TSS200", subset_16$UCSC_RefGene_Group), ]
subset_16 <- subset(subset_16, subset_16$Enhancer == "TRUE")
subset_16p <- subset(subset_16, subset_16$pos<35335801)
subset_16q <- subset(subset_16, subset_16$pos>38335801)

subset_17 <- subset(annot, annot$chr == "chr17")
#subset_17 <- subset_17[grep("TSS200", subset_17$UCSC_RefGene_Group), ]
subset_17 <- subset(subset_17, subset_17$Enhancer == "TRUE")
subset_17p <- subset(subset_17, subset_17$pos<22263006)
subset_17q <- subset(subset_17, subset_17$pos>25263006)

subset_18 <- subset(annot, annot$chr == "chr18")
#subset_18 <- subset_18[grep("TSS200", subset_18$UCSC_RefGene_Group), ]
subset_18 <- subset(subset_18, subset_18$Enhancer == "TRUE")
subset_18p <- subset(subset_18, subset_18$pos<15460898)
subset_18q <- subset(subset_18, subset_18$pos>18460898)

subset_19 <- subset(annot, annot$chr == "chr19")
#subset_19 <- subset_19[grep("TSS200", subset_19$UCSC_RefGene_Group), ]
subset_19 <- subset(subset_19, subset_19$Enhancer == "TRUE")
subset_19p <- subset(subset_19, subset_19$pos<24681782)
subset_19q <- subset(subset_19, subset_19$pos>27681782)

subset_20 <- subset(annot, annot$chr == "chr20")
#subset_20 <- subset_20[grep("TSS200", subset_20$UCSC_RefGene_Group), ]
subset_20 <- subset(subset_20, subset_20$Enhancer == "TRUE")
subset_20p <- subset(subset_20, subset_20$pos<26369569)
subset_20q <- subset(subset_20, subset_20$pos>29369569)

subset_21 <- subset(annot, annot$chr == "chr21")
#subset_21 <- subset_21[grep("TSS200", subset_21$UCSC_RefGene_Group), ]
subset_21 <- subset(subset_21, subset_21$Enhancer == "TRUE")
subset_21p <- subset(subset_21, subset_21$pos<11288129)
subset_21q <- subset(subset_21, subset_21$pos>14288129)

subset_22 <- subset(annot, annot$chr == "chr22")
#subset_22 <- subset_22[grep("TSS200", subset_22$UCSC_RefGene_Group), ]
subset_22 <- subset(subset_22, subset_22$Enhancer == "TRUE")
subset_22p <- subset(subset_22, subset_22$pos<13000000)
subset_22q <- subset(subset_22, subset_22$pos>16000000)



subset_1p <- as.data.frame(subset_1p)
subset_1q <- as.data.frame(subset_1q)

subset_2p <- as.data.frame(subset_2p)
subset_2q <- as.data.frame(subset_2q)

subset_3p <- as.data.frame(subset_3p)
subset_3q <- as.data.frame(subset_3q)

subset_4p <- as.data.frame(subset_4p)
subset_4q <- as.data.frame(subset_4q)

subset_5p <- as.data.frame(subset_5p)
subset_5q <- as.data.frame(subset_5q)

subset_6p <- as.data.frame(subset_6p)
subset_6q <- as.data.frame(subset_6q)

subset_7p <- as.data.frame(subset_7p)
subset_7q <- as.data.frame(subset_7q)

subset_8p <- as.data.frame(subset_8p)
subset_8q <- as.data.frame(subset_8q)

subset_9p <- as.data.frame(subset_9p)
subset_9q <- as.data.frame(subset_9q)

subset_10p <- as.data.frame(subset_10p)
subset_10q <- as.data.frame(subset_10q)

subset_11p <- as.data.frame(subset_11p)
subset_11q <- as.data.frame(subset_11q)

subset_12p <- as.data.frame(subset_12p)
subset_12q <- as.data.frame(subset_12q)

subset_13p <- as.data.frame(subset_13p)
subset_13q <- as.data.frame(subset_13q)

subset_14p <- as.data.frame(subset_14p)
subset_14q <- as.data.frame(subset_14q)

subset_15p <- as.data.frame(subset_15p)
subset_15q <- as.data.frame(subset_15q)

subset_16p <- as.data.frame(subset_16p)
subset_16q <- as.data.frame(subset_16q)

subset_17p <- as.data.frame(subset_17p)
subset_17q <- as.data.frame(subset_17q)

subset_18p <- as.data.frame(subset_18p)
subset_18q <- as.data.frame(subset_18q)

subset_19p <- as.data.frame(subset_19p)
subset_19q <- as.data.frame(subset_19q)

subset_20p <- as.data.frame(subset_20p)
subset_20q <- as.data.frame(subset_20q)

subset_21p <- as.data.frame(subset_21p)
subset_21q <- as.data.frame(subset_21q)

subset_22p <- as.data.frame(subset_22p)
subset_22q <- as.data.frame(subset_22q)


df <-data.frame(
  dim(subset_1p),
  dim(subset_1q),
  dim(subset_2p),
  dim(subset_2q),
  dim(subset_3p),
  dim(subset_3q),
  dim(subset_4p),
  dim(subset_4q),
  dim(subset_5p),
  dim(subset_5q),
  dim(subset_6p),
  dim(subset_6q),
  dim(subset_7p),
  dim(subset_7q),
  dim(subset_8p),
  dim(subset_8q),
  dim(subset_9p),
  dim(subset_9q),
  dim(subset_10p),
  dim(subset_10q),
  dim(subset_11p),
  dim(subset_11q),
  dim(subset_12p),
  dim(subset_12q),
  dim(subset_13p),
  dim(subset_13q),
  dim(subset_14p),
  dim(subset_14q),
  dim(subset_15p),
  dim(subset_15q),
  dim(subset_16p),
  dim(subset_16q),
  dim(subset_17p),
  dim(subset_17q),
  dim(subset_18p),
  dim(subset_18q),
  dim(subset_19p),
  dim(subset_19q),
  dim(subset_20p),
  dim(subset_20q),
  dim(subset_21p),
  dim(subset_21q),
  dim(subset_22p),
  dim(subset_22q))

write.excel <- function(x,row.names=FALSE,col.names=TRUE,...) {
  write.table(x,"clipboard",sep="\t",row.names=row.names,col.names=col.names,...)
}

write.excel(df)