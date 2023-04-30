################################################################################
# Master 2 - DSDM
# Code 02.01 - Get beta values mean for Mayne CpGs
################################################################################
# Author: Orlane LE QUELLENNEC
# Date: 30/08/2022
################################################################################
# Main objectives:
#------------------------------------------------------------------------------#
# For Mayne CpGs we compute mean beta values 
# Real mean beta values were requested, but no answer until now 
################################################################################
################################################################################
# R packages
#------------------------------------------------------------------------------#

# To get access to public dataset used by Mayne
library(GEOquery) 
# Memento GEOquery data 
# pData(gse) ## print the sample information
# fData(gse) ## print the gene annotation
# exprs(gse)[1:2,1:8] ## print the expression data

library(dplyr)

################################################################################
# Partie I : Get dataset data
################################################################################
# Datasets used for Mayne's clock development are : 
# - GSE31781
# - GSE36829
# - GSE74738
# - GSE44667

# 1) GSE31781 : used all subjects 
gse <- getGEO("GSE31781")
gse <- gse[[1]]
mean_31781 <- as.data.frame(apply(exprs(gse),1,mean))

# 2) GSE36829 : used all subjects 
gse <- getGEO("GSE36829")
gse <- gse[[1]]
mean_36829 <- as.data.frame(apply(exprs(gse),1,mean))

# 3) GSE74738 : used 28 subjects among 79
# inclusion criterion were : healthy placental samples aged between [36; 42] week
gse <- getGEO("GSE74738")
gse <- gse[[1]]

# We select only placental sample 
pData(gse) = pData(gse)[which(pData(gse)[, "sample tissue:ch1"] %in% c("placental trophoblast", 
                                               "placental decidua", 
                                               "placental mesenchyme", 
                                               "placental chorionic villi")),]

# We select gestational age btw [36; 42]
pData(gse) = pData(gse)[which(!(pData(gse)[, "fetal gestational age (weeks)/trimester:ch1"] %in% c("1st_trimester", 
                                                                   "2nd_trimester", 
                                                                   "<10w", "8w 48 days DA", 
                                                                   "8w 37d DA", "9w 26-30d DA", 
                                                                   "12w 35-38d", "15w by dates 9w size 53d"))),]

pData(gse) = pData(gse)[which(as.numeric(pData(gse)[, "fetal gestational age (weeks)/trimester:ch1"])>=36),]

# We select healthy placenta 
pData(gse) = pData(gse)[which(pData(gse)[,"status/group:ch1"]=="control"),]

exp_gene = exprs(gse)[,which(names(as.data.frame(exprs(gse))) %in% rownames(pData(gse)))]
mean_74738 <- as.data.frame(apply(exp_gene,1,mean))

# 4) GSE44667 : used 20 subjects among 40
# inclusion criterion were : non-preeclamptic placental samples 
gse <- getGEO("GSE44667")
gse <- gse[[1]]

# those without preeclampsia 
pData(gse) = pData(gse)[which(pData(gse)[, "condition:ch1"]=="Control"),]  

exp_gene = as.data.frame(exprs(gse)) %>% select(rownames(pData(gse)))
mean_44667 <- as.data.frame(apply(exp_gene,1,mean))

################################################################################
# Partie II : Get mean values for Mayne's CpGs
################################################################################
# Only shared CpGs between all datasets were selected
tmp_1 = merge(mean_31781, mean_36829, by="row.names")
names(tmp_1) <- c("CpG", "mean.x", "mean.y")

tmp_2 = merge(mean_74738, mean_44667, by="row.names")
names(tmp_2) <- c("CpG", "mean.x", "mean.y")

tmp = merge(tmp_1, tmp_2, by="CpG")

dim(tmp)

# Mean were computed for each CpG
global.mean = (tmp[,2]+tmp[,3]+tmp[,4]+tmp[,5])/4
mean_CpG = cbind(tmp[,1], global.mean)

################################################################################
# Partie III : Save the dataset
################################################################################
# This file miss the Mayne's dataset with 22 samples 
setwd("~/data/")
saveRDS(mean_CpG,  file = "mean_value_mayne_clock.rds")

################################################################################
rm(list = ls())
################################################################################


