################################################################################
# Epigenetic gestational age estimations with existing clocks
################################################################################
# Author: Orlane Rossini
# Date: 30/08/2022
################################################################################
# Main objectives:
#------------------------------------------------------------------------------#
# For EDEN (450K / EPIC) PELAGIE and SEPAGES 
# 1 - Get the mean of missing CpGs
# 2 - Get epigenetic gestational age estimations
# 3 - Compute residuals with imputed epigenetic GA and non imputed
################################################################################
################################################################################
# Paths 
#------------------------------------------------------------------------------#
setwd("~/data/")

# all db
File <- "private_db/db_pooled_imp.rds"

# Methylation data
File_cpgs = "methylation_data/ssnoob_processed.rds"

# Mean values
File_mean_Mayne = "methylation_data/mean_value_mayne_clock.rds"
File_mean_RPC = "methylation_data/mean_value_lee_rpc_clock.rds"
File_mean_CPC = "methylation_data/mean_value_lee_cpc_clock.rds"
File_mean_RRPC = "methylation_data/mean_value_lee_rrpc_clock.rds"
File_Mayne_CpG = "clock_coefs/Maynes_clock_coefs.csv"
File_Lee_CpG = "clock_coefs/Lees_clock_coefs.csv"

################################################################################
# R packages
#------------------------------------------------------------------------------#

# INSTALLATION 
# install.packages("ps")
# library(ps)
# install.packages("backports"); install.packages("devtools")
# library(devtools)
# if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
# BiocManager::install("methylclockData",force=T)
# BiocManager::install("methylclock",force=T)

library(methylclockData)
library(methylclock) # for Mayne and Lee estimations
library(dplyr) 

################################################################################
# Import data into R session
#------------------------------------------------------------------------------#
# pooled db from Emie 
db_pooled <- readRDS(File)

# cpgs file
cpgs <- readRDS(File_cpgs)

# Mean values
db_mean_mayne <- readRDS(File_mean_Mayne)
db_mean_rpc <- readRDS(File_mean_RPC)
db_mean_cpc <- readRDS(File_mean_CPC)
db_mean_rrpc <- readRDS(File_mean_RRPC)
db_mayne_cpg <- read.csv(File_Mayne_CpG, skip = 3, header = T, sep = ";")
db_lee_cpg <- read.csv(File_Lee_CpG, header = T, sep = ",")

################################################################################
# Partie I : Get DNAm data 
################################################################################
# I.1) Get nb of missing CpGs for each experiment / each clock
list_mayne = db_mayne_cpg[,1]
length(list_mayne) # 62 CpGs
list_rpc = db_lee_cpg[which(db_lee_cpg$Coefficient_RPC!=0),1]; list_rpc = list_rpc[-1]
length(list_rpc) # 558 CpGs
list_cpc = db_lee_cpg[which(db_lee_cpg$Coefficient_CPC!=0),1]; list_cpc = list_cpc[-1]
length(list_cpc) # 546 CpGs
list_rrpc = db_lee_cpg[which(db_lee_cpg$Coefficient_refined_RPC!=0),1]; list_rrpc = list_rrpc[-1]
length(list_rrpc) # 395 CpGs

# for (c in list(list_mayne, list_rpc, list_cpc, list_rrpc))
# {
  # print(length(setdiff(c, rownames(cpgs))))
  # 9 cpgs / 43 cpgs / 49 cpgs / 32 cpgs
# }

# I.2) Add CpG mean value if the CpG is missing 
# We keep only necessary CpGs
cpgs.missing.GA = checkClocksGA(cpgs[1:5,])
Mayne_coefs = commonClockCpgs(cpgs.missing.GA, "Mayne")
Lee_coefs = commonClockCpgs(cpgs.missing.GA, "Lee")
length(Mayne_coefs) #  62 cpgs manquants 
length(Lee_coefs) #  1125 cpgs manquants

list_CpG = unique(c(Mayne_coefs, Lee_coefs, rownames(cpgs[1:5,]))); 
length(list_CpG) # 1,191 CpGs are necessary for the all 4 clocks

dim(cpgs) # 333,229 CpGs in total 

cpgs = cpgs[which(rownames(cpgs)%in%list_CpG), ]
dim(cpgs) # 1,068 CpGs among the 1,191 required

# create a commun db for mean values of cpgs 
db_mean = as.data.frame(rbind(db_mean_mayne, db_mean_rpc, db_mean_cpc, db_mean_rrpc))
db_mean$global.mean = as.numeric(db_mean$global.mean)

db_mean = aggregate(global.mean~V1,mean, data = db_mean) # take the mean for common CpGs
dim(db_mean) # All mean values for CpGs used in Lee's clocks and Mayne's clock

# Get all missing CpGS 
cpgs_imp = cpgs
missing = setdiff(list_CpG, rownames(cpgs_imp))
tmp = db_mean[which(db_mean[,1] %in% missing),] # Mean values for missing cpgs
print(dim(tmp)) # 123 CpGs are missing 

# Add the mean for missing CpGs
for (i in 1:nrow(tmp))
{
  print(i)
  cpgs_imp = rbind(cpgs_imp, rep(as.numeric(tmp[i,2]), ncol(cpgs_imp)))
  rownames(cpgs_imp)[nrow(cpgs_imp)] <- tmp[i,1]
}

dim(cpgs_imp) # 1191 CpGs among the 1191 required
dim(cpgs) # 1068 CpGs
# No more missing CpGs

################################################################################
# Partie II : eGA estimation
################################################################################
# II.1) For non-imputed cpgs
# We keep only Mayne and Lee's clocks
cpgs_eGA <- DNAmGA(cpgs); cpgs_eGA = cpgs_eGA[,c(1,4:7)] 
dim(cpgs_eGA) # 1539 obs and 5 observations

names(cpgs_eGA) = c("Sample_Id","ega_mayne", "ega_rpc", "ega_cpc", "ega_rrpc")

# II.2) For imputed cpgs
# We keep only Mayne and Lee's clocks
cpgs_eGA_imp <- DNAmGA(cpgs_imp); cpgs_eGA_imp = cpgs_eGA_imp[,c(1,4:7)] 
dim(cpgs_eGA_imp) # 1539 obs and 5 observations

names(cpgs_eGA_imp) = c("Sample_Id","ega_mayne_imp", "ega_rpc_imp", "ega_cpc_imp", "ega_rrpc_imp")

# II.3) We merge cohorte by cohorte 
db_pooled = merge(db_pooled, cpgs_eGA, by = "Sample_Id")
db_pooled = merge(db_pooled, cpgs_eGA_imp, by = "Sample_Id")
dim(db_pooled) 
# 1525 obs and 63 variables 

################################################################################
# Partie II : RAA = residuals lm (eGA ~ GA)
################################################################################
# For non imputed cpgs
reg.mayne = lm(ega_mayne~clim_gest_duration, data = db_pooled)
db_pooled$raa.mayne= rstudent(reg.mayne)

reg.rpc = lm(ega_rpc~clim_gest_duration, data = db_pooled)
db_pooled$raa.rpc = rstudent(reg.rpc)

reg.cpc = lm(ega_cpc~clim_gest_duration, data = db_pooled)
db_pooled$raa.cpc= rstudent(reg.cpc)

reg.rrpc = lm(ega_rrpc~clim_gest_duration, data = db_pooled)
db_pooled$raa.rrpc = rstudent(reg.rrpc)

dim(db_pooled) # 1525 obs and 67 variables

# For imputed cpgs
reg.mayne_imp = lm(ega_mayne_imp~clim_gest_duration, data = db_pooled)
db_pooled$raa.mayne_imp= rstudent(reg.mayne_imp)

reg.rpc_imp = lm(ega_rpc_imp~clim_gest_duration, data = db_pooled)
db_pooled$raa.rpc_imp = rstudent(reg.rpc_imp)

reg.cpc_imp = lm(ega_cpc_imp~clim_gest_duration, data = db_pooled)
db_pooled$raa.cpc_imp= rstudent(reg.cpc_imp)

reg.rrpc_imp = lm(ega_rrpc_imp~clim_gest_duration, data = db_pooled)
db_pooled$raa.rrpc_imp = rstudent(reg.rrpc_imp)

dim(db_pooled) # 1525 obs and 71 variables

################################################################################
# Partie V: Save db
################################################################################

saveRDS(db_pooled, "private_db/db_pooled_eGA.rds")

################################################################################
rm(list = ls())
################################################################################
