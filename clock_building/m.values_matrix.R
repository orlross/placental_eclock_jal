################################################################################
# M-values matrix
################################################################################
# Author: Orlane Rossini
# Date: 08/09/2022
################################################################################
# Main objectives:
#------------------------------------------------------------------------------#
# For DNAm data, transform beta values into M values
# (M-values are a logit transformation of Beta values)
################################################################################
################################################################################
# Paths 
#------------------------------------------------------------------------------#
setwd("~/data/")

# all db_pooled splitted into training and testing sets
File <- "~/private_db/db_pooled_split.rds"

# DNAm data
File_cpgs = "methylation_data/ssnoob_processed.rds"

# annotation for sexual chr
File_annot <- "methylation_data/annotation.csv"

################################################################################
# R packages
#------------------------------------------------------------------------------#

# duplicates 
library(dplyr) 
# stratified function
library(splitstackshape) 

################################################################################
# Import data into R session
#------------------------------------------------------------------------------#
# imputed db 
db <- readRDS(File)

# cpgs file
cpgs <- readRDS(File_cpgs)

# annotation for sexual chr
annotation <- read.csv(File_annot, header = T, sep = ";")

################################################################################
# Partie I : Keep only necessary subject
################################################################################
dim(cpgs) #333 229 CpGs and 1539 observations 
cpgs = cpgs[,db$Sample_Id]
dim(cpgs) #333 229 CpGs and 1517 observations 

################################################################################
# Partie II : Remove sexual chr
################################################################################
`%notin%` <- Negate(`%in%`)
cpgs_sexuals = annotation$Name[which(annotation$Chromosome_36%in%c("X","Y"))]
length(cpgs_sexuals) # 11 646 CpGs

# cpgs[1:3,1:3]

cpgs = cpgs[which(rownames(cpgs)%notin%cpgs_sexuals),]
dim(cpgs) # 1517 obs 333 229 cpgs
# sexual Chr were already remove by Lucile B. 

################################################################################################################
# Partie III : Remove non-variable CpGs (IQR > 0.05)
################################################################################
# compute IQR for each CpGs 
cpg_iqr = apply(cpgs, 1, IQR)
range(cpg_iqr) # 0.001727618 - 0.847155231

cpg_variability = names(cpg_iqr)[which(cpg_iqr>=0.05)]
length(cpg_variability) # 172 743 cpgs

cpgs = cpgs[which(rownames(cpgs)%in%cpg_variability),]
dim(cpgs) # 172 743 cpgs and 1 517 obs

################################################################################################################
# Partie IV : logit transformation 
################################################################################
M = log(cpgs/(1-cpgs))
range(M) # -5.924878  5.808687

################################################################################################################
# Partie V : Save M matrix
################################################################################
saveRDS(M , "private_db/M_values.rds")

################################################################################
rm(list = ls())
################################################################################
