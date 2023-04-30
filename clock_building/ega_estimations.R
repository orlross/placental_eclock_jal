################################################################################
# Master 2 - DSDM
# Code 09 - eGA estimations
################################################################################
# Author: Orlane LE QUELLENNEC
# Date: 13/09/2022
################################################################################
# Main objectives:
#------------------------------------------------------------------------------#
# For the all models created 
# we compute epigenetic gestational age estimations 
# on all database
################################################################################
################################################################################
# Paths 
#------------------------------------------------------------------------------#
setwd("~/data/")

# all db_pooled
File <- "private_db/db_pooled_split.rds"

# CpGs files 
File_cpgs <- "private_db/M_ewas.rds"
File_cpgs_as <- "private_db/as_M_ewas.rds"

# Models
file_eclock <- "models/eclock.rds"
file_eclock_as <- "models/as_eclock.rds"

################################################################################
# R packages
#------------------------------------------------------------------------------#

library(glmnet) # regression elastic net
library(glmnetUtils) # cva.glmnet
library(dplyr) # select
library(caret) # create folds
library(naniar)# missing values
library(car)

################################################################################
# Import data into R session
#------------------------------------------------------------------------------#
# imputed db 
db <- readRDS(File)

# CpGs file
M_values <- readRDS(File_cpgs)
as_M_values <- readRDS(File_cpgs_as)

# Models
eclock_model <- readRDS(file_eclock)
as_eclock_model <- readRDS(file_eclock_as)

################################################################################
# Partie I : Arrange the dataset
################################################################################
# M_values
indices <- match(db$Sample_Id, rownames(M_values))
M_values = M_values[indices,]

sum(db$Sample_Id==rownames(M_values)) # 1,517

M_values = as.data.frame(M_values)
M_values$Sample_Id = rownames(M_values)

# as_M_values 

indices <- match(db$Sample_Id, rownames(as_M_values))
as_M_values = as_M_values[indices,]

sum(db$Sample_Id==rownames(as_M_values)) # 1,517

as_M_values = as.data.frame(as_M_values)
as_M_values$Sample_Id = rownames(as_M_values)

################################################################################
# Partie II : Epigenetic gestational age estimation
################################################################################
# 1) Regular model
p = predict(eclock_model, M_values, s =  as.numeric(eclock_model$bestTune[2]))
tmp = data.frame(Sample_Id = names(p), ega_2_2 = p)

db = merge(db, tmp, by = "Sample_Id")
range(db$ega_2_2) # 33.87282 40.67211

# 2) Sensitivity analysis 
p = predict(as_eclock_model, as_M_values, s =  as.numeric(as_eclock_model$bestTune[2]))
tmp = data.frame(Sample_Id = names(p), ega_2_as = p)

db = merge(db, tmp, by = "Sample_Id")
range(db$ega_2_as) # 35.56877 42.12551

################################################################################
# Partie III : Residuals 
################################################################################
reg.lasso.re = lm(ega_2_as~po_gd, data = db)
db$raa.2_as = rstudent(reg.lasso.re)

reg.lasso.re = lm(ega_2_2~clim_gest_duration, data = db)
db$raa.2_2 = rstudent(reg.lasso.re)

################################################################################
# Partie V : Save db
################################################################################
saveRDS(db , file = "private_db/db_pooled_eGa_end.rds")

################################################################################
rm(list = ls())
################################################################################
