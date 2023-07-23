################################################################################
# Epigenetic gestational age estimations 
################################################################################
# Author: Orlane Rossini
# Date: 23/07/2023
################################################################################
# Main objectives:
#------------------------------------------------------------------------------#
# For the all created models
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
file_eclock_lmp <- "models/eclock_lmp.rds"
file_eclock_ultrasound <- "models/eclock_ultrasound.rds"

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
eclock_model_lmp <- readRDS(file_eclock_lmp)
eclock_model_ultrasound <- readRDS(file_eclock_ultrasound)

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
p = predict(eclock_model_lmp, M_values, s =  as.numeric(eclock_model_lmp$bestTune[2]))
tmp = data.frame(Sample_Id = names(p), ega_2_2 = p)

db = merge(db, tmp, by = "Sample_Id")
range(db$ega_lmp) # 33.87282 40.67211

# 2) Sensitivity analysis 
p = predict(eclock_model_ultrasound, as_M_values, s =  as.numeric(eclock_model_ultrasound$bestTune[2]))
tmp = data.frame(Sample_Id = names(p), ega_2_as = p)

db = merge(db, tmp, by = "Sample_Id")
range(db$ega_ultrasound) # 35.56877 42.12551

################################################################################
# Partie III : Residuals 
################################################################################
reg.lasso.re = lm(ega_ultrasound~po_gd, data = db)
db$raa.ultrasound = rstudent(reg.lasso.re)

reg.lasso.re = lm(ega_lmp~clim_gest_duration, data = db)
db$raa.lmp = rstudent(reg.lasso.re)

################################################################################
# Partie V : Save db
################################################################################
saveRDS(db , file = "private_db/db_pooled_eGa_end.rds")

################################################################################
rm(list = ls())
################################################################################
