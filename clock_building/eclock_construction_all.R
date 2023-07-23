################################################################################
# Epigenetic clock construction - Build on all subjects
################################################################################
# Author: Orlane Rossini
# Date: 23/07/2023
################################################################################
# Main objectives:
#------------------------------------------------------------------------------#
# For the all dataset (on training and validation)
# 1 - We identify the best hyperparameters
# 2 - We create 1000 bootstrap models 
# 2 - Thanks to validation, we select the optimal number of CpGs
################################################################################
################################################################################
# Paths 
#------------------------------------------------------------------------------#
setwd("~/data/")

# all db_pooled
File <- "private_db/db_pooled_split.rds"

# CpGs files 
File_cpgs <- "private_db/M_ewas.rds"

################################################################################
# R packages
#------------------------------------------------------------------------------#

library(glmnet) # lasso regression 
library(dplyr)
library(caret) # create folds
library(car)

################################################################################
# Import data into R session
#------------------------------------------------------------------------------#
# imputed db 
db <- readRDS(File)

# CpGs file
M <- readRDS(File_cpgs)

################################################################################
# Partie I : Select only necessary cpgs 
################################################################################
M[1:3,1:3]

# I.1) We merge select CpGs selected in EWAS
dim(M) # # 15,240 cpgs

# I.2) We get the variability of each CpGs
cpg_iqr = apply(M, 1, IQR)
range(cpg_iqr) #  1.154182 2.002832

# I.3) Re-order 
indices <- match(db$Sample_Id, rownames(M))
M = M[indices,]

sum(db$Sample_Id==rownames(M)) # 1,517 obs

################################################################################
# Partie II : Build Lasso model
################################################################################
# Model Building :Lasso Regression
set.seed(99)
control = trainControl(method ="cv", number = 5)
Grid_la_reg = expand.grid(alpha = 1,
                          lambda = seq(0.001, 0.1, by = 0.0002))

# Training lasso regression model
lasso_model = train(x = M,
                    y = db$clim_gest_duration,
                    method = "glmnet",
                    trControl = control,
                    tuneGrid = Grid_la_reg
)
lasso_model # 0.0116

# mean validation score
mean(lasso_model$resample$RMSE) # 1.036939

# Plot
plot(lasso_model, main = "Lasso Regression")

# Structure

coef_Mat = predict(lasso_model$finalModel, type = "coef", 
                   mode = "fraction", s = as.numeric(lasso_model$bestTune[2]))
class(coef_Mat) # matrix
dim(coef_Mat) # 15,241 cpgs 1 var 

coef_Mat_02 = coef_Mat[which(coef_Mat[,1]>0 | coef_Mat[,1]<0),]
length(coef_Mat_02[-1]) # 985 CpGs
range(coef_Mat_02[-1]) # -0.4885118  0.4098394
mean(coef_Mat_02[-1]) # 0.001360334

################################################################################
# Partie III : Save final models
################################################################################

saveRDS(lasso_model , file = "models/eclock_all.rds")

################################################################################
rm(list = ls())
################################################################################
