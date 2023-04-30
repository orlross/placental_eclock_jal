################################################################################
# Master 2 - DSDM
# Code 07.01 - eClock construction rlm.m - Lasso
################################################################################
# Author: Orlane LE QUELLENNEC
# Date: 13/09/2022
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
File_cpgs <- "private_db/as_M_ewas.rds"

################################################################################
# R packages
#------------------------------------------------------------------------------#

library(glmnet) # regression elastic net
library(dplyr) # select
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
dim(M) # # 18,071 cpgs

# I.2) We get the variability of each CpGs
cpg_iqr = apply(M, 1, IQR)
range(cpg_iqr) # 1.183394 2.006298

# I.3) Split training / validation / testing
#/!\ Train and Validation are combined bcs in Lasso we don't need validation set
id_train = db$Sample_Id[which(db$set2%in%c("Training"))]; # , "Validation"
length(id_train) # 456 ind

M_train = M[id_train,]

db_train = db[which(db$set2%in%c("Training")),] # , "Validation"

# I.4) Re-order 
indices <- match(db_train$Sample_Id, rownames(M_train))
M_train = M_train[indices,]

sum(db_train$Sample_Id==rownames(M_train)) # 456 obs

################################################################################
# Partie II : Build Lasso model
################################################################################
# Model Building :Lasso Regression
set.seed(99)
control = trainControl(method ="cv", number = 5)
Grid_la_reg = expand.grid(alpha = 1,
                          lambda = seq(0.001, 0.1, by = 0.0002))

# Training lasso regression model
lasso_model = train(x = M_train,
                    y = db_train$po_gd,
                    method = "glmnet",
                    trControl = control,
                    tuneGrid = Grid_la_reg
)
lasso_model # 0.039

# mean validation score
mean(lasso_model$resample$RMSE) # 1.148033

# Plot
plot(lasso_model, main = "Lasso Regression")

# Structure

coef_Mat = predict(lasso_model$finalModel, type = "coef", 
                   mode = "fraction", s = as.numeric(lasso_model$bestTune[2]))
class(coef_Mat) # matrix
dim(coef_Mat) # 4,617 cpgs 1 var 

coef_Mat_02 = coef_Mat[which(coef_Mat[,1]>0 | coef_Mat[,1]<0),]
length(coef_Mat_02[-1]) # 286 CpGs
range(coef_Mat_02[-1]) # -0.5544920  0.5972249
mean(coef_Mat_02[-1]) # 0.004371185

################################################################################
# Partie III : Save final models
################################################################################

saveRDS(lasso_model , file = "models/as_eclock.rds")

################################################################################
rm(list = ls())
################################################################################

