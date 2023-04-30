# Master 2 - DSDM
# Code EWAS
################################################################################
# Author: Orlane LE QUELLENNEC
# Date: 08/09/2022
################################################################################
# Main objectives:
#------------------------------------------------------------------------------#
# Shortlisting CpGs to build our clock via a EWAS 
# Here we use robust linear regressions and M-values 
################################################################################
################################################################################
# Paths 
#------------------------------------------------------------------------------#
setwd("~/data/")

# all db_pooled
# File <- "Data V2/db_pooled_split_p2.rds"
File <- "private_db/db_pooled_split.rds"

# DNAm data
# File_M = "methylation_data/M_values.rds"
File_M = "private_db/M_values.rds"

################################################################################
# R packages
#------------------------------------------------------------------------------#
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("bacon")

library(dplyr)
library(bacon)
library(lme4)

################################################################################
# Import data into R session
#------------------------------------------------------------------------------#
# all db_pooled
db <- readRDS(File)

# M-values
M <- readRDS(File_M)

################################################################################
# Partie I : Mixted model for each XP adjusted EWAS
################################################################################
# /!\ We conduct our EWAS on the all dataset (Training + Validation + Testing)
# Because, if we don't, we won't have enough individuals included in PELAGIE /!\

dim(M) # 328,826 cpgs and 1,517 obs
# We transpose M matrix to get CpGs as variable and participant as observation 
M = t(M)

# We merge M_values with adjustement variables
M = as.data.frame(M)
M$Sample_Id = rownames(M)
tmp = merge(M, db[,c(1,28:38,44:46)], by = "Sample_Id")

# check order of subject
sum(tmp$Sample_Id==rownames(M)) # 1,517 obs

# The function to check the association of a cpg againt gestational duration 
r_val_lmer <- function(cpgs) {
  mod = try(
    { 
      reg = lme4::lmer(formula = as.formula(paste(cpgs, "clim_gest_duration + clim_mat_edu + 
      clim_mat_age + clim_parity_c + clim_delivery_mode + clim_labour + clim_sex + clim_season_conception + 
      clim_mat_pre_preg_bmi_2 + clim_smoke_2 + clim_area_2 + (1|Batch) + (1|Row) + (1|Column)", sep = "~")),
                       data = tmp)
    }
  )
  if( !("try-error" %in% class(mod)) ){
    return(coef(summary(reg))["clim_gest_duration","t value"])
  }else{
    return(coef(summary(reg))["clim_gest_duration","t value"])
  }
}

# parallelisation of the task 
cores = detectCores(all.tests = FALSE, logical = TRUE)-1 # 7 cores 

T1<-Sys.time()
result <- pbmcapply::pbmcmapply(r_val_lmer, colnames(M)[1:333196], mc.cores = cores)
T2<-Sys.time()
difftime(T2, T1) # 9 hours 

# store results into dataframe form 
results = data.frame(cpgs = colnames(M)[1:333196],
                     t.test = unlist(result))

dim(results) # 149,007 CpGs  

# get the adjusted p-values for multiple testing (FDR method) 
results = results%>%mutate(p.value = bacon::pval(bacon::bacon(results$t.test)))
results = results%>%mutate(p.value.adj =  p.adjust(results$p.value, method = "BH"))
range(results$p.value) #  1.580086e-19 9.999934e-01
range(results$p.value.adj) # 5.264783e-14 9.999934e-01
length(which(results$p.value.adj<0.05)) # 15,240 cpgs

# list of significant cpgs with a FDR values lower than 0.05
cpgs_ewas = results$cpgs[which(results$p.value.adj<0.05)]
length(cpgs_ewas) # 15,240 cpgs

# keep only significant FDR values (lower than 0.05)
M_ewas = M[,cpgs_ewas]
dim(M_ewas)# 15,240 cpgs

################################################################################
# Partie II : Save results
################################################################################

saveRDS(M_ewas , "private_db/M_ewas.rds")
saveRDS(cpgs_ewas , "private_db/ewas_list.rds")


################################################################################
# Partie III : AS - gestational age measured by ultrasound 
################################################################################
# /!\ We conduct our EWAS on the all dataset (Training + Validation + Testing)
# Because, if we don't, we won't have enough individuals included in PELAGIE /!\

# The function to check the association of a cpg againt gestational duration 
as_r_val_lmer <- function(cpgs) {
  mod = try(
    { 
      reg = lme4::lmer(formula = as.formula(paste(cpgs, "po_gd + clim_mat_edu + 
      clim_mat_age + clim_parity_c + clim_delivery_mode + clim_labour + clim_sex + clim_season_conception + 
      clim_mat_pre_preg_bmi_2 + clim_smoke_2 + clim_area_2 + (1|Batch) + (1|Row) + (1|Column)", sep = "~")),
                       data = tmp)
    }
  )
  if( !("try-error" %in% class(mod)) ){
    return(coef(summary(reg))["po_gd","t value"])
  }else{
    return(coef(summary(reg))["po_gd","t value"])
  }
}

# parallelisation of the task 
cores = detectCores(all.tests = FALSE, logical = TRUE)-1 # 7 cores 

T1<-Sys.time()
result <- pbmcapply::pbmcmapply(as_r_val_lmer, colnames(M)[1:333196], mc.cores = cores)
T2<-Sys.time()
difftime(T2, T1) # 9 hours 

# store results into dataframe form 
results = data.frame(cpgs = colnames(M)[1:333196],
                     t.test = unlist(result))

dim(results) # 149,007 CpGs  

# get the adjusted p-values for multiple testing (FDR method) 
results = results%>%mutate(p.value = bacon::pval(bacon::bacon(results$t.test)))
results = results%>%mutate(p.value.adj =  p.adjust(results$p.value, method = "BH"))
range(results$p.value) #  1.580086e-19 9.999934e-01
range(results$p.value.adj) # 5.264783e-14 9.999934e-01
length(which(results$p.value.adj<0.05)) # 15,240 cpgs

# list of significant cpgs with a FDR values lower than 0.05
as_cpgs_ewas = results$cpgs[which(results$p.value.adj<0.05)]
length(cpgs_ewas) # 15,240 cpgs

# keep only significant FDR values (lower than 0.05)
as_M_ewas = M[,cpgs_ewas]
dim(M_ewas)# 15,240 cpgs

################################################################################
# Partie IV : Save results
################################################################################

saveRDS(as_M_ewas , "private_db/as_M_ewas.rds")
saveRDS(as_cpgs_ewas , "private_db/as_ewas_list.rds")

################################################################################
# Partie V : EWAS comparison 
################################################################################

length(cpgs_ewas) 

length(intersect(cpgs_ewas, as_cpgs_ewas))

################################################################################
rm(list = ls())
################################################################################
