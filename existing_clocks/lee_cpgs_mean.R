################################################################################
# Get beta values mean for Lee CpGs
################################################################################
# Author: Orlane Rossini
# Date: 30/08/2022
################################################################################
# Main objectives:
#------------------------------------------------------------------------------#
# For Lee CpGs we compute mean beta values 
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
# Partie I : Get mean values for Lee's CpGs
################################################################################

# For RPC we can take every subject 
# For CPC we take sain subject 
# For RRPC we take sain subject with GA >= 36

# All geo dataset used in Lee's clock development 
# GSE71678
# GSE75248
# GSE74738
# GSE44667

# 1) Get mean values for each geo public dataset used in Lee's clock
# I-a) 343 subjects
gse <- getGEO("GSE71678")
gse <- gse[[1]]

# for rpc we take all subject
mean_rpc_71678 <- as.data.frame(apply(exprs(gse),1,mean))

# for cpc we take uncomplicated pregnancy
# 304 subjects
pData(gse) = pData(gse)[which(pData(gse)[, "gestational_diabetes:ch1"]=="No"),]

exp_gene = as.data.frame(exprs(gse)) %>% select(rownames(pData(gse)))
mean_cpc_71678 <- as.data.frame(apply(exp_gene,1,mean))

# for rrpc we take uncomplicated pregnancy with GA >= 36weeks
# 293 subjects
pData(gse) = pData(gse)[which(pData(gse)[, "gestational age:ch1"]>=36),]

exp_gene = as.data.frame(exprs(gse)) %>% select(rownames(pData(gse)))
mean_rrpc_71678 <- as.data.frame(apply(exp_gene,1,mean))

# I-b) 334 among 335 
gse <- getGEO("GSE75248")
gse <- gse[[1]]

# for rpc we take all subject
mean_rpc_75248 <- as.data.frame(apply(exprs(gse),1,mean))

# for cpc we take uncomplicated pregnancy
# 293 subjects
pData(gse) = pData(gse)[which(pData(gse)[, "diabetes_gestational_pregnancy:ch1"]=="No"),]

exp_gene = as.data.frame(exprs(gse)) %>% select(rownames(pData(gse)))
mean_cpc_75248 <- as.data.frame(apply(exp_gene,1,mean))

# for rrpc we take uncomplicated pregnancy with GA >= 36weeks
# 293 subjects
pData(gse) = pData(gse)[which(pData(gse)[, "gestational age:ch1"]>=36),]

exp_gene = as.data.frame(exprs(gse)) %>% select(rownames(pData(gse)))
mean_rrpc_75248 <- as.data.frame(apply(exp_gene,1,mean))

# I-c) 44 among 46 
# gse <- getGEO("GSE71719")
# gse <- gse[[1]]
# 
# # for rpc we take all subject
# moyenne_rpc_71719 <- as.data.frame(apply(exprs(gse),1,mean))
# 
# # for cpc we take uncomplicated pregnancy
# moyenne_cpc_71719 <- as.data.frame(apply(exprs(gse),1,mean))
# 
# # for rrpc we take uncomplicated pregnancy with GA >= 36weeks
# moyenne_rrpc_71719 <- as.data.frame(apply(exprs(gse),1,mean))

# I-d)  86 among 88 
gse <- getGEO("GSE120250")
gse <- gse[[1]]

# for rpc we take all subject
mean_rpc_120250 <- as.data.frame(apply(exprs(gse),1,mean))

# for cpc we take uncomplicated pregnancy
mean_cpc_120250 <- as.data.frame(apply(exprs(gse),1,mean))

# for rrpc we take uncomplicated pregnancy with GA >= 36weeks
mean_rrpc_120250 <- as.data.frame(apply(exprs(gse),1,mean))

# II - Combine all dataset 
# II - a) for rpc clock
tmp_1 = merge(mean_rpc_71678, mean_rpc_75248, by="row.names")
names(tmp_1) <- c("CpG", "mean.x", "mean.y")

mean_rpc_120250$CpG = rownames(mean_rpc_120250)
names(mean_rpc_120250)<-c("mean", "CpG")

tmp_2 = merge(tmp_1, mean_rpc_120250, by="CpG")

global.mean = (tmp_2[,2]+tmp_2[,3]+tmp_2[,4])/3
mean_rpc_CpG = cbind(tmp_2[,1], global.mean)

# II - b) for cpc clock
tmp_1 = merge(mean_cpc_71678, mean_cpc_75248, by="row.names")
names(tmp_1) <- c("CpG", "mean.x", "mean.y")

mean_cpc_120250$CpG = rownames(mean_cpc_120250)
names(mean_cpc_120250)<-c("mean", "CpG")

tmp_2 = merge(tmp_1, mean_cpc_120250, by="CpG")

global.mean = (tmp_2[,2]+tmp_2[,3]+tmp_2[,4])/3
mean_cpc_CpG = cbind(tmp_2[,1], global.mean)

# II - c) for rrpc clock
tmp_1 = merge(mean_rrpc_71678, mean_rrpc_75248, by="row.names")
names(tmp_1) <- c("CpG", "mean.x", "mean.y")

mean_rrpc_120250$CpG = rownames(mean_rrpc_120250)
names(mean_rrpc_120250)<-c("mean", "CpG")

tmp_2 = merge(tmp_1, mean_rrpc_120250, by="CpG")

global.mean = (tmp_2[,2]+tmp_2[,3]+tmp_2[,4])/3
mean_rrpc_CpG = cbind(tmp_2[,1], global.mean)

################################################################################
# Partie III : Save all datasets 
################################################################################
saveRDS(mean_rpc_CpG, 
        file = "mean_value_lee_rpc_clock.rds")

saveRDS(mean_cpc_CpG, 
        file = "mean_value_lee_cpc_clock.rds")

saveRDS(mean_rrpc_CpG, 
        file = "mean_value_lee_rrpc_clock.rds")

################################################################################
rm(list = ls())
################################################################################
