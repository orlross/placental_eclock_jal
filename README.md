# The placental epigenetic clock development and association study
R code for the creation of the JAL placental epigenetic clock and results of association study. 

## Introduction 
This work was realised during a final intership at IAB, Grenoble, under the supervision of Johanna Lepeule, Aurélie Nakamura, Lucile Broseus and François Septier. 
The file main.R is used to run all programs in the right order 

## Available data 

The data folder contains some results. 
It only contains public results such as Mayne's and Lee's clock coefficients or our epigenetic clock coefficient. 
It does not contains any information about cohorts used in the paper. 

## Existing placental epigenetic clocks 

There exist four different placental epigenetic clocks : 
- The Mayne's clock developed by Mayne 
- The RPC clock develop by Lee 
- The CPC clock develop by Lee 
- The RRPC clock develop by Lee 

Compared to our data, all of those epigenetic clocks contains missing cpgs. 
One idea to deal with this trouble was to imputed missing cpgs with the mean of cpgs included in original dataset. 

## The construction of the epigenetic clock 

First all beta values were transformed into M-values. 
Then an EWAS was conducted to shortlisted CpGs associated with gestational age. 
Then a Lasso regression was performed. 
A sensitivity analysis (as) was conduction with gestational age estimated by ultrasound.  


