rm(list=ls())

## source the R script
# source("library_v2.R") 
# source("library_v3.R") 
source("adaptive/library_v2.R") 
source("adaptive/library_v3.R") 
## load the data
# load("df_grouped.RData")
load("adaptive/df_grouped.RData")
dat0=df_AFEP[,-1]
patient_num=1:dim(dat0)[1] 

# !!label should be in the 1st column
dat0=dat0[,c("Label", setdiff(colnames(dat0), "Label"))]

### Adaptive lasso, 10 fold cross validation, repeat 5 times
fit.alasso=convert(ROC.FUN.ALASSO.revise.mycv(data=dat0, yes.CV=T,yes.seed=T,rep=5,K=10,regularize=T))

# Cross-validated AUC
fit.alasso$ROC.cv.auc
# prediction
fit.alasso$Y.hat
# beta coefficients
fit.alasso$beta

