rm(list=ls())

setwd("~/Dropbox/Collaboration/MS CLIMB Paper2/Step3_BuildModel/Adaptive")
## source the R script
source("Library.R")
source("Fun.R")
## load the data
load("df_grouped.RData")
dat0=df_AFEP[,-1]
patient_num=1:dim(dat0)[1] 

# !!label should be in the 1st column
dat0=dat0[,c("Label", setdiff(colnames(dat0), "Label"))]

### Adaptive lasso, BIC selection
bet = Est.ALASSO.GLM.Approx(dat0, fam0 = "binomial")

pred = g.logit(pmin(as.matrix(cbind(rep(1, nrow(dat0[,-1])), dat0[,-1])) %*% bet,100))



set.seed(1023)
n = 2000
p = 150

X = mvrnorm(n, mu = rnorm(p, mean = 0, sd = 1), Sigma = diag(x = 1, p,p))
beta =rep(0, p)
beta[sample(1:p, 40, replace = FALSE)] = rnorm(40, mean = 0, sd = 1)
z = X %*% beta + rnorm(n)
p = 1/(1+exp(-z))
y = rbinom(n, 1, p)
dat0 = cbind(y,X)


start_time <- Sys.time()
bet = Est.ALASSO.GLM.Approx(dat0, fam0 = "binomial")
end_time <- Sys.time()
