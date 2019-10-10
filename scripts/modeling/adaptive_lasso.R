# Adaptive LASSO models
# Library/functions -------------------------------------------------------
source("adaptive_bic/Library.R") 
source("adaptive_bic/Fun.R") 

library(pROC)
library(tidyverse)

# Load --------------------------------------------------------------------
CC_comb <- readRDS('modeling_data/train24_3.rds')
CC_comb_val <- readRDS('modeling_data/test24_3.rds')


# Format data -------------------------------------------------------------

# Use regular data
dat0_train <- cbind(CC_comb$CC, CC_comb[,-c(1:9)])
dat0_val <- cbind(CC_comb_val$CC, CC_comb_val[,-c(1:9)])
names(dat0_train)[1] <- 'CC'
names(dat0_val)[1] <- 'CC'

# Run ---------------------------------------------------------------------
### Adaptive lasso, BIC selection

get.test.result = function(pt_num) {
  bet = Est.ALASSO.GLM.Approx(dat0_train, fam0 = "binomial", adap = T, N_groups = pt_num)
  pred = g.logit(pmin(as.matrix(cbind(rep(1, nrow(dat0_val[,-1])), dat0_val[,-1])) %*% bet,100))
  tmp = roc(as.numeric(CC_comb_val$CC), as.numeric(pred))
  c(tmp$auc, sum(bet != 0))
}

result = sapply(seq(1000, 5000, by = 1000), function(i) get.test.result(i))
result2 = sapply(seq(1000, 4000, by = 50), function(i) get.test.result(i))
result3 = sapply(seq(1000, 4000, by = 100), function(i) get.test.result(i))

df = data.frame(pt_num = seq(1000, 5000, by = 1000),
                auc = result[1,],
                num_vars = result[2,])
df2 = data.frame(pt_num = seq(1000, 4000, by = 50),
                auc = result2[1,],
                num_vars = result2[2,])
df3 = data.frame(pt_num = seq(1000, 4000, by = 100),
                 auc = result3[1,],
                 num_vars = result[2,])


# Need to estimate in-patient correlation, rho, for effective df N*(1-rho)
# Try different values for N_group = sequence of 20,000, 10,000, 5,000, 1,000
# 1000 - 26657

# Get metrics -------------------------------------------------------------
# AUC
pred = g.logit(pmin(as.matrix(cbind(rep(1, nrow(dat0_val[,-1])), dat0_val[,-1])) %*% bet,100))
roc(as.numeric(dat0_val$CC), as.numeric(pred))

# Get number of selected features (subtract intercept)
(bet %>% length) - ((bet == 0) %>% sum) - 1 

