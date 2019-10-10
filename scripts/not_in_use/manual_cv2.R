# New analysis 2019-09-21 
# Library -----------------------------------------------------------------
library(pROC)
library(tidyverse)
library(glmnet)

# Load --------------------------------------------------------------------
CC_comb <- readRDS('intermediate_data/CC_comb_24_3_clin.rds')
CC_comb_val <- readRDS('intermediate_data/CC_comb_val_24_3_clin.rds')

# Apply filters -----------------------------------------------------------
CC_comb <- CC_comb %>% filter(StartDate >= as.Date('2006-01-01'))
CC_comb_val <- CC_comb_val %>% filter(StartDate >= as.Date('2006-01-01'))
CC_comb$RACE2 <- NULL; CC_comb_val$RACE2 <- NULL
CC_comb$ETHNICITY2 <- NULL; CC_comb_val$ETHNICITY2 <- NULL

# cv.glmnet(type = ‘auc’) -------------------------------------------------
glmnet.fit = cv.glmnet(as.matrix(CC_comb[,-c(1:9)]), unlist(CC_comb$CC),
                       family = "binomial", type.measure = "auc")
# Predict
preds <- predict(glmnet.fit, 
                 newx = as.matrix(CC_comb_val[,-c(1:9)]),
                 type = "response", s = glmnet.fit$lambda.1se)
# Compute AUC
roc(as.numeric(CC_comb_val$CC), as.numeric(preds))
# Get nonzero coefficients
coefs <- coef(glmnet.fit, s = "lambda.1se") %>% as.matrix()
coefs <- data.frame(term = row.names(coefs), beta = coefs[,1]); row.names(coefs) <- NULL
nonzero_coefs <- coefs %>% filter(beta != 0)
nonzero_coefs %>% nrow - 1 # Subtract 1 for intercept



# Manual 10-fold CV -------------------------------------------------------

lambdas <- c(0.0001, 0.001, 0.0015, 0.002, 0.003, 0.004, 0.005, 0.01, 0.05, 0.1, 0.5, 1) %>% 
  as.list() %>% set_names()

set.seed(12345)
# Split data into N=10 groups
nreps <- 3
N = 10

mean_auc_list <- map(lambdas, function(lambda) {
  
  for (j in 1:nreps) {
    # Split data
    split_list <- split(CC_comb, sample(1:N, nrow(CC_comb), replace=T))
    aucs <- rep(NA, N)
    num_coefs <- rep(NA, N)
    for (i in 1:N) {
      test_set <- split_list[[i]]
      tmp_list <- split_list
      tmp_list[[1]] <- NULL
      train_set <- do.call(rbind, tmp_list)
      glmnet.fit <- glmnet(as.matrix(train_set[, -c(1:9)]), unlist(train_set$CC),
                           family = "binomial", lambda = lambda, alpha = 1)
      preds <- predict(glmnet.fit, newx = as.matrix(test_set[,-c(1:9)]), type = "response")
      # Store
      aucs[i] <- roc(as.numeric(test_set$CC), as.numeric(preds))$auc
    }
  }
  
  return(mean(aucs))
})



# accuracies <- map(lambdas, function(lambda) {
#   glmnet.fit <- cv.glmnet(as.matrix(CC_comb[,-c(1:9)]), unlist(CC_comb$CC), 
#                        family = "binomial", lambda = lambda, alpha = 1)
#   preds <- predict(glmnet.fit, 
#                    newx = as.matrix(CC_comb_val[,-c(1:9)]),
#                    type = "response")
#   # Get nonzero coefficients
#   coefs <- coef(glmnet.fit, s = "lambda.1se") %>% as.matrix()
#   coefs <- data.frame(term = row.names(coefs), beta = coefs[,1]); row.names(coefs) <- NULL
#   
#   nonzero_coefs <- coefs %>% filter(beta != 0)
#   
#   # print(roc(as.numeric(CC_comb_val$CC), as.numeric(preds))$auc)
#   # print(nonzero_coefs %>% nrow - 1)
#   return(nonzero_coefs)
# })



