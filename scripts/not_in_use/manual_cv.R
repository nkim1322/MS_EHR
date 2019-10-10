# Sanity check for CV.glmnet
# Library -----------------------------------------------------------------
library(glmnet)
library(tidyverse)
library(pROC)
library(bsts)


# Load --------------------------------------------------------------------
CC_comb <- readRDS('intermediate_data/CC_comb_24_3_clin_clean.rds') # Training data
CC_comb_val <- readRDS('intermediate_data/CC_comb_val_24_3_clin_clean.rds') # Validation data


# Filter ------------------------------------------------------------------
CC_comb <- CC_comb %>% filter(StartDate >= as.Date('2006-01-01'))
CC_comb_val <- CC_comb_val %>% filter(StartDate >= as.Date('2006-01-01'))
CC_comb$RACE2 <- NULL; CC_comb$ETHNICITY2 <- NULL
CC_comb_val$RACE2 <- NULL; CC_comb_val$ETHNICITY2 <- NULL 


# Create double variables -------------------------------------------------
# PERIOD: indicator for whether date is after 1/1/2010
CC_comb$PERIOD <- ifelse(CC_comb$StartDate > as.Date('2010-01-01'), 2, 1)
CC_comb_val$PERIOD <- ifelse(CC_comb_val$StartDate > as.Date('2010-01-01'), 2, 1)

## Create interaction terms between PERIOD and code
# Variables with which to multiply by PERIOD
# Training data
vars <- names(CC_comb)[-c(1:9, 206)]
attach(CC_comb)
for (var in vars) {
  col_name <- paste0('PERIOD:', var)
  CC_comb[, col_name] <- PERIOD*get(var)
}
# Validation data
vars <- names(CC_comb_val)[-c(1:9, 206)]
attach(CC_comb_val)
for (var in vars) {
  col_name <- paste0('PERIOD:', var)
  CC_comb_val[, col_name] <- PERIOD*get(var)
}


# Model and get AUC -------------------------------------------------------

glmnet.fit = cv.glmnet(as.matrix(CC_comb[,-c(1:9)]), unlist(CC_comb$CC),
                       family = "binomial", type.measure = "deviance")

df <- CC_comb_val[,1:6]
df$Prob <- predict(glmnet.fit, 
                   newx = as.matrix(CC_comb_val[,-c(1:9)]),
                   type = "response", s = glmnet.fit$lambda.1se)
# Compute AUC
roc(as.numeric(CC_comb_val$CC), as.numeric(df$Prob))

# Get nonzero coefficients
coefs <- coef(glmnet.fit, s = "lambda.1se") %>% as.matrix()
coefs <- data.frame(term = row.names(coefs), beta = coefs[,1]); row.names(coefs) <- NULL

nonzero_coefs <- coefs %>% filter(beta != 0)
nonzero_coefs %>% nrow - 1 # Subtract 1 for intercept

# Set lambda --------------------------------------------------------------

lambdas <- c(0.0001, 0.001, 0.005, 0.006, 0.007, 0.008, 0.009, 0.01, 0.02, 0.03, 0.04, 0.05, 0.1, 0.5, 1, 2, 4) %>% 
  as.list() %>% set_names()

# Test for predictive accuracy
accuracies <- map(lambdas, function(lambda) {
  glmnet.fit <- glmnet(as.matrix(CC_comb[,-c(1:9)]), unlist(CC_comb$CC), 
                       family = "binomial", lambda = lambda, alpha = 1)
  preds <- predict(glmnet.fit, 
                   newx = as.matrix(CC_comb_val[,-c(1:9)]),
                   type = "response", s = lambda)
  
  # Get nonzero coefficients
  coefs <- coef(glmnet.fit, s = "lambda.1se") %>% as.matrix()
  coefs <- data.frame(term = row.names(coefs), beta = coefs[,1]); row.names(coefs) <- NULL
  
  nonzero_coefs <- coefs %>% filter(beta != 0)
  
  # print(roc(as.numeric(CC_comb_val$CC), as.numeric(preds))$auc)
  # print(nonzero_coefs %>% nrow - 1)
  return(nonzero_coefs)
})


# Use averaged deviance to select lambda ----------------------------------



