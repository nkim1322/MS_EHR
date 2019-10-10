# Change stopping criterion of glmnet
# Library -----------------------------------------------------------------

library(tidyverse)
library(glmnet)
library(pROC)

# Load --------------------------------------------------------------------
CC_comb <- readRDS('intermediate_data/CC_comb_24_3_clin_clean.rds') # Training data
CC_comb_val <- readRDS('intermediate_data/CC_comb_val_24_3_clin_clean.rds') # Validation data


# Filter ------------------------------------------------------------------
# Using data only after 2006
CC_comb <- CC_comb %>% filter(StartDate >= as.Date('2006-01-01'))
CC_comb_val <- CC_comb_val %>% filter(StartDate >= as.Date('2006-01-01'))

# Using only race/ethnicity combo variable
CC_comb$RACE2 <- NULL; CC_comb$ETHNICITY2 <- NULL
CC_comb_val$RACE2 <- NULL; CC_comb_val$ETHNICITY2 <- NULL 

# Model 1 -----------------------------------------------------------------
# Set thresh = 1E-7/Number of encounters
start_time <- Sys.time()
glmnet.fit = cv.glmnet(as.matrix(CC_comb[,-c(1:9)]), unlist(CC_comb$CC),
                       family = "binomial", type.measure = "deviance",
                       thresh = 1e-07/nrow(CC_comb))
end_time <- Sys.time(); end_time - start_time
probs <- predict(glmnet.fit, 
                   newx = as.matrix(CC_comb_val[,-c(1:9)]),
                   type = "response", s = glmnet.fit$lambda.1se)

# Compute AUC
roc(as.numeric(CC_comb_val$CC), as.numeric(probs))
# Get nonzero coefficients
coefs <- coef(glmnet.fit, s = "lambda.1se") %>% as.matrix()
coefs <- data.frame(term = row.names(coefs), beta = coefs[,1]); row.names(coefs) <- NULL

nonzero_coefs <- coefs %>% filter(beta != 0)
nonzero_coefs %>% nrow - 1 # Subtract 1 for intercept

# Model 2 -----------------------------------------------------------------
# Set thres = 1E-7/Number of patients 

start_time <- Sys.time()
glmnet.fit = cv.glmnet(as.matrix(CC_comb[,-c(1:9)]), unlist(CC_comb$CC),
                       family = "binomial", type.measure = "deviance",
                       thresh = 1e-07/length(CC_comb$PatientID %>% unique))
end_time <- Sys.time(); end_time - start_time
probs <- predict(glmnet.fit, 
                 newx = as.matrix(CC_comb_val[,-c(1:9)]),
                 type = "response", s = glmnet.fit$lambda.1se)

# Compute AUC
roc(as.numeric(CC_comb_val$CC), as.numeric(probs))
# Get nonzero coefficients
coefs <- coef(glmnet.fit, s = "lambda.1se") %>% as.matrix()
coefs <- data.frame(term = row.names(coefs), beta = coefs[,1]); row.names(coefs) <- NULL

nonzero_coefs <- coefs %>% filter(beta != 0)
nonzero_coefs %>% nrow - 1 # Subtract 1 for intercept







