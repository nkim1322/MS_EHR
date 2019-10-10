# Fit three-part epoch-stratified model 
## This was done before we decided to remove all data before 2006
# Library -----------------------------------------------------------------
library(tidyverse)
library(glmnet)
library(pROC)

# Load --------------------------------------------------------------------
# CC_comb <- readRDS()
# CC_comb_val <- readRDS()


# Three part model --------------------------------------------------------

glmnet.list <- list()
# 2000-2005
tmp <- CC_comb %>% filter(StartDate < as.Date('2006-01-01'))
glmnet.list[['2000-2005']] <- cv.glmnet(as.matrix(tmp[,-c(1:9)]), unlist(tmp$CC),
                                        family = "binomial", type.measure = "deviance")
# 2005-2010
tmp <- CC_comb %>% filter(StartDate >= as.Date('2006-01-01') & StartDate < as.Date('2011-01-01'))
glmnet.list[['2006-2010']] <- cv.glmnet(as.matrix(tmp[,-c(1:9)]), unlist(tmp$CC),
                                        family = "binomial", type.measure = "deviance")
# 2010-2016
tmp <- CC_comb %>% filter(StartDate >= as.Date('2011-01-01'))
glmnet.list[['2011-2016']] <- cv.glmnet(as.matrix(tmp[,-c(1:9)]), unlist(tmp$CC),
                                        family = "binomial", type.measure = "deviance")

# Predict
df <- CC_comb_val[,1:6]
## Stratified model
df$Prob[CC_comb_val$StartDate < as.Date('2006-01-01')] <- 
  predict(glmnet.list[['2000-2005']], 
          newx = as.matrix(CC_comb_val %>% 
                             filter(StartDate < as.Date('2006-01-01')) %>% 
                             select(-c(1:9))),
          type = "response", s = glmnet.list[['2000-2005']]$lambda.1se)
# 2005-2010
df$Prob[CC_comb_val$StartDate >= as.Date('2006-01-01') &
          CC_comb_val$StartDate < as.Date('2011-01-01')] <- 
  predict(glmnet.list[['2006-2010']], 
          newx = as.matrix(CC_comb_val %>% 
                             filter(StartDate >= as.Date('2006-01-01') & 
                                      StartDate < as.Date('2011-01-01')) %>% 
                             select(-c(1:9))),
          type = "response", s = glmnet.list[['2006-2010']]$lambda.1se)
# 2010-2016
df$Prob[CC_comb_val$StartDate >= as.Date('2011-01-01')] <- 
  predict(glmnet.list[['2011-2016']], 
          newx = as.matrix(CC_comb_val %>% 
                             filter(StartDate >= as.Date('2011-01-01')) %>% 
                             select(-c(1:9))),
          type = "response", s = glmnet.list[['2011-2016']]$lambda.1se)

# Compute AUC

# plot(roc(as.numeric(CC_comb_val$CC), as.numeric(df$Prob)), legacy.axes = TRUE)
# Overall
roc(as.numeric(CC_comb_val$CC), as.numeric(df$Prob))

# ROC plot
plot(roc(as.numeric(CC_comb_val$CC), as.numeric(df$Prob)), 
     col = 2, main = "ROC, Time Period = 1 month", legacy.axes = TRUE)

## Epoch-specific
# 2000-2005
roc(as.numeric(CC_comb_val %>% filter(StartDate < as.Date('2006-01-01')) %>% select(CC) %>% unlist()), 
    as.numeric(df %>% filter(StartDate < as.Date('2006-01-01')) %>% select(Prob) %>% unlist())) 
# 2006-2010
roc(as.numeric(CC_comb_val %>% filter(StartDate >= as.Date('2006-01-01') & 
                                        StartDate < as.Date('2011-01-01')) %>% 
                 select(CC) %>% unlist()), 
    as.numeric(df %>% filter(StartDate >= as.Date('2006-01-01') & 
                               StartDate < as.Date('2011-01-01')) %>% select(Prob) %>% unlist()))
# 2011-2016
roc(as.numeric(CC_comb_val %>% filter(StartDate > as.Date('2011-01-01')) %>% select(CC) %>% unlist()), 
    as.numeric(df %>% filter(StartDate > as.Date('2011-01-01')) %>% select(Prob) %>% unlist())) 

# Get coefs
# 2000-2005
coefs <- coef(glmnet.list[['2000-2005']], s = "lambda.1se") %>% as.matrix()
coefs <- data.frame(term = row.names(coefs), beta = coefs[,1]); row.names(coefs) <- NULL
nonzero_coefs <- coefs %>% filter(beta != 0)
nonzero_coefs %>% nrow - 1
# 2006-2010
coefs <- coef(glmnet.list[['2006-2010']], s = "lambda.1se") %>% as.matrix()
coefs <- data.frame(term = row.names(coefs), beta = coefs[,1]); row.names(coefs) <- NULL
nonzero_coefs <- coefs %>% filter(beta != 0)
nonzero_coefs %>% nrow - 1
# 2011-2016
coefs <- coef(glmnet.list[['2011-2016']], s = "lambda.1se") %>% as.matrix()
coefs <- data.frame(term = row.names(coefs), beta = coefs[,1]); row.names(coefs) <- NULL
nonzero_coefs <- coefs %>% filter(beta != 0)
nonzero_coefs %>% nrow - 1



# Save --------------------------------------------------------------------
save(glmnet.list, file = 'models/glmnet.list.rda')