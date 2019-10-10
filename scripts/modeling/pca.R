# PCA
# Library -----------------------------------------------------------------
library(pROC)
library(tidyverse)
library(glmnet)
library(factoextra) # for eigenvalue plots

# Load --------------------------------------------------------------------
CC_comb <- readRDS('modeling_data/train24_3.rds')
CC_comb_val <- readRDS('modeling_data/test24_3.rds')

# Combine data ------------------------------------------------------------
CC_all <- rbind(CC_comb, CC_comb_val)

# Preliminary PCA run -----------------------------------------------------

pca_fit <- prcomp(as.matrix(CC_comb[,-c(1:9)], scale. = T))
# Loadings
loadings <- pca_fit$rotation %>% as.data.frame()

write.csv(loadings, 'model_output/pca_loadings.csv')

# Functions ---------------------------------------------------------------
## **_data is assumed to have the same column structure as CC_comb/CC_comb_val

# PCA
fit_pca <- function(train_data, test_data) {
  # Run PCA only on training data, project test data onto orthog basis
  pca_fit <- prcomp(as.matrix(train_data[, -c(1:9)]), scale. = TRUE)
  train_PC <- pca_fit$x
  test_PC <- predict(pca_fit, newdata = as.matrix(test_data[, -c(1:9)]))
  
  # Returns train and test data as PCs in the form of a labeled list
  return(list(train_PC = train_PC, test_PC = test_PC))
}


pca_predict <- function(train_data, test_data, train_PC, test_PC, n_PCs = 2) {
  # Note: n_PCs > 1 (cv.glmnet will not work for 1 predictor)
  # Select only first n_PCs
  X_pc_train <- train_PC[,c(1:n_PCs)]; X_pc_test <- test_PC[,c(1:n_PCs)]
  
  # Ordinary logistic regression
  tmp_df <- cbind(unlist(train_data$CC), data.frame(X_pc_train)); names(tmp_df)[1] <- 'CC'
  glm.fit <- glm(CC ~ ., tmp_df, family = binomial)
  preds <- predict(glm.fit, newdata = data.frame(X_pc_test), type = 'response')
  auc_logreg <- roc(as.numeric(test_data$CC), as.numeric(preds))$auc
  
  # LASSO, type.measure = 'deviance'
  set.seed(12345)
  glmnet.fit <- cv.glmnet(as.matrix(X_pc_train), unlist(train_data$CC),
                         family = "binomial", type.measure = "deviance", alpha = 1)
  preds <- predict(glmnet.fit, newx = as.matrix(X_pc_test), 
                   type = "response", s = glmnet.fit$lambda.1se)
  auc_lasso_dev <- roc(as.numeric(test_data$CC), as.numeric(preds))$auc
  # LASSO, type.measure = 'auc'
  set.seed(12345)
  glmnet.fit <- cv.glmnet(as.matrix(X_pc_train), unlist(train_data$CC),
                          family = "binomial", type.measure = "auc", alpha = 1)
  preds <- predict(glmnet.fit, newx = as.matrix(X_pc_test), 
                   type = "response", s = glmnet.fit$lambda.1se)
  auc_lasso_auc <- roc(as.numeric(test_data$CC), as.numeric(preds))$auc
  # Ridge, type.measure = 'deviance'
  set.seed(12345)
  glmnet.fit <- cv.glmnet(as.matrix(X_pc_train), unlist(train_data$CC),
                          family = "binomial", type.measure = "deviance", alpha = 0)
  preds <- predict(glmnet.fit, newx = as.matrix(X_pc_test), 
                   type = "response", s = glmnet.fit$lambda.1se)
  auc_ridge_dev <- roc(as.numeric(test_data$CC), as.numeric(preds))$auc
  # Ridge, type.measure = 'auc'
  set.seed(12345)
  glmnet.fit <- cv.glmnet(as.matrix(X_pc_train), unlist(train_data$CC),
                          family = "binomial", type.measure = "auc", alpha = 0)
  preds <- predict(glmnet.fit, newx = as.matrix(X_pc_test), 
                   type = "response", s = glmnet.fit$lambda.1se)
  auc_ridge_auc <- roc(as.numeric(test_data$CC), as.numeric(preds))$auc
  
  return(data.frame(num_pc = n_PCs,
                    lasso_dev = auc_lasso_dev, lasso_auc = auc_lasso_auc,
                    ridge_dev = auc_ridge_dev, ridge_auc = auc_ridge_auc,
                    ord_logreg = auc_logreg))
  
}


# Run ---------------------------------------------------------------------

# New splits
n_splits <- 6
set.seed(12345)

# Run loop
all_results <- list()
for (j in 1:n_splits) {
  print(paste0('Split #', j))
  
  # Resplit and subset data
  nval = floor(length(unique(CC_all$PatientID))*0.3); val  = sample(unique(CC_all$PatientID),nval)
  CC_val = CC_all[CC_all$PatientID%in%val,]; CC_train = CC_all[!CC_all$PatientID%in%val,]
  pca_data <- fit_pca(CC_train, CC_val)
  
  pca_results_list <- list()
  pcs_vec <- c(2,3,10,20,30,40,50,60,70,80,120)
  for (i in pcs_vec) {
    pca_results_list[[i]] <- pca_predict(CC_train, CC_val, pca_data[[1]], pca_data[[2]], i)
  }
  
  all_results[[j]] <- do.call(rbind, pca_results_list)
}
print('All done!')

# Average results
pca_results <- (all_results[[1]]+all_results[[2]]+all_results[[3]]+all_results[[4]]+all_results[[5]]+all_results[[6]])/6

# Write results to csv ----------------------------------------------------
write_csv(pca_results, 'model_output/pca_auc.csv')



# Plot eigenvalues --------------------------------------------------------
# See: http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/118-principal-component-analysis-in-r-prcomp-vs-princomp/#access-to-the-pca-results
# eigs <- get_eig(ms_pca)
# 
# fviz_eig(ms_pca, choice = 'eigenvalue', geom = 'line', ncp = 30,
#          addlabels = T, main = 'Eigenvalues of PCs')
# plot(ms_pca)

# Other fns: fviz_pca_ind, fviz_pca_biplot, fviz_pca_var

