# Sparse PCA
# Library -----------------------------------------------------------------
library(sparsepca)
library(pROC)
library(tidyverse)
library(glmnet)
library(factoextra) # for eigenvalue plots
library(PMA)
library(impute)

# Load --------------------------------------------------------------------
CC_comb <- readRDS('modeling_data/train24_3.rds')
CC_comb_val <- readRDS('modeling_data/test24_3.rds')

# Combine data ------------------------------------------------------------
CC_all <- rbind(CC_comb, CC_comb_val)

# Preliminary sparse PCA and get loadings --------------------------------------------------
spca_fit <- sparsepca::spca(CC_comb[,-c(1:9)], k = 20, scale = T, alpha = 1e-02)

loadings <- spca_fit$loadings %>% data.frame()

# Add descriptions to loadings 
colnames(loadings) <- paste0('PC', 1:20)
rownames(loadings) <- names(CC_comb[-c(1:9)])
loadings$Desc <- NA

# Add descriptions
# Load keys
wkpath  = "raw_data/Box/Boston/MS CLIMB Data/"
ICDPheCode   = read.csv(paste0(wkpath,"EHR/MS_AllEncounters_ICD_Data_03282019.csv"),
                        stringsAsFactors = FALSE);ICDPheCode$phecode[ICDPheCode$concept_cd == "LPA268"] = "335_"
CUIdictAll = read.xlsx(paste0(wkpath,"EHR/AllCUI_Database.xlsx"),
                       sheet = 1); colnames(CUIdictAll) = c("ConceptCd","Desc")
firstup <- function(x) {substr(x, 1, 1) <- toupper(substr(x, 1, 1)) 
  return(x)}

# CPT Codes
tmp <- strtrim(rownames(loadings), 3) == 'CPT'
loadings$Desc[tmp] <- str_split(rownames(loadings)[tmp], '\\.', simplify = T)[,2] %>% 
  as.character()  %>% firstup()
# ICD PheCodes
tmp <- strtrim(rownames(loadings), 3) == 'Phe'
tmp2 <- str_split(rownames(loadings)[tmp], '\\.', simplify = T)[,2]
tmp3 <- character(length(tmp2))
for (i in 1:length(tmp2)) {
  tmp3[i] <- ICDPheCode[ICDPheCode$phecode == tmp2[i], 'phecode_description'] %>% unique()
}
loadings$Desc[tmp] <- tmp3
# CUI
# CUIs
tmp <- strtrim(rownames(loadings), 3) == 'CUI'
tmp2 <- str_split(rownames(loadings)[tmp], '\\.', simplify = T)[,2]
tmp3 <- character(length(tmp2))
for (i in 1:length(tmp2)) {
  tmp3[i] <- CUIdictAll[CUIdictAll$ConceptCd == tmp2[i], 'Desc'] %>% unique()
}
loadings$Desc[tmp] <- tmp3 
# Clinical vars
tmp <- !(strtrim(rownames(loadings), 3) == 'CPT' | strtrim(rownames(loadings), 3) == 'Phe' | strtrim(rownames(loadings), 3) == 'CUI')
loadings$Desc[tmp] <- rownames(loadings)[tmp]

# Truncate rownames for readability
row.names(loadings) <- strtrim(row.names(loadings), 25)


write.csv(loadings, 'model_output/spca_loadings_tuned.csv') # Loadings using spca

# write.csv(t(summary(spca_fit)), 'coefs/spca_summary.csv')
# Functions ---------------------------------------------------------------
## **_data is assumed to have the same column structure as CC_comb/CC_comb_val

fit_spca <- function(train_data, test_data, k = 20, alpha = 1e-02) {
  # Run sparse PCA only on training data, project test data onto orthog basis
  spca_fit <- sparsepca::spca(train_data[,-c(1:9)], k = k, scale = T, alpha = alpha)
  train_PC <- spca_fit$scores
  test_PC <- scale(test_data[,-c(1:9)], 
                   center = spca_fit$center, 
                   scale = spca_fit$scale) %*% spca_fit$transform
  # Returns train and test data as PCs in the form of a labeled list
  return(list(train_PC = train_PC, test_PC = test_PC))
}

fit_spc <- function(train_data, test_data, k = 20) {
  # Run PCA only on training data, project test data onto orthog basis
  # First scale the data
  scaled_train <- scale(train_data[,-c(1:9)])
  mean_train <- attr(scaled_train,"scaled:center"); scale_train <- attr(scaled_train,"scaled:scale")
  # Tune sparsity parameter
  set.seed(12345)
  cv_obj <- PMA::SPC.cv(as.matrix(scaled_train), center = FALSE)
  # Fit new sparse PCA using tuned parameter
  spc_obj <- PMA::SPC(as.matrix(scaled_train), sumabsv = cv_obj$bestsumabsv, K = k,
                      cnames = paste0('PC', 1:k))
  scaled_test <- scale(test_data[,-c(1:9)], center = mean_train, scale = scale_train)
  test_PC <- as.matrix(scaled_test) %*% as.matrix(spc_obj$v)
  
  train_PC <- spc_obj$u
  # Returns train and test data as PCs in the form of a labeled list
  return(list(train_PC = train_PC, test_PC = test_PC))
}

## 'pca_predict' function from pca.R script
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



# Manually tune sparsity parameter alpha using spca -----------------------
# Use original train-test split and k = 20
get.test.result = function(alpha) {
  spca_data <- fit_spca(CC_comb, CC_comb_val, k = 20, alpha = alpha)
  # Fit LASSO using deviance
  set.seed(12345)
  glmnet.fit <- cv.glmnet(as.matrix(spca_data[[1]]), unlist(CC_comb$CC),
                          family = "binomial", type.measure = "deviance", alpha = 1)
  preds <- predict(glmnet.fit, newx = as.matrix(spca_data[[2]]), 
                   type = "response", s = glmnet.fit$lambda.1se)
  return(roc(as.numeric(CC_comb_val$CC), as.numeric(preds))$auc)
}

result = sapply(10^(seq(-4, -1, by = 1)), function(i) get.test.result(i))
df_result = data.frame(alpha = 10^(seq(-4, -1, by = 1)), auc = result)


# Run spca ----------------------------------------------------------------
n_splits <- 6
set.seed(12345)
# Option 1
# k <- 10; pcs_vec <- c(2,3,10)
# Option 2
k <- 20; pcs_vec <- c(2,3,10,20)
alph <- 1e-02 # found from tuning above

# Run loop
all_results <- list()
for (j in 1:n_splits) {
  print(paste0('Split #', j))
  # Resplit and subset data
  nval = floor(length(unique(CC_all$PatientID))*0.3); val  = sample(unique(CC_all$PatientID),nval)
  CC_val = CC_all[CC_all$PatientID%in%val,]; CC_train = CC_all[!CC_all$PatientID%in%val,]
 
  # Fit Sparse PCA and project validation data
  spca_data <- fit_spca(CC_train, CC_val, k = k, alpha = alph)
  
  # Run models
  pca_results_list <- list()
  for (i in pcs_vec) {
    pca_results_list[[i]] <- pca_predict(CC_train, CC_val, spca_data[[1]], spca_data[[2]], i)
  }
  all_results[[j]] <- do.call(rbind, pca_results_list)
}
print('All done!')

# Average results
spca_auc <- (all_results[[1]]+all_results[[2]]+all_results[[3]]+all_results[[4]]+all_results[[5]]+all_results[[6]])/6

write_csv(spca_auc, 'model_output/spca_auc.csv')



# Run SPC.cv/SPC ----------------------------------------------------------
# Fit tuning parameter sumabsvs (sum of absolute values of elements of v)
# Scale first
scaled_train <- scale(CC_comb[,-c(1:9)])
mean_train <- attr(scaled_train,"scaled:center")
scale_train <- attr(scaled_train,"scaled:scale")
cv_obj <- PMA::SPC.cv(as.matrix(scaled_train), center = FALSE)

# Now fit
k <- 20
spc_obj <- PMA::SPC(as.matrix(scaled_train), sumabsv = cv_obj$bestsumabsv, K = k,
                    cnames = paste0('PC', 1:k))
# loadings <- spc_obj$v %>% data.frame()
# write.csv(loadings, 'model_output/spca_loadings_tuned_sumabsv.csv') # Loadings using SPC.cv and SPC
# Center test data and project to orthogonal axes
scaled_test <- scale(CC_comb_val[,-c(1:9)], center = mean_train, scale = scale_train)
test_PC <- as.matrix(scaled_test) %*% as.matrix(spc_obj$v)

# Fit model
glmnet.fit <- cv.glmnet(as.matrix(spc_obj$u), unlist(CC_comb$CC), family = "binomial", type.measure = "deviance", alpha = 1)
preds <- predict(glmnet.fit, newx = as.matrix(test_PC), type = "response", s = glmnet.fit$lambda.1se)
roc(as.numeric(CC_comb_val$CC), as.numeric(preds))$auc


k <- 20; pcs_vec <- c(2,3,10,20)
n_splits <- 6
# Run loop
all_results <- list()
for (j in 1:n_splits) {
  print(paste0('Split #', j))
  # Resplit and subset data
  nval = floor(length(unique(CC_all$PatientID))*0.3); val  = sample(unique(CC_all$PatientID),nval)
  CC_val = CC_all[CC_all$PatientID%in%val,]; CC_train = CC_all[!CC_all$PatientID%in%val,]
  
  # Fit Sparse PCA and project validation data
  spca_data <- fit_spc(CC_train, CC_val, k = k)
  
  # Run models
  pca_results_list <- list()
  for (i in pcs_vec) {
    pca_results_list[[i]] <- pca_predict(CC_train, CC_val, spca_data[[1]], spca_data[[2]], i)
  }
  all_results[[j]] <- do.call(rbind, pca_results_list)
}
print('All done!')


# Average results
spc_auc <- (all_results[[1]]+all_results[[2]]+all_results[[3]]+all_results[[4]]+all_results[[5]]+all_results[[6]])/6

write_csv(spc_auc, 'model_output/spc_cv_auc.csv')
  
  







