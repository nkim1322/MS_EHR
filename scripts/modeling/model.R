# Fit Lasso model and examine coefficients
# Library -----------------------------------------------------------------
library(tidyverse)
library(glmnet)
library(pROC)

# Load --------------------------------------------------------------------
# Use the tw = 24, tp = 3 data
CC_comb <- readRDS('modeling_data/train24_3.rds')
CC_comb_val <- readRDS('modeling_data/test24_3.rds')

write_out = FALSE # write out selected coefficients and code descriptions

# Model -------------------------------------------------------------------
set.seed(12345)
glmnet.fit = cv.glmnet(as.matrix(CC_comb[,-c(1:9)]), unlist(CC_comb$CC),
                       family = "binomial", type.measure = "deviance")
preds <- predict(glmnet.fit, newx = as.matrix(CC_comb_val[,-c(1:9)]),
                   type = "response", s = glmnet.fit$lambda.1se)

# Compute AUC
roc(as.numeric(CC_comb_val$CC), as.numeric(preds))


# Write coefs to csv with code descriptions -------------------------------
if (write_out) {
  
  # Examine coefficients ----------------------------------------------------
  coefs <- coef(glmnet.fit, s = "lambda.1se") %>% as.matrix()
  coefs <- data.frame(term = row.names(coefs), beta = coefs[,1])
  row.names(coefs) <- NULL
  
  nonzero_coefs <- coefs %>% filter(beta != 0)
  nonzero_coefs %>% nrow - 1 
  
  # Add code descriptions ---------------------------------------------------
  ## Load references
  wkpath  = "raw_data/Box/Boston/MS CLIMB Data/"
  ICDPheCode   = read.csv(paste0(wkpath,"EHR/MS_AllEncounters_ICD_Data_03282019.csv"), stringsAsFactors = FALSE)
  ICDPheCode$phecode[ICDPheCode$concept_cd == "LPA268"] = "335_"
  CUIdictAll = read.xlsx(paste0(wkpath,"EHR/AllCUI_Database.xlsx"),sheet = 1)
  colnames(CUIdictAll) = c("ConceptCd","Desc")
  
  nonzero_coefs$desc <- NA
  ## CUIs
  tmp <- str_detect(nonzero_coefs$term, 'CUI')
  tmp2 <- str_split(nonzero_coefs$term[tmp], '\\.', simplify = T)[,2]
  tmp3 <- character(length(tmp2))
  for (i in 1:length(tmp2)) {
    tmp3[i] <- unique(CUIdictAll[CUIdictAll$ConceptCd == tmp2[i], 'Desc'])[1]
  }
  nonzero_coefs$desc[tmp] <- tmp3 
  ## PheCodes
  tmp <- str_detect(nonzero_coefs$term, 'PheCode')
  tmp2 <- str_split(nonzero_coefs$term[tmp], '\\.', simplify = T)[,2]
  tmp3 <- character(length(tmp2))
  for (i in 1:length(tmp2)) {
    tmp3[i] <- ICDPheCode[ICDPheCode$phecode == tmp2[i], 'phecode_description'] %>% unique()
  }
  nonzero_coefs$desc[tmp] <- tmp3 
  
  ## No need for CPT groups, description is self-evident
  
  # Write to csv
  write_csv(nonzero_coefs, 'model_output/cv_glmnet_lasso.csv')
}



rm(write_out)
# Ordinary logistic regression --------------------------------------------
# glm.fit <- glm(CC_comb$CC ~ ., CC_comb[, -c(1:9)], family = binomial)
# preds1 <- predict(glm.fit, newdata = data.frame(CC_comb_val[,-c(1:9)]), type = 'response')
# roc(as.numeric(CC_comb_val$CC), as.numeric(preds1))


