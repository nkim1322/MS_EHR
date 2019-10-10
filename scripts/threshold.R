# Threshold Analysis

# Library -----------------------------------------------------------------

library(tidyverse)

# Load --------------------------------------------------------------------
df <- readRDS('intermediate_data/3partclin_probs.rds') # Predicted probabilities of validation data

# try using confusion matrix function?

# Set cutoffs -------------------------------------------------------------

cutoffs <- seq(0,1, by = 0.001)

# Calculate statistics ----------------------------------------------------
stats_df <- data.frame(cutoffs = cutoffs)
cutoffs <- as.list(cutoffs) %>% set_names()
# tpr = tp / (tp + fn) 
stats_df$TPR <- map(cutoffs, function(co) {
  tp <- sum(df$CC == 1 & df$Prob > co)
  fn <- sum(df$CC == 1 & df$Prob < co)
  return(tp/(tp + fn))
}) %>% unlist(use.names = F)

# fpr = fp / (fp + tn) 
stats_df$FPR <- map(cutoffs, function(co) {
  fp <- sum(df$CC == 0 & df$Prob > co)
  tn <- sum(df$CC == 0 & df$Prob < co)
  return(fp/(fp + tn))
}) %>% unlist(use.names = F)

# specificity: 1 - FPR = TN / (FP + TN)
stats_df$specificity <- 1 - stats_df$FPR

# ppv = tp / (tp + fp)
stats_df$PPV <- map(cutoffs, function(co) {
  tp <- sum(df$CC == 1 & df$Prob > co)
  fp <- sum(df$CC == 0 & df$Prob > co)
  return(tp/(tp + fp))
}) %>% unlist(use.names = F)

# npv = tn / (tn + fn)
stats_df$NPV <- map(cutoffs, function(co) {
  tn <- sum(df$CC == 0 & df$Prob < co) 
  fn <- sum(df$CC == 1 & df$Prob < co)
  return(tn/(tn + fn))
}) %>% unlist(use.names = F)



# Save --------------------------------------------------------------------

write.csv(stats_df, 'intermediate_data/stats_df.csv')

