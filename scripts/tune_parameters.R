# Tune tw (time window) and tp (time period) parameters
# Library -----------------------------------------------------------------
library(tidyverse)
library(glmnet)
library(ggplot2)
library(pROC)

# Load --------------------------------------------------------------------
tw_list     = rep(c(6, 12, 24), each = 3) # time window in months
tp_list     = rep(c(1, 3, 6), 3) # tw/4, time period for sample cases and controls in months

# Model and save predictions ----------------------------------------------
models.list <- list()
roc.list <- list()

for (i in 1:length(tw_list)) {
  tw = tw_list[i]; tp = tp_list[i]
  print(paste('tw:', tw, 'tp:',tp))
  
  CC_comb <- readRDS(paste0('modeling_data/train', tw, '_', tp, '.rds'))
  CC_comb_val <- readRDS(paste0('modeling_data/test', tw, '_', tp, '.rds'))
  
  # Filter 2006
  CC_comb <- CC_comb %>% filter(StartDate > as.Date('2006-01-01'))
  CC_comb_val <- CC_comb_val %>% filter(StartDate > as.Date('2006-01-01'))
  
  glmnet.fit = cv.glmnet(as.matrix(CC_comb[,-c(1:9)]), unlist(CC_comb$CC),
                         family = "binomial", type.measure = "deviance")
  preds <- predict(glmnet.fit, 
                   newx = as.matrix(CC_comb_val[,-c(1:9)]),
                   type = "response", 
                   s = glmnet.fit$lambda.1se)
  models.list[[paste0(tw, "_", tp)]] <- glmnet.fit
  roc.list[[paste0(tw, "_", tp)]] <- roc(unlist(CC_comb_val[,4]), as.numeric(preds))
}



# Plot and Save -----------------------------------------------------------
pdf('plots/paper2_final/roc_tune.pdf')

# Plot: tp = 1 month
plot(roc.list[[1]], col = 2, main = "ROC, Time Period = 1 month", legacy.axes = TRUE)
plot(roc.list[[4]], col = 3, add = TRUE)
plot(roc.list[[7]], col = 4, add = TRUE)
legend("bottomright", 
       legend = c(paste("6 mons: AUC = ", round(roc.list[[1]]$auc, 4)), 
                  paste("12 mons: AUC = ", round(roc.list[[4]]$auc, 4)),
                  paste("24 mons: AUC = ", round(roc.list[[7]]$auc, 4))),
       col = c(2, 3, 4), lty = 1, cex = 0.8, 
       title = "Time window")

# Plot: tp = 3 months
plot(roc.list[[2]], col = 2, main = "ROC, Time Period = 3 months", legacy.axes = TRUE)
plot(roc.list[[5]], col = 3, add = TRUE)
plot(roc.list[[8]], col = 4, add = TRUE)
legend("bottomright", 
       legend = c(paste("6 mons: AUC = ", round(roc.list[[2]]$auc, 4)), 
                  paste("12 mons: AUC = ", round(roc.list[[5]]$auc, 4)),
                  paste("24 mons: AUC = ", round(roc.list[[8]]$auc, 4))),
       col = c(2, 3, 4), lty = 1, cex = 0.8,
       title = "Time window")

# Plot: tp = 6 months
plot(roc.list[[3]], col = 2, main = "ROC, Time Period = 6 months", legacy.axes = TRUE)
plot(roc.list[[6]], col = 3, add = TRUE)
plot(roc.list[[9]], col = 4, add = TRUE)
legend("bottomright", 
       legend = c(paste("6 mons: AUC = ", round(roc.list[[3]]$auc, 4)), 
                  paste("12 mons: AUC = ", round(roc.list[[6]]$auc, 4)),
                  paste("24 mons: AUC = ", round(roc.list[[9]]$auc, 4))),
       col = c(2, 3, 4), lty = 1, cex = 0.8,
       title = "Time window")

dev.off()



