# Random Forest
# Library -----------------------------------------------------------------
library(tidyverse)
library(randomForest)

# Load --------------------------------------------------------------------
CC_comb <- readRDS('intermediate_data/CC_comb_24_3_clin.rds')
CC_comb_val <- readRDS('intermediate_data/CC_comb_val_24_3_clin.rds')


# Filter ------------------------------------------------------------------
CC_comb <- CC_comb %>% filter(StartDate >= as.Date('2006-01-01'))
CC_comb_val <- CC_comb_val %>% filter(StartDate >= as.Date('2006-01-01'))

CC_comb$RACE2 <- NULL; CC_comb_val$RACE2 <- NULL
CC_comb$ETHNICITY2 <- NULL; CC_comb_val$ETHNICITY2 <- NULL


# Fit ---------------------------------------------------------------------
set.seed(12345)
RF.fit  = randomForest(x = CC_comb[,-c(1:9)], 
                       y = factor(CC_comb$CC,levels = c(0,1)),
                       mtry = 14,
                       importance = TRUE)
save(RF.fit, file = 'models/rf.rda')

# Predict -----------------------------------------------------------------
RF.prob = predict(RF.fit, newdata = CC_comb_val[,-c(1:9)], type="prob")
RF.roc  = roc(unlist(CC_comb_val[,4]), RF.prob[,2])


# Analyze importance ------------------------------------------------------
impt <- RF.fit$importance

df_imp <- data.frame(term = row.names(impt),
                     MeanDecreaseAccuracy = impt[,'MeanDecreaseAccuracy'],
                     MeanDecreaseGini = impt[,'MeanDecreaseGini'])
row.names(df_imp) <- NULL


df_imp$term <- strtrim(df_imp$term, width = 20)
df_imp[order(df_imp$MeanDecreaseGini, decreasing = T),] %>% head




