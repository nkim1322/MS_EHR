# Probability plots over time per patient, stratified 3 part model
# Library -----------------------------------------------------------------

library(tidyverse)
library(glmnet)
library(pROC)

# Load --------------------------------------------------------------------

CC_comb <- readRDS('intermediate_data/CC_comb_24_3_clin.rds') # Training data
CC_comb_val <- readRDS('intermediate_data/CC_comb_val_24_3_clin.rds') # Validation data
MS_attack <- read.csv("raw_data/Box/Boston/MS CLIMB Paper2 Code&Result/Step3_BuildModel/IntermediateDataResult/Cleaned_MS_attack.csv",
                      stringsAsFactors = FALSE,
                      colClasses = c("character", rep("integer",6), "Date", rep("numeric",3),
                                     rep("character",12), rep("Date",2), "integer",
                                     rep("Date",5), rep("integer",2)))

# Set parameters ----------------------------------------------------------
## WORKING WITH 24 MONTHS/3 MONTHS MODEL
tw = 24; tp = 3
readin = TRUE  # Read in 3 part model list

# Fit 3 models: 2000-2005, 2006-2010, 2011-2016 ---------------------------

load('models/glmnet.list.rda')
# load('models/glmnet.list.noclin.rda'); glmnet.list <- glmnet.list.noclin

# Predict probability of relapse ------------------------------------------
df <- CC_comb_val[,1:6]
df$Prob <- NA

# 2000-2005
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

# saveRDS(df, 'intermediate_data/3partclin_probs.rds')
saveRDS(df, 'intermediate_data/3part_probs.rds')

# Plot --------------------------------------------------------------------
# df should be a dataframe including predicted probability

plot_per_pt <- function(df) {
  df_tibble <- as_tibble(df)
  df_tibble$StartDate <- as.Date(df_tibble$StartDate)
  
  # Graph validation probs per pt
  pts <- unique(df$PatientID)
  graphs_list <- list()
  for (p_id in pts) {
    df_test = df_tibble[df_tibble$PatientID == p_id,]
    
    graphs_list[[p_id]] <- 
      if (df_test %>% filter(CC == 1) %>% nrow > 0) {
        ggplot() + 
          geom_vline(aes(xintercept = onset, color = 'red'),
                     data = MS_attack %>% filter(PatientID == df_test$PatientID[1])) +
          geom_text(data = MS_attack %>% filter(PatientID == df_test$PatientID[1])
                    , mapping = aes(x = onset, y = 0.75, label = onset), 
                    size = 4, angle = 90, vjust = -0.4, hjust = 0, col = 'red') +
          geom_smooth(data = df_test, aes(StartDate + 365, Prob), color = 'black') +
          ylim(0, 1) +
          ggtitle(paste0("Patient ID: ", p_id)) + 
          theme(plot.title = element_text(hjust = 0.5),
                axis.text.x = element_text(angle = 90),
                legend.position = "none") + 
          xlab('Visit Date') + ylab('Probability of Relapse')
      } else {
        ggplot() + 
          geom_smooth(data = df_test, aes(StartDate + 365, Prob), color = 'black') +
          ylim(0, 1) +
          ggtitle(paste0("Patient ID: ", p_id)) + 
          theme(plot.title = element_text(hjust = 0.5),
                axis.text.x = element_text(angle = 90),
                legend.position = "none") + 
          xlab('Visit Date') + ylab('Probability of Relapse')
      }
  }
  
  return(graphs_list)
}


my_graphs_list <- plot_per_pt(df)

pdf('plots/24mons3mons_probplots_3mdl_clin.pdf')
# pdf('plots/24mons3mons_probplots_3mdl.pdf')
my_graphs_list
dev.off()

roc(unlist(CC_comb_val[,4]), as.numeric(df$Prob))


# Plot ROC to compare clinical and no clinical vars -----------------------
# 
load('models/glmnet.list.rda')
load('models/glmnet.list.noclin.rda')

df <- readRDS('intermediate_data/3partclin_probs.rds')
CC_comb_val_noclin <- readRDS('intermediate_data/CC_comb_val_24_3.rds')

df2 <- CC_comb_val_noclin[,1:6]
df2$Prob <- NA

# 2000-2005
df2$Prob[CC_comb_val_noclin$StartDate < as.Date('2006-01-01')] <-
  predict(glmnet.list.noclin[['2000-2005']],
          newx = as.matrix(CC_comb_val_noclin %>%
                             filter(StartDate < as.Date('2006-01-01')) %>%
                             select(-c(1:9))),
          type = "response", s = glmnet.list.noclin[['2000-2005']]$lambda.1se)
# 2005-2010
df2$Prob[CC_comb_val_noclin$StartDate >= as.Date('2006-01-01') &
           CC_comb_val_noclin$StartDate < as.Date('2011-01-01')] <-
  predict(glmnet.list.noclin[['2006-2010']],
          newx = as.matrix(CC_comb_val_noclin %>%
                             filter(StartDate >= as.Date('2006-01-01') &
                                      StartDate < as.Date('2011-01-01')) %>%
                             select(-c(1:9))),
          type = "response", s = glmnet.list.noclin[['2006-2010']]$lambda.1se)
# 2010-2016
df2$Prob[CC_comb_val_noclin$StartDate >= as.Date('2011-01-01')] <-
  predict(glmnet.list.noclin[['2011-2016']],
          newx = as.matrix(CC_comb_val_noclin %>%
                             filter(StartDate >= as.Date('2011-01-01')) %>%
                             select(-c(1:9))),
          type = "response", s = glmnet.list.noclin[['2011-2016']]$lambda.1se)

# Plot ROC curve
roc1 <- roc(unlist(CC_comb_val[,4]), as.numeric(df$Prob))
roc2 <- roc(unlist(CC_comb_val_noclin[,4]), as.numeric(df2$Prob))

pdf('plots/ROC_clin_vs_noclin.pdf')
plot(roc1, col = 2, main = "ROC, 24mons/3mons", legacy.axes = TRUE)
plot(roc2, col = 3, add = TRUE)
legend("bottomright",
       legend = c(paste("Clinical, AUC = ", round(roc1$auc, 4)),
                  paste("No clinical, AUC = ", round(roc2$auc, 4))),
       col = c(2, 3), lty = 1, cex = 0.8)
dev.off()



# Analyze individual model AUC --------------------------------------------

# No clinical variables
# 2000-2005:
roc(CC_comb_val %>% filter(StartDate < as.Date('2006-01-01')) %>% select(CC) %>% unlist(use.names = F),
    df2 %>% filter(StartDate < as.Date('2006-01-01')) %>% select(Prob) %>% unlist(use.names = F))
# 2006-2010
roc(CC_comb_val %>% filter(StartDate >= as.Date('2006-01-01') &
                             StartDate < as.Date('2011-01-01')) %>% select(CC) %>% unlist(use.names = F),
    df2 %>% filter(StartDate >= as.Date('2006-01-01') &
                     StartDate < as.Date('2011-01-01')) %>% select(Prob) %>% unlist(use.names = F))
# 2011-2016
roc(CC_comb_val %>% filter(StartDate >= as.Date('2011-01-01')) %>% select(CC) %>% unlist(use.names = F),
    df2 %>% filter(StartDate >= as.Date('2011-01-01')) %>% select(Prob) %>% unlist(use.names = F))

# With clinical variables
# 2000-2005:
roc(CC_comb_val %>% filter(StartDate < as.Date('2006-01-01')) %>% select(CC) %>% unlist(use.names = F),
    df %>% filter(StartDate < as.Date('2006-01-01')) %>% select(Prob) %>% unlist(use.names = F))
# 2006-2010
roc(CC_comb_val %>% filter(StartDate >= as.Date('2006-01-01') &
                             StartDate < as.Date('2011-01-01')) %>% select(CC) %>% unlist(use.names = F),
    df %>% filter(StartDate >= as.Date('2006-01-01') &
                    StartDate < as.Date('2011-01-01')) %>% select(Prob) %>% unlist(use.names = F))
# 2011-2016
roc(CC_comb_val %>% filter(StartDate >= as.Date('2011-01-01')) %>% select(CC) %>% unlist(use.names = F),
    df %>% filter(StartDate >= as.Date('2011-01-01')) %>% select(Prob) %>% unlist(use.names = F))




