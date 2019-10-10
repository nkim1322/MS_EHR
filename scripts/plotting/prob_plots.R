# Probability plots over time per patient, single model
# Library -----------------------------------------------------------------
library(tidyverse)
library(glmnet)
# Load --------------------------------------------------------------------
CC_comb <- readRDS('modeling_data/train24_3.rds')
CC_comb_val <- readRDS('modeling_data/test24_3.rds')

# Relapse data
MS_attack <- read.csv("raw_data/Box/Boston/MS CLIMB Paper2 Code&Result/Step3_BuildModel/IntermediateDataResult/Cleaned_MS_attack.csv",
                      stringsAsFactors = FALSE,
                      colClasses = c("character", rep("integer",6), "Date", rep("numeric",3),
                                     rep("character",12), rep("Date",2), "integer",
                                     rep("Date",5), rep("integer",2)))

# Fit model ---------------------------------------------------------------
source('scripts/modeling/model.R')

# Predict -----------------------------------------------------------------
df <- CC_comb_val[,1:6]
df$Prob <- predict(glmnet.fit, newx = as.matrix(CC_comb_val[,-c(1:9)]),
                   type = "response", s = glmnet.fit$lambda.1se)

# Plot function -----------------------------------------------------------
# df should be a dataframe including predicted probability
plot_per_pt <- function(df) {
  df_tibble <- as_tibble(df)
  df_tibble$StartDate <- as.Date(df_tibble$StartDate)
  
  # Graph validation probs per pt
  pts <- unique(df$PatientID)
  graphs_list <- list()
  for (p_id in pts) {
    df_test = df_tibble[df_tibble$PatientID == p_id,]
    # first_date = df_test
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


# Plot and save -----------------------------------------------------------

my_graphs_list <- plot_per_pt(df)
pdf('plots/paper2_final/24mons3mons_probplots.pdf')
my_graphs_list
dev.off()