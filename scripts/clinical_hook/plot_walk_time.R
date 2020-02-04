# Library -----------------------------------------------------------------
library(tidyverse)
library(openxlsx)

# Load --------------------------------------------------------------------
wkpath2 = "raw_data/Box/Boston/MS CLIMB Paper2 Code&Result/Step3_BuildModel/IntermediateDataResult/"
MS_attack = read.csv(paste0('thesis_tmp/',"Cleaned_MS_attack.csv"),stringsAsFactors = FALSE,
                     colClasses = c("character",rep("integer",6),"Date",rep("numeric",3),
                                    rep("character",12),rep("Date",2),"integer",
                                    rep("Date",5),rep("integer",2)))
MS_attack_subs <- MS_attack %>% select(PatientID, onset, val_start, val_stop)

# T25-FW data
walk_data <- read.xlsx("raw_data/Box/Boston/MS CLIMB Data/CLIMB Cohort Add'l Data/25ftWalk121615.xlsx",
                       detectDates = F)
walk_data$visit_date <- convertToDate(walk_data$visit_date, origin = "1900-01-01")


# Filter ------------------------------------------------------------------

# Exclude those with time less than 0 or greater than 180 seconds
walk_data <- walk_data %>% filter(`25ft_walking_time` <= 180 & `25ft_walking_time` >= 0)
# Remove provider info
walk_data$provider <- NULL
# Rename columns
colnames(walk_data) <- c('PatientID', 'WalkDate', 'WalkTime')
# Log transform WalkTime due to right skewedness
# walk_data$WalkTime <- log(walk_data$WalkTime)

# Get validated patients


# Investigate -------------------------------------------------------------
val_pts <- MS_attack_subs %>% filter(!is.na(val_start)) %>% select(PatientID) %>% unlist(use.names = F) %>% unique()
unique_pts <- unique(MS_attack_subs$PatientID)

graphs_list <- list()
for (p_id in unique_pts) {
  walk_subset <- walk_data %>% filter(PatientID == p_id)
  MS_subset <- MS_attack_subs %>% filter(PatientID == p_id)
  if (nrow(MS_subset) != 0) {
    graphs_list[[p_id]] <- ggplot() + 
      geom_point(data = walk_subset, aes(WalkDate, log(WalkTime))) + 
      geom_smooth(data = walk_subset, aes(WalkDate, log(WalkTime))) + 
      geom_vline(data = MS_subset, aes(xintercept = onset, color = 'red')) +
      ggtitle(paste('PatientID: ',p_id)) + 
      ylim(-1, 6) + 
      theme(plot.title = element_text(hjust = 0.5),
            legend.position = "none") 
  } else {
    graphs_list[[p_id]] <- ggplot() + 
      geom_point(data = walk_subset, aes(WalkDate, log(WalkTime))) + 
      geom_smooth(data = walk_subset, aes(WalkDate, log(WalkTime))) + 
      ggtitle(paste('PatientID: ',p_id)) + 
      ylim(-1,6)
      theme(plot.title = element_text(hjust = 0.5),
            legend.position = "none")
  }
}

pdf('thesis_tmp/walk_data.pdf')
# pdf('thesis_tmp/walk_data_val.pdf')
graphs_list
dev.off()

rm(list = ls())

