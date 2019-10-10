# Find clinical hook association
# Library -----------------------------------------------------------------
library(tidyverse)
library(glmnet)
library(openxlsx)

# Load --------------------------------------------------------------------
# Use the tw = 24, tp = 3 data
CC_comb <- readRDS('modeling_data/train24_3.rds')
CC_comb_val <- readRDS('modeling_data/test24_3.rds')

# T25-FW data
walk_data <- read.xlsx("raw_data/Box/Boston/MS CLIMB Data/CLIMB Cohort Add'l Data/25ftWalk121615.xlsx",
                       detectDates = F)
walk_data$visit_date <- convertToDate(walk_data$visit_date, origin = "1900-01-01")

# Clean walk data ---------------------------------------------------------
# Exclude those with time = 999.99 
walk_data <- walk_data %>% filter(`25ft_walking_time` != 999.99)
# Remove provider info
walk_data$provider <- NULL
# Rename columns
colnames(walk_data) <- c('PatientID', 'StartDate', 'WalkTime')






