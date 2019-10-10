## Data screening
# Library -----------------------------------------------------------------
library(tidyverse)

# Load --------------------------------------------------------------------
CC_comb <- readRDS('intermediate_data/CC_comb_24_3_clin.rds')
CC_comb_val <- readRDS('intermediate_data/CC_comb_val_24_3_clin.rds')

wkpath  = "raw_data/Box/Boston/MS CLIMB Data/"
wkpath2 = "raw_data/Box/Boston/MS CLIMB Paper2 Code&Result/Step3_BuildModel/IntermediateDataResult/"
ICDPheCode   = read.csv(paste0(wkpath,"EHR/MS_AllEncounters_ICD_Data_03282019.csv"), stringsAsFactors = FALSE)
ICD_CPT_CUI_Comb = read.csv(paste0(wkpath2,"ICD_CPT_CUI_Comb.csv"), stringsAsFactors = FALSE)
MS_cohort = read.csv(paste0(wkpath2,"Cleaned_MS_cohort.csv"),
                     stringsAsFactors = FALSE,
                     colClasses = c(rep("character",4), rep("Date",4), rep("integer",6), rep("numeric",3),
                                    "integer", rep("Date",4), rep("numeric",3), "integer"))
MS_map  = read.xlsx(paste0(wkpath,"CLIMB Cohort/Spec95_i2b2_Mapping2017.xlsx"),
                    sheet = 1); colnames(MS_map) = c("PatientNum","PatientID")

# Set the first encounter to be  the first ICD9 code encounter (335) --------
# In the data set, 'PheCode.335_' is the code for Multiple sclerosis

# All ICD data
library(lubridate)
full_icd_converted <- read.csv('raw_data/Box/Boston/MS CLIMB Data/EHR/full_icd_data_converted.csv',
                               sep = '|')
full_icd_converted$icd10_to_icd9 <- NULL
full_icd_converted$start.date <- as.character(full_icd_converted$start.date)
full_icd_converted$start.date <- as.Date(full_icd_converted$start.date)

full_icd_converted_cohort <- full_icd_converted[full_icd_converted$patient.num%in%MS_map$PatientNum,]
full_icd_converted_cohort_MS <- full_icd_converted_cohort %>% 
  filter(code == 340 | code == 0340 | code == 'ICD10:G35')


# Get first ICD MS code for training and validation patients
## Training data
CC_comb$FIRST_MS_CODE <- NA
for (p_num in unique(CC_comb$PatientNum)) {
  p_subset <- full_icd_converted_cohort_MS %>% filter(patient.num == p_num)
  CC_comb[CC_comb$PatientNum == p_num, 'FIRST_MS_CODE'] <- min(p_subset$start.date)
}
CC_comb$FIRST_MS_CODE <- as.Date(CC_comb$FIRST_MS_CODE, origin = '1970-01-01')
# Filter encounters after First MS Code
data_list <- list()
for (p_num in unique(CC_comb$PatientNum)) {
  p_subset <- CC_comb %>% filter(PatientNum == p_num)
  first_MS_code <- p_subset$FIRST_MS_CODE[1]
  data_list[[p_num]] <- p_subset %>% filter(StartDate >= first_MS_code)
}
CC_comb <- do.call(rbind, data_list)
CC_comb_val$FIRST_MS_CODE <- NULL

## Validation data
CC_comb_val$FIRST_MS_CODE <- NA
for (p_num in unique(CC_comb_val$PatientNum)) {
  p_subset <- full_icd_converted_cohort_MS %>% filter(patient.num == p_num)
  CC_comb_val[CC_comb_val$PatientNum == p_num, 'FIRST_MS_CODE'] <- min(p_subset$start.date)
}
CC_comb_val$FIRST_MS_CODE <- as.Date(CC_comb_val$FIRST_MS_CODE, origin = '1970-01-01')
# Filter encounters after First MS Code
data_list2 <- list()
for (p_num in unique(CC_comb_val$PatientNum)) {
  p_subset <- CC_comb_val %>% filter(PatientNum == p_num)
  first_MS_code <- p_subset$FIRST_MS_CODE[1]
  data_list2[[p_num]] <- p_subset %>% filter(StartDate >= first_MS_code)
}
CC_comb_val <- do.call(rbind, data_list2)
CC_comb_val$FIRST_MS_CODE <- NULL

# Screening Option #1
# Filter out those with no subsequent encounter w/i 1, 2, and 3 years --------
yrs <- c(1,2,3)
screen_yrs_data <- list()
for (yr in yrs) {
  # Training data
  df_train <- CC_comb
  pts_train <- unique(df_train$PatientID)
  for (p_id in pts_train) {
    p_subset <- df_train %>% filter(PatientID == p_id)
    # Remove patients with only 1 encounter
    if (nrow(p_subset) == 1) {
      df_train <- df_train[!df_train$PatientID == p_id, ]
    } else { 
      for (i in 1:(nrow(p_subset) - 1)) {
        enctr_date <- p_subset[i,'StartDate']
        
        if (difftime(p_subset[i + 1,'StartDate'], enctr_date, units = 'days') > yr*365) {
          df_train <- df_train[!(df_train$PatientID == p_id & df_train$StartDate == enctr_date), ]
        }
      }
    }
  }
  
  # Validation data
  df_val <- CC_comb_val
  pts_val <- unique(df_val$PatientID)
  for (p_id in pts_val) {
    p_subset <- df_val %>% filter(PatientID == p_id)
    # Remove patients with only 1 encounter
    if (nrow(p_subset) == 1) {
      df_val <- df_val[!df_val$PatientID == p_id, ]
    } else {
      for (i in 1:(nrow(p_subset) - 1)) {
        enctr_date <- p_subset[i,'StartDate']
        if (difftime(p_subset[i + 1,'StartDate'], enctr_date, units = 'days') > yr*365) {
          df_val <- df_val[!(df_val$PatientID == p_id & df_val$StartDate == enctr_date), ]
        }
      }
    }
  }
  
  
  screen_yrs_data[[yr]] <- list(train = df_train,
                                val = df_val)
}
save(screen_yrs_data, file = 'intermediate_data/screen_yrs_data.rda')

# Screening Option #2
# Keep data only 2006 & after  -----------------------------------------------
# Can do this in modeling script***


