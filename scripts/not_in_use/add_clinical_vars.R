# Add clinical variables to 24mons/3mons aggregated data
# Packages ----------------------------------------------------------------
library(tidyverse)
library(openxlsx)
library(lubridate)

# Parameters --------------------------------------------------------------
wkpath  = "raw_data/Box/Boston/MS CLIMB Data/"
wkpath2 = "raw_data/Box/Boston/MS CLIMB Paper2 Code&Result/Step3_BuildModel/IntermediateDataResult/"
tw = 24; tp = 3

# Load --------------------------------------------------------------------
# MS data
MS_cohort <- read.csv(paste0(wkpath2,"Cleaned_MS_cohort.csv"),
                      stringsAsFactors = FALSE,
                      colClasses = c(rep("character",4), rep("Date",4), rep("integer",6),
                                     rep("numeric",3), "integer", rep("Date",4),
                                     rep("numeric",3), "integer"))
# 24mons/3mons data
CC_comb <- readRDS('intermediate_data/CC_comb_24_3.rds') # Original training data
CC_comb_val <- readRDS('intermediate_data/CC_comb_val_24_3.rds') # Original validation data


# Subset into relevant data
cohort_data <- MS_cohort[, c('PatientID', 'SEX', 'RACE_DESC', 'ETHNICITY_DESC', 'DURA',
                             'AGE_AT_LASTVISIT', 'LAST_VISIT_DATE', 'AGE_AT_FIRSTSYMPTOM')]
# 'DURA' variable is follow-up duration (in years)

# Convert Encounter Date to date type -------------------------------------
CC_comb$StartDate <- as.Date(CC_comb$StartDate)
CC_comb_val$StartDate <- as.Date(CC_comb_val$StartDate)

# Sex ---------------------------------------------------------------------
# Create binary indicator for being female
tmp <- cohort_data[cohort_data$SEX == 'F','PatientID']
CC_comb$FEMALE <- ifelse(CC_comb$PatientID %in% tmp, 1, 0)
CC_comb_val$FEMALE <- ifelse(CC_comb_val$PatientID %in% tmp, 1, 0)


# Race --------------------------------------------------------------------
CC_comb$RACE <- NA; CC_comb_val$RACE <- NA

races <- cohort_data$RACE_DESC %>% unique()
for (race in races) {
  tmp <- cohort_data[cohort_data$RACE_DESC == race,'PatientID']
  CC_comb$RACE[CC_comb$PatientID %in% tmp] <- race
  CC_comb_val$RACE[CC_comb_val$PatientID %in% tmp] <- race
}

# Ethnicity ---------------------------------------------------------------
CC_comb$ETHNICITY <- NA; CC_comb_val$ETHNICITY <- NA

ethnicities <- cohort_data$ETHNICITY_DESC %>% unique()
for (ethn in ethnicities) {
  tmp <- cohort_data[cohort_data$ETHNICITY_DESC == ethn,'PatientID']
  CC_comb$ETHNICITY[CC_comb$PatientID %in% tmp] <- ethn
  CC_comb_val$ETHNICITY[CC_comb_val$PatientID %in% tmp] <- ethn
}


# Create "non-Hispanic European descent" variable (Option 1) --------------
# RACE1 == 1 if "non-Hispanic European descent"

tmp <- cohort_data[cohort_data$RACE == 'White' & cohort_data$ETHNICITY == 'Not hispanic or latino',
                   'PatientID']
CC_comb$RACE1 <- ifelse(CC_comb$PatientID %in% tmp, 1, 0)
CC_comb_val$RACE1 <- ifelse(CC_comb_val$PatientID %in% tmp, 1, 0)


# Preserve race and ethnicity (Option 2) ----------------------------------
# RACE2 == 1 if "White", ETHNICITY2 == 1 if "Not hispanic or latino"

tmp <- cohort_data[cohort_data$RACE == 'White', 'PatientID']
CC_comb$RACE2 <- ifelse(CC_comb$PatientID %in% tmp, 1, 0)
CC_comb_val$RACE2 <- ifelse(CC_comb_val$PatientID %in% tmp, 1, 0)

tmp <- cohort_data[cohort_data$ETHNICITY == 'Not hispanic or latino', 'PatientID']
CC_comb$ETHNICITY2 <- ifelse(CC_comb$PatientID %in% tmp, 1, 0)
CC_comb_val$ETHNICITY2 <- ifelse(CC_comb_val$PatientID %in% tmp, 1, 0)

# Remove Race and Ethnicity categorical variables
CC_comb$RACE <- NULL; CC_comb_val$RACE <- NULL
CC_comb$ETHNICITY <- NULL; CC_comb_val$ETHNICITY <- NULL


# Disease duration --------------------------------------------------------

# Age at first symptom (Find a more efficient way to do this**)
CC_comb$AGE_AT_FIRSTSYMPTOM <- NA ; CC_comb_val$AGE_AT_FIRSTSYMPTOM <- NA
for (pt in unique(CC_comb$PatientID)) {
  CC_comb[CC_comb$PatientID == pt, 'AGE_AT_FIRSTSYMPTOM'] <- 
    cohort_data[cohort_data$PatientID == pt, 'AGE_AT_FIRSTSYMPTOM']
}
for (pt in unique(CC_comb_val$PatientID)) {
  CC_comb_val[CC_comb_val$PatientID == pt, 'AGE_AT_FIRSTSYMPTOM'] <- 
    cohort_data[cohort_data$PatientID == pt, 'AGE_AT_FIRSTSYMPTOM']
}

# Create birthday variable
cohort_data$BIRTHDAY <- with(cohort_data, LAST_VISIT_DATE - years(round(AGE_AT_LASTVISIT)))

CC_comb$BIRTHDAY <- NA ; CC_comb_val$BIRTHDAY <- NA
for (pt in unique(CC_comb$PatientID)) {
  CC_comb[CC_comb$PatientID == pt, 'BIRTHDAY'] <- 
    cohort_data[cohort_data$PatientID == pt, 'BIRTHDAY']
}
for (pt in unique(CC_comb_val$PatientID)) {
  CC_comb_val[CC_comb_val$PatientID == pt, 'BIRTHDAY'] <- 
    cohort_data[cohort_data$PatientID == pt, 'BIRTHDAY']
}
CC_comb$BIRTHDAY <- as.Date(CC_comb$BIRTHDAY, origin = '1970-01-01')
CC_comb_val$BIRTHDAY <- as.Date(CC_comb_val$BIRTHDAY, origin = '1970-01-01')

# Create disease duration variable (in years)
CC_comb$DISEASE_DURA <- with(CC_comb, 
                             difftime(StartDate, BIRTHDAY, units = 'days')/365 - AGE_AT_FIRSTSYMPTOM) %>% 
  as.numeric()
CC_comb_val$DISEASE_DURA <- with(CC_comb_val, 
                                 difftime(StartDate, BIRTHDAY, units = 'days')/365 - AGE_AT_FIRSTSYMPTOM) %>% 
  as.numeric()

# Remove variables BIRTHDAY and AGE_AT_FIRSTSYMPTOM (not needed as predictors)
CC_comb$BIRTHDAY <- NULL; CC_comb_val$BIRTHDAY <- NULL
CC_comb$AGE_AT_FIRSTSYMPTOM <- NULL; CC_comb_val$AGE_AT_FIRSTSYMPTOM <- NULL

# Follow-up duration ------------------------------------------------------
# Taken from the original variable DURA created by Liang (see data cleaning in tune_parameters.Rmd)
CC_comb$FOLLOWUP_DURA <- NA ; CC_comb_val$FOLLOWUP_DURA <- NA
for (pt in unique(CC_comb$PatientID)) {
  CC_comb[CC_comb$PatientID == pt, 'FOLLOWUP_DURA'] <- 
    cohort_data[cohort_data$PatientID == pt, 'DURA']
}
for (pt in unique(CC_comb_val$PatientID)) {
  CC_comb_val[CC_comb_val$PatientID == pt, 'FOLLOWUP_DURA'] <- 
    cohort_data[cohort_data$PatientID == pt, 'DURA']
}


# Remove NA's -------------------------------------------------------------

CC_comb <- CC_comb[!is.na(CC_comb$RACE1) & !is.na(CC_comb$RACE2) & !is.na(CC_comb$ETHNICITY2) & 
                     !is.na(CC_comb$DISEASE_DURA) & !is.na(CC_comb$FOLLOWUP_DURA),]
CC_comb_val <- CC_comb_val[!is.na(CC_comb_val$RACE1) & !is.na(CC_comb_val$RACE2) & !is.na(CC_comb_val$ETHNICITY2) & 
                             !is.na(CC_comb_val$DISEASE_DURA) & !is.na(CC_comb_val$FOLLOWUP_DURA),]


# Save --------------------------------------------------------------------

saveRDS(CC_comb, 'intermediate_data/CC_comb_24_3_clin.rds')
saveRDS(CC_comb_val, 'intermediate_data/CC_comb_val_24_3_clin.rds')






