# Preprocess data for DMT causal inference (must be run after preprocess)
# Library -----------------------------------------------------------------
library(tidyverse)
library(openxlsx)
library(lubridate)

# Load --------------------------------------------------------------------
readin = TRUE 
wkpath  = "raw_data/Box/Boston/MS CLIMB Data/"
wkpath2 = "raw_data/Box/Boston/MS CLIMB Paper2 Code&Result/Step3_BuildModel/IntermediateDataResult/"
year   = 2000:2016

MS_cohort = read.csv(paste0(wkpath2,"Cleaned_MS_cohort.csv"),stringsAsFactors = FALSE,
                     colClasses = c(rep("character",4), rep("Date",4), rep("integer",6),rep("numeric",3),
                                    "integer",rep("Date",4),rep("numeric",3),"integer"))
MS_attack = read.csv(paste0(wkpath2,"Cleaned_MS_attack.csv"),stringsAsFactors = FALSE,
                     colClasses = c("character",rep("integer",6),"Date",rep("numeric",3),
                                    rep("character",12),rep("Date",2),"integer",
                                    rep("Date",5),rep("integer",2)))
MS_trt    = read.csv(paste0(wkpath2,"Cleaned_MS_treatment.csv"),stringsAsFactors = FALSE)



# Group DMTs --------------------------------------------------------------

# Dimethyl fumarate
dmf = c('TECFIDERA (BG-12)')
# Fingolimod
fingo = c('GILENYA')
# Glatiramer acetate
glat_ac = c('COPAXONE')
# Interferon-beta 1
intb1 = c('AVONEX','BETASERON','EXTAVIA','PLEGRIDY','REBIF')
# Natalizumab
nata = c('TYSABRI')
# Teriflunomide
terif = c('AUBAGIO')


MS_trt %>% filter(medication_desc %in% dmf) %>% 
  select(PatientID) %>% 
  unlist(use.names = F) %>% 
  unique %>% 
  length


tmp = MS_trt %>% filter(medication_desc %in% dmf) 
tmp <- tmp %>% filter(!is.na(val_start) & !is.na(val_stop))
tmp$treat_dura <- as.numeric(difftime(tmp$val_stop, tmp$val_start, units = 'weeks')/4)
tmp$treat_dura %>% summary


MS_Enct_Uniq = read.csv("raw_data/HSPH/Data From Partners MS Center/MS_Encounter_Date_PatientID.csv", 
                        stringsAsFactors = FALSE)
MS_Enct_Uniq$StartDate <- as.Date(MS_Enct_Uniq$StartDate, '%m/%d/%y')
pt_nums <- MS_Enct_Uniq$PatientNum %>% unique
btw_visit_list <- list()
for (j in 1:length(pt_nums)) {
  pt_subset <- MS_Enct_Uniq %>% filter(PatientNum == pt_nums[j])
  pt_subset <- pt_subset[order(pt_subset$StartDate),]
  if (nrow(pt_subset) == 1) {
    btw_visit_list[[j]] <- 0
  }
  else {
    store <- rep(NA, nrow(pt_subset) - 1)
    for (i in 1:(nrow(pt_subset) - 1)) {
      store[i] <- difftime(pt_subset[i + 1,'StartDate'],pt_subset[i,'StartDate'], units = 'days')
    }
    btw_visit_list[[j]] <- store
  }
}

vec <- do.call(c,btw_visit_list)
vec <- vec[vec != 0]
mean(vec)
sd(vec)



