# Remove parent and subcodes per Aaron Sonabend
# Library -----------------------------------------------------------------
library(tidyverse)

# Load --------------------------------------------------------------------
CC_comb <- readRDS('intermediate_data/CC_comb_24_3_clin.rds') # Training data
CC_comb_val <- readRDS('intermediate_data/CC_comb_val_24_3_clin.rds') # Validation data

codes_master <- read.csv('raw_data/clean_RPDRml_features_all (1).csv')
codes_master$X <- NULL

wkpath  = "raw_data/Box/Boston/MS CLIMB Data/"
ICDPheCode   = read.csv(paste0(wkpath,"EHR/MS_AllEncounters_ICD_Data_03282019.csv"),
                        stringsAsFactors = FALSE);ICDPheCode$phecode[ICDPheCode$concept_cd == "LPA268"] = "335_"

# Work with data ----------------------------------------------------------
master_phecodes <- codes_master[str_detect(codes_master$feature_id, 'PheCode'),]
master_phecodes$code <- str_split(master_phecodes$feature_id, '\\:', simplify = T)[,2]
master_phecodes$code <- gsub("\\.","\\_", master_phecodes$code)

phecodes <- names(CC_comb)[str_detect(names(CC_comb), 'PheCode')]
phecodes <- str_split(phecodes, '\\.', simplify = T)[,2]

# Manually compare --------------------------------------------------------
# Aaron's list
master_phecodes$code %>% unique
# PheCodes in data set 
phecodes

# Remove ------------------------------------------------------------------
# Ultimately opt to remove codes 1005_, 1010_, 773_, 563_ 
CC_comb <- CC_comb %>% select(-PheCode.1005_, -PheCode.1010_, -PheCode.773_, -PheCode.563_)
CC_comb_val <- CC_comb_val %>% select(-PheCode.1005_, -PheCode.1010_, -PheCode.773_, -PheCode.563_)

# Save --------------------------------------------------------------------
saveRDS(CC_comb, 'intermediate_data/CC_comb_24_3_clin_clean.rds')
saveRDS(CC_comb_val, 'intermediate_data/CC_comb_val_24_3_clin_clean.rds')








