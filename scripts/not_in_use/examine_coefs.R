# Examine coefficients
# Library -----------------------------------------------------------------
library(tidyverse)
library(glmnet)
library(openxlsx)
library(lubridate)

# Load --------------------------------------------------------------------

CC_comb <- readRDS('intermediate_data/CC_comb_24_3_clin_clean.rds') # Training data
CC_comb_val <- readRDS('intermediate_data/CC_comb_val_24_3_clin_clean.rds') # Validation data
MS_attack <- read.csv("raw_data/Box/Boston/MS CLIMB Paper2 Code&Result/Step3_BuildModel/IntermediateDataResult/Cleaned_MS_attack.csv",
                      stringsAsFactors = FALSE,
                      colClasses = c("character", rep("integer",6), "Date", rep("numeric",3),
                                     rep("character",12), rep("Date",2), "integer",
                                     rep("Date",5), rep("integer",2)))
df <- readRDS('intermediate_data/3partclin_probs.rds')
load('models/glmnet.list.rda')
# Load Code databases
wkpath  = "raw_data/Box/Boston/MS CLIMB Data/"
wkpath2 = "raw_data/Box/Boston/MS CLIMB Paper2 Code&Result/Step3_BuildModel/IntermediateDataResult/"
CUIdictAll = read.xlsx(paste0(wkpath,"EHR/AllCUI_Database.xlsx"),
                       sheet = 1); colnames(CUIdictAll) = c("ConceptCd","Desc")
ICDPheCode   = read.csv(paste0(wkpath,"EHR/MS_AllEncounters_ICD_Data_03282019.csv"),
                        stringsAsFactors = FALSE);ICDPheCode$phecode[ICDPheCode$concept_cd == "LPA268"] = "335_"

# Convert relevant columns to Date format
CC_comb$StartDate <- as.Date(CC_comb$StartDate)
CC_comb_val$StartDate <- as.Date(CC_comb_val$StartDate)

# Set parameters 
tw = 24; tp = 3

# Determine Nonzero Betas -------------------------------------------------

coefs.list <- map(glmnet.list, function(mdl) {
  coefs <- coef(mdl, s = "lambda.1se") %>% as.matrix()
  coefs <- data.frame(term = row.names(coefs),
                      beta = coefs[,1])
  row.names(coefs) <- NULL
  return(coefs)
})
nonzero.coefs.list <- map(coefs.list, function(tbl) {tbl %>% filter(beta != 0)})


# Obtain intersections ----------------------------------------------------
# 2000-2005 to 2006-2010 model
int1 <- intersect(nonzero.coefs.list[[1]]$term, nonzero.coefs.list[[2]]$term)
dropped1 <- nonzero.coefs.list[[2]]$term[!(nonzero.coefs.list[[2]]$term %in% int1)] %>% as.character()

# Get code descriptions
tmp1 <- str_detect(dropped1, 'CUI')
data.frame(term = dropped1, 
           tmp1 = str_split(dropped1, '\\.'))

tmp1 <- str_split(dropped1, '\\.', simplify = T) %>% as.data.frame()
cuis <- tmp1[tmp1$V1 == 'CUI', 'V2'] %>% as.character() %>% as.vector()
desc_cui <- rep(NA, length(cuis))
for (i in 1:length(cuis)) {
  desc_cui[i] <- CUIdictAll[CUIdictAll$ConceptCd == cuis[i], 'Desc'][1]
}
phes <- tmp1[tmp1$V1 == 'PheCode', 'V2'] %>% as.character() %>% as.vector()
desc_phes <- rep(NA, length(phes))
for (i in 1:length(phes)) {
  desc_phes[i] <- ICDPheCode[ICDPheCode$phecode == phes[i], 'phecode_description'] %>% unique
}
tmp1$V3 <- c(rep(NA, 8), desc_cui, desc_phes)
colnames(tmp1) <- c('Group', 'Code', 'Desc')


# 2006-2010 to 2011-2016 model
int2 <- intersect(nonzero.coefs.list[[2]]$term, nonzero.coefs.list[[3]]$term)
dropped2 <- nonzero.coefs.list[[3]]$term[!(nonzero.coefs.list[[3]]$term %in% int2)] %>% as.character()

# Get code descriptions
tmp1 <- str_split(dropped2, '\\.', simplify = T) %>% as.data.frame()
# Remove FEMALE and RACE for now
tmp1 <- tmp1[1:34,]
# Resume 
cuis <- tmp1[tmp1$V1 == 'CUI', 'V2'] %>% as.character() %>% as.vector()
desc_cui <- rep(NA, length(cuis))
for (i in 1:length(cuis)) {
  desc_cui[i] <- CUIdictAll[CUIdictAll$ConceptCd == cuis[i], 'Desc'][1]
}
phes <- tmp1[tmp1$V1 == 'PheCode', 'V2'] %>% as.character() %>% as.vector()
desc_phes <- rep(NA, length(phes))
for (i in 1:length(phes)) {
  desc_phes[i] <- ICDPheCode[ICDPheCode$phecode == phes[i], 'phecode_description'] %>% unique
}
tmp1$V3 <- c(rep(NA, 7), desc_cui, desc_phes)
colnames(tmp1) <- c('Group', 'Code', 'Desc')
tmp1 <- rbind(tmp1, c('FEMALE', NA, NA), c('RACE', NA, NA))

# Standardize to SD = 1 ---------------------------------------------------

scaled_X <- scale(CC_comb_val[,-c(1:9)]) %>% as.data.frame
scaled_X$PatientID <- CC_comb_val$PatientID
scaled_X$StartDate <- CC_comb_val$StartDate




