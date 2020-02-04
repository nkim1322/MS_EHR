# DMT data
# Library -----------------------------------------------------------------
library(tidyverse)
library(zoo)
library(openxlsx)
library(lubridate)
library(glmnet)

# Load --------------------------------------------------------------------
wkpath  = "raw_data/Box/Boston/MS CLIMB Data/"
wkpath2 = "raw_data/Box/Boston/MS CLIMB Paper2 Code&Result/Step3_BuildModel/IntermediateDataResult/"

MS_cohort = read.csv(paste0(wkpath2,"Cleaned_MS_cohort.csv"), stringsAsFactors = FALSE,
                     colClasses = c(rep("character",4),rep("Date",4),rep("integer",6),
                                    rep("numeric",3),"integer",rep("Date",4),
                                    rep("numeric",3),"integer"))
MS_map  = read.xlsx(paste0(wkpath,"CLIMB Cohort/Spec95_i2b2_Mapping2017.xlsx"), sheet = 1)
colnames(MS_map) = c("PatientNum","PatientID")
MS_EHR    = read.xlsx(paste0(wkpath,"EHR/MS_first_last_date_whole_list.xlsx"),
                      sheet = 1); colnames(MS_EHR)[1] = "PatientNum"

MS_EHR$First_date <- convertToDate(MS_EHR$First_date)
MS_EHR$Last_date <- convertToDate(MS_EHR$Last_date)
MS_EHR$PatientID = MS_map$PatientID[match(MS_EHR$PatientNum,MS_map$PatientNum)]

RXNORM_cohort <- read.xlsx(paste0(wkpath,"MS DMT/MS_RxNorm_data.xlsx"), sheet = 4)
RXNORM_cohort$start_date <- convertToDate(RXNORM_cohort$start_date)
RXNORM_EHR <- read.xlsx(paste0(wkpath,"MS DMT/MS_RxNorm_data.xlsx"), sheet = 3)
RXNORM_EHR$start_date <- convertToDate(RXNORM_EHR$start_date)


# RX Norm to Drug mapping -------------------------------------------------
rx_drug_map <- data.frame(drug = c("rituximab","natalizumab","dimethyl_fumarate","fingolimod"),
                          RX_Norm_code = c(121191, 354770, 1373478, 1012892))
# Create new columns
RXNORM_cohort$drug_name <- sapply(RXNORM_cohort$RX_Norm_code, function(i) {
  rx_drug_map$drug[rx_drug_map$RX_Norm_code == i]
})
RXNORM_EHR$drug_name <- sapply(RXNORM_EHR$RX_Norm_code, function(i) {
  rx_drug_map$drug[rx_drug_map$RX_Norm_code == i]
})



# Codes data --------------------------------------------------------------
ICD_CPT_CUI_Comb = read.csv(paste0(wkpath2,"ICD_CPT_CUI_Comb.csv"), stringsAsFactors = FALSE)
ICD_CPT_CUI_Comb$StartDate = as.Date(ICD_CPT_CUI_Comb$StartDate)
ICD_CPT_CUI_Comb$ConceptCd = as.factor(ICD_CPT_CUI_Comb$ConceptCd)

ICDPheCode   = read.csv(paste0(wkpath,"EHR/MS_AllEncounter_ICDCodes_08072019.csv"), stringsAsFactors = FALSE)
ICDPheCode$phecode[ICDPheCode$concept_cd == "LPA268"] = "335_"
ICDPheCode$start_date  = as.Date(ICDPheCode$start_date,format = "%m/%d/%y")
CUIdictAll = read.xlsx(paste0(wkpath,"EHR/AllCUI_Database.xlsx"), sheet = 1)
colnames(CUIdictAll) = c("ConceptCd","Desc")

# Subset for EHR patients only
ICD_CPT_CUI_Comb <- ICD_CPT_CUI_Comb[ICD_CPT_CUI_Comb$PatientNum%in%MS_EHR$PatientNum, ]
# Further subset CLIMB vs non-CLIMB
ICD_CPT_CUI_Comb_cohort = ICD_CPT_CUI_Comb[ICD_CPT_CUI_Comb$PatientID%in%MS_cohort$PatientID,]
ICD_CPT_CUI_Comb_UL     = ICD_CPT_CUI_Comb[!ICD_CPT_CUI_Comb$PatientID%in%MS_cohort$PatientID,]

# Use first codified Rx of meds as start date -----------------------------
trt_patients <- unique(RXNORM_EHR$`patient_num/i2b2`[RXNORM_EHR$drug_name %in% c("rituximab","natalizumab")])

# Get codes for DMTs
# ritux_cui <- c(paste0("CUI.", CUIdictAll$ConceptCd[str_detect(CUIdictAll$Desc, "rituximab")]),
#                paste0("CUI.", CUIdictAll$ConceptCd[str_detect(CUIdictAll$Desc, "Rituxan")]))
# nataliz_cui <- c(paste0("CUI.", CUIdictAll$ConceptCd[str_detect(CUIdictAll$Desc, "natalizumab")]),
#                  paste0("CUI.", CUIdictAll$ConceptCd[str_detect(CUIdictAll$Desc, "Tysabri")]))
# 
# # Get patients with either of these codes: 
# trt_patients <- unique(ICD_CPT_CUI_Comb_UL$PatientNum[ICD_CPT_CUI_Comb_UL$ConceptCd %in% c(ritux_cui, nataliz_cui)])
# ritux_pts <- unique(ICD_CPT_CUI_Comb_UL$PatientNum[ICD_CPT_CUI_Comb_UL$ConceptCd %in% ritux_cui])
# nataliz_pts <- unique(ICD_CPT_CUI_Comb_UL$PatientNum[ICD_CPT_CUI_Comb_UL$ConceptCd %in% nataliz_cui])

# overlap <- intersect(ritux_pts, nataliz_pts)
# overlap %>% length

trt_subset <- data.frame(PatientNum = trt_patients)
ICD_CPT_CUI_Comb_UL_trt <- ICD_CPT_CUI_Comb %>% filter(PatientNum %in% trt_patients)

# Get treatment assignment and start date ---------------------------------
trt_subset$medication_desc <- ""
trt_subset$start_date <- ""

for (p_num in trt_subset$PatientNum) {
  p_subset <- RXNORM_EHR %>% filter(`patient_num/i2b2` == p_num)
  p_subset <- p_subset[order(p_subset$start_date),]
  p_ritux <- p_subset %>% filter(drug_name == "rituximab")
  p_nataliz <- p_subset %>% filter(drug_name == "natalizumab")
  # For patients only on one drug
  if (nrow(p_nataliz) == 0) {
    trt_subset[trt_subset$PatientNum == p_num, c("medication_desc","start_date")] <- 
      c("RITUXAN", as.character(p_ritux[1,"start_date"]))
  } else if (nrow(p_ritux) == 0) {
    trt_subset[trt_subset$PatientNum == p_num, c("medication_desc","start_date")] <- 
      c("TYSABRI", as.character(p_nataliz[1,"start_date"]))
  } else {
    # For patients ever on both drugs: assign treatment based on whichever drug comes first
    if (as.Date(p_ritux[1, "start_date"]) >= as.Date(p_nataliz[1, "start_date"])) {
      trt_subset[trt_subset$PatientNum == p_num, c("medication_desc","start_date")] <- 
        c("RITUXAN", as.character(p_ritux[1,"start_date"]))
    } else {
      trt_subset[trt_subset$PatientNum == p_num, c("medication_desc","start_date")] <- 
        c("TYSABRI", as.character(p_nataliz[1,"start_date"]))
    }
  }
}
trt_subset$start_date <- as.Date(trt_subset$start_date)


# Not in use: use CUIs to define treatment start -----------------------

# for (p_num in trt_subset$PatientNum) {
#   p_subset <- ICD_CPT_CUI_Comb_UL_trt %>% filter(PatientNum == p_num)
#   p_ritux <- p_subset %>% filter(ConceptCd %in% ritux_cui)
#   p_ritux <- p_ritux[order(p_ritux$StartDate),]
#   p_nataliz <- p_subset %>% filter(ConceptCd %in% nataliz_cui)
#   p_nataliz <- p_nataliz[order(p_nataliz$StartDate),]
#   # For patients only on one drug
#   if (nrow(p_nataliz) == 0) {
#     trt_subset[trt_subset$PatientNum == p_num, "medication_desc"] <- "RITUXAN"
#     trt_subset[trt_subset$PatientNum == p_num, "start_date"] <- p_ritux[1, "StartDate"] %>% as.character()
#   } else if (nrow(p_ritux) == 0) {
#     trt_subset[trt_subset$PatientNum == p_num, "medication_desc"] <- "TYSABRI"
#     trt_subset[trt_subset$PatientNum == p_num, "start_date"] <- p_nataliz[1, "StartDate"] %>% as.character()
#   } else {
#     # For patients ever on both drugs: assign treatment based on whichever drug comes first
#     if (as.Date(p_ritux[1, "StartDate"]) >= as.Date(p_nataliz[1, "StartDate"])) {
#       trt_subset[trt_subset$PatientNum == p_num, "medication_desc"] <- "RITUXAN"
#       trt_subset[trt_subset$PatientNum == p_num, "start_date"] <- p_ritux[1, "StartDate"] %>% as.character()
#     } else {
#       trt_subset[trt_subset$PatientNum == p_num, "medication_desc"] <- "RITUXAN"
#       trt_subset[trt_subset$PatientNum == p_num, "start_date"] <- p_nataliz[1, "StartDate"] %>% as.character()
#     }
#   }
# }
# trt_subset$start_date <- as.Date(trt_subset$start_date)


# Add clinical vars -------------------------------------------------------
# Add PatientID for later use
trt_subset$PatientID <- sapply(trt_subset$PatientNum, function(p_num) {
  MS_map$PatientID[MS_map$PatientNum == p_num][1]
})
trt_subset <- trt_subset %>% select(PatientID, PatientNum, everything())

# Add sex and race
MS_demographics <- read.xlsx(paste0(wkpath, "EHR/ms_demographics.xlsx"))
MS_demographics <- MS_demographics %>% 
  dplyr::select(PATIENT_NUM, BIRTH_DATE, SEX_CD, RACE_CD)
MS_demographics$BIRTH_DATE <- convertToDate(MS_demographics$BIRTH_DATE)

trt_subset$FEMALE <- ifelse(trt_subset$PatientNum %in% MS_demographics$PATIENT_NUM[MS_demographics$SEX_CD == 'F'], 1, 0)
trt_subset$RACE <- ifelse(trt_subset$PatientNum %in% MS_demographics$PATIENT_NUM[MS_demographics$RACE_CD == 'WHITE'], 1, 0)

# Age at first MS code (MS ICD: 'PheCode.335_')
trt_subset$FIRSTMSICD_DATE <- sapply(trt_subset$PatientNum, function(pnum) {
  tmp <- ICDPheCode %>% filter(patient_num == pnum & phecode == '335_')
  tmp <- tmp[order(tmp$start_date), ]
  return(tmp$start_date[1])
})
trt_subset$FIRSTMSICD_DATE <- as.Date(trt_subset$FIRSTMSICD_DATE)

# Exclude those with missing FIRSTMSICD_DATE (who also have no codes/encounters data)
trt_subset <- trt_subset %>% filter(!is.na(FIRSTMSICD_DATE))

trt_subset$BIRTHDAY <- sapply(trt_subset$PatientNum, function(p_num) {
  MS_demographics$BIRTH_DATE[MS_demographics$PATIENT_NUM == p_num] %>% as.character()
})
trt_subset$BIRTHDAY <- as.Date(trt_subset$BIRTHDAY)
trt_subset$AGE_AT_FIRSTMSICD <- with(trt_subset, difftime(FIRSTMSICD_DATE,BIRTHDAY,units='weeks')/52.25) %>% as.numeric()

# Follow-up duration: (in years)
## time period between date of occurrence of first of any ICD code in the EHR and time of treatment initiation
trt_subset$FIRSTICD_DATE <- sapply(trt_subset$PatientNum, function(pnum) {
  tmp <- ICDPheCode %>% filter(patient_num == pnum)
  tmp <- tmp[order(tmp$start_date),]
  return(tmp$start_date[1])
})
trt_subset$FIRSTICD_DATE <- as.Date(trt_subset$FIRSTICD_DATE)
trt_subset$FOLLOWUP_DURA <- with(trt_subset, difftime(start_date, FIRSTICD_DATE, units='weeks')/52.25) %>% as.numeric()

# Healthcare utilization, overall: 
## total number of all ICD codes during entire f/u duration
trt_subset$HCUTIL_OVERALL <- sapply(trt_subset$PatientNum, function(pnum) {
  tmp <- ICDPheCode %>% filter(patient_num == pnum)
  tmp <- tmp[order(tmp$start_date),]
  tmp <- tmp %>% filter(start_date >= trt_subset$FIRSTICD_DATE[trt_subset$PatientNum == pnum] & 
                          start_date <= trt_subset$start_date[trt_subset$PatientNum == pnum])
  return(nrow(tmp))
})

# Healthcare utilization, in the 3* months prior to DMT start
trt_subset$HCUTIIL_3MONS <- sapply(trt_subset$PatientNum, function(pnum) {
  tmp <- ICDPheCode %>% filter(patient_num == pnum)
  tmp <- tmp[order(tmp$start_date),]
  tmp <- tmp %>% 
    filter(start_date >= trt_subset$start_date[trt_subset$PatientNum == pnum] - months(3) & 
             start_date <= trt_subset$start_date[trt_subset$PatientNum == pnum])
  return(nrow(tmp))
})

# Adjusted frequency of ICD code for MS, in 3* months prior to DMT start: 
## number of ICD code for MS divided by healthcare utilization
trt_subset$MSICDFREQ_ADJ_3MONS <- sapply(trt_subset$PatientNum, function(pnum) {
  tmp <- ICDPheCode %>% filter(patient_num == pnum)
  tmp <- tmp[order(tmp$start_date),]
  tmp <- tmp %>% 
    filter(start_date >= trt_subset$start_date[trt_subset$PatientNum == pnum] - months(3) & 
             start_date <= trt_subset$start_date[trt_subset$PatientNum == pnum])
  if (nrow(tmp) != 0) {
    x <- (tmp %>% filter(phecode == '335_') %>% nrow)/nrow(tmp)
    return(x)
  } else {
    return(0)
  }
})

# Adjusted frequency of CUI code for “multiple sclerosis” (MS: 'CUI.C0026769'), in the 3* months prior to DMT start: 
## number of CUI code for MS divided by healthcare utilization
CUISelected  = read.csv(paste0(wkpath,"EHR/MS_AllEncounters_CUIs_08072019.csv"), stringsAsFactors = FALSE, 
                        col.names = c('patient_num', 'encounter_num', 'start_date', 'concept_cd'))
CUISelected$start_date <- format(as.POSIXct(CUISelected$start_date, format = '%Y-%m-%d %H:%M:%S'), format = '%Y-%m-%d')
trt_subset$MSCUIFREQ_ADJ_3MONS <- sapply(trt_subset$PatientNum, function(pnum) {
  tmp <- CUISelected %>% filter(patient_num == pnum)
  tmp <- tmp[order(tmp$start_date),]
  tmp <- tmp %>% 
    filter(start_date >= trt_subset$start_date[trt_subset$PatientNum == pnum] - months(3) & 
             start_date <= trt_subset$start_date[trt_subset$PatientNum == pnum])
  if (nrow(tmp) != 0) {
    x <- (tmp %>% filter(concept_cd == 'C0026769') %>% nrow)/nrow(tmp)
    return(x)
  } else {
    return(0)
  }
  
})


# Duration of prior DMT treatment (in months)
## (Is this a valid way to get this info for EHR patients?)
dmt_cui <- c("C0058218", "C3556178", "C1699926", "C2938762", "C0717787", "C0289884", "C0528175", 
             "C4027077", "C0015980", "C0244713", "C0254119", "C0284968", "C0594372", "C0752980", 
             "C2719461", "C3848580", "C3848664", "C1718383", "C3497721")
dmt_cui <- paste0("CUI.", dmt_cui)

trt_subset$PRIORDMT_DURA <- sapply(trt_subset$PatientNum, function(pnum) {
  tmp <- ICD_CPT_CUI_Comb_UL_trt %>% 
    filter(PatientNum == pnum) %>% 
    filter(ConceptCd %in% dmt_cui)
  tmp <- tmp[order(tmp$StartDate),]
  if (nrow(tmp) == 0) {
    return(0)
  } else {
    # Calculate difference between first and last DMT code mention
    tmp2 <- difftime(tmp$StartDate[length(tmp$StartDate)], tmp$StartDate[1], units = 'weeks')/4
    return(as.numeric(tmp2))
  }
})

# Disease duration: (in years)
## age at the time of DMT start - age at first MS ICD code
trt_subset$DISEASE_DURA <- with(trt_subset, as.numeric(difftime(start_date, BIRTHDAY)/365)) - 
  trt_subset$AGE_AT_FIRSTMSICD


# Response ----------------------------------------------------------------
aggregate = FALSE
model = FALSE

# For preprocessing imputation data for EHR patients
if (aggregate) {
  tw = 24; tp = 3
  ### Aggregate codes
  # MS_Enct_Uniq = read.csv(paste0("raw_data/HSPH/Data From Partners MS Center/MS_Encounter_Date_PatientID.csv"), stringsAsFactors = FALSE)
  # MS_Enct_Uniq = MS_Enct_Uniq[MS_Enct_Uniq$PatientNum%in%trt_patients,]
  # MS_Enct_Uniq$StartDate = as.Date(MS_Enct_Uniq$StartDate,format = "%m/%d/%y")
  # tmp = plyr::count(MS_Enct_Uniq[,c("PatientNum","PatientID","StartDate")])[,-4]
  tmp = ICD_CPT_CUI_Comb_UL_trt[,c("PatientNum",  "PatientID", "StartDate")]
  tmp = plyr::count(tmp)[,-4]
  tmp = tmp[order(tmp$PatientNum,tmp$StartDate),]
  tmp2 = sapply(1:nrow(tmp),function(i){
    if (i %% 100 == 0) {
      print(i)
    }
    aa = ICD_CPT_CUI_Comb_UL_trt[ICD_CPT_CUI_Comb_UL_trt$PatientNum == tmp$PatientNum[i],]
    bb = difftime(aa$StartDate,tmp$StartDate[i],units = "days")/30
    aa = aa[bb < tw & bb >= 0,]
    aa = table(aa$ConceptCd);names(aa) = NULL
    aa
  })
  tmp2 = t(tmp2)
  colnames(tmp2) = levels(ICD_CPT_CUI_Comb_UL_trt$ConceptCd)
  EHR_comb = cbind(tmp,tmp2)
  EHR_comb$StartDate <- as.Date(EHR_comb$StartDate)
  # order by PatientID & StartDate
  EHR_comb = EHR_comb[order(EHR_comb$PatientID,EHR_comb$StartDate),]
  # Add First and last date as first and last encounter in EHR
  EHR_comb$FIRST_DATE <- NA; EHR_comb$LAST_DATE <- NA
  for (pnum in unique(EHR_comb$PatientNum)) {
    tmp <- EHR_comb %>% filter(PatientNum == pnum)
    EHR_comb[EHR_comb$PatientNum == pnum, "FIRST_DATE"] <- as.character(tmp$StartDate[1])
    EHR_comb[EHR_comb$PatientNum == pnum, "LAST_DATE"] <- as.character(tmp$StartDate[length(tmp$StartDate)])
  }
  EHR_comb$FIRST_DATE <- as.Date(EHR_comb$FIRST_DATE)
  EHR_comb$LAST_DATE <- as.Date(EHR_comb$LAST_DATE)
  EHR_comb$Period = as.numeric(difftime(EHR_comb$StartDate,EHR_comb$FIRST_DATE,units="days")/30/tp)
  EHR_comb = EHR_comb[EHR_comb$Period>=0,]
  EHR_comb$Period = ceiling(EHR_comb$Period)
  EHR_comb <- EHR_comb %>% 
    dplyr::select(PatientNum, PatientID, StartDate, FIRST_DATE, LAST_DATE, Period, everything())
  
  # Save intermediate data
  write.csv(EHR_comb,
            paste0("intermediate_data/EHR_ICD_CPT_CUI_Count_",tw,"mons_",tp,"mons.csv"),
            row.names = FALSE)
  # EHR_comb <- read.csv(paste0("intermediate_data/EHR_ICD_CPT_CUI_Count_",tw,"mons_",tp,"mons.csv"))
  
  
  
  ### Order CPT groups by number, deal with discrepency due to special characters 
  CPTgroup = read.csv(paste0(wkpath2,"MS_ConceptCode_Mapping.csv"),stringsAsFactors = FALSE)
  tmp = grep("CPTGroup.", names(EHR_comb))
  aa = names(EHR_comb)[tmp]
  aa = substr(aa, 10, nchar(aa))
  aa = gsub("\\."," ",aa)
  aa = gsub("  "," ",aa)
  bb = CPTgroup$cpt_group
  bb = gsub("[-/(),]"," ",bb)
  bb = gsub("  "," ",bb)
  names(EHR_comb)[tmp] = paste0("CPTGroup.", CPTgroup$cpt_group[match(aa,bb)])
  
  
  ## Add clinical variables
  # Sex
  tmp <- trt_subset[trt_subset$FEMALE == 1,'PatientID']
  EHR_comb$FEMALE <- ifelse(EHR_comb$PatientNum %in% tmp, 1, 0)
  # Race
  tmp <- trt_subset[trt_subset$RACE == 1,'PatientID']
  EHR_comb$RACE <- ifelse(EHR_comb$PatientNum %in% tmp, 1, 0)
  # Follow-up duration
  EHR_comb$FOLLOWUP_DURA <- NA 
  for (pnum in unique(EHR_comb$PatientNum)) {
    tmp <- trt_subset[trt_subset$PatientNum == pnum, 'FOLLOWUP_DURA']
    if (length(tmp) == 0) {
      print("No followup duration")
    } else {
      EHR_comb[EHR_comb$PatientNum == pnum, 'FOLLOWUP_DURA'] <- tmp
    }
  }
  # Disease duration proxy: duration between encounter date and 1st MS ICD Code (in years)
  EHR_comb$DISEASE_DURA <- NA
  for (i in 1:nrow(EHR_comb)) {
    first_ms_icd_date <- trt_subset[trt_subset$PatientNum == EHR_comb$PatientNum[i], "FIRSTMSICD_DATE"]
    tmp <- as.numeric(difftime(EHR_comb$StartDate[i], first_ms_icd_date)/365)
    if (length(tmp) == 0) {
      print("No disease duration")
    } else {
      EHR_comb$DISEASE_DURA[i] <- tmp
    }
  }
  
  ### Read in CLIMB data for code mapping
  CC_comb <- readRDS(paste0("intermediate_data/CC_comb_", tw, "_", tp, ".rds"))
  # Save codes only in CLIMB dataset
  EHR_comb <- EHR_comb[,intersect(names(CC_comb), names(EHR_comb))]
  EHR_comb$StartDate <- as.Date(EHR_comb$StartDate)
  
  # Remove rows with missing data
  # EHR_comb <- EHR_comb[!is.na(EHR_comb$DISEASE_DURA) & !is.na(EHR_comb$FOLLOWUP_DURA),]
  
  ## Data processing
  # Log transform
  # Do not transform ID columns and clinical vars 
  no_transform <- c(1:6, (length(names(EHR_comb)) - 3):length(names(EHR_comb)))
  EHR_comb[,-no_transform] = log(1 + EHR_comb[,-no_transform])
  
  # Remove data before 2006?
  # EHR_comb <- EHR_comb %>% filter(StartDate >= as.Date('2006-01-01'))
  
  # Save
  saveRDS(EHR_comb, "modeling_data/causalRN_impute_24_3_EHR.rds")
} else {
  EHR_comb <- readRDS("modeling_data/causalRN_impute_24_3_EHR.rds")
}

if (model) {
  # Fit model using CLIMB
  train <- readRDS("modeling_data/train24_3.rds")
  test <- readRDS("modeling_data/test24_3.rds")
  comb <- rbind(train, test)
  # Remove all data before 2006??
  comb <- comb %>% filter(StartDate >= as.Date('2006-01-01'))
  # Fit Lasso model
  lasso <- cv.glmnet(as.matrix(all[,-c(1:9)]), all$CC, family = 'binomial', type.measure = 'auc')
  save(lasso, file = 'models/lasso_imputation.rda')
  rm(train, test, all)
} else {
  load("models/lasso_imputation.rda")
}

# Only keep patients for which we have EHR data
# trt_subset <- trt_subset[trt_subset$PatientNum %in% unique(EHR_comb$PatientNum), ]

# Predict probability of relapse
trt_subset$PROB_RELAPSE <- sapply(trt_subset$PatientNum, function(pnum) {
  start_date <- trt_subset$start_date[trt_subset$PatientNum == pnum]
  pred_data <- EHR_comb %>% filter(PatientNum == pnum & StartDate == start_date)
  if (nrow(pred_data) == 0) {
    pred_data <- EHR_comb %>% 
      filter(PatientNum == pnum) %>% 
      filter(StartDate >= start_date - months(1) & StartDate <= start_date + months(1))
  }
  if (nrow(pred_data) == 0) {
    pred_data <- EHR_comb %>% 
      filter(PatientNum == pnum) %>% 
      filter(StartDate >= start_date - months(3) & StartDate <= start_date + months(3))
  }
  
  pred <- predict(lasso, newx = as.matrix(pred_data[1,-c(1:6)]),
                  type = "response", s = lasso$lambda.1se)
  return(pred)
})

# Remove patients without predicted relapse probability
trt_subset <- trt_subset[!is.na(trt_subset$PROB_RELAPSE),]


# Test --------------------------------------------------------------------
tmp <- sapply(trt_subset$PatientNum, function(pnum) {
  start_date <- trt_subset$start_date[trt_subset$PatientNum == pnum]
  pred_data <- EHR_comb %>% filter(PatientNum == pnum & StartDate == start_date)
  ctr = 0
  if (nrow(pred_data) == 0) {
    pred_data <- EHR_comb %>%
      filter(PatientNum == pnum) %>%
      filter(StartDate >= start_date - months(1) & StartDate <= start_date + months(1))
    # filter(abs(difftime(StartDate, start_date, units = "weeks")/4) < 1)
    ctr = ctr + 1
  }

  if (nrow(pred_data) == 0) {
    pred_data <- EHR_comb %>%
      filter(PatientNum == pnum) %>%
      filter(StartDate >= start_date - months(3) & StartDate <= start_date + months(3))
    # filter(abs(difftime(StartDate, start_date, units = "weeks")/4) < 1)
    ctr = ctr + 1
  }
  if (nrow(pred_data) == 0) {
    return(3)
  } else {
    return(ctr)
  }
})

# Save --------------------------------------------------------------------
## Move all non-modeling variables to first few columns of dataframe
trt_subset <- trt_subset %>% select(PatientID, PatientNum, medication_desc, start_date, 
                                    BIRTHDAY, FIRSTICD_DATE, FIRSTMSICD_DATE, PROB_RELAPSE,
                                    everything())

# First 7 columns are non-modeling, 8th column is outcome
saveRDS(trt_subset, 'modeling_data/causalRN_3mons_EHR.rds')

rm(list = ls())


# ms_demographics ---------------------------------------------------------
# MS_demographics <- read.xlsx(paste0(wkpath, "EHR/ms_demographics.xlsx"), sheet = 1)
# MS_demographics <- MS_demographics %>% select(PATIENT_NUM, BIRTH_DATE); colnames(MS_demographics)[1] <- 'PatientNum'
# MS_demographics$BIRTH_DATE <- convertToDate(MS_demographics$BIRTH_DATE)
# MS_demographics$PatientID <- sapply(MS_demographics$PatientNum, function(pnum) {
#   MS_map$PatientID[MS_map$PatientNum == pnum][1]
# })
# trt_subset$BIRTHDAY <- sapply(trt_subset$PatientNum, function(pnum) {
#   MS_demographics$BIRTH_DATE[MS_demographics$PatientNum == pnum]
# })













