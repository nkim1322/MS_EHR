# DMT data
# Library -----------------------------------------------------------------
library(tidyverse)
library(openxlsx)
library(zoo)
library(lubridate)
library(glmnet)

# Load --------------------------------------------------------------------
wkpath  = "raw_data/Box/Boston/MS CLIMB Data/"
wkpath2 = "raw_data/Box/Boston/MS CLIMB Paper2 Code&Result/Step3_BuildModel/IntermediateDataResult/"
MS_cohort = read.csv(paste0(wkpath2,"Cleaned_MS_cohort.csv"), stringsAsFactors = FALSE,
                     colClasses = c(rep("character",4),rep("Date",4),rep("integer",6),
                                    rep("numeric",3),"integer",rep("Date",4),
                                    rep("numeric",3),"integer"))
MS_attack = read.csv(paste0(wkpath2,"Cleaned_MS_attack.csv"), stringsAsFactors = FALSE,
                     colClasses = c("character",rep("integer",6),"Date",rep("numeric",3),
                                    rep("character",12),rep("Date",2),"integer",
                                    rep("Date",5),rep("integer",2)))
MS_trt    = read.csv(paste0(wkpath2,"Cleaned_MS_treatment.csv"), stringsAsFactors = FALSE)
MS_trt$start_date <- as.Date(with(MS_trt, paste(start_year, start_month, start_day,sep="-")), "%Y-%m-%d")
MS_trt$stop_date <- as.Date(with(MS_trt, paste(stop_year, stop_month, stop_day,sep="-")), "%Y-%m-%d")

MS_map  = read.xlsx(paste0(wkpath,"CLIMB Cohort/Spec95_i2b2_Mapping2017.xlsx"), sheet = 1)
colnames(MS_map) = c("PatientNum","PatientID")



# First line drugs --------------------------------------------------------
inf_beta <- c('AVONEX','BETASERON','PLEGRIDY','REBIF')
glat_acet <- c('COPAXONE')
# first_line_pts <- MS_trt %>% filter(medication_desc %in% c(inf_beta, glat_acet)) %>% 
#   select(PatientID) %>% unlist(use.names = F) %>% unique()

trt_subset <- MS_trt
# Begin tidy --------------------------------------------------------------
# Coerce dates to format
trt_subset$start_date <- as.Date(with(trt_subset, paste(start_year, start_month, start_day,sep="-")), "%Y-%m-%d")
trt_subset$stop_date <- as.Date(with(trt_subset, paste(stop_year, stop_month, stop_day,sep="-")), "%Y-%m-%d")
trt_subset <- trt_subset[,-c(3:8)]
trt_subset <- trt_subset %>% dplyr::select(PatientID, medication_desc, start_date, stop_date, everything())

# Add treatment duration in months
trt_subset$tx_dura <- with(trt_subset, difftime(stop_date,start_date,units = 'weeks')/4) %>% as.numeric()
# Remove patients with treatment duration = 0
trt_subset <- trt_subset %>% filter(tx_dura != 0)

MS_trt_subset <- trt_subset

# Dimethyl fumarate vs fingolimod -----------------------------------------
# dmf_pts <- MS_trt %>% filter(medication_desc == 'TECFIDERA (BG-12)') %>%
#   select(PatientID) %>% unlist(use.names = F) %>% unique()
# 
# fingo_pts <- MS_trt %>% filter(medication_desc == 'GILENYA') %>%
#   select(PatientID) %>% unlist(use.names = F) %>% unique()
# 
# intersect(dmf_pts, fingo_pts) %>% length

# Rituximab vs natalizumab ------------------------------------------------
ritux_pts <- trt_subset %>% filter(medication_desc == 'RITUXAN') %>%
  dplyr::select(PatientID) %>% unlist(use.names = F) %>% unique()

nataliz_pts <- trt_subset %>% filter(medication_desc == 'TYSABRI') %>%
  dplyr::select(PatientID) %>% unlist(use.names = F) %>% unique()

# overlap <- intersect(nataliz_pts, ritux_pts)
# overlap %>% length

# Define treatment status -------------------------------------------------
# patients <- unique(trt_subset$PatientID)
patients <- union(ritux_pts, nataliz_pts)

options <- rep(NA, length(patients)) # keep track of what options are applying to patients
drug_assign <- data.frame(matrix(ncol = 5, nrow = length(patients)))
keep_cols <- c('PatientID', 'medication_desc', 'start_date', 'stop_date', 'tx_dura')
# Reformat start and stop dates later
trt_subset$start_date <- as.character(trt_subset$start_date)
trt_subset$stop_date <- as.character(trt_subset$stop_date)
colnames(drug_assign) <- keep_cols
for (i in 1:length(patients)) {
  p_id <- patients[i]
  p_subset <- trt_subset %>% filter(PatientID == p_id) %>% plyr::arrange(start_date)
  first_line <- p_subset %>% filter(medication_desc %in% c(inf_beta, glat_acet))
  
  # If patient was never on a first line drug, assign patient to first of 2 drugs patient was on
  if (nrow(first_line) == 0) {
    drug_assign[i,] <- (p_subset %>% filter(medication_desc %in% c('RITUXAN','TYSABRI')))[1,keep_cols]
    options[i] <- 'Option 1'
  } else if (p_subset %>% 
             filter(medication_desc %in% c('RITUXAN', 'TYSABRI')) %>% 
             filter(start_date > first_line[1,'start_date']) %>% 
             nrow != 0) {
    # If patient was on first line drug AND eventually switched to one of 2 drugs, 
    # assign patient to first of 2 drugs that patient switched to
    drug_assign[i,] <- (p_subset %>% filter(medication_desc %in% c('RITUXAN','TYSABRI')) %>%
                          filter(start_date > first_line[1,'start_date']))[1,keep_cols]
    options[i] <- 'Option 2'
  } else {
    # If patient was on first line drug BUT did not switch to one of 2 drugs afterwards
    # assign patient to first of 2 drugs patient was on
    drug_assign[i,] <- (p_subset %>% filter(medication_desc %in% c('RITUXAN','TYSABRI')))[1,keep_cols]
    options[i] <- 'Option 3'
  }
}
drug_assign$start_date <- as.Date(drug_assign$start_date)
drug_assign$stop_date <- as.Date(drug_assign$stop_date)


# ritux_data <- drug_assign %>% filter(medication_desc == 'RITUXAN')
# nataliz_data <- drug_assign %>% filter(medication_desc == 'TYSABRI')


# Examine treatment duration and overlap ----------------------------------
# ritux_pts <- trt_subset %>% filter(medication_desc == 'RITUXAN') %>% 
#   select(PatientID) %>% unlist(use.names = F) %>% unique()
# nataliz_pts <- trt_subset %>% filter(medication_desc == 'TYSABRI') %>% 
#   select(PatientID) %>% unlist(use.names = F) %>% unique()
# 
# 
# tmp1 <- trt_subset %>% filter(PatientID %in% ritux_pts)
# tmp1$drug <- 'Rituximab'
# tmp2 <- trt_subset %>% filter(PatientID %in% nataliz_pts)
# tmp2$drug <- 'Natalizumab'
# plot_df <- rbind(tmp1, tmp2); rm(tmp1, tmp2)
# 
# ggplot(plot_df) +
#   geom_density(aes(tx_dura, fill = drug), alpha = 0.5) +
#   ggtitle('Rituximab vs. Natalizumab Treatment Duration') + theme(plot.title = element_text(hjust = 0.5)) + 
#   xlab('Treatment Duration (months)') + ylab('Density')
# 
# ggplot(plot_df) + 
#   geom_boxplot(aes(x = drug, y= tx_dura, color = drug), alpha = 0.5) + 
#   ggtitle('Rituximab vs. Natalizumab Treatment Duration') + 
#   theme(plot.title = element_text(hjust = 0.5), legend.position = "none") +
#   xlab('Drug') + ylab('Treatment Duration (months)')
#   
# 
# # Examine overlap patients
# overlap <- intersect(nataliz_pts, ritux_pts)
# ritux_seq <- rep(0, length(overlap)); nataliz_seq <- rep(0, length(overlap))
# for (i in 1:length(overlap)) {
#   p_id <- overlap[i]
#   tmp <- trt_subset %>% filter(PatientID == p_id) %>% filter(medication_desc == 'RITUXAN')
#   ritux_seq[i] <- mean(tmp$tx_dura)
#   tmp <- trt_subset %>% filter(PatientID == p_id) %>% filter(medication_desc == 'TYSABRI')
#   nataliz_seq[i] <- mean(tmp$tx_dura)
# }
# 
# ggplot() +
#   geom_point(aes(ritux_seq, nataliz_seq, size = )) + 
#   geom_smooth(aes(1:length(overlap), 1:length(overlap)), se = F) + 
#   geom_hline(aes(yintercept = 24), color = 'red') +
#   geom_text(aes(-2,26,label = 24), color = 'red') +
#   geom_text(aes(26,-2,label = 24), color = 'red') +
#   annotate(geom = "text", x = 70, y = 72, label = "y = x", color = "blue", angle = 35) + 
#   geom_vline(aes(xintercept = 24), color = 'red') + 
#   xlab('Rituximab (months)') + ylab('Natalizumab (months)') + 
#   ggtitle('Mean Rituximab vs. Natalizumab Treatment Duration for Overlapping Patients') + 
#   theme(plot.title = element_text(hjust = 0.5))
#   
# ritux_info <- list()
# for (p_id in ritux_pts) {
#   ritux_info[[p_id]] <- trt_subset %>% filter(PatientID == p_id & medication_desc == 'RITUXAN')
# }
# nataliz_info <- list()
# for (p_id in nataliz_pts) {
#   nataliz_info[[p_id]] <- trt_subset %>% filter(PatientID == p_id & medication_desc == 'TYSABRI')
# }


# Add clinical vars -------------------------------------------------------
trt_subset <- drug_assign
# Add PatientNum for later use
trt_subset$PatientNum <- sapply(trt_subset$PatientID, function(pid) {
  MS_map$PatientNum[MS_map$PatientID == pid][1]
})
trt_subset <- trt_subset %>% dplyr::select(PatientID, PatientNum, everything())


# Add sex and race
trt_subset$FEMALE <- ifelse(trt_subset$PatientID %in% MS_cohort$PatientID[MS_cohort$SEX == 'F'], 1, 0)
trt_subset$race <- ifelse(trt_subset$PatientID %in% MS_cohort$PatientID[MS_cohort$RACE_DESC == 'White'], 1, 0)
trt_subset$ethn <- ifelse(trt_subset$PatientID %in% MS_cohort$PatientID[MS_cohort$ETHNICITY_DESC == 'Not hispanic or latino'], 1, 0)
trt_subset$RACE <- with(trt_subset, race*ethn); trt_subset$race <- NULL; trt_subset$ethn <- NULL

# Age at first MS code (MS ICD: 'PheCode.335_')
ICDPheCode   = read.csv(paste0(wkpath,"EHR/MS_AllEncounter_ICDCodes_08072019.csv"), stringsAsFactors = FALSE)
ICDPheCode$phecode[ICDPheCode$concept_cd == "LPA268"] = "335_"
ICDPheCode$start_date  = as.Date(ICDPheCode$start_date,format = "%m/%d/%y")

trt_subset$FIRSTMSICD_DATE <- sapply(trt_subset$PatientNum, function(pnum) {
  tmp <- ICDPheCode %>% filter(patient_num == pnum & phecode == '335_')
  tmp <- tmp[order(tmp$start_date), ]
  return(tmp$start_date[1])
})
trt_subset$FIRSTMSICD_DATE <- as.Date(trt_subset$FIRSTMSICD_DATE)


# Exclude those with missing FIRSTMSICD_DATE (who also have no codes/encounters data)
# Doing this removes 138 patients total
trt_subset <- trt_subset %>% filter(!is.na(FIRSTMSICD_DATE))

MS_cohort$FIRSTSYMP_DATE <- with(MS_cohort, as.Date(paste(FIRSTSYMP_YEAR,FIRSTSYMP_MONTH,FIRSTSYMP_DAY,sep="-"),format="%Y-%m-%d"))
MS_cohort$BIRTHDAY <- with(MS_cohort, FIRSTSYMP_DATE - years(round(AGE_AT_FIRSTSYMPTOM)))
trt_subset$BIRTHDAY <- sapply(trt_subset$PatientID, function(pid) {
  MS_cohort$BIRTHDAY[MS_cohort$PatientID == pid]
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

# Disease duration proxy: Age at the time of DMT start - age at first MS ICD code (in years)
trt_subset$DISEASE_DURA <- with(trt_subset, as.numeric(difftime(start_date, BIRTHDAY)/365)) - 
  trt_subset$AGE_AT_FIRSTMSICD

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
trt_subset$PRIORDMT_DURA <- sapply(trt_subset$PatientID, function(pid) {
  tmp <- MS_trt_subset %>% 
    filter(PatientID == pid) %>% 
    filter(start_date < trt_subset$start_date[trt_subset$PatientID == pid])
  if (nrow(tmp) == 0) {
    return(0)
  } else {
    return(sum(tmp$tx_dura))
  }
})


# Response ----------------------------------------------------------------
# Add indicator for relapse within 2 years 
trt_subset$RELAPSE <- sapply(trt_subset$PatientID, function(pid) {
  tmp <- MS_attack %>% 
    filter(PatientID == pid) %>% 
    filter(onset >= trt_subset$start_date[trt_subset$PatientID == pid] &
             onset <= trt_subset$start_date[trt_subset$PatientID == pid] + years(2))
  if (nrow(tmp) != 0) {
    return(1)
  } else {
    return(0)
  }
})



# Predict Relapse Probability for CLIMB -----------------------------------
aggregate = FALSE
model = FALSE

if (aggregate) {
  tw = 24; tp = 3
  
  # Load imputation data
  ICD_CPT_CUI_Comb = read.csv(paste0(wkpath2,"ICD_CPT_CUI_Comb.csv"), stringsAsFactors = FALSE)
  ICD_CPT_CUI_Comb$StartDate = base::as.Date(ICD_CPT_CUI_Comb$StartDate)
  ICD_CPT_CUI_Comb$ConceptCd = as.factor(ICD_CPT_CUI_Comb$ConceptCd)
  ICDPheCode   = read.csv(paste0(wkpath,"EHR/MS_AllEncounter_ICDCodes_08072019.csv"), stringsAsFactors = FALSE)
  ICDPheCode$phecode[ICDPheCode$concept_cd == "LPA268"] = "335_"
  
  ICD_CPT_CUI_Comb_trt <- ICD_CPT_CUI_Comb %>% filter(PatientNum %in% trt_subset$PatientNum)
  
  tmp = ICD_CPT_CUI_Comb_trt[,c("PatientNum",  "PatientID", "StartDate")]
  tmp = plyr::count(tmp)[,-4]
  tmp = tmp[order(tmp$PatientNum,tmp$StartDate),]
  tmp2 = sapply(1:nrow(tmp),function(i){
    if (i %% 100 == 0) {
      print(i)
    }
    aa = ICD_CPT_CUI_Comb_trt[ICD_CPT_CUI_Comb_trt$PatientNum == tmp$PatientNum[i],]
    bb = difftime(aa$StartDate,tmp$StartDate[i],units = "days")/30
    aa = aa[bb < tw & bb >= 0,]
    aa = table(aa$ConceptCd);names(aa) = NULL
    aa
  })
  tmp2 = t(tmp2)
  colnames(tmp2) = levels(ICD_CPT_CUI_Comb_trt$ConceptCd)
  CC_comb = cbind(tmp,tmp2)
  CC_comb$StartDate <- as.Date(CC_comb$StartDate)
  # order by PatientID & StartDate
  CC_comb = CC_comb[order(CC_comb$PatientID,CC_comb$StartDate),]
  # Save intermediate data
  write.csv(CC_comb,
            paste0("intermediate_data/CLIMB_ICD_CPT_CUI_Count_",tw,"mons_",tp,"mons.csv"),
            row.names = FALSE)
  # CC_comb <- read.csv(paste0("intermediate_data/CLIMB_ICD_CPT_CUI_Count_",tw,"mons_",tp,"mons.csv"))
  
  CPTgroup = read.csv(paste0(wkpath2,"MS_ConceptCode_Mapping.csv"),stringsAsFactors = FALSE)
  tmp = grep("CPTGroup.", names(CC_comb))
  aa = names(CC_comb)[tmp]
  aa = substr(aa, 10, nchar(aa))
  aa = gsub("\\."," ",aa)
  aa = gsub("  "," ",aa)
  bb = CPTgroup$cpt_group
  bb = gsub("[-/(),]"," ",bb)
  bb = gsub("  "," ",bb)
  names(CC_comb)[tmp] = paste0("CPTGroup.", CPTgroup$cpt_group[match(aa,bb)])
  
  tmp <- readRDS("modeling_data/test24_3.rds")
  # Keep only variables in final modeling data
  CC_comb <- CC_comb[, intersect(names(tmp), names(CC_comb))]
  
  ## Add clinical variables
  # Sex
  tmp <- trt_subset[trt_subset$FEMALE == 1,'PatientID']
  CC_comb$FEMALE <- ifelse(CC_comb$PatientNum %in% tmp, 1, 0)
  # Race
  tmp <- trt_subset[trt_subset$RACE == 1,'PatientID']
  CC_comb$RACE <- ifelse(CC_comb$PatientNum %in% tmp, 1, 0)
  # Follow-up duration
  CC_comb$FOLLOWUP_DURA <- NA 
  for (pnum in unique(CC_comb$PatientNum)) {
    tmp <- trt_subset[trt_subset$PatientNum == pnum, 'FOLLOWUP_DURA']
    if (length(tmp) == 0) {
      print("No followup duration")
    } else {
      CC_comb[CC_comb$PatientNum == pnum, 'FOLLOWUP_DURA'] <- tmp
    }
  }
  # Disease duration proxy: duration between encounter date and 1st MS ICD Code (in years)
  CC_comb$DISEASE_DURA <- NA
  for (i in 1:nrow(CC_comb)) {
    first_ms_icd_date <- trt_subset[trt_subset$PatientNum == CC_comb$PatientNum[i], "FIRSTMSICD_DATE"]
    tmp <- as.numeric(difftime(CC_comb$StartDate[i], first_ms_icd_date)/365)
    if (length(tmp) == 0) {
      print("No disease duration")
    } else {
      CC_comb$DISEASE_DURA[i] <- tmp
    }
  }
  # Remove 1 problematic patient
  CC_comb$StartDate <- as.Date(CC_comb$StartDate)
  ## Data processing
  # Log transform
  # Do not transform ID columns and clinical vars 
  no_transform <- c(1:3, (length(names(CC_comb)) - 3):length(names(CC_comb)))
  CC_comb[,-no_transform] = log(1 + CC_comb[,-no_transform])
  
  # Remove data before 2006?
  # CC_comb <- CC_comb %>% filter(StartDate >= as.Date('2006-01-01'))
  
  saveRDS(CC_comb, "modeling_data/causalRN_impute_24_3_CLIMB.rds")
} else {
  CC_comb <- readRDS("modeling_data/causalRN_impute_24_3_CLIMB.rds")
}

# Model
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

# Predict probability of relapse 
trt_subset$PROB_RELAPSE <- sapply(trt_subset$PatientNum, function(pnum) {
  start_date <- trt_subset$start_date[trt_subset$PatientNum == pnum]
  pred_data <- CC_comb %>% filter(PatientNum == pnum & StartDate == start_date)
  if (nrow(pred_data) == 0) {
    pred_data <- CC_comb %>% 
      filter(PatientNum == pnum) %>% 
      filter(StartDate >= start_date - months(1) & StartDate <= start_date + months(1))
  }
  
  if (nrow(pred_data) == 0) {
    pred_data <- CC_comb %>% 
      filter(PatientNum == pnum) %>% 
      filter(StartDate >= start_date - months(3) & StartDate <= start_date + months(3))
  }
  # print(ctr)
  pred <- predict(lasso, newx = as.matrix(pred_data[1,-c(1:3)]),
                  type = "response", s = lasso$lambda.1se)
  return(pred)
})


# Remove patients without predicted relapse probability
trt_subset <- trt_subset[!is.na(trt_subset$PROB_RELAPSE),]

# Test --------------------------------------------------------------------
# tmp <- sapply(trt_subset$PatientNum, function(pnum) {
#   start_date <- trt_subset$start_date[trt_subset$PatientNum == pnum]
#   pred_data <- CC_comb %>% filter(PatientNum == pnum & StartDate == start_date)
#   ctr = 0
#   if (nrow(pred_data) == 0) {
#     pred_data <- CC_comb %>% 
#       filter(PatientNum == pnum) %>% 
#       filter(StartDate >= start_date - months(1) & StartDate <= start_date + months(1))
#     # filter(abs(difftime(StartDate, start_date, units = "weeks")/4) < 1)
#     ctr = ctr + 1
#   }
#   
#   if (nrow(pred_data) == 0) {
#     pred_data <- CC_comb %>% 
#       filter(PatientNum == pnum) %>% 
#       filter(StartDate >= start_date - months(3) & StartDate <= start_date + months(3))
#     # filter(abs(difftime(StartDate, start_date, units = "weeks")/4) < 1)
#     ctr = ctr + 1
#   }
#   if (nrow(pred_data) == 0) {
#     return(3)
#   } else {
#     return(ctr)
#   }
# })




# Save --------------------------------------------------------------------
## Move all non-modeling variables to first few columns of dataframe
trt_subset <- trt_subset %>% select(PatientID, PatientNum, medication_desc, start_date, stop_date, 
                                    tx_dura, BIRTHDAY, FIRSTICD_DATE, FIRSTMSICD_DATE, RELAPSE,
                                    PROB_RELAPSE, everything())

# First 9 columns are non-modeling, 10th column is outcome
saveRDS(trt_subset, 'modeling_data/causalRN_3mons_CLIMB.rds')

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













