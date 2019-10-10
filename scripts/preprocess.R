# Data cleaning and preprocessing
## This script takes raw data for cleaning and aggregation according tuning params
# Library -----------------------------------------------------------------
library(tidyverse)
library(openxlsx)
library(lubridate)

# Set parameters and paths ------------------------------------------------
readin = TRUE 
wkpath  = "raw_data/Box/Boston/MS CLIMB Data/"
wkpath2 = "raw_data/Box/Boston/MS CLIMB Paper2 Code&Result/Step3_BuildModel/IntermediateDataResult/"
year   = 2000:2016
rarethresh = 0.05 # threshold for rare codes; slightly better than 0.02

# Cohort, Treatment, Attack Data ------------------------------------------
MS_map  = read.xlsx(paste0(wkpath,"CLIMB Cohort/Spec95_i2b2_Mapping2017.xlsx"), sheet = 1)
colnames(MS_map) = c("PatientNum","PatientID")
# Read in and clean Cohort/Treatment/Attack data
if (readin) {
  MS_cohort = read.csv(paste0(wkpath2,"Cleaned_MS_cohort.csv"),stringsAsFactors = FALSE,
                       colClasses = c(rep("character",4), rep("Date",4), rep("integer",6),rep("numeric",3),
                                      "integer",rep("Date",4),rep("numeric",3),"integer"))
  MS_attack = read.csv(paste0(wkpath2,"Cleaned_MS_attack.csv"),stringsAsFactors = FALSE,
                       colClasses = c("character",rep("integer",6),"Date",rep("numeric",3),
                                      rep("character",12),rep("Date",2),"integer",
                                      rep("Date",5),rep("integer",2)))
  MS_trt    = read.csv(paste0(wkpath2,"Cleaned_MS_treatment.csv"),stringsAsFactors = FALSE)
  
} else {
  ## cohort and relapse
  MS_cohort = read.xlsx(paste0(wkpath,"CLIMB Cohort/ClimbCohortReviewSent081717.xlsx"), sheet = 1)
  MS_trt    = read.xlsx(paste0(wkpath,"CLIMB Cohort/ClimbDMTAttackZongqi062817Sent.xlsx"), sheet = 1)
  MS_attack = read.xlsx(paste0(wkpath,"CLIMB Cohort/ClimbDMTAttackZongqi062817Sent.xlsx"), sheet = 2)
  colnames(MS_trt)[1] = colnames(MS_attack)[1] = colnames(MS_cohort)[1] = "PatientID"
  ### remove duplicates
  MS_cohort = MS_cohort[!duplicated(MS_cohort$PatientID),]
  # There is one patient in cohort data with missing enrollment date. Further checking shows all the information of this patients is missing. 
  MS_cohort = MS_cohort[!is.na(MS_cohort$BASELINE_DATE),]
  ## remove duplicates based on onset_year,onset_month,onset_day
  MS_attack = MS_attack[!duplicated(MS_attack[,c("PatientID","onset_year","onset_month","onset_day")]),]
  ## further remove 28 duplicates based on onset_year & onset_month
  # MS_attack = MS_attack[order(MS_attack$PatientID,MS_attack$onset_year,MS_attack$onset_month,MS_attack$onset_day),]
  MS_attack = MS_attack[!duplicated(MS_attack[,c("PatientID","onset_year","onset_month")]),]
  
  MS_EHR    = read.xlsx(paste0(wkpath,"EHR/MS_first_last_date_whole_list.xlsx"),
                        sheet = 1); colnames(MS_EHR)[1] = "PatientNum"
  MS_EHR$PatientID = MS_map$PatientID[match(MS_EHR$PatientNum,MS_map$PatientNum)]  
  
  ### Belongs to EHR or not
  MS_cohort$EHR = as.numeric(MS_cohort$PatientID%in%MS_EHR$PatientID)
  MS_attack$EHR = as.numeric(MS_attack$PatientID%in%MS_EHR$PatientID)
  
  ### add last code date and last date = max(last code date, last visit date) to both cohort and attack data
  MS_cohort$LAST_CODE   = MS_EHR$Last_date[match(MS_cohort$PatientID,MS_EHR$PatientID)]
  MS_cohort$FIRST_CODE  = MS_EHR$First_date[match(MS_cohort$PatientID,MS_EHR$PatientID)]
  MS_attack$LAST_CODE   = MS_EHR$Last_date[match(MS_attack$PatientID,MS_EHR$PatientID)]
  MS_attack$FIRST_CODE  = MS_EHR$First_date[match(MS_attack$PatientID,MS_EHR$PatientID)]
  
  ### Impute first diagnosis date (cohort)
  ## impute missing month/day
  tmp                  = is.na(MS_cohort$FIRSTDIAG_MONTH)
  MS_cohort$FIRSTDIAG_MONTH[tmp] = 7; MS_cohort$FIRSTDIAG_DAY[tmp]   = 1
  tmp                  = is.na(MS_cohort$FIRSTDIAG_DAY)
  MS_cohort$FIRSTDIAG_DAY[tmp]   = 15
  a                    = lapply(1:nrow(MS_cohort),function(i){
    as.Date(paste0(c(as.character(MS_cohort$FIRSTDIAG_YEAR[i]),as.character(MS_cohort$FIRSTDIAG_MONTH[i]),as.character(MS_cohort$FIRSTDIAG_DAY[i])),collapse ="-"),format="%Y-%m-%d")
  })
  a                    = do.call(c,a)
  ## get first/last date
  MS_cohort$LAST_DATE  = pmax(as.Date(MS_cohort$LAST_VISIT_DATE,
                                      format = "%Y-%m-%d", origin = "1900-01-01"),
                              as.Date(MS_cohort$LAST_CODE,
                                      format = "%Y-%m-%d", origin = "1900-01-01"),
                              na.rm = TRUE) # max(last vist, last code)
  MS_cohort$FIRST_DATE = pmin(a, # first diagnosis date
                              as.Date(MS_cohort$FIRST_CODE, 
                                      format = "%Y-%m-%d", origin="1900-01-01"),
                              as.Date(MS_cohort$BASELINE_DATE, 
                                      format = "%Y-%m-%d",origin="1900-01-01"),
                              na.rm = TRUE)
  ## delete patients with first date>last date
  tmp                  = with(MS_cohort,FIRST_DATE>LAST_DATE)
  tmp[is.na(tmp)]      = FALSE
  MS_cohort            = MS_cohort[!tmp,]
  a                    = a[!tmp]
  MS_attack            = MS_attack[MS_attack$PatientID %in% MS_cohort$PatientID &
                                     MS_attack$onset_year%in%year,]
  MS_trt               = MS_trt[MS_trt$PatientID%in%MS_cohort$PatientID,]
  # Match patients in attack and cohort dataframes
  MS_attack$LAST_DATE  = MS_cohort$LAST_DATE[match(MS_attack$PatientID,MS_cohort$PatientID)]
  MS_attack$FIRST_DATE = MS_cohort$FIRST_DATE[match(MS_attack$PatientID,MS_cohort$PatientID)]
  
  
  ### Impute first symptom date (cohort)
  ## impute missing month/day
  tmp2                 = is.na(MS_cohort$FIRSTSYMP_MONTH)
  MS_cohort$FIRSTSYMP_MONTH[tmp2] = 7; MS_cohort$FIRSTSYMP_DAY[tmp2]   = 1
  tmp2                 = is.na(MS_cohort$FIRSTSYMP_DAY)
  MS_cohort$FIRSTSYMP_DAY[tmp2]   = 15
  
  ## Impute onset date (attack)
  tmp                  = is.na(MS_attack$onset_day)
  MS_attack$onset_day[tmp] = 15
  MS_attack$onset = with(MS_attack,as.Date(paste0(onset_year,"-",onset_month,"-",onset_day),
                                           format = "%Y-%m-%d"))
  
  ### excludes patients with problematic entries (trt)
  ## exclude entries with missing stop year and history == stopped 
  ## as well as entries with missing start year but history == ever
  tmp  = with(MS_trt,(is.na(stop_year) & history == "stopped"))
  tmp2 = with(MS_trt,(is.na(start_year) & history == "ever"))
  tmp3 = with(MS_trt,start_year > stop_year)
  tmp4 = with(MS_trt,start_year == stop_year & start_month > stop_month)
  tmp5 = with(MS_trt,start_year == stop_year & start_month ==stop_month & start_day > stop_day)
  tmp3[is.na(tmp3)] = FALSE; tmp4[is.na(tmp4)] = FALSE; tmp5[is.na(tmp5)] = FALSE
  MS_trt     = MS_trt[!(tmp|tmp2|tmp3|tmp4|tmp5),]
  ## for patient with a missing stop year but history == current, impute by 2017/7/1
  tmp        = is.na(MS_trt$stop_year)
  MS_trt$stop_year[tmp]   = 2017; MS_trt$stop_month[tmp]  = 7; MS_trt$stop_day[tmp]    = 1
  ## for patient with a missing stop month, impute by 7/1
  tmp        = is.na(MS_trt$stop_month)
  MS_trt$stop_month[tmp]  = 7; MS_trt$stop_day[tmp]    = 1
  ## for patient with a missing stop day, impute by 15
  tmp        = is.na(MS_trt$stop_day)
  MS_trt$stop_day[tmp]    = 15
  ## for patient with a missing start month, impute by 7/1
  tmp        = is.na(MS_trt$start_month)
  MS_trt$start_month[tmp] = 7; MS_trt$start_day[tmp]   = 1
  ## for patient with a missing start day, impute by 15
  tmp        = is.na(MS_trt$start_day)
  MS_trt$start_day[tmp]   = 15
  ## check consistency
  tmp3   = with(MS_trt,start_year>stop_year)
  tmp4   = with(MS_trt,start_year==stop_year & start_month>stop_month)
  tmp5   = with(MS_trt,start_year==stop_year & start_month==stop_month & start_day > stop_day)
  MS_trt = MS_trt[!(tmp3|tmp4|tmp5),]
  
  # Convert relevant columns to Date format
  MS_cohort$LAST_VISIT_DATE = as.Date(MS_cohort$LAST_VISIT_DATE,format = "%Y-%m-%d", origin = "1900-01-01")
  MS_cohort$FIRST_CODE = as.Date(MS_cohort$FIRST_CODE, format = "%Y-%m-%d", origin = "1900-01-01")
  MS_cohort$FIRST_DATE = as.Date(MS_cohort$FIRST_DATE, format = "%Y-%m-%d", origin = "1900-01-01")
  
  ## Age at first code
  MS_cohort$AGE_AT_FIRSTCODE = with(MS_cohort,AGE_AT_LASTVISIT - as.numeric(difftime(LAST_VISIT_DATE, FIRST_CODE,units="days"))/365)
  
  ## Age at sympton onset
  MS_cohort$AGE_AT_FIRSTSYMPTOM = with(MS_cohort,AGE_AT_LASTVISIT-as.numeric(difftime(LAST_VISIT_DATE,                                               as.Date(paste(FIRSTSYMP_YEAR,FIRSTSYMP_MONTH,FIRSTSYMP_DAY,sep="-"),format="%Y-%m-%d"),units="days"))/365)
  
  ## Age at first date
  MS_cohort$AGE_AT_FIRSTDATE = with(MS_cohort,AGE_AT_LASTVISIT-as.numeric(difftime(LAST_VISIT_DATE,FIRST_DATE,units="days"))/365)
  
  ## Duration of follow-up (last date - first date)
  MS_cohort$DURA = with(MS_cohort,as.numeric(difftime(LAST_DATE, FIRST_DATE,
                                                      units = "days"))/365)
  
  ## Number of treatments
  tmp    = unique(MS_trt$PatientID)
  tmp2   = with(MS_trt,sapply(tmp,function(x){ length(unique(medication_desc[PatientID == x]))}))
  MS_cohort$NumTrt = tmp2[match(MS_cohort$PatientID,tmp)]
  tmp    = is.na(MS_cohort$NumTrt)
  MS_cohort$NumTrt[tmp] = 0
  
  ## clinical relapse and radiographic relapse
  clinical   = names(MS_attack)[12:23]
  radiograph = clinical[(1:12%%3) == 2]; clinical   = clinical[(1:12%%3) != 2]
  MS_attack$clinical    = sapply(1:nrow(MS_attack),function(x){
    as.numeric(sum(MS_attack[x,clinical] == "Y", na.rm = TRUE) >= 1)
  })
  MS_attack$radiographic = sapply(1:nrow(MS_attack),function(x){
    as.numeric(sum(MS_attack[x,radiograph] == "Y",na.rm = TRUE) >= 1)
  }) 
  
  
  write.csv(MS_cohort,paste0(wkpath2,"Cleaned_MS_cohort.csv"),row.names = FALSE)
  write.csv(MS_attack,paste0(wkpath2,"Cleaned_MS_attack.csv"),row.names = FALSE)
  write.csv(MS_trt,paste0(wkpath2,"Cleaned_MS_treatment.csv"),row.names = FALSE)
}




# ICD, CPT, CUI Codes -----------------------------------------------------
CUIdictAll = read.xlsx(paste0(wkpath,"EHR/AllCUI_Database.xlsx"), sheet = 1)
colnames(CUIdictAll) = c("ConceptCd","Desc")
if (readin) {
  ICD_CPT_CUI_Comb = read.csv(paste0(wkpath2,"ICD_CPT_CUI_Comb.csv"), stringsAsFactors = FALSE)
  ICD_CPT_CUI_Comb$StartDate = as.Date(ICD_CPT_CUI_Comb$StartDate)
  ICD_CPT_CUI_Comb$ConceptCd = as.factor(ICD_CPT_CUI_Comb$ConceptCd)
  ICDPheCode   = read.csv(paste0(wkpath,"EHR/MS_AllEncounter_ICDCodes_08072019.csv"), stringsAsFactors = FALSE)
  ICDPheCode$phecode[ICDPheCode$concept_cd == "LPA268"] = "335_"
} else {
  ICDPheCode   = read.csv(paste0(wkpath,"EHR/MS_AllEncounter_ICDCodes_08072019.csv"),stringsAsFactors = FALSE)
  tmp = ICDPheCode$concept_cd == "LPA268"; ICDPheCode$phecode[tmp] = "335_"
  ICDPheCode$phecode_description[tmp] = "Multiple sclerosis"
  tmp = ICDPheCode[ICDPheCode$concept_cd%in%c("V700","V762"),-6]
  ICDPheCode  = ICDPheCode[!ICDPheCode$concept_cd%in%c("V700","V762"),]
  CPTGrouped  = read.csv(paste0(wkpath,"EHR/MS_AllEncounters_CPTCodes_08072019.csv"),stringsAsFactors = FALSE)
  colnames(tmp) = colnames(CPTGrouped)
  CPTGrouped  = rbind(tmp, CPTGrouped)
  
  # ***Upddated CUI sheet no longer has header***
  CUISelected  = read.csv(paste0(wkpath,"EHR/MS_AllEncounters_CUIs_08072019.csv"), stringsAsFactors = FALSE, 
                          col.names = c('patient_num', 'encounter_num', 'start_date', 'concept_cd'))
  colnames(ICDPheCode)  = c("PatientNum","EncounterNum","StartDate","ConceptCd","PheCode","Desc")
  colnames(CPTGrouped)  = c("PatientNum","EncounterNum","StartDate","ConceptCd","CPTGroup")
  colnames(CUISelected) = c("PatientNum","EncounterNum","StartDate","ConceptCd")
  ICDPheCode$StartDate  = as.Date(ICDPheCode$StartDate,format = "%m/%d/%y")
  CPTGrouped$StartDate  = as.Date(CPTGrouped$StartDate,format = "%m/%d/%y")
  CUISelected$StartDate = format(as.POSIXct(CUISelected$StartDate, format = '%Y-%m-%d %H:%M:%S'), format = '%Y-%m-%d')
  # Add Code descriptions
  CUISelected$Desc = CUIdictAll$Desc[match(CUISelected$ConceptCd,CUIdictAll$ConceptCd)]
  # Format
  CPTGrouped$ConceptCd = toupper(CPTGrouped$ConceptCd)
  CPTGrouped$CPTGroup  = substr(CPTGrouped$CPTGroup, 1, nchar(CPTGrouped$CPTGroup))
  # Alter specific codes
  CPTGrouped$CPTGroup[CPTGrouped$ConceptCd%in%c("V700","CG0463")] = "outpatient clinic visit"
  CPTGrouped$CPTGroup[CPTGrouped$ConceptCd%in%c("V762","CG0101","CG0202","CQ0091")] ="cancer screening"
  CPTGrouped$CPTGroup[CPTGrouped$ConceptCd%in%c("C90765","C90766","C90780","C90781","CC8950","CC8951")] = "hydration, therapeutic, prophylactic, diagnostic injections and infusions, and chemotherapy and other highly complex drug or highly complex biologic agent administration"
  CPTGrouped$CPTGroup[CPTGrouped$ConceptCd%in%c("C70540","C70542","C70543")] = "MRI orbit"
  CPTGrouped$CPTGroup[CPTGrouped$ConceptCd%in%c("C70551","C70552","C70553")] = "MRI brain"
  CPTGrouped$CPTGroup[CPTGrouped$ConceptCd%in%c("C72141","C72142","C72146",
                                                "C72147","C72156","C72157")] = "MRI spine (cervical, thoracic)"
  
  ICDPheCode$ConceptCd  = paste0("PheCode.", ICDPheCode$PheCode)
  CPTGrouped$ConceptCd  = paste0("CPTGroup.", CPTGrouped$CPTGroup)
  CUISelected$ConceptCd = paste0("CUI.", CUISelected$ConceptCd)
  
  # Combine dataframes
  ICD_CPT_CUI_Comb = rbind(ICDPheCode[,c("PatientNum","EncounterNum","StartDate","ConceptCd")],
                           CPTGrouped[,c("PatientNum","EncounterNum","StartDate","ConceptCd")],
                           CUISelected[,c("PatientNum","EncounterNum","StartDate","ConceptCd")])
  ICD_CPT_CUI_Comb = count(ICD_CPT_CUI_Comb)[,c("PatientNum", "StartDate", "ConceptCd")]
  ICD_CPT_CUI_Comb$PatientID = MS_map$PatientID[match(ICD_CPT_CUI_Comb$PatientNum, MS_map$PatientNum)] 
  ICD_CPT_CUI_Comb$ConceptCd = as.factor(ICD_CPT_CUI_Comb$ConceptCd)
  
  write.csv(ICD_CPT_CUI_Comb, paste0(wkpath2,"ICD_CPT_CUI_Comb.csv"),row.names = FALSE)
}




# Cohort data vs Unlabeled EMR: divide codes data -------------------------
ICD_CPT_CUI_Comb_cohort = ICD_CPT_CUI_Comb[ICD_CPT_CUI_Comb$PatientID%in%MS_cohort$PatientID,]
ICD_CPT_CUI_Comb_UL     = ICD_CPT_CUI_Comb[!ICD_CPT_CUI_Comb$PatientID%in%MS_cohort$PatientID,]


# Aggregate codes ---------------------------------------------------------
tw_list     = rep(c(6, 12, 24), each = 3) # time window in months
tp_list     = rep(c(1, 3, 6), 3) # tw/4, time period for sample cases and controls in months


if (!readin) {
  for (i in 1:length(tw_list)) {
    tw = tw_list[i]; tp = tp_list[i]
    print(paste('tw:', tw, 'tp:',tp))
    
    MS_Enct_Uniq = read.csv(paste0("raw_data/HSPH/Data From Partners MS Center/MS_Encounter_Date_PatientID.csv"), stringsAsFactors = FALSE)
    MS_Enct_Uniq = MS_Enct_Uniq[MS_Enct_Uniq$PatientID%in%MS_cohort$PatientID,]
    MS_Enct_Uniq$StartDate = as.Date(MS_Enct_Uniq$StartDate,format = "%m/%d/%y")
    tmp = plyr::count(MS_Enct_Uniq[,c("PatientNum","PatientID","StartDate")])[,-4]
    tmp2 = sapply(1:nrow(tmp),function(i){
      aa = MS_attack[MS_attack$PatientID == tmp$PatientID[i],]
      CC = abs(difftime(tmp$StartDate[i], aa$onset, units = "days")/30) < tw # 
      CR = c(any(aa$clinical[CC] >= 1),any(aa$radiographic[CC] >= 1))  # clinical & radiographic relapse
      CC = any(CC) # case-control status
      aa = ICD_CPT_CUI_Comb_cohort[ICD_CPT_CUI_Comb_cohort$PatientID == tmp$PatientID[i],]
      bb = difftime(aa$StartDate,tmp$StartDate[i],units = "days")/30
      aa = aa[bb < tw & bb >= 0,]
      aa = c(CC,CR,table(aa$ConceptCd));names(aa) = NULL
      aa
    })
    tmp2 = t(tmp2)
    colnames(tmp2) = c("CC","Clinical","Radiographic",levels(ICD_CPT_CUI_Comb$ConceptCd)) 
    CC_comb = cbind(tmp,tmp2)
    CC_comb$LAST_DATE = MS_cohort$LAST_DATE[match(CC_comb$PatientID,MS_cohort$PatientID)]
    CC_comb$FIRST_DATE = MS_cohort$FIRST_DATE[match(CC_comb$PatientID,MS_cohort$PatientID)]
    # order by PatientID & StartDate
    tmpp    = order(CC_comb$PatientID,CC_comb$StartDate)
    CC_comb = CC_comb[tmpp,]
    # add time period to min(First date, First encounter) by tp months
    tmpp = tapply(CC_comb$StartDate,CC_comb$PatientID,min)
    tmpp = tmpp[match(CC_comb$PatientID,names(tmpp))]
    tmpp = pmax(CC_comb$FIRST_DATE, as.Date(tmpp,origin = "1970-01-01"))
    tmpp = data.frame(CC_comb[,c(1:4)],FIRST_DATE2 = tmpp,stringsAsFactors = FALSE)
    tmpp$StartDate   = as.Date(tmpp$StartDate)
    tmpp$FIRST_DATE2 = as.Date(tmpp$FIRST_DATE2)
    CC_comb$Period = as.numeric(difftime(tmpp$StartDate,tmpp$FIRST_DATE2,units="days")/30/tp)
    CC_comb = CC_comb[CC_comb$Period>=0,]
    CC_comb$Period = ceiling(CC_comb$Period)
    # organize columns
    CC_comb = CC_comb[,c(1:6,229:231,7:228)]
    
    CPTgroup = read.csv(paste0(wkpath2,"MS_ConceptCode_Mapping.csv"), stringsAsFactors = FALSE)
    
    tmp  = grep("CPTGroup.",names(CC_comb))
    aa   = names(CC_comb)[tmp]
    aa   = substr(aa,10,nchar(aa))
    indx = 1:ncol(CC_comb)
    indx[tmp] = sort(match(aa,CPTgroup$cpt_group),index.return=TRUE)$ix + 9
    CC_comb = CC_comb[,indx]
    
    write.csv(CC_comb,
              paste0("intermediate_data/CC_ICD_CPT_CUI_Count_",tw,"mons_",tp,"mons.csv"),
              row.names = FALSE)
  }
}


# Format and add clinical variables ---------------------------------------

for (i in 1:length(tw_list)) {
  tw = tw_list[i]; tp = tp_list[i]
  print(paste('tw:', tw, 'tp:',tp))
  CC_comb <- read.csv(paste0("intermediate_data/CC_ICD_CPT_CUI_Count_",tw,"mons_",tp,"mons.csv"))
  # Order CPT groups by number, deal with discrepency due to special characters 
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
  
  # Remove rare codes
  CodeFreq = data.frame(Cases = apply(CC_comb[CC_comb$CC == 1, -c(1:9)], 2, mean),
                        Controls = apply(CC_comb[CC_comb$CC == 0, -c(1:9)], 2, mean),
                        Overall = apply(CC_comb[, -c(1:9)], 2, mean))
  CodeFreq = round(CodeFreq,4)
  tmp = with(CodeFreq, Cases > rarethresh | Controls > rarethresh)
  CC_comb = CC_comb[,c(rep(TRUE,9),tmp)]
  
  ## Remove two-digit phecode if one-digit exists
  tmp = names(CC_comb)[grep("PheCode.*_[0-9][0-9]",names(CC_comb))]
  tmp = tmp[substr(tmp,1,nchar(tmp)-1)%in%names(CC_comb)]
  CC_comb = CC_comb[,!names(CC_comb)%in%tmp]
  ## remove zero-digit phecode if one-digit exists
  tmp = names(CC_comb)[grep("PheCode.*_[0-9]",names(CC_comb))]
  tmp = tmp[substr(tmp,1,nchar(tmp)-1)%in%names(CC_comb)]
  tmp = substr(tmp,1,nchar(tmp)-1)
  CC_comb = CC_comb[,!names(CC_comb)%in%tmp]
  
  # Add clinical variables
  # Subset relevant data, 'DURA' is follow-up duration (in years)
  cohort_data <- MS_cohort[, c('PatientID', 'SEX', 'RACE_DESC', 'ETHNICITY_DESC', 'DURA',
                               'AGE_AT_LASTVISIT', 'LAST_VISIT_DATE', 'AGE_AT_FIRSTSYMPTOM')]
  CC_comb$StartDate <- as.Date(CC_comb$StartDate)
  # Sex
  tmp <- cohort_data[cohort_data$SEX == 'F','PatientID']
  CC_comb$FEMALE <- ifelse(CC_comb$PatientID %in% tmp, 1, 0)
  # Race
  CC_comb$RACE <- NA
  races <- cohort_data$RACE_DESC %>% unique()
  for (race in races) {
    tmp <- cohort_data[cohort_data$RACE_DESC == race,'PatientID']
    CC_comb$RACE[CC_comb$PatientID %in% tmp] <- race
  }
  # Ethnicity
  CC_comb$ETHNICITY <- NA
  ethnicities <- cohort_data$ETHNICITY_DESC %>% unique()
  for (ethn in ethnicities) {
    tmp <- cohort_data[cohort_data$ETHNICITY_DESC == ethn,'PatientID']
    CC_comb$ETHNICITY[CC_comb$PatientID %in% tmp] <- ethn
  }
  # RACE == 1 if "non-Hispanic European descent"
  tmp <- cohort_data[cohort_data$RACE == 'White' & cohort_data$ETHNICITY == 'Not hispanic or latino','PatientID']
  CC_comb$RACE <- ifelse(CC_comb$PatientID %in% tmp, 1, 0)
  CC_comb$ETHNICITY <- NULL
  # **Disease duration = age at encounter - age at 1st symptom
  cohort_data$FIRSTSYMP_DATE <- with(MS_cohort, as.Date(paste(FIRSTSYMP_YEAR,FIRSTSYMP_MONTH,FIRSTSYMP_DAY,sep="-"),format="%Y-%m-%d"))
  cohort_data$BIRTHDAY <- with(cohort_data, FIRSTSYMP_DATE - years(round(AGE_AT_FIRSTSYMPTOM)))
  
  for (pt in unique(CC_comb$PatientID)) {
    CC_comb[CC_comb$PatientID == pt, 'BIRTHDAY'] <- cohort_data[cohort_data$PatientID == pt, 'BIRTHDAY']
    CC_comb[CC_comb$PatientID == pt, 'AGE_AT_FIRSTSYMPTOM'] <- cohort_data[cohort_data$PatientID == pt, 'AGE_AT_FIRSTSYMPTOM']
  }
  
  CC_comb$AGE_AT_VISIT <- as.numeric(with(CC_comb, difftime(StartDate, BIRTHDAY, unit = 'weeks')/52.25))
  CC_comb$DISEASE_DURA <- with(CC_comb, AGE_AT_VISIT - AGE_AT_FIRSTSYMPTOM)
  
  # CC_comb$AGE_AT_FIRSTSYMPTOM <- NA; CC_comb$BIRTHDAY <- NA
  # cohort_data$BIRTHDAY <- with(cohort_data, LAST_VISIT_DATE - years(round(AGE_AT_LASTVISIT)))
  # for (pt in unique(CC_comb$PatientID)) {
  #   CC_comb[CC_comb$PatientID == pt, 'AGE_AT_FIRSTSYMPTOM'] <- cohort_data[cohort_data$PatientID == pt, 'AGE_AT_FIRSTSYMPTOM']
  #   CC_comb[CC_comb$PatientID == pt, 'BIRTHDAY'] <- cohort_data[cohort_data$PatientID == pt, 'BIRTHDAY']
  #   CC_comb$BIRTHDAY <- as.Date(CC_comb$BIRTHDAY, origin = '1970-01-01')
  # }
  # 
  # CC_comb$DISEASE_DURA <- with(CC_comb, difftime(StartDate, BIRTHDAY, units = 'days')/365 - AGE_AT_FIRSTSYMPTOM) %>% as.numeric()
  CC_comb$BIRTHDAY <- NULL; CC_comb$AGE_AT_FIRSTSYMPTOM <- NULL; CC_comb$AGE_AT_VISIT <- NULL
  CC_comb$FOLLOWUP_DURA <- NA 
  for (pt in unique(CC_comb$PatientID)) {
    # DURA variable is equivalent to the follow-up duration taht we want
    CC_comb[CC_comb$PatientID == pt, 'FOLLOWUP_DURA'] <- cohort_data[cohort_data$PatientID == pt, 'DURA']
  }
  # Remove visits with negative disease duration: this does not affect number of patients
  # CC_comb <- CC_comb %>% filter(DISEASE_DURA > 0)
  
  # Remove NAs
  CC_comb <- CC_comb[!is.na(CC_comb$RACE) & !is.na(CC_comb$DISEASE_DURA) & !is.na(CC_comb$FOLLOWUP_DURA),]
  
  # Remove ICD Parent
  # Ultimately opt to remove codes 1005_, 1010_, 773_, 563_ 
  drops <- c('PheCode.1005_', 'PheCode.1010_', 'PheCode.773_', 'PheCode.563_')
  CC_comb <- CC_comb[ , !(names(CC_comb) %in% drops)]
  
  ## Save final data
  saveRDS(CC_comb, paste0('intermediate_data/CC_comb_', tw, '_', tp, '.rds'))
}


# Transform codes and split -----------------------------------------------
# Counts to log(1+counts)
for (i in 1:length(tw_list)) {
  tw = tw_list[i]; tp = tp_list[i]
  print(paste('tw:', tw, 'tp:',tp))
  
  CC_comb <- readRDS(paste0('intermediate_data/CC_comb_', tw, '_', tp, '.rds'))
  # Do not transform logistic columns and clinical vars 
  no_transform <- c(1:9, (length(names(CC_comb)) - 3):length(names(CC_comb)))
  CC_comb[,-no_transform] = log(1 + CC_comb[,-no_transform])
  
  # Sample (Yuri Edits: )
  group   = paste(CC_comb$PatientID,CC_comb$Period,sep=".")
  indices = tapply(1:nrow(CC_comb), group, function(x){ifelse(length(x) == 1, x, sample(x,1))})
  counts = tapply(1:nrow(CC_comb), group, function(x){length(x)})
  CC_comb = CC_comb[indices,]
  # CC_comb$periodCounts <- counts
  
  # Remove data before 2006
  CC_comb <- CC_comb %>% filter(StartDate >= as.Date('2006-01-01'))
  
  ## Split into training and validation
  set.seed(1000)
  nval = floor(length(unique(CC_comb$PatientID))*0.3)
  val  = sample(unique(CC_comb$PatientID),nval)
  CC_comb_val = CC_comb[CC_comb$PatientID%in%val,]
  CC_comb     = CC_comb[!CC_comb$PatientID%in%val,]
  ## Subsample case-control data; one obs per time period 
  # group   = paste(CC_comb$PatientID, CC_comb$Period, sep = ".")
  # indices = tapply(1:nrow(CC_comb), group, sample, size = 1)
  # CC_comb = CC_comb[indices, ]
  # print(sub_rate = round(mean(CC_comb$CC),3))
  
  
  
  # Save
  saveRDS(CC_comb, paste0('modeling_data/train', tw, '_', tp, '.rds'))
  saveRDS(CC_comb_val, paste0('modeling_data/test', tw, '_', tp, '.rds'))
}
