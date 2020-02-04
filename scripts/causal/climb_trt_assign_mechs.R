# CLIMB: test different treatment assignment mechanisms

# Load --------------------------------------------------------------------
wkpath  = "raw_data/Box/Boston/MS CLIMB Data/"
wkpath2 = "raw_data/Box/Boston/MS CLIMB Paper2 Code&Result/Step3_BuildModel/IntermediateDataResult/"

# CLIMB DMT data
MS_trt    = read.csv(paste0(wkpath2,"Cleaned_MS_treatment.csv"), stringsAsFactors = FALSE)
MS_trt$start_date <- as.Date(with(MS_trt, paste(start_year, start_month, start_day,sep="-")), "%Y-%m-%d")
MS_trt$stop_date <- as.Date(with(MS_trt, paste(stop_year, stop_month, stop_day,sep="-")), "%Y-%m-%d")

# RXNORM Data
RXNORM_cohort <- read.xlsx(paste0(wkpath,"MS DMT/MS_RxNorm_data_withRXGenericName.xlsx"), sheet = 6)
RXNORM_cohort$start_date <- convertToDate(RXNORM_cohort$start_date)
# Interferon beta 
RXNORM_interferon_cohort <- read.csv(paste0(wkpath,"EHR/interferon_beta_20180723.csv"),
                                     header=TRUE,sep = "|", stringsAsFactors = FALSE)
RXNORM_interferon_cohort <- RXNORM_interferon_cohort[,-c(5:6)]
colnames(RXNORM_interferon_cohort) = c("patient_num","concept_cd","concept_name", "start_date")
RXNORM_interferon_cohort$start_date = as.Date(RXNORM_interferon_cohort$start_date,
                                              "%Y-%m-%d %H:%M:%S")
RXNORM_interferon_cohort$RX_Generic_Name <- "Interferon-beta"
RXNORM_interferon_cohort <- RXNORM_interferon_cohort %>% dplyr::select(patient_num, start_date, everything())
# Additional drugs
RXNORM_addl_cohort <- read.xlsx(paste0(wkpath,"MS DMT/MS_AdditionalRXNorm_data.xlsx"), sheet = 6)
RXNORM_addl_cohort$start_date <- convertToDate(RXNORM_addl_cohort$start_date)
# Combine
RXNORM_cohort <- RXNORM_cohort[, names(RXNORM_cohort)[names(RXNORM_cohort) %in% names(RXNORM_interferon_cohort)]]
RXNORM_addl_cohort <- RXNORM_addl_cohort[, names(RXNORM_addl_cohort)[names(RXNORM_addl_cohort) %in% names(RXNORM_interferon_cohort)]]
RXNORM_cohort <- rbind(RXNORM_cohort, RXNORM_addl_cohort, RXNORM_interferon_cohort)

# Using special CLIMB DMT data --------------------------------------------
trt_subset <- MS_trt

trt_subset <- trt_subset[,-c(3:8)]
trt_subset <- trt_subset %>% dplyr::select(PatientID, medication_desc, start_date, stop_date, everything())
# Add treatment duration in months
trt_subset$tx_dura <- with(trt_subset, difftime(stop_date,start_date,units = 'weeks')/4) %>% as.numeric()
# Remove patients with treatment duration = 0
trt_subset <- trt_subset %>% filter(tx_dura != 0)
MS_trt_subset <- trt_subset

ritux_pts <- trt_subset %>% filter(medication_desc == 'RITUXAN') %>%
  dplyr::select(PatientID) %>% unlist(use.names = F) %>% unique()
nataliz_pts <- trt_subset %>% filter(medication_desc == 'TYSABRI') %>%
  dplyr::select(PatientID) %>% unlist(use.names = F) %>% unique()

# First line drugs
inf_beta <- c('AVONEX','BETASERON','PLEGRIDY','REBIF')
glat_acet <- c('COPAXONE')

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
    drug_assign[i,] <- (p_subset %>% filter(medication_desc %in% c('RITUXAN','TYSABRI')))[1, keep_cols]
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
    drug_assign[i,] <- (p_subset %>% filter(medication_desc %in% c('RITUXAN','TYSABRI')))[1, keep_cols]
    options[i] <- 'Option 3'
  }
}
drug_assign$start_date <- as.Date(drug_assign$start_date)
drug_assign$stop_date <- as.Date(drug_assign$stop_date)


# First RXNORM appearance of either drug ----------------------------------
trt_patients <- unique(RXNORM_cohort$patient_num[RXNORM_cohort$RX_Generic_Name %in% c("Rituximab","Natalizumab")])
trt_subset <- data.frame(PatientNum = trt_patients)

trt_subset$medication_desc <- ""
trt_subset$start_date <- ""
# Use first codified Rx of meds as start date
for (p_num in trt_subset$PatientNum) {
  p_subset <- RXNORM_cohort %>% filter(patient_num == p_num)
  p_subset <- p_subset[order(p_subset$start_date),]
  p_ritux <- p_subset %>% filter(RX_Generic_Name == "Rituximab")
  p_nataliz <- p_subset %>% filter(RX_Generic_Name == "Natalizumab")
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


# High-efficacy criteria using RXNORM data --------------------------------

ritux_pts <- RXNORM_cohort %>% filter(RX_Generic_Name == 'Rituximab') %>%
  dplyr::select(patient_num) %>% unlist(use.names = F) %>% unique()
nataliz_pts <- RXNORM_cohort %>% filter(RX_Generic_Name == 'Natalizumab') %>%
  dplyr::select(patient_num) %>% unlist(use.names = F) %>% unique()

patients <- union(ritux_pts, nataliz_pts)
options <- rep(NA, length(patients)) # keep track of what options are applying to patients

trt_subset <- data.frame(PatientNum = patients)
trt_subset$medication_desc <- ""
trt_subset$start_date <- ""
keep_cols <- c("patient_num","RX_Generic_Name", "start_date")


RXNORM_cohort$start_date <- as.character(RXNORM_cohort$start_date)
for (i in 1:length(patients)) {
  pnum <- patients[i]
  p_subset <- RXNORM_cohort %>% filter(patient_num == pnum) %>% plyr::arrange(start_date)
  first_line <- p_subset %>% filter(RX_Generic_Name %in% c("Interferon-beta", "Glatiramer Acetate"))
  
  # If patient was never on a first line drug, assign patient to first of 2 drugs patient was on
  if (nrow(first_line) == 0) {
    trt_subset[i,] <- (p_subset %>% filter(RX_Generic_Name %in% c('Rituximab','Natalizumab')))[1, keep_cols]
    options[i] <- 'Option 1'
  } else if (p_subset %>%
             filter(RX_Generic_Name %in% c('Rituximab','Natalizumab')) %>%
             filter(start_date > first_line[1,'start_date']) %>%
             nrow != 0) {
    # If patient was on first line drug AND eventually switched to one of 2 drugs,
    # assign patient to first of 2 drugs that patient switched to
    trt_subset[i,] <- (p_subset %>% filter(RX_Generic_Name %in% c('Rituximab','Natalizumab')) %>%
                          filter(start_date > first_line[1,'start_date']))[1,keep_cols]
    options[i] <- 'Option 2'
  } else {
    # If patient was on first line drug BUT did not switch to one of 2 drugs afterwards
    # assign patient to first of 2 drugs patient was on
    trt_subset[i,] <- (p_subset %>% filter(RX_Generic_Name %in% c('Rituximab','Natalizumab')))[1, keep_cols]
    options[i] <- 'Option 3'
  }
}
RXNORM_cohort$start_date <- as.Date(RXNORM_cohort$start_date)
trt_subset$start_date <- as.Date(trt_subset$start_date)
trt_subset$stop_date <- as.Date(trt_subset$stop_date)
# Filter pre-2006 
trt_subset <- trt_subset %>% filter(start_date >= as.Date("2006-01-01"))











