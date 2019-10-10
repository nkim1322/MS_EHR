# Heatmap
# Library -----------------------------------------------------------------
library(tidyverse)
library(glmnet)
library(openxlsx)
library(lubridate)

# Load --------------------------------------------------------------------
CC_comb <- readRDS('modeling_data/train24_3.rds')
CC_comb_val <- readRDS('modeling_data/test24_3.rds')

MS_attack <- read.csv("raw_data/Box/Boston/MS CLIMB Paper2 Code&Result/Step3_BuildModel/IntermediateDataResult/Cleaned_MS_attack.csv",
                      stringsAsFactors = FALSE,
                      colClasses = c("character", rep("integer",6), "Date", rep("numeric",3),
                                     rep("character",12), rep("Date",2), "integer",
                                     rep("Date",5), rep("integer",2)))
# Load Code databases
wkpath  = "raw_data/Box/Boston/MS CLIMB Data/"
wkpath2 = "raw_data/Box/Boston/MS CLIMB Paper2 Code&Result/Step3_BuildModel/IntermediateDataResult/"
ICDPheCode   = read.csv(paste0(wkpath,"EHR/MS_AllEncounters_ICD_Data_03282019.csv"),
                        stringsAsFactors = FALSE);ICDPheCode$phecode[ICDPheCode$concept_cd == "LPA268"] = "335_"
CUIdictAll = read.xlsx(paste0(wkpath,"EHR/AllCUI_Database.xlsx"),
                       sheet = 1); colnames(CUIdictAll) = c("ConceptCd","Desc")
# Set parameters 
# tw = 24; tp = 3; rarethresh = 0.05


# Fit model  --------------------------------------------------------------
source('scripts/modeling/model.R')

# Functions ---------------------------------------------------------------
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  return(x)
}

# Keep only nonzero coefficients ------------------------------------------
coefs <- coef(glmnet.fit, s = "lambda.1se") %>% as.matrix()
coefs <- data.frame(term = row.names(coefs),
                    beta = coefs[,1])
row.names(coefs) <- NULL
nonzero_coefs <- coefs %>% filter(beta != 0)
nonzero_coefs_names <- nonzero_coefs$term
# Remove intercept term
nonzero_coefs_names <- nonzero_coefs_names[-1]
nonzero_coefs_names <- nonzero_coefs_names %>% as.character()

# Convert to matrix form
coefs <- coefs[coefs$term %in% nonzero_coefs_names,'beta'] %>% as.matrix()


# Standardize to SD = 1 ---------------------------------------------------
scaled_X <- scale(CC_comb_val[, names(CC_comb_val) %in% nonzero_coefs_names], center = F) %>% as.data.frame
scaled_X$PatientID <- CC_comb_val$PatientID
scaled_X$StartDate <- CC_comb_val$StartDate

# Select patients ---------------------------------------------------------
## This was done manually
lowprob_highrelapse <- list('BWH-492479', 'BWH-482798', 'BWH-489265') %>% set_names()
highprob_lowrelapse <- list('BWH-487063', 'BWH-486782', 'BWH-487706') %>% set_names()

pts_maxprob_norelapse <- c("BWH-489888", "BWH-493354", "BWH-492416", "BWH-490352", "BWH-486174",
                           "BWH-486150", "BWH-482462", "BWH-492315", "BWH-486782", "BWH-493418") %>% 
  set_names()

pts_minprob_withrelapse <- c("BWH-492013", "BWH-491072", "BWH-487099", "BWH-490848", "BWH-491854",
                             "BWH-480939", "BWH-492479", "BWH-480822", "BWH-486920", "BWH-489994") %>% 
  set_names()


# Sample encounters -------------------------------------------------------
# High prob, low relapse
X_highprob_lowrelapse <- map(pts_maxprob_norelapse, function(p_id) {
  p_data <- scaled_X %>% filter(PatientID == p_id)
  sample_list <- list()
  for (i in 1:nrow(p_data)) {
    enctr_date <- p_data[i, 'StartDate']
    # For two-year follow-up after given encounter date
    p_subset <- p_data %>% filter(StartDate > enctr_date & StartDate <= enctr_date + years(2))
    if (nrow(p_subset) > 0) {
      # Divide follow-up period into 3 month periods, only loop over as many as 3 month periods as there are
      for (j in 1:ceiling(((p_subset[nrow(p_subset), 'StartDate'] - enctr_date) %>% as.double)/(30*tp))) {
        p_subsubset <- p_subset %>%
          filter(StartDate > enctr_date + (j-1)*months(3) & StartDate <= enctr_date + j*months(3))
        # If more than 1 encounter within 3 month period, sample 1 visit and save data
        if (nrow(p_subsubset) > 0) {
          set.seed(12345)
          sample_list[[paste0(i,j)]] <- sample_n(p_subsubset, size = 1)
        }}}}
  # Bind to create dataframe
  p_df <- do.call(rbind, sample_list)
  # Remove duplicate sampled encounters (?)
  p_df <- distinct(p_df)
})

# Low prob, high relapse
X_lowprob_highrelapse <- map(pts_minprob_withrelapse, function(p_id) {
  p_data <- scaled_X %>% filter(PatientID == p_id)
  sample_list <- list()
  for (i in 1:nrow(p_data)) {
    enctr_date <- p_data[i, 'StartDate']
    # For two-year follow-up after given encounter date
    p_subset <- p_data %>% filter(StartDate > enctr_date & StartDate <= enctr_date + years(2))
    if (nrow(p_subset) > 0) {
      # Divide follow-up period into 3 month periods, only loop over as many as 3 month periods as there are
      for (j in 1:ceiling(((p_subset[nrow(p_subset), 'StartDate'] - enctr_date) %>% as.double)/(30*tp))) {
        p_subsubset <- p_subset %>%
          filter(StartDate > enctr_date + (j-1)*months(3) & StartDate <= enctr_date + j*months(3))
        # If more than 1 encounter within 3 month period, sample 1 visit and save data
        if (nrow(p_subsubset) > 0) {
          set.seed(12345)
          sample_list[[paste0(i,j)]] <- sample_n(p_subsubset, size = 1)
        }}}}
  # Bind to create dataframe
  p_df <- do.call(rbind, sample_list)
  # Remove duplicate sampled encounters
  p_df <- distinct(p_df)
})

# Multiply coefs & X's -----------------------------------------------------
# High prob, low relapse
data_highprob_lowrelapse <- map(X_highprob_lowrelapse, function(tmp) {
  p_id <- tmp$PatientID %>% unique()
  tmp_dates <- tmp$StartDate
  tmp$StartDate = tmp$PatientID = NULL
  coefs <- t(coefs)
  
  X_beta <- tmp
  for (i in 1:length(coefs)) {X_beta[,i] <- coefs[i]*tmp[,i]}
  X_beta <- X_beta %>% as.matrix()
  # Reinsert encounter date
  rownames(X_beta) <- tmp_dates
  
  X_beta_df <- X_beta %>% as.data.frame()
  X_beta_df$StartDate <- tmp_dates
  X_beta_df <- X_beta_df %>% select(StartDate, everything())
  test <- X_beta_df %>% 
    gather(names(X_beta_df)[-1], key = "code", value = "value")
  # Convert encounter date to string for graphing purposes
  test$StartDate <- as.character(test$StartDate)
  # Replace phecode with phenotype string for interpretability
  # CPT Codes
  tmp2 <- strtrim(test$code, 3) == 'CPT'
  test$Code[tmp2] <- str_split(test$code[tmp2], '\\.', simplify = T)[,2] %>% 
    as.character() %>% strtrim(width = 20) %>% firstup()
  # ICD PheCodes
  tmp2 <- strtrim(test$code, 3) == 'Phe'
  tmp3 <- str_split(test$code[tmp2], '\\.', simplify = T)[,2]
  tmp4 <- character(length(tmp3))
  for (i in 1:length(tmp3)) {
    tmp4[i] <- ICDPheCode[ICDPheCode$phecode == tmp3[i], 'phecode_description'] %>% unique()
  }
  test$Code[tmp2] <- tmp4 %>% strtrim(width = 20)
  # CUIs
  tmp2 <- strtrim(test$code, 3) == 'CUI'
  tmp3 <- str_split(test$code[tmp2], '\\.', simplify = T)[,2]
  tmp4 <- character(length(tmp3))
  for (i in 1:length(tmp3)) {
    tmp4[i] <- unique(CUIdictAll[CUIdictAll$ConceptCd == tmp3[i], 'Desc'])[1]
  }
  test$Code[tmp2] <- tmp4 %>% strtrim(width = 20)
  # Clinical variables
  tmp2 <- !(strtrim(test$code, 3) == 'CPT' | strtrim(test$code, 3) == 'Phe')
  test$Code[tmp2] <- test$code[tmp2]
  # Add back in PatientID
  test$PatientID <- rep(p_id, nrow(test))
  return(test)
})

# Low prob, high relapse
data_lowprob_highrelapse <- map(X_lowprob_highrelapse, function(tmp) {
  p_id <- tmp$PatientID %>% unique()
  tmp_dates <- tmp$StartDate
  tmp$StartDate = tmp$PatientID = NULL
  
  X_beta <- tmp
  for (i in 1:length(coefs)) { X_beta[,i] <- coefs[i]*tmp[,i]}
  X_beta <- X_beta %>% as.matrix()
  # Reinsert encounter date
  rownames(X_beta) <- tmp_dates
  
  X_beta_df <- X_beta %>% as.data.frame()
  X_beta_df$StartDate <- tmp_dates
  X_beta_df <- X_beta_df %>% select(StartDate, everything())
  test <- X_beta_df %>% 
    gather(names(X_beta_df)[-1], key = "code", value = "value")
  # Convert encounter date to string for graphing purposes
  test$StartDate <- as.character(test$StartDate)
  # Replace phecode with phenotype string for interpretability
  # CPT Codes
  tmp2 <- strtrim(test$code, 3) == 'CPT'
  test$Code[tmp2] <- str_split(test$code[tmp2], '\\.', simplify = T)[,2] %>% 
    as.character() %>% strtrim(width = 20) %>% firstup()
  # ICD PheCodes
  tmp2 <- strtrim(test$code, 3) == 'Phe'
  tmp3 <- str_split(test$code[tmp2], '\\.', simplify = T)[,2]
  tmp4 <- character(length(tmp3))
  for (i in 1:length(tmp3)) {
    tmp4[i] <- ICDPheCode[ICDPheCode$phecode == tmp3[i], 'phecode_description'] %>% unique()
  }
  test$Code[tmp2] <- tmp4 %>% strtrim(width = 20)
  # CUIs
  tmp2 <- strtrim(test$code, 3) == 'CUI'
  tmp3 <- str_split(test$code[tmp2], '\\.', simplify = T)[,2]
  tmp4 <- character(length(tmp3))
  for (i in 1:length(tmp3)) {
    tmp4[i] <- unique(CUIdictAll[CUIdictAll$ConceptCd == tmp3[i], 'Desc'])[1]
  }
  test$Code[tmp2] <- tmp4 %>% strtrim(width = 20)
  # Clinical variables
  tmp2 <- !(strtrim(test$code, 3) == 'CPT' | strtrim(test$code, 3) == 'Phe')
  test$Code[tmp2] <- test$code[tmp2]
  # Add back in PatientID
  test$PatientID <- rep(p_id, nrow(test))
  return(test)
})

# Fix CUIs
data_highprob_lowrelapse <- map(data_highprob_lowrelapse, function(df){
  tmp <- strtrim(df$Code, 3) == 'CUI'
  tmp2 <- str_split(df$Code[tmp], '\\.', simplify = T)[,2]
  tmp3 <- character(length(tmp2))
  for (i in 1:length(tmp2)) {
    tmp3[i] <- unique(CUIdictAll[CUIdictAll$ConceptCd == tmp2[i], 'Desc'])[1]
  }
  df$Code[tmp] <- tmp3
  return(df)
})

data_lowprob_highrelapse <- map(data_lowprob_highrelapse, function(df){
  tmp <- strtrim(df$Code, 3) == 'CUI'
  tmp2 <- str_split(df$Code[tmp], '\\.', simplify = T)[,2]
  tmp3 <- character(length(tmp2))
  for (i in 1:length(tmp2)) {
    tmp3[i] <- unique(CUIdictAll[CUIdictAll$ConceptCd == tmp2[i], 'Desc'])[1]
  }
  df$Code[tmp] <- tmp3
  return(df)
})

# Plot --------------------------------------------------------------------
graphs_highprob_lowrelapse <- map(data_highprob_lowrelapse, function(tmp) {
  g <- ggplot(tmp, aes(x = StartDate, y = Code, fill = value)) +
    geom_tile() + 
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.y=element_text(size=4),
          axis.text.x = element_text(angle = 90, size = 4)) + 
    ggtitle(paste0("Patient ID: ", tmp$PatientID[1])) + 
    scale_fill_gradient2()
  return(g)
})

graphs_lowprob_highrelapse <- map(data_lowprob_highrelapse, function(tmp) {
  g <- ggplot(tmp, aes(x = StartDate, y = Code, fill = value)) +
    geom_tile() + 
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.y = element_text(size = 4),
          axis.text.x = element_text(angle = 90, size = 4)) +
    ggtitle(paste0("Patient ID: ", tmp$PatientID[1])) + 
    scale_fill_gradient2()
  return(g)
})


# Output ------------------------------------------------------------------
pdf('plots/heatmap_maxprob_norelapse.pdf')
graphs_highprob_lowrelapse
dev.off()

pdf('plots/heatmap_minprob_withrelapse.pdf')
graphs_lowprob_highrelapse
dev.off()




