# Clinical hook: examining T25-FW vs actual & predicted 2 year relapse rate
# Library -----------------------------------------------------------------
library(tidyverse)
library(openxlsx)
library(gee)
library(glmnet)

# Load --------------------------------------------------------------------
# Relapse data
wkpath2 = "raw_data/Box/Boston/MS CLIMB Paper2 Code&Result/Step3_BuildModel/IntermediateDataResult/"
MS_attack = read.csv(paste0(wkpath2,"Cleaned_MS_attack.csv"),stringsAsFactors = FALSE,
                     colClasses = c("character",rep("integer",6),"Date",rep("numeric",3),
                                    rep("character",12),rep("Date",2),"integer",
                                    rep("Date",5),rep("integer",2)))

# T25-FW data
walk_data <- read.xlsx("raw_data/Box/Boston/MS CLIMB Data/CLIMB Cohort Add'l Data/25ftWalk121615.xlsx",
                       detectDates = F)
walk_data$visit_date <- convertToDate(walk_data$visit_date, origin = "1900-01-01")

# Clean walk data ---------------------------------------------------------
# Exclude those with time less than 0 or greater than 180 seconds
walk_data <- walk_data %>% filter(`25ft_walking_time` <= 180 & `25ft_walking_time` >= 0)
# Remove provider info
walk_data$provider <- NULL
# Rename columns
colnames(walk_data) <- c('PatientID', 'WalkDate', 'WalkTime')
# Log transform WalkTime due to right skewedness
walk_data$WalkTime <- log(walk_data$WalkTime)


# Subset relapse data -----------------------------------------------------
MS_attack_subs <- MS_attack %>% select(PatientID, onset, val_start, val_stop)


# Correlation: subset with validated relapse only -------------------------
# Select only validated patients and the encounters within validated time periods
val_pts <- MS_attack_subs %>% filter(!is.na(val_start)) %>% select(PatientID) %>% unlist(use.names = F) %>% unique()
MS_attack_val <- MS_attack_subs %>% filter(!is.na(val_start))
walk_data_val <- walk_data %>% filter(PatientID %in% val_pts)

walk_data_val$CC_1yr <- NA
walk_data_val$CC_2yr <- NA
for (i in 1:nrow(walk_data_val)) {
  p_id <- walk_data_val$PatientID[i]
  curr_date <- walk_data_val$WalkDate[i]
  attack_subset <- MS_attack_val %>% filter(PatientID == p_id)
  walk_data_val$CC_1yr[i] <- ifelse(nrow(attack_subset %>% filter(onset >= curr_date & onset <= curr_date + years(1))) > 0, 1, 0)
  walk_data_val$CC_2yr[i] <- ifelse(nrow(attack_subset %>% filter(onset >= curr_date & onset <= curr_date + years(2))) > 0, 1, 0)
}

logreg_val_1yr <- glm(CC_1yr ~ WalkTime, walk_data_val, family = binomial)
logreg_val_2yr <- glm(CC_2yr ~ WalkTime, walk_data_val, family = binomial)

# Correlation: all data with relapse info ---------------------------------
walk_data <- walk_data %>% filter(PatientID %in% unique(MS_attack$PatientID))

walk_data$CC_1yr <- NA
walk_data$CC_2yr <- NA
for (i in 1:nrow(walk_data)) {
  p_id <- walk_data$PatientID[i]
  curr_date <- walk_data$WalkDate[i]
  attack_subset <- MS_attack_subs %>% filter(PatientID == p_id)
  walk_data$CC_1yr[i] <- ifelse(nrow(attack_subset %>% filter(onset >= curr_date & onset <= curr_date + years(1))) > 0, 1, 0)
  walk_data$CC_2yr[i] <- ifelse(nrow(attack_subset %>% filter(onset >= curr_date & onset <= curr_date + years(2))) > 0, 1, 0)
}

logreg_all_1yr <- glm(CC_1yr ~ WalkTime, walk_data, family = binomial)
logreg_all_2yr <- glm(CC_2yr ~ WalkTime, walk_data, family = binomial)

# Correlation: predicted relapse ------------------------------------------
# GEE with quasibinomial link

# Tuning parameters
tw_list     = rep(c(6, 12, 24), each = 3) # time window in months
tp_list     = rep(c(1, 3, 6), 3) # tw/4, time period for sample cases and controls in months

mdls_list <- list()
for (i in 1:length(tw_list)) {
  tw = tw_list[i]; tp = tp_list[i]
  print(paste0('tw: ',tw,', tp:',tp))
  # Read in appropriate data
  CC_comb <- readRDS(paste0('modeling_data/train', tw, '_', tp, '.rds'))
  CC_comb_val <- readRDS(paste0('modeling_data/test', tw, '_', tp, '.rds'))
  CC_all <- rbind(CC_comb, CC_comb_val)
  # Fit 
  set.seed(12345)
  glmnet.fit = cv.glmnet(as.matrix(CC_comb[,-c(1:9)]), unlist(CC_comb$CC),
                         family = "binomial", type.measure = "deviance")
  CC_subset <- CC_all[,c(1:4)]
  CC_subset$RelapsePred <- predict(glmnet.fit, newx = as.matrix(CC_all[,-c(1:9)]),
                                type = "response", s = glmnet.fit$lambda.1se) %>% as.numeric()

  walk_data$CC_pred <- NA
  # Add WalkTime 
  for (j in 1:nrow(walk_data)) {
    p_id <- walk_data$PatientID[j]
    curr_date <- walk_data$WalkDate[j]
    preds_subset <- CC_subset %>% filter(PatientID == p_id) %>% 
      filter(StartDate >= curr_date) %>% filter(StartDate <= curr_date + months(tw))
    if (nrow(preds_subset) > 0) {
      walk_data$CC_pred[j] <- mean(preds_subset$RelapsePred)
    }
  }
  # Fit GEE 
  tmp_data <- walk_data[,c('PatientID','WalkTime', 'CC_pred')] %>% filter(!is.na(CC_pred))
  print(paste0("N = ", nrow(tmp_data)))
  tmp_data$PatientID <- factor(tmp_data$PatientID)
  gee.fit <- gee(formula = CC_pred ~ WalkTime, 
                 id = PatientID, 
                 data = tmp_data, 
                 na.action = na.omit, maxiter = 100,
                 family = quasi)
  print(tmp <- summary(gee.fit)$coefficients['WalkTime', c('Estimate', 'Robust z')])
  mdls_list[[paste0('tw',tw,'tp',tp)]] <- tmp
}


# Add relapse data --------------------------------------------------------
# comb_data <- walk_data
# # Number of relapses within 2 year follow-up of a given walk date
# comb_data$RelapseRate <- NA
# 
# for (i in 1:nrow(comb_data)) {
#   p_id <- comb_data[i, 'PatientID']
#   walkdate <- comb_data[i, 'WalkDate']
#   attacks <- MS_attack_subs %>% filter(PatientID == p_id) %>% filter(onset >= walkdate & onset <= walkdate + 365*2)
#   comb_data[i, 'RelapseRate'] <- nrow(attacks)
# }






