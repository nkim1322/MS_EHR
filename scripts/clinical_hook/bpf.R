# Brain parenchymal fraction (BPF)
# Library -----------------------------------------------------------------
library(tidyverse)
library(openxlsx)
library(gee)

# Load --------------------------------------------------------------------
# Relapse data
wkpath2 = "raw_data/Box/Boston/MS CLIMB Paper2 Code&Result/Step3_BuildModel/IntermediateDataResult/"
MS_attack = read.csv(paste0(wkpath2,"Cleaned_MS_attack.csv"),stringsAsFactors = FALSE,
                     colClasses = c("character",rep("integer",6),"Date",rep("numeric",3),
                                    rep("character",12),rep("Date",2),"integer",
                                    rep("Date",5),rep("integer",2)))

# BPF Data
bpf <- read.xlsx("raw_data/Box/Boston/MS CLIMB Data/CLIMB Cohort Add'l Data/MS Database Query 2014-03-14 i2b2All_MS_moreMShx.xlsx",
                  sheet = 5)
bpf$MRI_DATE <- as.Date(with(bpf, paste(MRI_YEAR, MRI_MONTH, MRI_DAY,sep="-")), "%Y-%m-%d")
## Select only relevant columns
bpf <- bpf %>% select(subject, MRI_DATE, BPF) %>% filter(!is.na(BPF))
# Remove negative values of BPF
bpf <- bpf %>% filter(BPF >= 0)

# Subset relapse data -----------------------------------------------------
MS_attack_subs <- MS_attack %>% select(PatientID, onset, val_start, val_stop)

# Correlation: subset with validated relapse only -------------------------
# Select only validated patients and the encounters within validated time periods
val_pts <- MS_attack_subs %>% filter(!is.na(val_start)) %>% select(PatientID) %>% unlist(use.names = F) %>% unique()
MS_attack_val <- MS_attack_subs %>% filter(!is.na(val_start))
bpf_val <- bpf %>% filter(subject %in% val_pts)

bpf_val$CC_1yr <- NA
bpf_val$CC_2yr <- NA
for (i in 1:nrow(bpf_val)) {
  p_id <- bpf_val$subject[i]
  curr_date <- bpf_val$MRI_DATE[i]
  attack_subset <- MS_attack_val %>% filter(PatientID == p_id)
  bpf_val$CC_1yr[i] <- ifelse(nrow(attack_subset %>% filter(onset >= curr_date & onset <= curr_date + years(1))) > 0, 1, 0)
  bpf_val$CC_2yr[i] <- ifelse(nrow(attack_subset %>% filter(onset >= curr_date & onset <= curr_date + years(2))) > 0, 1, 0)
}

logreg_val_1yr <- glm(CC_1yr ~ BPF, bpf_val, family = binomial)
logreg_val_2yr <- glm(CC_2yr ~ BPF, bpf_val, family = binomial)
summary(logreg_val_1yr); summary(logreg_val_2yr)



# Correlation: all data with relapse info ---------------------------------
bpf <- bpf %>% filter(subject %in% unique(MS_attack$PatientID))

bpf$CC_1yr <- NA
bpf$CC_2yr <- NA
for (i in 1:nrow(bpf)) {
  p_id <- bpf$subject[i]
  curr_date <- bpf$MRI_DATE[i]
  attack_subset <- MS_attack_subs %>% filter(PatientID == p_id)
  bpf$CC_1yr[i] <- ifelse(nrow(attack_subset %>% filter(onset >= curr_date & onset <= curr_date + years(1))) > 0, 1, 0)
  bpf$CC_2yr[i] <- ifelse(nrow(attack_subset %>% filter(onset >= curr_date & onset <= curr_date + years(2))) > 0, 1, 0)
}

logreg_all_1yr <- glm(CC_1yr ~ BPF, bpf, family = binomial)
logreg_all_2yr <- glm(CC_2yr ~ BPF, bpf, family = binomial)
summary(logreg_all_1yr); summary(logreg_all_2yr)


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
  
  bpf$CC_pred <- NA
  # Add WalkTime 
  for (j in 1:nrow(bpf)) {
    p_id <- bpf$subject[j]
    curr_date <- bpf$MRI_DATE[j]
    preds_subset <- CC_subset %>% filter(PatientID == p_id) %>% 
      filter(StartDate >= curr_date) %>% filter(StartDate <= curr_date + months(tw))
    if (nrow(preds_subset) > 0) {
      bpf$CC_pred[j] <- mean(preds_subset$RelapsePred)
    }
  }
  # Fit GEE 
  tmp_data <- bpf[,c('subject','BPF', 'CC_pred')] %>% filter(!is.na(CC_pred))
  print(paste0("N = ", nrow(tmp_data)))
  tmp_data$subject <- factor(tmp_data$subject)
  gee.fit <- gee(formula = CC_pred ~ BPF, 
                 id = subject, 
                 data = tmp_data, 
                 na.action = na.omit, maxiter = 100,
                 family = quasi)
  print(tmp <- summary(gee.fit)$coefficients['BPF', c('Estimate', 'Robust z')])
  mdls_list[[paste0('tw',tw,'tp',tp)]] <- tmp
}

