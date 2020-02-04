# Run causal inference on all patients
# Library -----------------------------------------------------------------
library(tidyverse)
source('scripts/causal/causalmethods.R')

# Load --------------------------------------------------------------------
EHR_trt <- readRDS('modeling_data/causalRN_3mons_EHR.rds')
CLIMB_trt <- readRDS("modeling_data/causalRN_3mons_CLIMB.rds")
# EHR_impute <- readRDS("modeling_data/causalRN_impute_24_3_EHR.rds")
# train <- readRDS("modeling_data/train24_3.rds"); test <- readRDS("modeling_data/test24_3.rds")
# CLIMB_impute <- rbind(train, test); rm(train,test)

# View distribution of predicted probabilities ----------------------------
png("plots/causal/hist_imputed_probs.png", units = "in", width = 7, height = 5, res = 300)
hist(EHR_trt$PROB_RELAPSE, 
     main = "Histogram of EHR Imputed 2-Year Relapse Probability",
     xlab = "Imputed Probability of Relapse")
dev.off()

# Combine data ------------------------------------------------------------
W_i <- c(CLIMB_trt$PROB_RELAPSE, EHR_trt$PROB_RELAPSE)

CLIMB_trt$PROB_RELAPSE <- NULL

names(EHR_trt)[names(EHR_trt) == "PROB_RELAPSE"] <- "RELAPSE"
EHR_trt$RELAPSE <- 1

CLIMB_trt$stop_date <- NULL; CLIMB_trt$tx_dura <- NULL


trt_subset <- rbind(CLIMB_trt, EHR_trt)


  
# Get estimate ------------------------------------------------------------
tx_vec <- ifelse(trt_subset$medication_desc == 'TYSABRI', 1, 0)
R_vec <- c(rep(1, nrow(CLIMB_trt)), rep(0, nrow(EHR_trt)))

n <- nrow(trt_subset)
pihat <- dips.ss(Yi = trt_subset$RELAPSE,
                 Ti = tx_vec,
                 Ri = R_vec,
                 Wi = W_i, 
                 Xi = as.matrix(trt_subset[,-c(1:8)]), 
                 Gi = rep(1, n),
                 fam = "binomial")

est <- ipw(Yi = trt_subset$RELAPSE, Ti = tx_vec, pihat, normalize = T)
est

# Bootstrap ---------------------------------------------------------------
bootsims <- 10000
boot_ests <- rep(NA, bootsims)

set.seed(1234)
for (i in 1:bootsims) {
  if (i %% 50 == 0) {
    print(i)
  }
  boot_idx <- sample(1:nrow(trt_subset), replace = T)
  boot_data <- trt_subset[boot_idx, ]
  tx_vec <- ifelse(boot_data$medication_desc == 'TYSABRI', 1, 0)
  
  W_boot <- W_i[boot_idx]
  R_boot <- R_vec[boot_idx]
  
  pihat <- dips.ss(Yi = trt_subset$RELAPSE,
                   Ti = tx_vec,
                   Ri = R_boot,
                   Wi = as.matrix(W_boot), 
                   Xi = as.matrix(trt_subset[boot_idx,-c(1:8)]), 
                   Gi = rep(1, n),
                   fam = "binomial") # Is this correct?

  boot_ests[i] <- ipw(Yi = boot_data$RELAPSE, Ti = tx_vec, pihat, normalize = T)
}

write.csv(boot_ests, "causal_output/boot_ests_RN_SS.csv")


# Check for significance --------------------------------------------------
ci_bounds <- quantile(boot_ests, probs = c(0.025, 0.975))
print(ci_bounds)
ci_bounds[1] > 0 | ci_bounds[2] < 0

# Get point estimate  -----------------------------------------------------

mean(boot_ests)

# mean(abs(boot_ests) >= abs(est)
# hist(boot_ests)
# abline(v = 0, col = 2)










