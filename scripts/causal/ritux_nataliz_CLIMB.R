# Run causal inference on CLIMB
# Library -----------------------------------------------------------------
# library(tidyverse)
source('scripts/causal/causalmethods.R')


# Load --------------------------------------------------------------------
trt_subset <- readRDS('modeling_data/causalRN_3mons_CLIMB.rds')

# Get estimate ------------------------------------------------------------
tx_vec <- ifelse(trt_subset$medication_desc == 'TYSABRI', 1, 0)

n <- nrow(trt_subset)
pihat <- dips(Yi = trt_subset$RELAPSE,
              Ti = tx_vec,
              Xi = as.matrix(trt_subset[,-c(1:10)]),
              Gi = rep(1, n),
              fam = "binomial")

est <- ipw(Yi = trt_subset$RELAPSE, Ti = tx_vec, pihat, normalize=T)
# aipw(Yi,Ti,pihat.glm,muhat1.glm,muhat0.glm)
est

# Bootstrap ---------------------------------------------------------------
bootsims <- 10000
boot_ests <- rep(NA, bootsims)

for (i in 1:bootsims) {
  boot_idx <- sample(1:nrow(trt_subset), replace = T)
  boot_data <- trt_subset[boot_idx, ]
  tx_vec <- ifelse(boot_data$medication_desc == 'TYSABRI', 1, 0)
  pihat <- dips(Yi = boot_data$RELAPSE,
                Ti = tx_vec,
                Xi = as.matrix(boot_data[,-c(1:10)]),
                Gi = rep(1,nrow(boot_data)),
                fam = "binomial")
  
  # ipw(Yi = boot_data$RELAPSE, Ti = tx_vec, pihat, normalize = T)
  boot_ests[i] <- ipw(Yi = boot_data$RELAPSE, Ti = tx_vec, pihat, normalize = T)
}

write.csv(boot_ests, "causal_output/boot_ests_RN_CLIMB.csv")


# Get point estimate ------------------------------------------------------
mean(boot_ests)


# Check for significance --------------------------------------------------
ci_bounds <- quantile(boot_ests, probs = c(0.025, 0.975))
print(ci_bounds)
# ci_bounds[1] > 0 | ci_bounds[2] < 0


# mean(abs(boot_ests) >= abs(est)
# hist(boot_ests)
# abline(v = 0, col = 2)










