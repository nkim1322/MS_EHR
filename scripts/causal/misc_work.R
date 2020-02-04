# Script for miscellaneous work

# Investigate divergent bootstrap values ----------------------------------

# Estimate
# Semi-supervised all data
set.seed(1234)
n <- nrow(trt_subset)
bootsims <- 10000
boot_indices <- list()
for (i in 1:bootsims) {
  boot_indices[[i]] <- sample(1:n, n, replace = T)
}

boot_ests <- rep(NA, bootsims)
for (i in 1:bootsims) {
  if (i %% 50 == 0) {
    print(i)
  }
  
  boot_data <- trt_subset[boot_indices[[i]], ]
  
  boot_ests[i] <- dips.ss(Yi = boot_data$RELAPSE,
                          Ti = ifelse(boot_data$medication_desc == 'Natalizumab', 1, 0),
                          Ri = boot_data$COHORT,
                          Wi = boot_data$PROB_RELAPSE, 
                          Xi = as.matrix(boot_data[,-c(1:10)]), 
                          Gi = rep(1, n),
                          fam = "binomial")
}

# Investigate outliers
outlier_idx <- which(abs(boot_ests) > 100)

graphs_list <- list()
for (j in 1:length(outlier_idx)) {
  tmp_data <- trt_subset[boot_indices[[outlier_idx[j]]], ]
  tmp_data <- data.frame(tmp_data)
  trt_subset <- data.frame(trt_subset)
  for (i in c(8:20)) {
    x <- data.frame(x = tmp_data[,i])
    y <- data.frame(x = trt_subset[,i])
    graphs_list[[paste0(as.character(j), ":", as.character(i))]] <- ggplot() + 
      geom_histogram(data = x, mapping = aes(x), fill = "red", alpha = 0.5) +
      geom_histogram(data = y, mapping = aes(x), fill = "blue", alpha = 0.5) +
      ggtitle(paste0(names(trt_subset)[i], ", Est: ", boot_ests[outlier_idx[j]]))
  }
}
pdf("plots/causal/hist_stability.pdf")
graphs_list
dev.off()



# Semi-supervised CLIMB only
set.seed(1234)
n <- nrow(trt_subsubset_val)
bootsims <- 10000
boot_indices_climb <- list()
for (i in 1:bootsims) {
  boot_indices_climb[[i]] <- sample(1:n, n, replace = T)
}

boot_ests_climb <- rep(NA, bootsims)
for (i in 1:bootsims) {
  if (i %% 50 == 0) {
    print(i)
  }
  
  boot_data <- trt_subsubset_val[boot_indices_climb[[i]], ]
  
  boot_ests_climb[i] <- dips.ss(Yi = boot_data$RELAPSE,
                          Ti = ifelse(boot_data$medication_desc == 'Natalizumab', 1, 0),
                          Ri = boot_data$COHORT,
                          Wi = boot_data$PROB_RELAPSE, 
                          Xi = as.matrix(boot_data[,-c(1:10)]), 
                          Gi = rep(1, n),
                          fam = "binomial")
}

outlier_idx <- which(abs(boot_ests_climb) > 100)

graphs_list <- list()
for (j in 1:length(outlier_idx)) {
  tmp_data <- trt_subsubset_val[boot_indices_climb[[outlier_idx[j]]], ]
  for (i in c(8:20)) {
    x <- data.frame(x = tmp_data[,i])
    y <- data.frame(x = trt_subsubset_val[,i])
    graphs_list[[paste0(as.character(j), ":", as.character(i))]] <- ggplot() + 
      geom_histogram(data = x, mapping = aes(x), fill = "red", alpha = 0.5) +
      geom_histogram(data = y, mapping = aes(x), fill = "blue", alpha = 0.5) +
      ggtitle(paste0(names(trt_subset)[i], ", Est: ", boot_ests_climb[outlier_idx[j]]))
  }
}
pdf("plots/causal/hist_stability_climb.pdf")
graphs_list
dev.off()


# Comparing CLIMB and non-CLIMB covariates --------------------------------
graphs_list <- list()
climb <- trt_subset %>% filter(COHORT == 1)
non_climb <- trt_subset %>% filter(COHORT == 0)
graph_cols <- names(trt_subset)
for (i in 8:20) {
  x <- data.frame(x = climb[,i])
  y <- data.frame(x = non_climb[,i])
  graphs_list[[as.character(i)]] <- ggplot() + 
    geom_histogram(data = x, mapping = aes(x), fill = "red", alpha = 0.5) +
    geom_histogram(data = y, mapping = aes(x), fill = "blue", alpha = 0.5) +
    # geom_density(data = x, mapping = aes(x, ..scaled..), fill = "red", alpha = 0.5) + 
    # geom_density(data = y, mapping = aes(x, ..scaled..), fill = "blue", alpha = 0.5) +
    ggtitle(graph_cols[i])
}

pdf("plots/causal/hist_subset_comparison.pdf")
graphs_list
dev.off()


# Compare kernel weights --------------------------------------------------

graphs_list <- list()
tmp_df <- data.frame(weight = pi.ric,
                     cohort = trt_subset$COHORT)
tmp_df$cohort <- ifelse(tmp_df$cohort == 1, "Labeled", "Unlabeled")
tmp_df$cohort <- factor(tmp_df$cohort)

colors = c("1" = "Labeled", "2" = "Unlabeled")

graphs_list[[1]] <- ggplot() + 
  geom_histogram(data = tmp_df, aes(x = weight, fill = cohort), alpha = 0.5) +
  labs(x = "weights", fill = "Legend") + 
  scale_color_manual(values = colors) + 
  ggtitle("Kernel Smoothed Weights of the joint OR and PS (All Data)") 

tmp_df <- data.frame(weight = pi.ric,
                     cohort = trt_subsubset_val$COHORT)
tmp_df$cohort <- ifelse(tmp_df$cohort == 1, "Labeled", "Unlabeled")
tmp_df$cohort <- factor(tmp_df$cohort)

graphs_list[[2]] <- ggplot() + 
  # geom_histogram(data = tmp_df, aes(x = weight, fill = cohort), alpha = 0.5) +
  geom_density(data = tmp_df, mapping = aes(weight, ..scaled.., fill = cohort), alpha = 0.5)+
  labs(x = "weights", fill = "Legend") + 
  scale_color_manual(values = colors) + 
  ggtitle("Kernel Smoothed Weights of the joint OR and PS (CLIMB only)") 

pdf("plots/causal/hist_kernel_smoothed_weights.pdf")
graphs_list
dev.off()

# Lasso-HMM data ----------------------------------------------------------

idx1 <- sample(c(1:nrow(CLIMB_comb)), floor(0.1*nrow(CLIMB_comb)))
t1 <- CLIMB_comb[idx1,]


# x <- runAll(t1, t2)
source("scripts/modeling/lasso_hmm.R")
source("scripts/modeling/lasso_hmm_v2.R")

t3 <- cbind(as.numeric(rbernoulli(nrow(t2), p = 0.2)), 1, 1, 1, 1, t2)
colnames(t3)[1:5] <- c("CC", "Dummy2", "Dummy3", "Dummy4", "Dummy5")
x <- runAll(CLIMB_comb, t3)
y <- runAll(CLIMB_comb, t3)

idx2 <- sample(c(1:nrow(CC_comb)), floor(0.2*nrow(CC_comb)))
t2 <- CC_comb[idx2,]
x <- runAllSS(CLIMB_comb, t2)


output <- runAll(train, test)
