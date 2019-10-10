# Sort patients by degree of misclassification

# Generate dataframe of predicted probabilities called df
# Library -----------------------------------------------------------------
library(tidyverse)
library(glmnet)
# Load --------------------------------------------------------------------
CC_comb <- readRDS('intermediate_data/CC_comb_24_3_clin_clean.rds') # Training data
CC_comb_val <- readRDS('intermediate_data/CC_comb_val_24_3_clin_clean.rds') # Validation data
MS_attack <- read.csv("raw_data/Box/Boston/MS CLIMB Paper2 Code&Result/Step3_BuildModel/IntermediateDataResult/Cleaned_MS_attack.csv",
                      stringsAsFactors = FALSE,
                      colClasses = c("character", rep("integer",6), "Date", rep("numeric",3),
                                     rep("character",12), rep("Date",2), "integer",
                                     rep("Date",5), rep("integer",2)))
# load('models/glmnet_cv_deviance_24mons_3mons_0.05rarethresh.rda')
# load('models/glmnet_cv_deviance_24mons_3mons_0.05rarethresh_clin.rda')
load('models/glmnet_cv_deviance_24mons_3mons_0.05rarethresh_clin_2006on_clean.rda')

# Predict -----------------------------------------------------------------
# Filter data for final model
CC_comb <- CC_comb %>% filter(StartDate >= as.Date('2006-01-01'))
CC_comb_val <- CC_comb_val %>% filter(StartDate >= as.Date('2006-01-01'))
CC_comb$RACE2 <- NULL; CC_comb$ETHNICITY2 <- NULL
CC_comb_val$RACE2 <- NULL; CC_comb_val$ETHNICITY2 <- NULL 

df <- CC_comb_val[,1:6]
df$Prob <- predict(glmnet.fit, 
                   newx = as.matrix(CC_comb_val[,-c(1:9)]),
                   type = "response", s = glmnet.fit$lambda.1se)


# Sort --------------------------------------------------------------------

# Sort patients into relapse or no relapse
no_relapse <- unique(CC_comb_val$PatientID)[!(unique(CC_comb_val$PatientID) %in% unique(MS_attack$PatientID))]
with_relapse <- unique(CC_comb_val$PatientID)[unique(CC_comb_val$PatientID) %in% unique(MS_attack$PatientID)]

# Create lists to use map function
no_relapse <- no_relapse %>% as.list %>% set_names()
with_relapse <- with_relapse %>% as.list %>% set_names()


# For patients with no relapse --------------------------------------------
# Obtain max predicted prob
maxprob_no_relapse <- map(no_relapse, function(p_id){
  p_subset <- df %>% filter(PatientID == p_id)
  return(max(p_subset$Prob))
})
# Create dataframe
df_maxprob_no_relapse <- data.frame(PatientID = names(maxprob_no_relapse),
                                    Max_prob = unlist(maxprob_no_relapse, use.names = F))
# Sort by predicted prob, decreasing
df_maxprob_no_relapse <- df_maxprob_no_relapse[order(df_maxprob_no_relapse$Max_prob, decreasing = T),]

# Get top ten
pts_maxprob_norelapse <- df_maxprob_no_relapse %>% head(n= 10) %>% select(PatientID) %>% unlist(use.names = F)

# pts_maxprob_norelapse <- c("BWH-489888", "BWH-493354", "BWH-492416", "BWH-490352", "BWH-486174",
#                            "BWH-486150", "BWH-482462", "BWH-492315", "BWH-486782", "BWH-493418")

# For patients with relapse -----------------------------------------------
# Obtain min predicted prob
minprob_with_relapse <- map(with_relapse, function(p_id){
  p_subset <- df %>% filter(PatientID == p_id)
  return(min(p_subset$Prob))
})

# Create dataframe
df_minprob_with_relapse <- data.frame(PatientID = names(minprob_with_relapse),
                                      Min_prob = unlist(minprob_with_relapse, use.names = F))
# Sort by predicted prob, increasing
df_minprob_with_relapse <- df_minprob_with_relapse[order(df_minprob_with_relapse$Min_prob),]

# Get top ten
pts_minprob_withrelapse <- df_minprob_with_relapse %>% head(n = 10) %>% select(PatientID) %>% unlist(use.names = F)
# pts_minprob_withrelapse <- c("BWH-492013", "BWH-491072", "BWH-487099", "BWH-490848", "BWH-491854",
#                            "BWH-480939", "BWH-492479", "BWH-480822", "BWH-486920", "BWH-489994")



# Plot probability plots --------------------------------------------------
# **Run plot_per_pt function in 'prob_plots.R' script first!!

# Max prob, no relapse
df_subset <- df %>% filter(PatientID %in% pts_maxprob_norelapse)

graphs <- plot_per_pt(df_subset)
pdf('plots/top10_maxprob_norelapse.pdf')
graphs
dev.off()

# Min prob, with relapse
df_subset <- df %>% filter(PatientID %in% pts_minprob_withrelapse)

graphs <- plot_per_pt(df_subset)
pdf('plots/top10_minprob_withrelapse.pdf')
graphs
dev.off()



# Get PatientNum from PatientID -------------------------------------------
wkpath  = "raw_data/Box/Boston/MS CLIMB Data/"
MS_map  = read.xlsx(paste0(wkpath,"CLIMB Cohort/Spec95_i2b2_Mapping2017.xlsx"),
                    sheet = 1);  colnames(MS_map) = c("PatientNum","PatientID")

ptnum_maxprob_norelapse <- map(as.list(pts_maxprob_norelapse) %>% set_names(), function(p_id){
  MS_map %>% filter(PatientID == p_id) %>% dplyr::select(PatientNum) %>% unlist(use.names = F)
})
ptnum_maxprob_norelapse <- unlist(ptnum_maxprob_norelapse, use.names = F)

ptnum_minprob_withrelapse <- map(as.list(pts_minprob_withrelapse) %>% set_names(), function(p_id){
  MS_map %>% filter(PatientID == p_id) %>% dplyr::select(PatientNum) %>% unlist(use.names = F)
})
ptnum_minprob_withrelapse <- unlist(ptnum_minprob_withrelapse, use.names = F)

# Create heatmaps ---------------------------------------------------------




# Chart Review Workspace --------------------------------------------------
# Load keys
# wkpath  = "raw_data/Box/Boston/MS CLIMB Data/"
# wkpath2 = "raw_data/Box/Boston/MS CLIMB Paper2 Code&Result/Step3_BuildModel/IntermediateDataResult/"
# ICDPheCode   = read.csv(paste0(wkpath,"EHR/MS_AllEncounters_ICD_Data_03282019.csv"),
#                         stringsAsFactors = FALSE);ICDPheCode$phecode[ICDPheCode$concept_cd == "LPA268"] = "335_"
# CUIdictAll = read.xlsx(paste0(wkpath,"EHR/AllCUI_Database.xlsx"),
#                        sheet = 1); colnames(CUIdictAll) = c("ConceptCd","Desc")
# # Define function
# firstup <- function(x) {
#   substr(x, 1, 1) <- toupper(substr(x, 1, 1))
#   return(x)
# }
# 
# ## Isolate patient data
# p_num <- 35631
# tmp <- CC_comb_val %>% filter(PatientNum == p_num)

## Visualize predictors
tmp2 <- tmp[,-c(1:2, 4:9)]
test <- tmp2 %>% gather(names(tmp2)[-1], key = "code", value = "value")
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
  tmp4[i] <- CUIdictAll[CUIdictAll$ConceptCd == tmp3[i], 'Desc'] %>% unique()
}
test$Code[tmp2] <- tmp4 %>% strtrim(width = 20)
# Clinical variables
tmp2 <- !(strtrim(test$code, 3) == 'CPT' | strtrim(test$code, 3) == 'Phe' | strtrim(test$code, 3) == 'CUI')
test$Code[tmp2] <- test$code[tmp2]

#Plot
pdf('plots/test_pt.pdf')
ggplot(test, aes(x = StartDate, y = Code, fill = value)) +
  geom_tile() + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.y=element_text(size=4),
        axis.text.x = element_text(angle = 90, size = 4)) + 
  ggtitle(paste0("PatientNum: ", p_num)) + 
  scale_fill_gradient2()
dev.off()


