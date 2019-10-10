# DMT data
# Library -----------------------------------------------------------------
library(tidyverse)

# Load --------------------------------------------------------------------
wkpath2 = "raw_data/Box/Boston/MS CLIMB Paper2 Code&Result/Step3_BuildModel/IntermediateDataResult/"
MS_cohort = read.csv(paste0(wkpath2,"Cleaned_MS_cohort.csv"),
                     stringsAsFactors = FALSE,
                     colClasses = c(rep("character",4),rep("Date",4),rep("integer",6),
                                    rep("numeric",3),"integer",rep("Date",4),
                                    rep("numeric",3),"integer"))
MS_attack = read.csv(paste0(wkpath2,"Cleaned_MS_attack.csv"),
                     stringsAsFactors = FALSE,
                     colClasses = c("character",rep("integer",6),"Date",rep("numeric",3),
                                    rep("character",12),rep("Date",2),"integer",
                                    rep("Date",5),rep("integer",2)))
MS_trt    = read.csv(paste0(wkpath2,"Cleaned_MS_treatment.csv"),
                     stringsAsFactors = FALSE)