# BPF
# Library -----------------------------------------------------------------
library(tidyverse)
library(openxlsx)

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


# Investigate -------------------------------------------------------------
val_pts <- MS_attack_subs %>% filter(!is.na(val_start)) %>% select(PatientID) %>% unlist(use.names = F) %>% unique()
unique_pts <- unique(MS_attack_subs$PatientID)

graphs_list <- list()
for (p_id in unique_pts) {
  bpf_subset <- bpf %>% filter(subject == p_id)
  MS_subset <- MS_attack_subs %>% filter(PatientID == p_id)
  if (nrow(MS_subset) != 0) {
    graphs_list[[p_id]] <- ggplot() + 
      geom_point(data = bpf_subset, aes(MRI_DATE, BPF)) + 
      geom_smooth(data = bpf_subset, aes(MRI_DATE, BPF)) + 
      geom_vline(data = MS_subset, aes(xintercept = onset, color = 'red')) +
      ggtitle(paste('PatientID: ',p_id)) + 
      ylim(0, 1.5) +
      theme(plot.title = element_text(hjust = 0.5),
            legend.position = "none")
  } else {
    graphs_list[[p_id]] <- ggplot() + 
      geom_point(data = bpf_subset, aes(MRI_DATE, BPF)) + 
      geom_smooth(data = bpf_subset, aes(MRI_DATE, BPF)) + 
      ggtitle(paste('PatientID: ',p_id)) + 
      ylim(0, 1.5) +
      theme(plot.title = element_text(hjust = 0.5),
            legend.position = "none")
  }
}

pdf('thesis_tmp/bpf.pdf')
# pdf('thesis_tmp/bpf_val.pdf')
graphs_list
dev.off()

rm(list = ls())

