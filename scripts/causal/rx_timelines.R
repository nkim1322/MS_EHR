# Plot EHR treatment data timelines per patient - evaluating accuracy of data use 
# Source: http://benalexkeen.com/creating-a-timeline-graphic-using-r-and-ggplot2/
# Library -----------------------------------------------------------------
library(tidyverse)
library(scales)
library(lubridate)
library(openxlsx)
library(zoo)

# Set parameters and paths ------------------------------------------------
wkpath  = "raw_data/Box/Boston/MS CLIMB Data/"
wkpath2 = "raw_data/Box/Boston/MS CLIMB Paper2 Code&Result/Step3_BuildModel/IntermediateDataResult/"

tw = 24
tp = 3

# Load --------------------------------------------------------------------
# PatientNum <-> PatientID mapping
MS_map  = read.xlsx(paste0(wkpath,"CLIMB Cohort/Spec95_i2b2_Mapping2017.xlsx"), sheet = 1)
colnames(MS_map) = c("PatientNum","PatientID")

# MS Cohort data for CLIMB inclusion
MS_cohort = read.csv(paste0(wkpath2,"Cleaned_MS_cohort.csv"),stringsAsFactors = FALSE,
                     colClasses = c(rep("character",4), rep("Date",4), rep("integer",6),rep("numeric",3),
                                    "integer",rep("Date",4),rep("numeric",3),"integer"))
# Add PatientNum 
MS_cohort$PatientNum <- sapply(MS_cohort$PatientID, function(pid) {
  MS_map$PatientNum[MS_map$PatientID == pid][1]
})

# CLIMB treatment data
MS_trt    = read.csv(paste0(wkpath2,"Cleaned_MS_treatment.csv"), stringsAsFactors = FALSE)
MS_trt$start_date <- as.Date(with(MS_trt, paste(start_year, start_month, start_day,sep="-")), "%Y-%m-%d")
MS_trt$stop_date <- as.Date(with(MS_trt, paste(stop_year, stop_month, stop_day,sep="-")), "%Y-%m-%d")

# RxNorm data for drugs of interest
RXNORM_comb <- read.xlsx(paste0(wkpath,"MS DMT/MS_RxNorm_data_withRXGenericName.xlsx"), sheet = 5)
RXNORM_comb$start_date <- convertToDate(RXNORM_comb$start_date)

# CUI data
CUISelected  = read.csv(paste0(wkpath,"EHR/MS_AllEncounters_CUIs_08072019.csv"),
                        stringsAsFactors = FALSE, 
                        col.names = c('patient_num', 'encounter_num', 'start_date', 'concept_cd'))
CUISelected$start_date <- format(as.POSIXct(CUISelected$start_date, format = '%Y-%m-%d %H:%M:%S'), format = '%Y-%m-%d')
CUISelected$start_date <- as.Date(CUISelected$start_date)

CUIdictAll = read.xlsx(paste0(wkpath,"EHR/AllCUI_Database.xlsx"), sheet = 1)
colnames(CUIdictAll) = c("ConceptCd","Desc")
# Cleaned treatment data
trt_subset <- readRDS(paste0('modeling_data/causalRN_trt_', tw,'_', tp, 'mons_.rds'))

# Subset data for rituximab/natalizumab patients --------------------------
MS_cohort <- MS_cohort %>% filter(PatientNum %in% trt_subset$PatientNum)

# MS_trt <- MS_trt %>% filter(medication_desc %in% c("RITUXAN", "TYSABRI"))
# Add PatientNum 
MS_trt$PatientNum <- sapply(MS_trt$PatientID, function(pid) {
  MS_map$PatientNum[MS_map$PatientID == pid][1]
})
MS_trt <- MS_trt %>% dplyr::select(PatientID, PatientNum, everything())
MS_trt <- MS_trt %>% filter(!is.na(PatientNum))

RXNORM_comb <- RXNORM_comb %>% 
  filter(patient_num %in% trt_subset$PatientNum) %>% 
  filter(RX_Generic_Name %in% c("Rituximab", "Natalizumab"))

CUISelected <- CUISelected %>% filter(patient_num %in% trt_subset$PatientNum)


# DMT CUI mapping ---------------------------------------------------------

ritux_cui <- CUIdictAll$ConceptCd[str_detect(CUIdictAll$Desc, "itux")][-1]
nataliz_cui <- c(CUIdictAll$ConceptCd[str_detect(CUIdictAll$Desc, "nataliz")], 
                 CUIdictAll$ConceptCd[str_detect(CUIdictAll$Desc, "Tysabri")])
dmt_mapping <- data.frame(CUI = c(ritux_cui, nataliz_cui),
                          Drug = c(rep("Rituximab", length(ritux_cui)),
                                   rep("Natalizumab", length(nataliz_cui))))

# Add to CUI dataframe for future reference
CUISelected <- CUISelected %>% filter(concept_cd %in% c(ritux_cui, nataliz_cui))
CUISelected$Desc <- sapply(CUISelected$concept_cd, function(cui) {
  return(dmt_mapping$Drug[dmt_mapping$CUI == cui])
})


# Legend colors -----------------------------------------------------------
Type_colors <- c("#0070C0", "#00B050", "#FFC000")
Type_levels <- c("RxNorm", "CUI", "CLIMB Validated")


# Plot --------------------------------------------------------------------
graphs_list <- list()
for (pnum in unique(trt_subset$PatientNum)) {
  RXNORM_pt <- RXNORM_comb %>% filter(patient_num == pnum)
  if (nrow(RXNORM_pt) != 0) {
    RXNORM_pt <- RXNORM_pt[order(RXNORM_pt$start_date),]
    RXNORM_pt$Label <- paste(RXNORM_pt$RX_Generic_Name, "RxNorm")
    RXNORM_pt$Type <- "RxNorm"
    tmp1 <- RXNORM_pt[,c("patient_num", "start_date", "Label", "Type")]
    tmp1 <- unique(tmp1)
    colnames(tmp1) <- c("PatientNum", "StartDate", "Label", "Type")
    plot_df <- tmp1
  }
  
  CUISelected_pt <- CUISelected %>% filter(patient_num == pnum)
  if (nrow(CUISelected_pt) != 0) {
    CUISelected_pt <- CUISelected_pt[order(CUISelected_pt$start_date),]
    CUISelected_pt$Label <- paste(CUISelected_pt$Desc, "CUI")
    CUISelected_pt$Type <- "CUI"
    tmp2 <- CUISelected_pt[,c("patient_num", "start_date", "Label", "Type")]
    tmp2 <- unique(tmp2)
    colnames(tmp2) <- c("PatientNum", "StartDate", "Label", "Type")
    
    if (nrow(tmp1) != 0) {
      plot_df <- rbind(plot_df, tmp2)
    } else {
      plot_df <- tmp2
    }
  }
  
  # For CLIMB patients, add validated treatment start
  if (pnum %in% MS_cohort$PatientNum) {
    MS_trt_pt <- MS_trt %>% filter(PatientNum == pnum)
    MS_trt_pt <- MS_trt_pt[order(MS_trt_pt$start_date),]
    if (nrow(MS_trt_pt) != 0) {
      MS_trt_pt <- MS_trt_pt[, c("PatientNum", "start_date", "medication_desc")]
      colnames(MS_trt_pt) <- c("PatientNum", "StartDate", "Label")
      MS_trt_pt$Label <- paste(MS_trt_pt$Label, "Validated")
      MS_trt_pt$Type <- "CLIMB Validated"
      plot_df <- rbind(plot_df, MS_trt_pt)
    }
  }
  
  
  ###################################################################################
  ################################   Tidy Data   ####################################
  ###################################################################################
  plot_df <- plot_df[order(plot_df$StartDate),]
  plot_df$Type <- factor(plot_df$Type)
  
  # Assign lines and heights for milestones within same months to be same 
  positions <- c(0.5, -0.5, 1.0, -1.0, 1.5, -1.5)
  directions <- c(1, -1)
  
  line_pos <- data.frame(
    "StartDate" = unique(plot_df$StartDate),
    "position" = rep(positions, length.out = length(unique(plot_df$StartDate))),
    "direction" = rep(directions, length.out = length(unique(plot_df$StartDate)))
  )
  
  plot_df <- merge(x = plot_df, y = line_pos, by = "StartDate", all = TRUE)
  plot_df <- plot_df[with(plot_df, order(StartDate)), ]
  
  # If there are multiple codes for a given month, slightly alter their positions
  text_offset <- 0.05
  
  plot_df$month_count <- ave(plot_df$StartDate == plot_df$StartDate, plot_df$StartDate, FUN=cumsum)
  plot_df$text_position <- (plot_df$month_count * text_offset * plot_df$direction) + plot_df$position
  
  plot_df$Type <- factor(plot_df$Type, levels = Type_levels, ordered = TRUE)
  
  # Create a data frame for months to display, beginning with +/- 1 month from 1st/last StartDate
  month_buffer <- 1
  
  month_date_range <- seq(from = min(plot_df$StartDate) - lubridate::month(month_buffer), 
                          to = max(plot_df$StartDate) + lubridate::month(month_buffer), 
                          by = 'month')
  month_format <- format(month_date_range, '%b')
  month_df <- data.frame(month_date_range, month_format)
  
  # Do the same for years to display, for which there is a December/January crossover
  year_date_range <- seq(min(plot_df$StartDate) - lubridate::month(month_buffer), 
                         max(plot_df$StartDate) + lubridate::month(month_buffer), 
                         by = 'year')
  year_date_range <- as.Date(
    intersect(
      ceiling_date(year_date_range, unit = "year"),
      floor_date(year_date_range, unit = "year")
    ),  origin = "1970-01-01"
  )
  year_format <- format(year_date_range, '%Y')
  year_df <- data.frame(year_date_range, year_format)
  
  ###################################################################################
  ################################ Plot timeline ####################################
  ###################################################################################
  timeline_plot <- ggplot(plot_df, aes(x = StartDate, y = 0, color = Type, label = Label))
  timeline_plot <- timeline_plot + 
    labs(col = "Code Type") + 
    scale_color_manual(labels = Type_levels,
                       values = c("RxNorm" = "#0070C0",
                                  "CUI" = "#00B050",
                                  "CLIMB Validated" = "#FFC000"), 
                       drop = FALSE) + 
    theme_classic()
  
  # Plot horizontal black line for timeline
  timeline_plot <- timeline_plot + geom_hline(yintercept = 0, color = "black", size = 0.3)
  
  # Plot vertical segment lines for code appearances
  timeline_plot <- timeline_plot + geom_segment(data = plot_df[plot_df$month_count == 1,], 
                                                aes(y = position, yend = 0, xend = StartDate), 
                                                color = 'black', size = 0.2)
  
  # Plot scatter points at zero and date
  timeline_plot <- timeline_plot + geom_point(aes(y = 0), size = 3)
  
  # Don't show axes, appropriately position legend
  timeline_plot <- timeline_plot + theme(axis.line.y = element_blank(),
                                         axis.text.y = element_blank(),
                                         axis.title.x = element_blank(),
                                         axis.title.y = element_blank(),
                                         axis.ticks.y = element_blank(),
                                         axis.text.x = element_blank(),
                                         axis.ticks.x = element_blank(),
                                         axis.line.x = element_blank(),
                                         legend.position = "bottom",
                                         plot.title = element_text(hjust = 0.5))
  
  # Show text for each month
  timeline_plot <- timeline_plot + 
    geom_text(data = month_df, 
              aes(x = month_date_range, y = -0.1, label = month_format),
              size = 1, vjust = 0.5, color = 'black', angle = 90)
  
  # Show text for each year
  if (nrow(year_df) != 0) {
    timeline_plot <- timeline_plot +
      geom_text(data = year_df, 
                aes(x = year_date_range, y = -0.2, label = year_format, fontface = "bold"),
                size = 2.5, color = 'black')
  }
  
  # Show text for each code appearance
  timeline_plot <- timeline_plot + 
    geom_text(data = plot_df, aes(y = text_position), size = 2)
  
  # Add plot title
  start_dte <- as.character(trt_subset$start_date[trt_subset$PatientNum == pnum])
  timeline_plot <- timeline_plot + 
    ggtitle(paste0("PatientNum: ", pnum, 
                  ", Assigned Treatment: ", trt_subset$medication_desc[trt_subset$PatientNum == pnum],
                  ", Assigned Start Date: ", start_dte))

  # Plot assigned treatment start date
  timeline_plot <- timeline_plot + 
    geom_vline(xintercept = as.Date(start_dte), 
               color = "red", size = 0.3,
               linetype = "dashed") + 
    geom_text(data = data.frame(start_date = start_dte), 
              mapping = aes(x = as.Date(start_date), 
                            y = max(plot_df$text_position) - 0.3, 
                            label = start_dte), 
              size = 2.5, angle = 90, vjust = -0.4, hjust = 0, col = 'red') 
  
  # Store plot
  graphs_list[[as.character(pnum)]] <- timeline_plot
}



# Plot and save -----------------------------------------------------------

pdf('plots/causal/rx_timelines.pdf', width = 10, height = 6)
graphs_list
dev.off()





