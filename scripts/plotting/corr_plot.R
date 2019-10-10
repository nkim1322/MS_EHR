# Correlation plots
# Library -----------------------------------------------------------------
library(tidyverse)
library(corrplot)

# Load --------------------------------------------------------------------
CC_comb <- readRDS('modeling_data/train24_3.rds')
CC_comb_val <- readRDS('modeling_data/test24_3.rds')
# Combine
CC_all <- rbind(CC_comb, CC_comb_val)
rm(CC_comb, CC_comb_val)

# ICD PheCode
wkpath  = "raw_data/Box/Boston/MS CLIMB Data/"
wkpath2 = "raw_data/Box/Boston/MS CLIMB Paper2 Code&Result/Step3_BuildModel/IntermediateDataResult/"
ICDPheCode   = read.csv(paste0(wkpath,"EHR/MS_AllEncounters_ICD_Data_03282019.csv"),stringsAsFactors = FALSE)
ICDPheCode$phecode[ICDPheCode$concept_cd == "LPA268"] = "335_"
# CUI
CUIdictAll = read.xlsx(paste0(wkpath,"EHR/AllCUI_Database.xlsx"),sheet = 1)
colnames(CUIdictAll) = c("ConceptCd","Desc")

# Calculate correlations --------------------------------------------------

corr   = cor(CC_all[,-c(1:9)])
# aa  = which(abs(corr - 1) < 1e-8)
# aa  = aa[!aa %in% seq(1, nrow(corr)^2, nrow(corr) + 1)]
# cc  = ceiling(aa/nrow(corr))
# rr  = aa - nrow(corr)*(cc - 1)
# bb  = cc < rr
# aa  = aa[bb]
# cc  = cc[bb]
# rr  = rr[bb]
# if (length(cc) > 0) {
#   for (i in length(cc)) {
#     cat("Perfect correlation between",row.names(corr)[c(cc,rr)],"\n")
#     cat(row.names(corr)[rr],"deleted.\n")
#   }
#   corr    = corr[-rr,-rr]
#   CC_all = CC_all[,-c(rr + 9)]
# }


# CPT codes ---------------------------------------------------------------

indx = grep("CPT",row.names(corr))
tmp  = corr[indx, indx]
row.names(tmp) = colnames(tmp) = substr(row.names(tmp),10, 58)
# Save to PDF 
pdf('plots/corrplot_CPT.pdf')
# corrplot(tmp, title = "CPT Group", tl.cex = 0.5, mar = c(0,0,1,0))
corrplot(tmp, tl.cex = 0.5, mar = c(0,0,1,0))
dev.off()
# Get stats
cat("Among CPT Groups:\n")
summary(c(tmp[upper.tri(tmp)]))
cat("CPT Group vs others:\n")
summary(c(corr[indx,-indx]))



# ICD Phecodes ------------------------------------------------------------

indx = grep("PheCode.", row.names(corr))
tmp  = corr[indx,indx]
phecodes <- substr(row.names(tmp),9,nchar(row.names(tmp)))
for (i in 1:length(phecodes)) {
  phecodes[i] <- (ICDPheCode[ICDPheCode$phecode == phecodes[i], 'phecode_description'] %>% unique())
}
row.names(tmp) = colnames(tmp) = phecodes %>% strtrim(width = 40)
# Save to PDF
pdf('plots/corrplot_PheCode.pdf')
# corrplot(tmp,title = "PheCode", tl.cex = 0.5, mar = c(0,0,1,0))
corrplot(tmp, tl.cex = 0.4, mar = c(0,0,1,0))
dev.off()
# Get stats
cat("Among PheCodes:\n")
summary(c(tmp[upper.tri(tmp)]))
cat("PheCode vs others:\n")
summary(c(corr[indx,-indx]))

# CUI ---------------------------------------------------------------------
indx = grep("CUI.", row.names(corr))
tmp  = corr[indx, indx]
cuis <- substr(row.names(tmp), 5, nchar(row.names(tmp)))
for (i in 1:length(cuis)) {
  cuis[i] <- (CUIdictAll[CUIdictAll$ConceptCd == cuis[i], 'Desc'] %>% unique())
}
# row.names(tmp) = colnames(tmp) = cuis %>% strtrim(width = 15)
row.names(tmp) = colnames(tmp) = cuis 
# aa   = which(tmp > 0.85)
# aa   = aa[!aa %in% seq(1,nrow(tmp)^2,nrow(tmp)+1)]
# cc   = ceiling(aa/nrow(tmp))
# rr   = aa - nrow(tmp)*(cc - 1)
# bb   = cc < rr
# aa   = aa[bb]
# cc   = cc[bb]
# rr   = rr[bb]
# aa   = data.frame(row.names(tmp)[cc],row.names(tmp)[rr],round(tmp[aa],3))
# colnames(aa) = c("CUI","CUI","Corr")
# bb   = sort(aa$Corr,index.return = TRUE)$ix
# aa   = aa[bb,];row.names(aa) = NULL

# Save to PDF
pdf('plots/corrplot_CUI.pdf')
# corrplot(tmp, title = "CUI", tl.cex = 0.5, mar = c(0,0,1,0))
corrplot(tmp, tl.cex = 0.35, mar = c(0,0,1,0))
dev.off()
# Get stats
# cat("Among CUIs:\n")
# summary(c(tmp[upper.tri(tmp)]))
# cat("CUI vs others:\n")
# summary(c(corr[indx,-indx]))


# All codes ---------------------------------------------------------------
## Be wary of running this, as readability is difficult with all 196 codes

# tmp  = corr
# indx1 = grep("CPT", row.names(corr))
# indx2 = grep("PheCode.", row.names(corr))
# # CPT codes get truncated
# row.names(tmp)[indx1] = colnames(tmp)[indx1] = substr(row.names(tmp)[indx1],10,19)
# # PheCodes substitute for phenotypic string (already done in previous section)
# row.names(tmp)[indx2] = colnames(tmp)[indx2] = phecodes %>% strtrim(width = 10)
# # Save to PDF 
# pdf('plots/corrplot_all.pdf')
# corrplot(tmp, title = "All Codes", tl.cex = 0.5, mar = c(0,0,1,0))
# dev.off()
# Get stats
# cat("Among CPT Groups:\n")
# summary(c(tmp[upper.tri(tmp)]))
# cat("CPT Group vs others:\n")
# summary(c(corr[indx,-indx]))

