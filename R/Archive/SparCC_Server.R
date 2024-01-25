#### SparCC - ServerScript ####

setwd("~/SparCC/")

file_ASV <- "Data/Complete/"

### Dependencies ###

library(tidyverse)
library(BiocParallel)

source("Utilities/Datalist_Wrangling_Functions.R")
source("Utilities/SparCC_Dataimport.R")
source("Utilities/SparCC_Wrapper.R")
source("Utilities/Import_Data.R")

options(MulticoreParam=MulticoreParam(workers=10))

# Epi LPA Euk

SparCC_total <- sparCC_dataimport(file_ASV, kingdom = "Euk") %>%

### Data Wrangling ###
  
  mutate_meta_datalist(Depth_Grp = ifelse(Depth <= DCM, "Epi", "Meso")) %>%

  sparCC_wrapper(., envir_filter = Depth_Grp == "Epi" & Size_Fraction == 8, n_boot = 99, frac = F)

### Save Results ###
  
write.csv(SparCC_total$cor, "Output/Complete/Cor_SparCC_Euk_LPA_Epi.csv")
write.csv(SparCC_total$pVal, "Output/Complete/Pval_SparCC_Euk_LPA_Epi.csv")

# Meso LPA Euk

SparCC_total <- sparCC_dataimport(file_ASV, kingdom = "Euk") %>%
  
### Data Wrangling ###
  
  mutate_meta_datalist(Depth_Grp = ifelse(Depth <= DCM, "Epi", "Meso")) %>%
  
  sparCC_wrapper(., envir_filter = Depth_Grp == "Meso" & Size_Fraction == 8, n_boot = 99, frac = F)

### Save Results ###

write.csv(SparCC_total$cor, "Output/Complete/Cor_SparCC_Euk_LPA_Meso.csv")
write.csv(SparCC_total$pVal, "Output/Complete/Pval_SparCC_Euk_LPA_Meso.csv")
