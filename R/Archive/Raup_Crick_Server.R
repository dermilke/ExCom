library(tidyverse)
library(doParallel)

bootstrap_reps = 9999
parallel_cores = 6

datalist <- import_data("data/", kingdom = "Prok", min_counts = 2000, abundance_filter = T) 

rc_FL <- datalist %>% 
  filter_station_datalist(colSums(select_if(datalist$Count_Data, is.numeric) > 0) >= 2) %>%
  filter_station_datalist(Size_Fraction == 0.22) %>%
  
  raup_crick_abundance(., classic_metric = F, split_ties = T, reps = bootstrap_reps, set_all_species_equal = F,
                       as.distance.matrix = F, report_similarity = F, use.cores = parallel_cores)

write.csv(rc_FL, "output/Raup_Crick_Prok_FL.csv")

rc_SPA <- datalist %>% 
  filter_station_datalist(colSums(select_if(datalist$Count_Data, is.numeric) > 0) >= 2) %>%
  filter_station_datalist(Size_Fraction == 3) %>%
  raup_crick_abundance(., classic_metric = F, split_ties = T, reps = bootstrap_reps, set_all_species_equal = F,
                       as.distance.matrix = F, report_similarity = F, use.cores = parallel_cores)

write.csv(rc_SPA, "output/Raup_Crick_Prok_SPA.csv")

rc_LPA <- datalist %>% 
  filter_station_datalist(colSums(select_if(datalist$Count_Data, is.numeric) > 0) >= 2) %>%
  filter_station_datalist(Size_Fraction == 8) %>%
  raup_crick_abundance(., classic_metric = F, split_ties = T, reps = bootstrap_reps, set_all_species_equal = F,
                       as.distance.matrix = F, report_similarity = F, use.cores = parallel_cores)

write.csv(rc_LPA, "output/Raup_Crick_Prok_LPA.csv")
