### Workflow Example ExCom Pipeline

# 0. Load dependencies
# 1. Process raw bioinformatic datalists
# 2. Read processed datalists
# 3. Create datatables from datalists
# 4. Create diversitytables from datalist (also included in datatables)
# 5. Example: Distance_Decay-Analysis
# 6. Example: Diversity Analysis 

##### 0. Load dependencies #####

rm(list = ls(all=T))

library(tidyverse)

source("datalistFunctions.R")
source("datatableFunctions.R")

##### 1. Transform raw tables #####

prepare_raw("Data/Example_Data/Pacific-Project/V4V5_Primerset/Pool_2/",
            "Data/Example_Data/Pacific-Project/Meta_Data/Meta_Data_Pacific_Sample.tsv",
            kingdom = "Prok", confidence_lvl = 0.8, fillEmpty = T)

##### 2. Read processed tables #####

datalist <- data_select("Data/Example_Data/Pacific-Project/V4V5_Primerset/Pool_2/", "Prok")

##### 3. Create tidy datatables #####

datatable <- datalist %>%
  create_datatable(., grpBy = Class, otherThreshold = 0.005)

##### 4. Create diversitytables only #####

diversitytable <- datalist %>%
  rarefy_datalist(., rare_lim = 2000, drop = T) %>%
  diversity_datatable(.)

##### 5. Example: Distance_Decay-Analysis #####

distance_decay(datalist, dist_measure = "Bray_Curtis") %>%
  
  ggplot(., aes(x = Distance, y = Bray_Curtis, col = From_Province)) +
    geom_point() +
    facet_wrap(~To_Depth)

##### 6. Example: Diversity Analysis #####

diversitytable %>%
  
  ggplot(., aes(x = Pot_Temperature, y = Richness)) +
    geom_point()

##### 7. Example: Get Proportion of Rhodobacteraceae and plot them against Province #####

datalist %>%
  make_proportion_datalist(.) %>%
  filter_by_taxa(., 5, "Rhodobacteraceae") %>%
  create_datatable(., grpBy = Genus) %>%
  
  ggplot(., aes(x = Group, y = Abundance, col = Province)) +
    geom_boxplot()
