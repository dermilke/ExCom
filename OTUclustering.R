##### Creating OTUs based on aligned sequence distances #####

rm(list = ls(all=T))

library(tidyverse)

source("datalistFunctions.R")
source("datatableFunctions.R")

# Number of Processors
nproc <- 4

# Use DECIPHER to read DNA sequences, align them and calculate Hamming-Distances between sequences. Finally cut-off tree
# at 3% Dissimilarity = 97% OTUs to form OTU cluster.
# Use highest abundance ASV for taxonomic classification of that OTU.
dna <- Biostrings::readDNAStringSet("~/PhD/Data_Storage/Fasta/Complete/V4V5_Primerset/Prok/Sequences_Pacific_Complete.fasta")
aln <- DECIPHER::AlignSeqs(dna, processors = nproc)
d <- DECIPHER::DistanceMatrix(aln, processors = nproc)
clusters_97 <- DECIPHER::IdClusters(
  d, 
  method = "complete",
  cutoff = 0.03,
  processors = nproc)

clusters_df <- tibble(OTU_ID = rownames(clusters_97), OTU_cluster = clusters_97$cluster)

# Safe OTU-ASV-Cluster Assignments for later.
write_csv(clusters_df, "~/PhD/SoftwareBuilds/ExCom/OTU_Cluster_IDs.csv")

# Transform ASV-Datatable into 97% OTU table using calculated OTU cluster
datalist <- data_select("~/PhD/Data_Storage/ASV_Tabs/Pacific/V4V5_Primerset/Complete/", "Prok")

datalist_97 <- datalist

datalist_97$Count_Data <- datalist_97$Count_Data %>%
  mutate(OTU_cluster = clusters_df$OTU_cluster[match(.$OTU_ID, clusters_df$OTU_ID)])

merged_taxonomy <- datalist_97$Count_Data %>%
  select_if(is.character) %>%
  mutate(OTU_cluster = datalist_97$Count_Data$OTU_cluster) %>%
  mutate(Abundance = rowSums(select_if(datalist_97$Count_Data, is.numeric))) %>%
  group_by(OTU_cluster) %>%
  filter(Abundance == max(Abundance)) %>%
  arrange(OTU_cluster, Abundance) %>%
  .[!duplicated(.$OTU_cluster),]

datalist_97$Count_Data <- datalist_97$Count_Data %>%
  group_by(OTU_cluster) %>%
  select_if(is.numeric) %>%
  summarize_all(sum) %>%
  left_join(merged_taxonomy, ., by = "OTU_cluster") %>%
  ungroup() %>%
  select(-OTU_cluster, -Abundance)

write_tsv(datalist_97$Count_Data, "~/PhD/Data_Storage/ASV_Tabs/Pacific/V4V5_Primerset/OTU_97pct/Complete/Processed/Prok/Full_Prok_Count.tsv")