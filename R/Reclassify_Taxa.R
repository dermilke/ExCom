reclassify_taxa <- function(datalist, database, sequences, taxonomy_filter,
                            pident_min = 95, mismatch_max = 10, length_min = 250, drop = F) {
  
  reblast_taxonomy <- function(datalist, database, sequences, taxonomy_filter) {
    
    taxonomy_filter <- enquo(taxonomy_filter)
    
    OTU_IDs <- datalist$Count_Data %>%
      filter(!!taxonomy_filter) %>%
      .$OTU_ID
    
    taxonomy_db <- read_tsv(paste0(database, "taxonomy.tsv")) 
    
    blastResults <- map(OTU_IDs, function(x) {
      query <- tibble(OTU_ID = x, 
                      Sequence = system(paste0("grep -A 1 ", x, " ", sequences), intern = T)[2]) %>%
        mutate(OTU_ID = paste0(">", OTU_ID))
      
      do.call(rbind, lapply(seq(nrow(query)), function(i) t(query[i, ]))) %>%
        write.table(., file = "~/PhD/tmp.fasta", row.names = FALSE, col.names = FALSE, quote = FALSE)
      
      system2("blastn",
              args = c("-db", paste0(database, "Blast_DB"),
                       "-query", "~/PhD/tmp.fasta",
                       "-outfmt", 6), wait = T, stdout = T) %>%
        tibble::enframe(., name = NULL) %>%
        separate(col = value, sep = "\t", convert = T,
                 into = c("OTU_ID", "Feature_ID", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart",
                          "send", "evalue", "bitscore")) %>%
        dplyr::slice(1:50) %>%
        left_join(., taxonomy_db, by = "Feature_ID")
      
    })
    
    grep("8069#AGNK01006210", read_file(sequences))
    
    return(blastResults)
    
  }
  
  taxonomy_filter <- enquo(taxonomy_filter)
  
  new_taxonomy <- reblast_taxonomy(datalist, database, sequences, !!taxonomy_filter) %>%
    map(., function(x) {
      x %>%
        filter(pident >= pident_min & mismatch <= mismatch_max & length >= length_min) %>%
        select(-c(Feature_ID, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore)) %>%
        dplyr::slice(1)
    }) %>%
    bind_rows()
  
  table_merge <- datalist$Count_Data %>%
    filter(OTU_ID %in% new_taxonomy$OTU_ID) %>%
    group_by(OTU_ID) %>%
    select_if(is.numeric) %>%
    ungroup() %>%
    with(., right_join(new_taxonomy, ., by = "OTU_ID")) %>%
    select_if(names(.) %in% names(datalist$Count_Data))
  
  table_result <- datalist$Count_Data %>%
    filter(!(OTU_ID %in% table_merge$OTU_ID)) %>%
    select_if(names(.) %in% names(table_merge)) %>%
    bind_rows(., table_merge)
  
  if (drop) {
    table_result <- table_result %>%
      filter(!(OTU_ID %in% {filter(datalist$Count_Data, !!taxonomy_filter) %>% filter(!(OTU_ID %in% new_taxonomy$OTU_ID)) %>% .$OTU_ID}))
  }
  
  return(list(Count_Data = table_result, Meta_Data = datalist$Meta_Data))
  
}