read_taxonomy <- function(file, kingdom, Prefix = "D_0__", DB = "SILVA", fillEmpty = T) {
  
  fill_empty <- function(x, entity) {
    
    for (i in (3:ncol(select_if(x, is.character)))) {
      
      x[is.na(x[,i]),i] <- x[is.na(x[,i]),i-1] %>%
        pull(.) %>%
        paste(entity, ., sep = " ")
      
      x[,i] <- gsub(paste(entity, " ", entity, " ", sep = ""), paste(entity, " ", sep = ""), pull(x[,i]))
      
    }
    
    return(x)
    
  }
  
  if (kingdom == "Prok") {
    
    taxonomy_prok <- suppressWarnings(
      suppressMessages(read_tsv(file)) %>%
        dplyr::slice(c(grep(paste0(Prefix, "Bacteria"),
                            .$taxonomy, value = F),
                       grep(paste0(Prefix, "Archaea"),
                            .$taxonomy, value = F))) %>%
        separate(., 2, c("Kingdom","Phylum","Class","Order","Family","Genus","Species"), sep = ";") %>%
        mutate_if(is.character, 
                  str_replace_all, pattern = "D_[0-9]*__",
                  replacement = "") %>%
        dplyr::rename(., 'OTU_ID' = '#OTUID') %>%
        with(., if (fillEmpty) {fill_empty(., "unknown")} else {.})
    )
    
    taxonomy_chloroplast <- suppressWarnings(
      suppressMessages(read_tsv(file)) %>%
        dplyr::slice(grep("Kingdom.Eukaryota", .$taxonomy, value = F)) %>%
        separate(., 2, c("Kingdom","Supergroup","Phylum","Class","Subclass","Order","Suborder", "Family", "Genus", "Species"), sep = ";") %>%
        mutate_if(is.character, 
                  str_replace_all, pattern = ".*\\.", replacement = "") %>%
        dplyr::rename(., 'OTU_ID' = '#OTUID') %>%
        with(., if (fillEmpty) {fill_empty(., "unknown")} else {.})
    )
    
    taxonomy <- list(Prok = taxonomy_prok, Chloroplast = taxonomy_chloroplast)
    
  } else if (kingdom == "Euk") {
    
    if (DB == "SILVA") {
      
      taxonomy <- suppressWarnings(
        suppressMessages(read_tsv(file)) %>%
          separate(., 2, c("Kingdom","Phylum","Class","Order","Family","Genus","Species"), sep = ";") %>%
          mutate_if(is.character, 
                    str_replace_all, pattern = "D_[0-9]*__", replacement = "") %>%
          dplyr::rename(., 'OTU_ID' = '#OTUID') %>%
          with(., if (fillEmpty) {fill_empty(., "unknown")} else {.})
      )
      
    } else if (DB == "PR2") {
      
      taxonomy <- suppressWarnings(
        suppressMessages(read_tsv(file)) %>%
          separate(., 2, c("Kingdom","Supergroup","Phylum","Class","Order", "Family", "Genus", "Species"), sep = ";") %>%
          mutate_if(is.character, 
                    str_replace_all, pattern = "[a-z]*_", replacement = "") %>%
          dplyr::rename(., 'OTU_ID' = '#OTUID') %>%
          with(., if (fillEmpty) {fill_empty(., "unknown")} else {.})
      )
      
    }
  }
  
  return(taxonomy)
  
}

prepare_raw <- function(file_ASV, file_Meta, file_Fasta = NULL, confidence_lvl = 0.8, kingdom = "Prok", DB = "SILVA132", fill = T) {
  
  if (kingdom == "Euk") {
    
    taxonomy <- read_taxonomy(paste(file_ASV, "Raw/", kingdom, "/taxonomy-", DB, ".tsv", sep = ""), kingdom, DB, fill)
    
    Count_Data <- suppressMessages(read_delim(paste(file_ASV, "Raw/", kingdom, "/all-18S-seqs.with-", DB, "-tax.tsv", sep = ""),
                                              delim = "\t", skip = 1)) %>%
      dplyr::rename(., 'OTU_ID' = '#OTU ID') %>%
      right_join(taxonomy, ., by = 'OTU_ID') %>%
      filter(confidence >= confidence_lvl) %>%                                                           
      select(-confidence) %>%
      select_if(!colnames(.) == "taxonomy") %>%
      filter(rowSums(select_if(.,is.numeric)) != 0)
    
    write_tsv(Count_Data, paste(file_ASV, "Processed/Euk/Full_Euk_Count.tsv", sep = ""))
    
    Meta_Data <- suppressMessages(read_delim(file_Meta, del = "\t")) %>% 
      dplyr::slice(match(names(select_if(Count_Data, is.numeric)), .$Sample_ID))
    
    write_tsv(Meta_Data, paste(file_ASV, "Meta_Data/Euk/Meta_Data.tsv", sep = ""))
    
    if (!is.null(file_Fasta)) {
      tmp_seqs <- readLines(file_Fasta)
      
      fasta_table <- tibble(Seq_ID = tmp_seqs[seq(1, length(tmp_seqs)-1, 2)],
                            Sequence = tmp_seqs[seq(2, length(tmp_seqs), 2)]) %>%
        mutate(Seq_ID = str_replace_all(Seq_ID, pattern = ">", replacement = "")) %>%
        dplyr::slice(match(Count_Data$OTU_ID, Seq_ID))
      
      if (nrow(fasta_table) != nrow(Count_Data)) {
        warning(paste0("Fasta file not complete: Expected ", nrow(Count_Data), 
                       " sequences, but ", nrow(fasta_table), " sequences matched to count-table."))
      }
      
      fasta_table %>%
        select(Seq_ID, Sequence) %>%
        mutate(Seq_ID = paste0(">", Seq_ID)) %>%
        with(., do.call(rbind, lapply(seq(nrow(.)), function(i) t(.[i, ])))) %>%
        write.table(., file = "Fasta/Euk/Full_Euk_Sequences.fasta", row.names = FALSE, col.names = FALSE, quote = FALSE)
    }
    
  } else {
    
    taxonomy <- read_taxonomy(paste(file_ASV, "Raw/Prok/taxonomy.tsv", sep = ""), kingdom, fillEmpty = fill)
    
    Count_Data_Prok <- suppressMessages(read_delim(paste(file_ASV, "Raw/", kingdom, "/all-16S-seqs.with-tax.tsv", sep = ""), 
                                              delim = "\t", skip = 1)) %>%
      dplyr::rename(., 'OTU_ID' = '#OTU ID') %>%
      left_join(taxonomy$Prok, ., by = 'OTU_ID') %>%
      filter(confidence >= confidence_lvl) %>%                                                           
      select(-confidence) %>%
      select_if(!colnames(.) == "taxonomy") %>%
      filter(rowSums(select_if(.,is.numeric)) != 0)
    
    write_tsv(Count_Data_Prok, paste(file_ASV, "Processed/Prok/Full_Prok_Count.tsv", sep = ""))
    
    Meta_Data_Prok <- suppressMessages(read_delim(file_Meta, del = "\t")) %>% 
      dplyr::slice(match(names(select_if(Count_Data_Prok, is.numeric)), .$Sample_ID))
    
    write_tsv(Meta_Data_Prok, paste(file_ASV, "Meta_Data/Prok/Meta_Data.tsv", sep = ""))
    
    Count_Data_Chlo <- suppressMessages(read_delim(paste(file_ASV, "Raw/", kingdom, "/all-16S-seqs.with-tax.tsv", sep = ""), 
                                              delim = "\t", skip = 1)) %>%
      dplyr::rename(., 'OTU_ID' = '#OTU ID') %>%
      left_join(taxonomy$Chloroplast, ., by = 'OTU_ID') %>%
      filter(confidence >= confidence_lvl) %>%                                                           
      select(-confidence) %>%
      select_if(!colnames(.) == "taxonomy") %>%
      filter(rowSums(select_if(.,is.numeric)) != 0)
    
    write_tsv(Count_Data_Chlo, paste(file_ASV, "Processed/Chloroplast/Full_Chloroplast_Count.tsv", sep = ""))
    
    Meta_Data_Chlo <- suppressMessages(read_delim(file_Meta, del = "\t")) %>% 
      dplyr::slice(match(names(select_if(Count_Data_Chlo, is.numeric)), .$Sample_ID))
    
    write_tsv(Meta_Data_Chlo, paste(file_ASV, "Meta_Data/Chloroplast/Meta_Data.tsv", sep = ""))
    
    if (!is.null(file_Fasta)) {
      tmp_seqs <- readLines(file_Fasta)
      
      fasta_table <- tibble(Seq_ID = tmp_seqs[seq(1, length(tmp_seqs)-1, 2)],
                            Sequence = tmp_seqs[seq(2, length(tmp_seqs), 2)]) %>%
        mutate(Seq_ID = str_replace_all(Seq_ID, pattern = ">", replacement = "")) 
      
      fasta_table_Prok <- fasta_table %>%
        dplyr::slice(match(Count_Data_Prok$OTU_ID, Seq_ID))
      
      fasta_table_Chlo <- fasta_table %>%
        dplyr::slice(match(Count_Data_Chlo$OTU_ID, Seq_ID))
      
      if ((nrow(fasta_table_Prok) != nrow(Count_Data_Prok)) | (nrow(fasta_table_Chlo) != nrow(Count_Data_Chlo))) {
        warning(paste0("Fasta file not complete: Expected ", nrow(Count_Data_Prok), 
                       " sequences in Prok and ", nrow(Count_Data_Chlo), " in Chloroplast data, but ", nrow(fasta_table_Prok), 
                       " sequences in Prok and ", nrow(fasta_table_Chlo), " in Chloroplast fasta matched to count-table."))
      }
      
      dir.create(paste0(file_ASV, "Fasta/Prok/"), recursive = T)
      dir.create(paste0(file_ASV, "Fasta/Chloroplast/"), recursive = T)
      
      fasta_table_Prok %>%
        select(Seq_ID, Sequence) %>%
        mutate(Seq_ID = paste0(">", Seq_ID)) %>%
        with(., do.call(rbind, lapply(seq(nrow(.)), function(i) t(.[i, ])))) %>%
        write.table(., file = paste0(file_ASV, "Fasta/Prok/Full_Prok_Sequences.fasta"), row.names = FALSE, col.names = FALSE, quote = FALSE)
      
      fasta_table_Chlo %>%
        select(Seq_ID, Sequence) %>%
        mutate(Seq_ID = paste0(">", Seq_ID)) %>%
        with(., do.call(rbind, lapply(seq(nrow(.)), function(i) t(.[i, ])))) %>%
        write.table(., file = paste0(file_ASV, "Fasta/Chloroplast/Full_Chloroplast_Sequences.fasta"), row.names = FALSE, col.names = FALSE, quote = FALSE)
    }
    
  }
  
  cat("Raw tables converted into ", paste(file_ASV, "\n", sep = ""))
  
}
