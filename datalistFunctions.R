### Functions for Handling Datalist Objects ###

remove_NAs <- function(Count_Data, replace = 0) {

  Count_Data[is.na(Count_Data)] <- replace

}

make_proportion <- function(x) {
  result <- x/sum(x, na.rm = T)
  return(result)
}

#make_proportion_datalist <- function(datalist) {
#  
#  datalist$Count_Data <- datalist$Count_Data %>%
#    mutate_if(is.numeric, make_proportion) 
#  
#  return(datalist)
#  
#}

fill_empty <- function(x, entity) {
  
  for (i in (3:ncol(select_if(x, is.character)))) {
    
    x[is.na(x[,i]),i] <- x[is.na(x[,i]),i-1] %>%
                           pull(.) %>%
                           paste(., entity)
    
    x[,i] <- gsub(paste(" ",entity, " ", entity, sep = ""), paste(" ",entity, sep = ""), pull(x[,i]))
  
  }
  
  return(x)
  
}

read_taxonomy <- function(file, kingdom, DB, fillEmpty = T) {
  
  if (kingdom == "Prok") {
    
    taxonomy_prok <- suppressWarnings(
      suppressMessages(read_tsv(file)) %>%
        slice(c(grep("D_0__Bacteria", .$taxonomy, value = F),
                grep("D_0__Archaea", .$taxonomy, value = F))) %>%
        separate(., 2, c("Kingdom","Phylum","Class","Order","Family","Genus","Species"), sep = ";") %>%
        mutate_if(is.character, 
                  str_replace_all, pattern = "D_[0-9]*__", replacement = "") %>%
        rename(., 'OTU_ID' = '#OTUID') %>%
        with(., if (fillEmpty) {fill_empty(., "bacterium")} else {.})
    )
    
    taxonomy_chloroplast <- suppressWarnings(
      suppressMessages(read_tsv(file)) %>%
        slice(grep("Kingdom.Eukaryota", .$taxonomy, value = F)) %>%
        separate(., 2, c("Kingdom","Supergroup","Phylum","Class","Subclass","Order","Suborder", "Family", "Genus", "Species"), sep = ";") %>%
        mutate_if(is.character, 
                  str_replace_all, pattern = ".*\\.", replacement = "") %>%
        rename(., 'OTU_ID' = '#OTUID') %>%
        with(., if (fillEmpty) {fill_empty(., "organism")} else {.})
    )
    
    taxonomy <- list(Prok = taxonomy_prok, Chloroplast = taxonomy_chloroplast)
    
  } else if (kingdom == "Euk") {
    
    if (DB == "SILVA132") {
      
      taxonomy <- suppressWarnings(
        suppressMessages(read_tsv(file)) %>%
          separate(., 2, c("Kingdom","Phylum","Class","Order","Family","Genus","Species"), sep = ";") %>%
          mutate_if(is.character, 
                    str_replace_all, pattern = "D_[0-9]*__", replacement = "") %>%
          rename(., 'OTU_ID' = '#OTUID') %>%
          with(., if (fillEmpty) {fill_empty(., "organism")} else {.})
      )
      
    } else if (DB == "PR2") {
      
      taxonomy <- suppressWarnings(
          suppressMessages(read_tsv(file)) %>%
          separate(., 2, c("Kingdom","Supergroup","Phylum","Class","Subclass","Order","Suborder", "Family", "Genus", "Species"), sep = ";") %>%
          mutate_if(is.character, 
                    str_replace_all, pattern = "[a-z]*_", replacement = "") %>%
          rename(., 'OTU_ID' = '#OTUID') %>%
            with(., if (fillEmpty) {fill_empty(., "organism")} else {.})
      )
      
    }
  }
  
  return(taxonomy)

}

prepare_raw <- function(file_ASV, file_Meta, confidence_lvl = 0.8, kingdom = "Prok", DB = "SILVA132", fillEmpty = T) {
  
  if (kingdom == "Euk") {
    
    taxonomy <- read_taxonomy(paste(file_ASV, "Raw/", kingdom, "/taxonomy-", DB, ".tsv", sep = ""), kingdom, DB, fillEmpty)
    
    Count_Data <- suppressMessages(read_delim(paste(file_ASV, "Raw/", kingdom, "/all-18S-seqs.with-", DB, "-tax.tsv", sep = ""),
                                              delim = "\t", skip = 1)) %>%
      rename(., 'OTU_ID' = '#OTU ID') %>%
      right_join(taxonomy, ., by = 'OTU_ID') %>%
      filter(confidence >= confidence_lvl) %>%                                                           
      select(-confidence) %>%
      select_if(!colnames(.) == "taxonomy") %>%
      filter(rowSums(select_if(.,is.numeric)) != 0)
    
    write_tsv(Count_Data, paste(file_ASV, "Processed/Euk/Full_Euk_Count.tsv", sep = ""))
    
    Meta_Data <- suppressMessages(read_delim(file_Meta, del = "\t")) %>% 
      slice(match(names(select_if(Count_Data, is.numeric)), .$Sample_ID))
    
    write_tsv(Meta_Data, paste(file_ASV, "Meta_Data/Euk/Meta_Data.tsv", sep = ""))

  } else {
    
    taxonomy <- read_taxonomy(paste(file_ASV, "Raw/", kingdom, "/taxonomy.tsv", sep = ""), kingdom, fillEmpty = fillEmpty)
    
    Count_Data <- suppressMessages(read_delim(paste(file_ASV, "Raw/", kingdom, "/all-16S-seqs.with-tax.tsv", sep = ""), 
                                              delim = "\t", skip = 1)) %>%
      rename(., 'OTU_ID' = '#OTU ID') %>%
      left_join(taxonomy$Prok, ., by = 'OTU_ID') %>%
      filter(confidence >= confidence_lvl) %>%                                                           
      select(-confidence) %>%
      select_if(!colnames(.) == "taxonomy") %>%
      filter(rowSums(select_if(.,is.numeric)) != 0)
    
    write_tsv(Count_Data, paste(file_ASV, "Processed/Prok/Full_Prok_Count.tsv", sep = ""))
    
    Meta_Data <- suppressMessages(read_delim(file_Meta, del = "\t")) %>% 
      slice(match(names(select_if(Count_Data, is.numeric)), .$Sample_ID))
    
    write_tsv(Meta_Data, paste(file_ASV, "Meta_Data/Prok/Meta_Data.tsv", sep = ""))
    
    Count_Data <- suppressMessages(read_delim(paste(file_ASV, "Raw/", kingdom, "/all-16S-seqs.with-tax.tsv", sep = ""), 
                                              delim = "\t", skip = 1)) %>%
      rename(., 'OTU_ID' = '#OTU ID') %>%
      left_join(taxonomy$Chloroplast, ., by = 'OTU_ID') %>%
      filter(confidence >= confidence_lvl) %>%                                                           
      select(-confidence) %>%
      select_if(!colnames(.) == "taxonomy") %>%
      filter(rowSums(select_if(.,is.numeric)) != 0)
    
    write_tsv(Count_Data, paste(file_ASV, "Processed/Chloroplast/Full_Chloroplast_Count.tsv", sep = ""))
    
    Meta_Data <- suppressMessages(read_delim(file_Meta, del = "\t")) %>% 
      slice(match(names(select_if(Count_Data, is.numeric)), .$Sample_ID))
    
    write_tsv(Meta_Data, paste(file_ASV, "Meta_Data/Chloroplast/Meta_Data.tsv", sep = ""))
    
  }
  
  cat("Raw tables converted into ", paste(file_ASV, "\n", sep = ""))
  
}

data_select <- function(file_ASV, kingdom = "Prok") {
  
  Count_Data <- suppressMessages(read_delim(paste(file_ASV, "Processed/", kingdom, "/Full_", kingdom,"_Count.tsv", sep = ""), 
                                            del = "\t")) %>%
    select(-grep("[Mm]ock", names(.))) %>%
    select(-grep("NC", names(.)))
  
  Meta_Data <- suppressMessages(read_delim(paste(file_ASV, "Meta_Data/", kingdom, "/Meta_Data.tsv", sep = ""), 
                                           del = "\t"))
  
  return(list(Meta_Data = Meta_Data, 
              Count_Data = Count_Data))
  
}

combine_data <- function(datalist_1, datalist_2) {
  
  Count_Data <- full_join(datalist_1$Count_Data, datalist_2$Count_Data, 
                          by = names(select_if(datalist_1$Count_Data, is.character))) %>%
    mutate_if(is.numeric,
              replace_na, replace = 0)
  
  replaceX <- grep("\\.x$",names(Count_Data)) 
  replaceY <- grep("\\.y$",names(Count_Data)) 
  
  replaceX <- replaceX[match(sub(".y","",names(Count_Data)[replaceY]), 
                             sub(".x","",names(Count_Data)[replaceX]))]
  
  if (!is_empty(replaceX)) {
    
    for (i in 1:length(replaceX)) {
      Count_Data[,replaceX[i]] <- Count_Data[,replaceX[i]] + Count_Data[,replaceY[i]]
    }
    
    Count_Data <- Count_Data %>%
      select(-replaceY)
    
    names(Count_Data)[replaceX] <- sub("\\.x$", "", names(Count_Data)[replaceX])
    
  }
  
  Count_Data <- Count_Data %>%
    select_if(is.numeric) %>%
    select_if(colSums(.) > 0) %>%
    cbind(select_if(Count_Data, is.character), .) %>%
    as_tibble()
  
  Meta_Data <- bind_rows(datalist_1$Meta_Data, datalist_2$Meta_Data) %>%
    distinct() %>%
    slice(match(names(select_if(Count_Data, is.numeric)), .$Sample_ID))
  
  return(list(Count_Data = Count_Data, Meta_Data = Meta_Data))
  
}

filter_abundance <- function(datalist) {
  
  counts <- datalist$Count_Data
  prop <- datalist$Count_Data %>%
    mutate_if(is.numeric,
              make_proportion)
  
  counts_filtered <- counts %>%
    filter(((rowSums(select_if(counts, is.numeric))/sum(rowSums(select_if(counts, is.numeric)))) > 0.00001) &
            (apply(select_if(prop, is.numeric),1,max) > 0.01) |
            ((apply(select_if(prop, is.numeric) > 0.001, 1, sum) / ncol(select_if(prop, is.numeric))) > 0.02) |
            ((apply(select_if(prop, is.numeric) > 0, 1, sum) / ncol(select_if(prop, is.numeric))) > 0.05))
  
  datalist$Count_Data <- counts_filtered
  
  return(datalist)
  
}

rarefy_datalist <- function(datalist, rare_lim, drop = F) {
  
  count_rared <- datalist$Count_Data %>%
    select_if(is.numeric) %>%
    select_if(colSums(.) >= ifelse(drop, rare_lim, 0)) %>%
    t() %>% 
    vegan::rrarefy(., rare_lim) %>%
    t() %>%
    as_tibble() %>%
    bind_cols(select_if(datalist$Count_Data, is.character), .) %>%
    filter(rowSums(select_if(., is.numeric)) > 0)
  
  meta_subset <- datalist$Meta_Data %>%
    slice(match(names(select_if(count_rared, is.numeric)), datalist$Meta_Data$Sample_ID))
  
  datalist$Count_Data <- count_rared
  datalist$Meta_Data <- meta_subset
  
  return(datalist)
  
}

summarize_by_taxa <- function(datalist, tax_lvl = 7) {
  
  datalist$Count_Data <- datalist$Count_Data %>%
    group_by_at(2:(tax_lvl+1)) %>%
    summarize_if(is.numeric, sum) %>%
    ungroup()
  
  return(datalist)
  
}

make_depth_groups <- function(datalist, depth_vek) {
  
  datalist$Meta_Data <- datalist$Meta_Data %>%
    mutate(Depth_Grp = depth_vek[findInterval(datalist$Meta_Data$Depth, depth_vek)])
  
  return(datalist)
  
}

filter_station_datalist <- function(.datalist, ...) {
  
  exprFilter <- enexprs(...)
  
  .datalist$Meta_Data <- .datalist$Meta_Data %>%
    filter(!!! exprFilter)
  
  .datalist$Count_Data <- .datalist$Count_Data %>%
    select_if(is.numeric) %>%
    select_if(names(.) %in% .datalist$Meta_Data$Sample_ID) %>%
    bind_cols(select_if(.datalist$Count_Data, is.character), .) %>%
    filter(select_if(., is.numeric) %>% rowSums(.) > 0)

  return(.datalist)
  
}

filter_taxa_datalist <- function(.datalist, ...) {
  
  exprFilter <- enquos(...)
  
  .datalist$Count_Data <- filter(.datalist$Count_Data, !!!(exprFilter))
    
  
  return(.datalist)
  
}

slice_station_datalist <- function(.datalist, ...) {

  exprSlice <- enexprs(...)
  
  .datalist$Count_Data <- .datalist$Count_Data %>%
      select_if(is.numeric) %>%
      select(!!! exprSlice) %>%
      bind_cols(select_if(.datalist$Count_Data, is.character), .) %>%
      filter(select_if(., is.numeric) %>% rowSums(.) > 0)
    
  .datalist$Meta_Data <- .datalist$Meta_Data %>%
      slice(!!! exprSlice)
    
  return(.datalist)
  
}

mutate_count_datalist <- function(.datalist, func, ...) {
  
  exprFunc <- rlang::enexpr(func)
  arguments <- rlang::list2(...)
  
  .datalist$Count_Data <- .datalist$Count_Data %>%
    select_if(is.numeric) %>%
    mutate_all(eval(exprFunc), !!! arguments) %>%
    cbind(select_if(.datalist$Count_Data, is.character), .) %>%
    as_tibble()
  
  return(.datalist)
  
}

mutate_meta_datalist <- function(.datalist, ...) {
  
  exprFunc <- enquos(...)
  
  .datalist$Meta_Data <- .datalist$Meta_Data %>%
    mutate(!!!exprFunc)
  
  return(.datalist)
  
}

NMDS_ordination_datalist <- function(datalist) {
  
  nmds_result <- datalist$Count_Data %>%
    select_if(is.numeric) %>%
    t(.) %>%
    vegan::vegdist(.) %>%
    as.matrix(.) %>%
    vegan::metaMDS(try = 30)
  
  return(nmds_result)
  
}

sizefraction_communities <- function(datalist) {
  
  datalist_022 <- datalist %>%
    filter_station_datalist(Size_Fraction == "0.22") %>%
    filter_taxa_datalist(!(OTU_ID %in% filter_station_datalist(datalist, Size_Fraction == "3")$Count_Data$OTU_ID | 
                           OTU_ID %in% filter_station_datalist(datalist, Size_Fraction == "8")$Count_Data$OTU_ID))
  
  datalist_3 <- datalist %>%
    filter_station_datalist(Size_Fraction == "3") %>%
    filter_taxa_datalist(!(OTU_ID %in% filter_station_datalist(datalist, Size_Fraction == "0.22")$Count_Data$OTU_ID | 
                           OTU_ID %in% filter_station_datalist(datalist, Size_Fraction == "8")$Count_Data$OTU_ID))
  
  datalist_8 <- datalist %>%
    filter_station_datalist(Size_Fraction == "8") %>%
    filter_taxa_datalist(!(OTU_ID %in% filter_station_datalist(datalist, Size_Fraction == "0.22")$Count_Data$OTU_ID | 
                           OTU_ID %in% filter_station_datalist(datalist, Size_Fraction == "3")$Count_Data$OTU_ID))
  
  return(list(Fraction_022 = datalist_022,
              Fraction_3 = datalist_3,
              Fraction_8 = datalist_8))
  
}

#decostand_datalist <- function(datalist, method) {
#  
#  datalist$Count_Data <- datalist$Count_Data %>%
#    select_if(is.numeric) %>%
#    vegan::decostand(method = method, MARGIN = 2) %>%
#    cbind(select_if(datalist$Count_Data, is.character), .) %>%
#    as_tibble()
#
#  return(datalist)
#  
#}

