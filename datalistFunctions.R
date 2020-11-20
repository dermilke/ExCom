### Functions for Handling Datalist Objects ###

remove_NAs <- function(Count_Data, replace = 0) {

  Count_Data[is.na(Count_Data)] <- replace

}

make_proportion <- function(x) {
  result <- x/sum(x, na.rm = T)
  return(result)
}

fill_empty <- function(x, entity) {
  
  for (i in (3:ncol(select_if(x, is.character)))) {
    
    x[is.na(x[,i]),i] <- x[is.na(x[,i]),i-1] %>%
                           pull(.) %>%
                           paste(entity, ., sep = " ")
    
    x[,i] <- gsub(paste(entity, " ", entity, " ", sep = ""), paste(entity, " ", sep = ""), pull(x[,i]))
  
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
        dplyr::rename(., 'OTU_ID' = '#OTUID') %>%
        with(., if (fillEmpty) {fill_empty(., "unknown")} else {.})
    )
    
    taxonomy_chloroplast <- suppressWarnings(
      suppressMessages(read_tsv(file)) %>%
        slice(grep("Kingdom.Eukaryota", .$taxonomy, value = F)) %>%
        separate(., 2, c("Kingdom","Supergroup","Phylum","Class","Subclass","Order","Suborder", "Family", "Genus", "Species"), sep = ";") %>%
        mutate_if(is.character, 
                  str_replace_all, pattern = ".*\\.", replacement = "") %>%
        dplyr::rename(., 'OTU_ID' = '#OTUID') %>%
        with(., if (fillEmpty) {fill_empty(., "unknown")} else {.})
    )
    
    taxonomy <- list(Prok = taxonomy_prok, Chloroplast = taxonomy_chloroplast)
    
  } else if (kingdom == "Euk") {
    
    if (DB == "SILVA132") {
      
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
          separate(., 2, c("Kingdom","Supergroup","Phylum","Class","Subclass","Order","Suborder", "Family", "Genus", "Species"), sep = ";") %>%
          mutate_if(is.character, 
                    str_replace_all, pattern = "[a-z]*_", replacement = "") %>%
          dplyr::rename(., 'OTU_ID' = '#OTUID') %>%
          with(., if (fillEmpty) {fill_empty(., "unknown")} else {.})
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
      dplyr::rename(., 'OTU_ID' = '#OTU ID') %>%
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
      dplyr::rename(., 'OTU_ID' = '#OTU ID') %>%
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
      dplyr::rename(., 'OTU_ID' = '#OTU ID') %>%
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
  
  Count_Data <- full_join(mutate_if(datalist_1$Count_Data, is.logical, as.character), 
                          mutate_if(datalist_2$Count_Data, is.logical, as.character), 
                          by = names(select_if(datalist_1$Count_Data, function(x) (is.character(x) | is.logical(x))))) %>%
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
    dplyr::slice(match(names(select_if(count_rared, is.numeric)), datalist$Meta_Data$Sample_ID))
  
  datalist$Count_Data <- count_rared
  datalist$Meta_Data <- meta_subset
  
  return(datalist)
  
}

summarize_by_taxa <- function(datalist, tax_lvl = 7) {
  
  tax_lvl <- enquo(tax_lvl)
  
  if (is.name(rlang::quo_get_expr(tax_lvl))) {
  
    datalist$Count_Data <- datalist$Count_Data %>%
      group_by(!! tax_lvl) %>%
      summarize_if(is.numeric, sum) %>%
      ungroup()

  } else {
    
    datalist$Count_Data <- datalist$Count_Data %>%
      group_by_at(2:(tax_lvl+1)) %>%
      summarize_if(is.numeric, sum) %>%
      ungroup()
    
  }
      
  return(datalist)
  
}

make_depth_groups <- function(datalist, depth_vek) {
  
  datalist$Meta_Data <- datalist$Meta_Data %>%
    mutate(Depth_Grp = depth_vek[findInterval(datalist$Meta_Data$Depth, depth_vek)])
  
  return(datalist)
  
}

filter_station_datalist <- function(.datalist, ...,  removeEmpty = T) {
  
  exprFilter <- enexprs(...)
  
  .datalist$Meta_Data <- .datalist$Meta_Data %>%
    filter(!!! exprFilter)
  
  .datalist$Count_Data <- .datalist$Count_Data %>%
    select_if(is.numeric) %>%
    select_if(names(.) %in% .datalist$Meta_Data$Sample_ID) %>%
    bind_cols(select_if(.datalist$Count_Data, is.character), .) %>%
    with(., if(removeEmpty) filter(., {select_if(., is.numeric) %>% rowSums(.)} > 0) else .)

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

select_meta_datalist <- function(.datalist, ...) {
  
  exprFunc <- enquos(...)
  
  .datalist$Meta_Data <- .datalist$Meta_Data %>%
    select(!!!exprFunc)
  
  return(.datalist)
  
}

mutate_meta_datalist <- function(.datalist, ...) {
  
  exprFunc <- enquos(...)
  
  .datalist$Meta_Data <- .datalist$Meta_Data %>%
    mutate(!!!exprFunc)
  
  return(.datalist)
  
}

NMDS_ordination_datalist <- function(datalist) {
  
  nmds <- datalist$Count_Data %>%
    select_if(is.numeric) %>%
    t(.) %>%
    vegan::vegdist(.) %>%
    as.matrix(.) %>%
    vegan::metaMDS(try = 30) 
  
  table <- nmds%>%
    .$points %>%
    as_tibble(rownames = "Sample_ID") %>%
    full_join(., datalist$Meta_Data, by = "Sample_ID")
  
  nmds_result <- list(table = table, stress = nmds$stress, nmds_obj = nmds)
  
  return(nmds_result)
  
}

sizefraction_communities <- function(datalist, unique_taxa = F) {
  
  datalist_022 <- datalist %>%
    filter_station_datalist(Size_Fraction == "0.22") %>%
    with(., if (unique_taxa) 
      filter_taxa_datalist(!(OTU_ID %in% filter_station_datalist(datalist, Size_Fraction == "3")$Count_Data$OTU_ID | 
                             OTU_ID %in% filter_station_datalist(datalist, Size_Fraction == "8")$Count_Data$OTU_ID))
            else .)
  
  datalist_3 <- datalist %>%
    filter_station_datalist(Size_Fraction == "3") %>%
    with(., if (unique_taxa) 
      filter_taxa_datalist(!(OTU_ID %in% filter_station_datalist(datalist, Size_Fraction == "0.22")$Count_Data$OTU_ID | 
                             OTU_ID %in% filter_station_datalist(datalist, Size_Fraction == "8")$Count_Data$OTU_ID))
      else .)
  
  datalist_8 <- datalist %>%
    filter_station_datalist(Size_Fraction == "8") %>%
    with(., if (unique_taxa) 
      filter_taxa_datalist(!(OTU_ID %in% filter_station_datalist(datalist, Size_Fraction == "0.22")$Count_Data$OTU_ID | 
                             OTU_ID %in% filter_station_datalist(datalist, Size_Fraction == "3")$Count_Data$OTU_ID))
      else .)
  
  return(list(Fraction_022 = datalist_022,
              Fraction_3 = datalist_3,
              Fraction_8 = datalist_8))
  
}

correct_ambiguous <- function(datalist, fromTaxLvl = 8) {
  
  replacer <- function(Count_Data, taxLvl, replaceLvl, pattern) {
    Count_Data[grep(pattern, x = as_vector(Count_Data[, taxLvl])), taxLvl] <- paste0("Unknown ", as_vector(Count_Data[grep(pattern, x = as_vector(Count_Data[, taxLvl])), replaceLvl]))
    return(Count_Data)
  }
  
  for (taxLvl in fromTaxLvl:1) {
  
    for (i in (taxLvl-1):1) {
      
      datalist$Count_Data <- replacer(datalist$Count_Data, taxLvl, replaceLvl = i, pattern = 'bacterium$') %>%
        replacer(., taxLvl, replaceLvl = i, pattern = "uncultured") %>%
        replacer(., taxLvl, replaceLvl = i, pattern = "metagenome") %>%
        replacer(., taxLvl, replaceLvl = i, pattern = "unidentified") %>%
        replacer(., taxLvl, replaceLvl = i, pattern = "Ambiguous_taxa") %>%
        replacer(., taxLvl, replaceLvl = i, pattern = "unknown") %>%
        replacer(., taxLvl, replaceLvl = i, pattern = "unidentified marine bacterioplankton") %>%
        replacer(., taxLvl, replaceLvl = i, pattern = "^Unknown.*bacterium$") %>%
        replacer(., taxLvl, replaceLvl = i, pattern = "^Uncultured.*bacterium$") 
      
    }
  
  }  
    
  return(datalist)
  
}

summarize_by_param <- function(datalist, ...) {
  
  param <- enquos(...)
  
  datalist <- datalist %>%
    mutate_count_datalist(make_proportion)
  
  tmp_final <- select_if(datalist$Count_Data, is.character)
    
  tmp_obj <- unique(select(datalist$Meta_Data, !!!param)) %>%
    mutate_if(is.ordered, as.character)
      
  for (j in 1:nrow(tmp_obj)) {
      
    tmp <- datalist 
      
    for (i in 1:ncol(tmp_obj)) {
        
          tmp <- tmp %>%
            filter_station_datalist(!!param[[i]] == !!deframe(tmp_obj[j,i]), removeEmpty = F)
          
    }
    
    tmp_final <- tmp_final %>%
      mutate(!! paste(tmp_obj[j,], collapse = "_") := rowMeans(select_if(tmp$Count_Data, is.numeric)))
  }
  
  Meta_Data_new <- datalist$Meta_Data %>%
    group_by(!!!param) %>%
    summarize_if(is.numeric, mean) %>%
    ungroup() %>%
    mutate(Sample_ID = deframe(unite(select(., !!! param), "")))
  
  datalist_return <- list(Count_Data = tmp_final, Meta_Data = Meta_Data_new) %>%
    mutate_count_datalist(make_proportion)
  
  return(datalist_return)
  
}

plot_VennDiagram <- function(datalist, param, file) {
  
  param <- enquo(param)
  
  VD_list <- list()
  
  for (i in 1:nrow(unique(select(datalist$Meta_Data, !! param)))) {
    obj_tmp <- deframe(unique(select(datalist$Meta_Data, !! param))[i,])
    VD_list[[i]] <- datalist$Count_Data$OTU_ID[rowSums(select_if(select_if(datalist$Count_Data, is.numeric), select(datalist$Meta_Data, !! param) == obj_tmp)) > 0]  
  }
  
  venn.diagram(x = VD_list, category.names = deframe(unique(select(datalist$Meta_Data, !! param))),
               filename = file,
               output = F, imagetype="png", height = 1000, width = 1000, resolution = 180, compression = "lzw",
               lwd = 2, lty = 'blank', fill = pal_aaas(palette = c("default"))(i), cex = .6, fontface = "bold", fontfamily = "sans",
               cat.cex = 0.6, cat.fontface = "bold",  
               cat.fontfamily = "sans", rotation.degree = 1,
               main = "Venn Diagram", main.fontfamily = "sans")
  
}


plot_heatmap <- function(datalist, taxa_filter = NULL, envir_filter = NULL, cluster_cols = F, row_label = NULL,
                         scale_method = make_proportion, col_label = NULL, title = NULL, cluster_order = NULL,
                         gaps_col = NULL) {
  
  taxa_filter <- enquo(taxa_filter) 
  envir_filter <- enquo(envir_filter)
  col_label <- enquo(col_label)
  row_label <- enquo(row_label)
  
  datalist_filtered <- datalist %>%
    with(., if (!rlang::quo_is_null(taxa_filter)) filter_taxa_datalist(., !!taxa_filter) else .) %>%
    with(., if(!rlang::quo_is_null(envir_filter)) filter_station_datalist(., !!envir_filter) else .)
  
  Count_Matrix <- datalist_filtered %>%
    .$Count_Data %>%
    select_if(is.numeric) %>%
    select(datalist_filtered$Meta_Data %>% 
             with(., 
               if (!rlang::quo_is_null(col_label)) {
                 if ("Size_Fraction" %in% names(datalist$Meta_Data)) {
                   arrange(., Size_Fraction, !!col_label)
                 } else {
                   arrange(., !!col_label)
                 }
               } else {
                 if ("Size_Fraction" %in% names(datalist$Meta_Data)) {
                   arrange(., Size_Fraction, Latitude)
                 } else {
                   arrange(., Latitude)
                 }
               }) %>% .$Sample_ID) %>%
    as.matrix() %>%
    apply(., 1, function(x) scale_method(x)) %>% 
    t()
  
  rownames(Count_Matrix) <- datalist_filtered$Count_Data %>%
    select(OTU_ID) %>%
    deframe()
  
  row_label <- if(!rlang::quo_is_null(row_label)) {
    datalist_filtered$Count_Data %>%
      select(!!row_label) %>%
      deframe()
  } else NULL
  
  col_label <- if (!rlang::quo_is_null(col_label)) {
    
    if ("Size_Fraction" %in% names(datalist$Meta_Data)) {
      datalist_filtered$Meta_Data %>%
        arrange(Size_Fraction, !!col_label) %>%
        select(!!col_label) %>%
        mutate(!!col_label := round(!!col_label, digits = 0)) %>%
        deframe() %>%
        as.character()
    } else {
      datalist_filtered$Meta_Data %>%
        arrange(!!col_label) %>%
        mutate(!!col_label := if (is.numeric(!!col_label)) round(!!col_label, digits = 0) else !!col_label) %>%
        select(!!col_label) %>%
        deframe() %>%
        as.character()
    }
    
  } else NULL
  
  clust_num <- factoextra::fviz_nbclust(Count_Matrix, cluster::pam, method = "silhouette", k.max = min(nrow(Count_Matrix)-1, 8))$data$y %>%
    which.max()
  
  cluster_obj <- Count_Matrix %>%
    vegan::vegdist(method = "euclidian") %>%
    hclust(., method = "ward.D2") %>%
    with(., if ("Latitude" %in% names(datalist_filtered$Meta_Data)) {
      reorder(., wts = abs(datalist_filtered$Meta_Data$Latitude), agglo.FUN = "uwmean")
      } else .)
  
  
  #%>%
  #  with(., reorder(., wts = apply(Count_Matrix, 1, function(x) abs(datalist_filtered$Meta_Data$Latitude[which.max(x)])*sd(x)),
  #                  agglo.FUN = "mean"))
  
  row_annotation <- data.frame(Cluster = paste("Cluster", cutree(cluster_obj, k = clust_num)), row.names = rownames(Count_Matrix))
  row_annotation_color <- list(Cluster = structure(ggsci::pal_rickandmorty(palette = c("schwifty"))(clust_num), names = unique(row_annotation$Cluster)))
  
  if (is.null(title)) {
    title <- paste(str_split_fixed(as_label(taxa_filter), pattern = "==", 2)[1] %>% str_replace(., " ", ""), ": ",
                   str_split_fixed(as_label(taxa_filter), pattern = "==", 2)[2] %>% str_replace_all(., "\\\"", "") %>% str_replace(., " ", ""), sep = "")
  }
  
  pheatmap(Count_Matrix, cluster_cols = cluster_cols, cluster_rows = cluster_obj, cutree_rows = clust_num, #clustering_method = "ward.D2", 
           main = title, border_color = "grey60",
           legend_breaks = c(min(Count_Matrix), max(Count_Matrix)),legend_labels = c("Min  ", "Max  "),
            angle_col = "90", annotation_legend = F, labels_col = col_label,
           labels_row = row_label, show_rownames = ifelse(!is.null(row_label),T,F),
           annotation_row = row_annotation, fontsize_row = 6, fontsize_col = 7,
           annotation_colors = row_annotation_color, fontsize = 7,
           gaps_col = gaps_col, 
           clustering_callback = callback)
  
}

phylogenetic_distance <- function(datalist, fasta_location = NULL) {
  
  my.read.lines2=function(fname) {
    s = file.info( fname )$size 
    buf = readChar( fname, s, useBytes=T)
    strsplit( buf,"\n",fixed=T,useBytes=T)[[1]]
  }
  
  library(msa)
  
  subset <- datalist$Count_Data$OTU_ID
  Sequence <- NULL
  
  for (i in 1:length(subset)) {
    tmp_index <- grep(subset[i], my.read.lines2(fasta_location))
    Sequence <- c(Sequence, readLines(fasta_location)[tmp_index + 1])
  }
  
  tibble(OTU_ID = subset, Sequence = Sequence) %>%
    mutate(OTU_ID = paste(">", OTU_ID, sep = "")) %>%
    with(., do.call(rbind, lapply(seq(nrow(.)), function(i) t(.[i, ])))) %>%
    write.table(., "~/PhD/tmp.fasta", row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  DNA_set <- Biostrings::readDNAStringSet("~/PhD/tmp.fasta", format = "fasta")
  DNA_align <- msa(DNA_set, method = "ClustalW") %>%
    msaConvert(type="seqinr::alignment")
  
  Phylo_Distances <- seqinr::dist.alignment(DNA_align)
  
  Count_Matrix <- datalist$Count_Data %>%
    dplyr::slice(match(attr(Phylo_Distances, "Labels"), datalist$Count_Data$OTU_ID)) %>%
    select_if(is.numeric) %>%
    data.frame(row.names = attr(Phylo_Distances, "Labels"), check.names = F) %>%
    as.matrix()
  
  Taxo_Distances <- Count_Matrix %>%
    apply(., 1, function(x) (x-mean(x))/sd(x)) %>%
    t() %>%
    vegan::vegdist(method = "euclidean")
  
  Tax_Matrix <-  datalist$Count_Data %>%
    dplyr::slice(match(attr(Phylo_Distances, "Labels"), datalist$Count_Data$OTU_ID)) %>%
    dplyr::rename("Domain" = "Kingdom") %>%
    select_if(is.character) %>%
    select(-OTU_ID) %>%
    data.frame(row.names = attr(Phylo_Distances, "Labels")) %>%
    as.matrix()
  
  Meta_Matrix <- datalist$Meta_Data %>%
    data.frame(row.names = .$Sample_ID)
  
  Unifrac_Distances <- phyloseq(otu_table(Count_Matrix, taxa_are_rows = T), 
                                tax_table(Tax_Matrix), 
                                sample_data(Meta_Matrix), 
                                refseq(DNA_set), 
                                ape::nj(Phylo_Distances)) %>%
    UniFrac(., weighted=F, normalized=TRUE, parallel=TRUE, fast=TRUE)
  
  return(list(Phylogenetic_Distance = Phylo_Distances,
              Unifrac_Distance = Unifrac_Distances))
  
}

unassigned_taxonomy <- function(datalist, tax_grp) {
  
  tax_grp <- enquo(tax_grp)
  
  datalist %>%
    mutate_count_datalist(function(x) ifelse(x > 0, 1, 0))  %>%
    .$Count_Data %>%
    mutate(!!tax_grp := ifelse(is.na(select(., !!tax_grp)), "Unassigned", "Assigned")) %>%
    group_by(!!tax_grp) %>%
    summarize_if(is.numeric, sum) %>%
    t() %>%
    .[-1,] %>%
    as_tibble(rownames = "Sample_ID") %>%
    dplyr::rename("Assigned" = "V1", "Unassigned" = "V2") %>%
    mutate(Assigned = as.numeric(Assigned)) %>%
    mutate(Unassigned = as.numeric(Unassigned)) %>%
    right_join(., datalist$Meta_Data) %>%
    mutate(Unassigned_Rel = Unassigned/(Unassigned + Assigned)) 
  
}
