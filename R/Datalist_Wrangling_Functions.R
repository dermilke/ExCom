#### Datalist Wrangling Functions ####

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
    dplyr::slice(match(names(select_if(Count_Data, is.numeric)), .$Sample_ID))
  
  datalist_return <- list(Count_Data = Count_Data, Meta_Data = Meta_Data)
  
  if (("Sequence" %in% names(datalist_1)) & ("Sequence" %in% names(datalist_2))) {
    datalist_return$Sequence <- bind_rows(datalist_1$Sequence, datalist_2$Sequence) %>%
      distinct() %>%
      slice(match(datalist_return$Count_Data$OTU_ID, Seq_ID))
  }
  
  return(datalist_return)
  
}

singleton_filter <- function(datalist, min_count = 5, min_station = 2) {
  
  Count_Data <- datalist$Count_Data %>%
    filter(select_if(., is.numeric) %>% with(., .>0) %>% rowSums() >= min_station) %>%
    filter(select_if(., is.numeric) %>% rowSums() >= min_count)
  
  datalist$Count_Data <- Count_Data
  
  if ("Sequences" %in% names(datalist)) {
    datalist$Sequences <- dplyr::slice(datalist$Sequences, match(datalist$Count_Data$OTU_ID, Seq_ID))
  }
  
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
  
  if ("Sequences" %in% names(.datalist)) {
    .datalist$Sequences <- dplyr::slice(.datalist$Sequences, match(.datalist$Count_Data$OTU_ID, Seq_ID))
  }
  
  return(.datalist)
  
}

filter_taxa_datalist <- function(.datalist, ...) {
  
  exprFilter <- enquos(...)
  
  .datalist$Count_Data <- filter(.datalist$Count_Data, !!!(exprFilter))
  
  if ("Sequence" %in% names(.datalist)) {
    .datalist$Sequence <- dplyr::slice(.datalist$Sequence, match(.datalist$Count_Data$OTU_ID, Seq_ID))
  }
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
    dplyr::slice(!!! exprSlice)
  
  if ("Sequences" %in% names(.datalist)) {
    .datalist$Sequences <- dplyr::slice(.datalist$Sequences, match(.datalist$Count_Data$OTU_ID, Seq_ID))
  }
  
  return(.datalist)
  
}

select_meta_datalist <- function(.datalist, ...) {
  
  exprFunc <- enquos(...)
  
  .datalist$Meta_Data <- .datalist$Meta_Data %>%
    select(!!!exprFunc)
  
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

summarize_by_taxa <- function(datalist, tax_lvl = Species) {
  
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
  
  if ("Sequences" %in% names(datalist)) {
    datalist <- datalist[names(datalist) != "Sequences"]
  }
  
  return(datalist)
  
}

arrange_datalist <- function(datalist, ...) {
  
  param <- enquos(...)
  
  datalist$Meta_Data <- datalist$Meta_Data %>%
    dplyr::arrange(!!!param)
  
  datalist$Count_Data <- bind_cols(select_if(datalist$Count_Data, is.character),
                                   select(datalist$Count_Data, datalist$Meta_Data$Sample_ID))
  
  return(datalist)
  
}

filter_abundance <- function(datalist) {
  
  counts <- datalist$Count_Data
  prop <- datalist$Count_Data %>%
    mutate_if(is.numeric,
              function(x) x/sum(x))
  
  counts_filtered <- counts %>%
    filter(((rowSums(select_if(counts, is.numeric))/sum(rowSums(select_if(counts, is.numeric)))) > 0.00001) &
             (apply(select_if(prop, is.numeric),1,max) > 0.01) |
             ((apply(select_if(prop, is.numeric) > 0.001, 1, sum) / ncol(select_if(prop, is.numeric))) > 0.02) |
             ((apply(select_if(prop, is.numeric) > 0, 1, sum) / ncol(select_if(prop, is.numeric))) > 0.05))
  
  datalist$Count_Data <- counts_filtered
  
  if ("Sequences" %in% names(.datalist)) {
    datalist$Sequences <- dplyr::slice(datalist$Sequences, match(datalist$Count_Data$OTU_ID, Seq_ID))
  }
  
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
  
  if ("Sequences" %in% names(.datalist)) {
    datalist$Sequences <- dplyr::slice(datalist$Sequences, match(datalist$Count_Data$OTU_ID, Seq_ID))
  }
  
  return(datalist)
  
}

summarize_by_param <- function(datalist, ...) {
  
  param <- enquos(...)
  
  datalist <- datalist %>%
    mutate_count_datalist(function(x) x/sum(x))
  
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
      mutate(!! paste(tmp_obj[j,], collapse = "_") := 
               rowSums(select_if(tmp$Count_Data, is.numeric)))
  } 
  
  Meta_Data_new <- datalist$Meta_Data %>%
    group_by(!!!param) %>%
    summarize_if(is.numeric, mean) %>%
    ungroup() %>%
    mutate(Sample_ID = deframe(unite(select(., !!! param), "_"))) %>%
    dplyr::slice(match(names(select_if(tmp_final, is.numeric)),Sample_ID))
  
  
  datalist_return <- list(Count_Data = tmp_final, Meta_Data = Meta_Data_new) %>%
    mutate_count_datalist(function(x) x/sum(x))
  
  if ("Sequences" %in% names(.datalist)) {
    datalist$Sequences <- dplyr::slice(datalist$Sequences, match(datalist$Count_Data$OTU_ID, Seq_ID))
  }
  
  return(datalist_return)
  
}

create_datatable <- function(datalist, grpBy = Family, upper_grp = Class, otherThreshold = 0.01, selectGrps = NULL, 
                             addColorScheme = F, addOthers = T) {
  
  make_proportion <- function(x) {
    result <- x/sum(x, na.rm = T)
    return(ifelse(is.finite(result), result, 0))
  }
  
  grpBy <- enquo(grpBy)
  upper_grp <- enquo(upper_grp)
  
  if (addColorScheme) {
    
    datalist$Count_Data <- mutate(datalist$Count_Data, !!quo_name(grpBy) := paste(!!upper_grp, !!grpBy, sep = " - "))
    
  }
  
  if (is.null(selectGrps)) {
    
    if (addOthers) {
      
      tmp_Upper <- datalist$Count_Data %>%
        group_by(!!upper_grp) %>%
        select_if(is.numeric) %>%
        summarize_all(sum, na.rm = T)
      
      tmp_Lower <- datalist$Count_Data %>%
        group_by(!!grpBy) %>%
        select_if(is.numeric) %>%
        summarize_all(sum, na.rm = T) %>%
        filter(apply(select(., -1), 2, make_proportion) %>% apply(., 1, mean, na.rm = T) > otherThreshold) 
      
      tmp_both <- tmp_Lower %>%
        mutate(!!quo_name(upper_grp) := pull(left_join(tmp_Lower, distinct(select(datalist$Count_Data, !!upper_grp, !!grpBy))), !!upper_grp), .before = 1) 
      
      tmp_Upper_red <- tmp_both %>%
        group_by(!!upper_grp) %>%
        select_if(is.numeric) %>%
        summarize_all(sum, na.rm = T)
      
      tmp <- (select_if(filter(tmp_Upper, !!upper_grp %in% pull(tmp_Upper_red, !!upper_grp)), is.numeric) - 
                select_if(tmp_Upper_red, is.numeric)) %>%
        as_tibble() %>%
        mutate(!!quo_name(grpBy) := paste0(pull(tmp_Upper_red, !!upper_grp), " - Others"), .before = 1) %>%
        rbind(., tmp_Lower) %>%
        mutate_at(1, ~replace(., is.na(.), "")) %>%
        mutate_at(2:ncol(.), as.numeric) %>%
        select_if(~sum(is.na(.)) == 0) %>%
        add_case({colSums(select_if(datalist$Count_Data, is.numeric)) -
            summarise(., across(where(is.numeric), sum))} %>%
                   mutate(!!quo_name(grpBy) := paste0("Other ", as_label(upper_grp)))) %>%
        filter(rowSums(select_if(., is.numeric)) > 10^(-16))
      
    } else {
      
      tmp <- datalist$Count_Data %>%
        group_by(!!grpBy) %>%
        select_if(is.numeric) %>%
        summarize_all(sum, na.rm = T) %>%
        filter(apply(select(., -1), 2, make_proportion) %>% apply(., 1, mean, na.rm = T) > otherThreshold) %>%
        rbind(., c("Others", apply(select_if(datalist$Count_Data, is.numeric), 2, sum) - apply(select(., -1), 2, sum))) %>%
        mutate_at(1, ~replace(., is.na(.), "")) %>%
        mutate_at(2:ncol(.), as.numeric) %>%
        select_if(~sum(is.na(.)) == 0)
      
    }
    
    
  } else {
    
    tmp <- datalist$Count_Data %>%
      group_by(!! grpBy) %>%
      select_if(is.numeric) %>%
      summarize_all(sum, na.rm = T) %>%
      filter(!! grpBy %in% selectGrps) %>%
      rbind(., c(paste0("Others (<", otherThreshold * 100,"%)"), apply(select_if(datalist$Count_Data, is.numeric), 2, sum) - apply(select(., -1), 2, sum))) %>%
      mutate_at(1, ~replace(., is.na(.), "")) %>%
      mutate_at(2:ncol(.), as.numeric) %>%
      select_if(~sum(is.na(.)) == 0)
    
  }
  
  if (addOthers) {
    
    group_label <- NA
    
    for (i in 1:length(pull(tmp_Upper_red, !!upper_grp))) {
      
      group_label <- c(group_label, pull(filter(tmp_both, !!upper_grp == pull(tmp_Upper_red, !!upper_grp)[i]), !!grpBy), 
                       paste0(pull(tmp_Upper_red, !!upper_grp)[i], " - Others"))
      
    }
    
    group_label <- c(group_label[-1], paste0("Other ", as_label(upper_grp)))
    
    datatable <- tmp %>%
      reshape2::melt(.) %>%
      dplyr::rename("Sample_ID" = "variable", "Abundance" = "value") %>%
      left_join(., distinct(datalist$Meta_Data), by = "Sample_ID") %>%
      as_tibble() %>%
      dplyr::rename("Group" = as_label(grpBy)) %>%
      mutate(Group = ordered(Group, levels = group_label)) %>%
      arrange(Group)
    
  } else {
    
    tmp_richness <- datalist$Count_Data %>%
      mutate_if(is.numeric, function(x) as.numeric(x > 0)) %>%
      group_by(!!grpBy) %>%
      select_if(is.numeric) %>%
      summarize_all(sum, na.rm = T) %>%
      filter(!!grpBy %in% deframe(select(tmp, as_label(grpBy)))) %>%
      rbind(., c("Others", apply(select_if(datalist$Count_Data, is.numeric) > 0, 2, sum) - apply(select(., -1), 2, sum))) %>%
      mutate_at(1, ~replace(., is.na(.), "")) %>%
      mutate_at(2:ncol(.), as.numeric) %>%
      select_if(~sum(is.na(.)) == 0) 
    
    datatable <- tmp %>%
      reshape2::melt(.) %>%
      dplyr::rename("Sample_ID" = "variable", "Abundance" = "value") %>%
      left_join(., {tmp_richness %>% reshape2::melt(.) %>% dplyr::rename("Sample_ID" = "variable", "Richness" = "value")}, by = c("Sample_ID", as_label(grpBy))) %>%
      left_join(., distinct(datalist$Meta_Data), by = "Sample_ID") %>%
      as_tibble() %>%
      dplyr::rename("Group" = as_label(grpBy)) %>%
      mutate(Group = ordered(Group, levels = c(unique(Group)[-which(Group == "Others")], "Others"))) %>%
      arrange(Group)
    
  }
  
  if (addColorScheme) {
    
    #colorTaxonomy <- read_csv("https://raw.githubusercontent.com/dermilke/ExCom/master/Data/Colors/Taxonomy_Colour.csv")
    colorTaxonomy <- read_csv("~/PhD/SoftwareBuilds/ExCom/Data/Colors/Taxonomy_Colour.csv")
    
    derivative_color = str_split_fixed(datatable$Group, " - ", 2) %>%
      unique() %>%
      .[,1] %>%
      table() %>%
      as_tibble() %>%
      dplyr::rename("Group" = ".") %>%
      mutate(Group = ordered(Group, levels = c(unique(Group[Group != paste0("Other ", as_label(upper_grp))]), 
                                               paste0("Other ", as_label(upper_grp))))) %>%
      arrange(Group) %>%
      mutate(main_color = colorTaxonomy$Colour[match(Group, colorTaxonomy$Group)] %>%
               ifelse(is.na(.), "grey", .)) %>%
      mutate(rowNumber = 1:n()) %>%
      with(., {
        derivative_color <- NULL
        for (i in rowNumber) {
          
          derivative_color <- c(derivative_color,
                                colorspace::lighten(main_color[i], 
                                                    amount = {if (n[i] == 1) 0
                                                      else if (n[i] > 6) (seq_len(n[i])/n[i])*1.6-0.8
                                                      else (seq_len(n[i])/n[i])-0.5}
                                )
          )
          
        }
        derivative_color
      })
    
    result <- list(table = datatable, color = derivative_color)
    
  } else {
    
    result <- datatable
    
  }
  
  return(result)
  
}
