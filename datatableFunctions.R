### Transformation and Statistics ###

create_datatable <- function(datalist, grpBy = Class, otherThreshold = 0.01) {

  grpBy <- enquo(grpBy)
  
  tmp <- datalist$Count_Data %>%
    group_by(!! grpBy) %>%
    select_if(is.numeric) %>%
    mutate(Richness_Groupmerge = 1) %>%
    summarize_all(sum) %>%
    filter(apply(select(., -1, -Richness_Groupmerge), 2, make_proportion) %>% apply(., 1, mean, na.rm = T) > otherThreshold) %>%
    rbind(., c("Others", c(apply(select_if(datalist$Count_Data, is.numeric), 2, sum), nrow(datalist$Count_Data)) - apply(select(., -1), 2, sum))) %>%
    mutate_at(1, ~replace(., is.na(.), "")) %>%
    mutate_at(2:ncol(.), as.numeric) %>%
    select_if(~sum(is.na(.)) == 0)
  
  result <- tmp %>%
    select(-Richness_Groupmerge) %>%
    reshape2::melt(.) %>%
    rename("Sample_ID" = "variable", "Abundance" = "value") %>%
    left_join(., distinct(datalist$Meta_Data), by = "Sample_ID") %>%
    as_tibble() %>%
    left_join(., select(tmp, 1, Richness_Groupmerge), by = as_label(grpBy)) %>%
    rename("Group" = as_label(grpBy)) 
  
  return(result)
  
}

distance_decay <- function(datalist, dist_measure = "Bray_Curtis") {
  
  if (!(dist_measure %in% c("Bray_Curtis", "Shared_ASVs"))) {
    stop("dist_measure must be either \"Bray_Curtis\" or \"Shared_ASVs\".")
  }
  
  geographic_distance <- function(long1, lat1, long2, lat2) {
    
    dlong <- long1 * pi/180 - long2 * pi/180
    dlat <- lat1 * pi/180 - lat2 * pi/180
    
    R <- 6371
    a <- sin(dlat/2)^2 + cos(lat1 * pi/180) * cos(lat2 * pi/180) * sin(dlong/2)^2
    c <- 2 * atan2(sqrt(a), sqrt(1-a))
    
    distance <- R * c
    
    return(distance)
    
  }
  
  shared_ASVs_distance <- function(Count_Data) {
    
    Count_Data <- Count_Data %>%
      select_if(is.numeric)
    
    shared <- matrix(data = NA, nrow = ncol(Count_Data), ncol = ncol(Count_Data))
    
    for (i in 1:ncol(Count_Data)) {
      for (j in i:ncol(Count_Data)) {
        shared[i,j] <- sum((Count_Data[,j] > 0) & (Count_Data[,i] > 0))
        shared[j,i] <- shared[i,j]
      }
      shared[i,] <- shared[i,] / shared[i,i]
    }
    
    result <- shared 
    
    colnames(result) <- colnames(Count_Data)
    rownames(result) <- colnames(Count_Data)
    
    return(result)
    
  }
  
  DD_datatable <- datalist$Count_Data %>%
    select_if(is.numeric) %>%
    with(., if (dist_measure == "Bray_Curtis") {t(.) %>% vegan::vegdist(.) %>% t(.) %>% as.matrix(.)} else {shared_ASVs_distance(.)}) %>%
    reshape2::melt() %>%
    rename("From" = "Var1", "Sample_ID" = "Var2", !!dist_measure := "value") %>%
    left_join(., datalist$Meta_Data %>%
                select(Sample_ID, Latitude, Longitude, Province, Depth), by = "Sample_ID") %>%
    rename("To" = "Sample_ID", "To_Latitude" = "Latitude", "To_Province" = "Province", 
           "To_Longitude" = "Longitude", "Sample_ID" = "From", "To_Depth" = "Depth") %>%
    left_join(., datalist$Meta_Data %>%
                select(Sample_ID, Latitude, Longitude, Province, Depth), by = "Sample_ID") %>%
    rename("From" = "Sample_ID", "From_Latitude" = "Latitude", "From_Longitude" = "Longitude", 
           "From_Province" = "Province", "From_Depth" = "Depth") %>%
    mutate(Distance = geographic_distance(To_Longitude, To_Latitude, From_Longitude, From_Latitude)) %>%
    as_tibble()
  
  return(DD_datatable)
  
}

diversity_datatable <- function(datalist) {
  
  diversity <- tibble(Richness = apply(select_if(datalist$Count_Data, is.numeric) > 0, 2, sum),
                      Shannon = apply(select_if(datalist$Count_Data, is.numeric), 2, vegan::diversity),
                      Evenness = Shannon/log(apply(select_if(datalist$Count_Data, is.numeric), 2, vegan::specnumber)),
                      Simpson = apply(select_if(datalist$Count_Data, is.numeric), 2, vegan::diversity, index = "simpson")) %>%
    bind_cols(datalist$Meta_Data, .)
  
  return(diversity)
  
} 

running_mean <- function(datatable, interv = 3) {
  # This function is the wrapper-Function for runningMeanCore(), which calculates the running mean for a defined 
  # window of samples (+/- interv). It requires a datatable object. 
  
  running_mean_core <- function(x, interv = 3) {
    # This function calculates the running mean for a predefined window (+/- interv) of samples. It is the core function
    # for the Wrapper-function runningMean()
    
    result <- x
    
    for (i in interv+1:(length(x)-interv)) {
      result[i] <- mean(x[(i-interv):(i+interv)], na.rm = T)
    }
    
    return(result)
    
  }
  
  datatable <- datatable %>%
    arrange(Latitude)
  
  Prop <- datatable$Abundance
  Size_Levels <- unique(datatable$Size_Fraction)
  Depth_Levels <- if ("Depth_Grp" %in% names(datatable)) {unique(datatable$Depth_Grp)} else {unique(datatable$Depth)}
  
  for (i in 1:length(Size_Levels)) {
    for (j in 1:length(Depth_Levels)) {
      Prop[datatable$Size_Fraction == Size_Levels[i] & datatable$Depth_Grp == Depth_Levels[j]] <- running_mean_core(datatable$Proportion[datatable$Size_Fraction == Size_Levels[i] & datatable$Depth_Grp == Depth_Levels[j]], interv)
    }
  }
  
  datatable <- datatable %>%
    mutate(Proportion_Smooth = Prop)
  
  return(data_table)
  
}

