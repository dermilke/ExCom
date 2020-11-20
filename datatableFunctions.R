### Transformation and Statistics ###

create_datatable <- function(datalist, grpBy = Class, otherThreshold = 0.01, selectGrps = NULL, addColorScheme = F) {

  grpBy <- enquo(grpBy)

  if (addColorScheme) {
    
    datalist$Count_Data <- mutate(datalist$Count_Data, !!quo_name(grpBy) := paste(Class, !!grpBy, sep = ";"))
    
  }
  
  if (is.null(selectGrps)) {
    
    tmp <- datalist$Count_Data %>%
      group_by(!!grpBy) %>%
      select_if(is.numeric) %>%
      summarize_all(sum, na.rm = T) %>%
      filter(apply(select(., -1), 2, make_proportion) %>% apply(., 1, mean, na.rm = T) > otherThreshold) %>%
      rbind(., c("Others", apply(select_if(datalist$Count_Data, is.numeric), 2, sum) - apply(select(., -1), 2, sum))) %>%
      mutate_at(1, ~replace(., is.na(.), "")) %>%
      mutate_at(2:ncol(.), as.numeric) %>%
      select_if(~sum(is.na(.)) == 0)
    
  } else {
    
    tmp <- datalist$Count_Data %>%
      group_by(!! grpBy) %>%
      select_if(is.numeric) %>%
      summarize_all(sum, na.rm = T) %>%
      filter(!! grpBy %in% selectGrps) %>%
      rbind(., c("Others", apply(select_if(datalist$Count_Data, is.numeric), 2, sum) - apply(select(., -1), 2, sum))) %>%
      mutate_at(1, ~replace(., is.na(.), "")) %>%
      mutate_at(2:ncol(.), as.numeric) %>%
      select_if(~sum(is.na(.)) == 0)
    
  }
  
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
  
  if (addColorScheme) {
    
    colorTaxonomy <- read_csv("https://raw.githubusercontent.com/dermilke/ExCom/master/Data/Colors/Taxonomy_Colour.csv")
    
    derivative_color = str_split_fixed(datatable$Group, ";", 2) %>%
      unique() %>%
      .[,1] %>%
      table() %>%
      as_tibble() %>%
      rename("Group" = ".") %>%
      mutate(Group = ordered(Group, levels = c(unique(Group)[-which(Group == "Others")], "Others"))) %>%
      arrange(Group) %>%
      mutate(main_color = colorTaxonomy$Colour[match(Group, colorTaxonomy$Group)]) %>%
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
        derivative_color[length(derivative_color)] <- "grey40"
        derivative_color
      })
    
    result <- list(table = datatable, color = derivative_color)
    
  } else {
    
    result <- datatable
    
  }
  
  return(result)
  
}

distance_decay <- function(datalist, dist_measure = "Bray_Curtis") {
  
  if (!("Depth_Grp" %in% names(datalist$Meta_Data))) {
    warning("Depth_Grps not found, using actual Depth instead.")
    datalist <- datalist %>%
      mutate_meta_datalist(Depth_Grp = Depth)
  }
  
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
                select(Sample_ID, Latitude, Longitude, Province, Depth_Grp), by = "Sample_ID") %>%
    rename("To" = "Sample_ID", "To_Latitude" = "Latitude", "To_Province" = "Province", 
           "To_Longitude" = "Longitude", "Sample_ID" = "From", "To_Depth" = "Depth_Grp") %>%
    left_join(., datalist$Meta_Data %>%
                select(Sample_ID, Latitude, Longitude, Province, Depth_Grp), by = "Sample_ID") %>%
    rename("From" = "Sample_ID", "From_Latitude" = "Latitude", "From_Longitude" = "Longitude", 
           "From_Province" = "Province", "From_Depth" = "Depth_Grp") %>%
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

profiling <- function(datatable, groupFilter) {
  
  tmp <- datatable %>%
    filter(Group == groupFilter)
    
  # 1. Size_Fraction:
  # Only best SF if there are significant differences between SFs
  
  SF_tmp <- group_by(tmp, Size_Fraction) %>% 
    summarize(Abundance = sum(Abundance)) 
  
  if (kruskal.test(tmp$Abundance, tmp$Size_Fraction)$p.value < 0.05) {
    best_SF <- SF_tmp %>%
      arrange(desc(Abundance)) %>%
      slice(1) %>%
      .$Size_Fraction
    # Check if best SF is significantly different to other SF
    # and if not: add the other SF that is not different to best SF
    a <- DescTools::DunnTest(tmp$Abundance, as.factor(tmp$Size_Fraction))[[1]][,2] %>%
      tibble::enframe() %>%
      separate(name, c("A", "B"), "-") %>%
      filter(A %in% best_SF | B %in% best_SF) %>%
      filter(value > 0.05)
    
    best_SF <- c(best_SF, a$A[a$A != best_SF], a$B[a$B != best_SF]) %>%
      unique() %>%
      as.character()
    
  } else if (max(SF_tmp$Abundance) / (min(SF_tmp$Abundance)+0.000001) > 10) {
    # Get best SF: If summed Abundance > 10 times smallest: Define best_SF 
    # Then check highest SF_Abundance and look which other SF is at least 33% of Abundance
    # of highest SF Abundance
    best_SF <- SF_tmp %>%
      filter(Abundance > 0) %>%
      filter((Abundance)/max(Abundance) > 0.33) %>%
      select(Size_Fraction) %>%
      as_vector() 
    names(best_SF) <- NULL
    
  } else {
    best_SF <- "All" # or "not 
  }
  
  best_SF <- if (length(best_SF) == 3) "All" else best_SF
  
  # 2. Oceans:
  # Select Ocean if ASV is 90% (or other value) more abundant in one Ocean
  best_Ocean <- if (length(unique(tmp$Cruise)) != 1) {
    if (kruskal.test(tmp$Abundance, tmp$Cruise)$p.value < 0.05) {
    
      group_by(tmp, Cruise) %>%
        summarize(Abundance = sum(Abundance)) %>%
        arrange(desc(Abundance)) %>%
        slice(1) %>%
        .$Cruise 
      
    } else {
      
      best_Ocean <- tmp %>%
        group_by(Cruise) %>%
        summarize(Abundance = sum(Abundance))
      
      if (sum(best_Ocean$Abundance == 0) == 1) {
        best_Ocean <- best_Ocean %>%
          filter(Abundance > 0) %>%
          .$Cruise
      } else "All"
      
    }
  } else best_Ocean = NA
  
  # 3. Depths:
  # Best Depth is its 90% Quantile
  best_Depth <- tmp %>%
    {if (best_SF[1] != "All") filter(., Size_Fraction %in% best_SF) else .} %>%
    {if (best_Ocean != "All" & !is.na(best_Ocean)) filter(., Cruise %in% best_Ocean) else .} %>%
    group_by(Depth) %>%
    summarize(Abundance = sum(Abundance)) %>%
    mutate(Abundance = Abundance / sum(Abundance)) %>%
    with(., quantile(rep(x = Depth, times = round(Abundance*1000)), c(0.1,0.90)))
  
  # 4. Abundance-Maxima:
  
  # First: Calculate a loess-model-fit and find the maxima of the fitted curve. Then take only those maxima
  # that are at least 75% of the highest maximum found. Because it is the fitted curve, outliers wont 
  # affect maximum height very much, thats why i think 15% should be sufficient.
  Model_Fit <- tmp %>%
    {if (best_SF[1] != "All") filter(., Size_Fraction %in% best_SF) else .} %>%
    {if (best_Ocean != "All" & !is.na(best_Ocean)) filter(., Cruise %in% best_Ocean) else .} %>%
    filter(Depth >= best_Depth[1] & Depth <= best_Depth[2])
  
  if (nrow(Model_Fit) >= 3) {
    
    Model_Fit <- Model_Fit %>%
      mutate(Abundance_Model = loess(Abundance ~ Latitude, data = ., degree = 2, span = 0.4)$fitted) %>%
      arrange(Latitude) %>%
      select(Latitude, Abundance_Model) %>%
      distinct()
    
    # A minimum can not be the neighbour of a maximum
    # And a maximum must be at least 20% of highest maximum abundance
    # And borders of interval can be maximum if they have higher abundance than 80% of lowest Maximum Abundance Peak
    # And all Abundances are the fitted model-Abundances
    best_Latitude <- Model_Fit %>%
      slice(which(diff(sign(diff(.$Abundance_Model)))==-2 &
                  c(diff(sign(diff(.$Abundance_Model)))[-1], -2) != 2 &
                  c(-2, diff(sign(diff(.$Abundance_Model)))[-(nrow(.)-2)]) != 2)+1) %>%
      filter(Abundance_Model >= 0.2*max(Abundance_Model)) %>%
      rbind(., filter(Model_Fit[c(1, nrow(Model_Fit)),], Abundance_Model > min(.$Abundance_Model)*0.8))
    
  } else {
    
    best_Latitude <- Model_Fit %>%
      filter(Abundance == max(Abundance)) %>%
      select(Latitude, Abundance) %>%
      rename("Abundance_Model" = "Abundance")
    
  }
  
  # 5. Hemisphere:
  # Simple approach: if there is a 10 times difference of abundances between hemispheres, assign a hemisphere.
  # If not: just "All"
  best_Hemisphere <- tmp %>%
    {if (best_SF[1] != "All") filter(., Size_Fraction %in% best_SF) else .} %>%
    {if (best_Ocean != "All" & !is.na(best_Ocean)) filter(., Cruise %in% best_Ocean) else .} %>%
    filter(Depth >= best_Depth[1] & Depth <= best_Depth[2]) %>%
    mutate(Hemisphere = ifelse(Latitude <= 0, "Southern", "Northern")) %>%
    group_by(Hemisphere) %>%
    summarize(Abundance = sum(Abundance)) %>%
    arrange(Hemisphere) %>%
    with(., ifelse(Abundance[1]/Abundance[2] >= 10, "Northern", 
                   ifelse(Abundance[1]/Abundance[2] <= 0.1, "Southern", "All")))
  
  best_Temp <- tmp %>%
    {if (best_SF[1] != "All") filter(., Size_Fraction %in% best_SF) else .} %>%
    {if (best_Ocean != "All" & !is.na(best_Ocean)) filter(., Cruise %in% best_Ocean) else .} %>%
    filter(Depth >= best_Depth[1] & Depth <= best_Depth[2]) %>%
    with(., rep(x = Pot_Temperature, times = Abundance*1000)) %>%
    quantile(.,  c(0.1,0.90))
  
  profile_data <- list(Size_Fraction = best_SF,
                       Cruise = best_Ocean,
                       Depth_Interval = best_Depth,
                       Hemisphere = best_Hemisphere,
                       Temp_Interval = best_Temp,
                       Latitude_Max = best_Latitude)
  
  profile_plotter(tmp, profile_data, groupFilter)
  
  return(profile_data)

}

profile_plotter <- function(datatable, profile_data, groupFilter) {
  
  layout(matrix(c(1,2,3,4,4,5,6,6,7,7), nrow = 2, byrow = T))
  
  boxplot(log(datatable$Abundance+1)~datatable$Size_Fraction, xlab = "", ylab = "log Abundance", 
          main = "Size Fractions")
  points(x = ifelse(profile_data$Size_Fraction == 0.22, 1, 
                    ifelse(profile_data$Size_Fraction == 3, 2, 
                           ifelse(profile_data$Size_Fraction == 8, 3, 0))), 
         y = rep(0, length(profile_data$Size_Fraction)), bg = "red", pch = 21, cex = 2.8)
  
  boxplot(log(datatable$Abundance+1)~datatable$Cruise, xlab = "", ylab = "log Abundance", main = "Cruise")
  points(x = ifelse(profile_data$Cruise == "ANT28-4/5", 1, 
                    ifelse(profile_data$Cruise == "SO248-254", 2, 0)), y = 0, bg = "red", 
         pch = 21, cex = 2.8)
  
  datatable %>%
    {if (profile_data$Size_Fraction[1] != "All") filter(., Size_Fraction %in% profile_data$Size_Fraction) else .} %>%
    {if (profile_data$Cruise != "All" & !is.na(profile_data$Cruise)) filter(., Cruise %in% profile_data$Cruise) else .} %>%
    with(., boxplot(log(Abundance+1)~ordered(Depth_Grp, levels = c(20, 40, 60, 100, 200, 300)), xlab = "", ylab = "log Abundance", main = "Depth"))
  
  datatable %>%
    {if (profile_data$Size_Fraction[1] != "All") filter(., Size_Fraction %in% profile_data$Size_Fraction) else .} %>%
    {if (profile_data$Cruise != "All" & !is.na(profile_data$Cruise)) filter(., Cruise %in% profile_data$Cruise) else .} %>%
    with(., lines(x = findInterval(profile_data$Depth_Interval, c(20, 40, 60, 100, 200, 300)) + c(-0.4, 0.4), 
                  y = c(min(log(Abundance+1)),min(log(Abundance+1))), col = "red", lwd = 5))
  
  datatable %>%
    {if (profile_data$Size_Fraction[1] != "All") filter(., Size_Fraction %in% profile_data$Size_Fraction) else .} %>%
    {if (profile_data$Cruise != "All" & !is.na(profile_data$Cruise)) filter(., Cruise %in% profile_data$Cruise) else .} %>%
    filter(Depth_Grp >= profile_data$Depth_Interval[1] & Depth_Grp <= profile_data$Depth_Interval[2]) %>%
    with(., plot(Pot_Temperature, Abundance))
  
  lines(profile_data$Temp_Interval, c(0,0), col = "red", lwd = 5)
  
  boxplot(log(datatable$Abundance+1)~as.factor(ifelse(datatable$Latitude >= 0, "Northern", "Southern")), ylab = "log Abundance",
          xlab = "", main = "Hemisphere")
  points(x = ifelse(profile_data$Hemisphere == "Northern", 1, 
                    ifelse(profile_data$Hemisphere == "Southern", 2, 0)), y = 0,
         bg = "red", pch = 21, cex = 2.8)
  
  datatable %>%
    {if (profile_data$Size_Fraction[1] != "All") filter(., Size_Fraction %in% profile_data$Size_Fraction) else .} %>%
    {if (profile_data$Cruise != "All" & !is.na(profile_data$Cruise)) filter(., Cruise %in% profile_data$Cruise) else .} %>%
    filter(Depth_Grp >= profile_data$Depth_Interval[1] & Depth_Grp <= profile_data$Depth_Interval[2]) %>%
    with(., plot(Latitude, Abundance, xlab = "", ylab = "Abundance", main = "Latitude", xlim = range(datatable$Latitude)))
  
  datatable %>%
    {if (profile_data$Size_Fraction[1] != "All") filter(., Size_Fraction %in% profile_data$Size_Fraction) else .} %>%
    {if (profile_data$Cruise != "All" & !is.na(profile_data$Cruise)) filter(., Cruise %in% profile_data$Cruise) else .} %>%
    filter(Depth_Grp >= profile_data$Depth_Interval[1] & Depth_Grp <= profile_data$Depth_Interval[2]) %>%
    mutate(Abundance_Model = loess(Abundance ~ Latitude, data = ., degree = 2, span = 0.5)$fitted) %>%
    with(., points(Latitude, Abundance_Model, col = "red"))
  
  points(x = profile_data$Latitude_Max$Latitude, y = rep(0, length(profile_data$Latitude_Max$Latitude)), pch = 21, bg = "red", cex = 2.8)
  
  mtext(groupFilter, 3, 3, cex = 0.6)
  
  datatable %>%
    {if (profile_data$Size_Fraction[1] != "All") filter(., Size_Fraction %in% profile_data$Size_Fraction) else .} %>%
    {if (profile_data$Cruise != "All" & !is.na(profile_data$Cruise)) filter(., Cruise %in% profile_data$Cruise) else .} %>%
    filter(Depth_Grp >= profile_data$Depth_Interval[1] & Depth_Grp <= profile_data$Depth_Interval[2]) %>%
    with(., plot(Latitude, Richness, xlab = "Latitude", ylab = "Richness", xlim = range(datatable$Latitude)))
  
  datatable %>%
    {if (profile_data$Size_Fraction[1] != "All") filter(., Size_Fraction %in% profile_data$Size_Fraction) else .} %>%
    {if (profile_data$Cruise != "All" & !is.na(profile_data$Cruise)) filter(., Cruise %in% profile_data$Cruise) else .} %>%
    filter(Depth_Grp >= profile_data$Depth_Interval[1] & Depth_Grp <= profile_data$Depth_Interval[2]) %>%
    mutate(Richness_Model = loess(Richness ~ Latitude, data = ., degree = 2, span = 0.5)$fitted) %>%
    with(., points(Latitude, Richness_Model, col = "red"))
  
}

