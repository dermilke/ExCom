diversity_datatable <- function(datalist) {
  
  diversity <- tibble(Richness = apply(select_if(datalist$Count_Data, is.numeric) > 0, 2, sum),
                      Shannon = apply(select_if(datalist$Count_Data, is.numeric), 2, vegan::diversity),
                      Evenness = Shannon/log(apply(select_if(datalist$Count_Data, is.numeric), 2, vegan::specnumber)),
                      Simpson = apply(select_if(datalist$Count_Data, is.numeric), 2, vegan::diversity, index = "simpson"),
                      ESN = apply(select_if(datalist$Count_Data, is.numeric), 2, vegan::diversity, index = "invsimpson")) %>%
    bind_cols(datalist$Meta_Data, .)
  
  return(diversity)
  
} 