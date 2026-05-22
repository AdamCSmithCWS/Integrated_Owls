# function to identify which rows of sf polygon layer are nearest and less than
# dist metres from route start lecations.
which_min_LT <- function(outside,
                         strata_map,
                         distance_to_strata = 2000){
  
  x <- sf::st_is_within_distance(outside,
                                 strata_map,
                                 dist = distance_to_strata,
                                 sparse = TRUE)
  
  miss_dist <- sf::st_distance(outside,strata_map)
  
  distance_to_strata <- units::set_units(distance_to_strata, "m")
  mtch <- vector(mode = "integer",length = nrow(x))
  for(j in 1:nrow(miss_dist)){
    tmp <- which.min(miss_dist[j,])
    mtch[j] <- ifelse(miss_dist[j,tmp] < distance_to_strata, tmp,NA)
  }
  
  strata_details <- strata_map |> 
    sf::st_drop_geometry()
  strata_details_join <- strata_details[mtch,]
  outside_ret <- dplyr::bind_cols(outside,
                                  strata_details_join)  
  return(outside_ret)
}
