

## altered version of generate_indices that can handle partial overlap between
## strata and regions. e.g., it can deal with degree blocks that overlap multiple
## regions and so some portion of the information on trajectory from that degree block
## needs to be included in multiple regions
## it also can then estimate broader regional trajectories based on summaries of intermediate
## regional trajectories - e.g., all testing was done with the bbs-strata estimate used
## to summarize up to national trajectories and provincial trajectories

generate_indices_shared <- function(
    model_fit = NULL,
    meta_strata = NULL, # df with columns stratum, strata_name
    meta_years = NULL,
    raw_data = NULL,
    quantiles = c(0.025, 0.05, 0.25, 0.75, 0.95, 0.975),
    regions = c("strata_level","bbs_strata"),
    regional_sub = c("strata_level","strata_level"),# same length as regions indicates which regions form the subunits of each level
    regions_index = NULL,
    alternate_n = "n",
    gam_smooths = FALSE,
    start_year = NULL,
    max_backcast = NULL,
    drop_exclude = FALSE,
    hpdi = FALSE,
    quiet = FALSE,
    fractional_weights = NULL, # optional dataframe that quantifies the fractional contribution each strata makes to a broader region (e.g. if not every strata is clearly allocated to a single region)
    rel_abundance_weights = NULL #list same length and names as regions, each element is a dataframe with columns for each week of the survey
                                # and rows for each subregion that forms the larger regions identified in in the list names
    ) {
  
  # Checks
  check_numeric(quantiles)
  check_numeric(start_year, max_backcast, allow_null = TRUE)
  check_logical(drop_exclude, quiet, hpdi)
  if(is.null(rel_abundance_weights)){
    stop("first element of rel_abundance_weights must indicate the relative abundance weights
         for each strata.")
  }
  if(names(rel_abundance_weights)[1] != "strata_level"){
    stop("first element of rel_abundance_weights must indicate the relative abundance weights
         for each strata and must be named strata_level.")
  }
  if(!"strata_level" %in% regions){
    warning("indices must be calculated for strata, 
            adding strata_level indices to list of regions")
    regions <- c("strata_level",regions)
  }
  
  if(length(rel_abundance_weights) == 1 &
     length(regions) > 1){
    warning(paste("rel_abundance_weights only include strata_level weights.",
                  "assuming that all regions including
                  ",
                  paste(regions[-which(regions == "strata_level")], collapse = " & "),
                  ", 
                  are based on weighted summaries of strata_level trajectories.
                  If some regions should be considered weighted summaries of
                  intermediate regions (e.g., if survey_wide should be weighted
                  summaries of bcr trajectories), then supply a dataframe that
                  lists which regions should be combined to form the higher-level
                  regional estimates to the argument region_hierarchy."))
    for(r in regions[-which(regions == "strata_level")]){
      rel_abundance_weights[[r]] <- rel_abundance_weights[["strata_level"]]
    }
  }
  
 # Get data
  

meta_strata <- meta_strata %>% 
  arrange(stratum)
 

  names(regional_sub) <- regions

  # Start years
  if(!is.null(start_year)){
    inity <- min(meta_years$year)-1
    
    if(inity > start_year){
      warning(
        "Value of ", start_year, " for `start_year` is earlier than the ",
        "earliest year of the data, using ",
        start_year <- min(meta_years$year),
        " instead", call. = FALSE)
    }
    
  } else{
    start_year <- min(meta_years$year)
  }
  end_year <- max(meta_years$year)
  
  raw_data_sel <- raw_data %>%
    # Set start year
    dplyr::group_by(.data[["stratum"]]) %>%
    dplyr::mutate(first_year = min(.data$year[.data$count > 0],
                                   na.rm = TRUE)) %>%
    dplyr::ungroup() %>%
    # Trim year range
    dplyr::filter(.data$year >= .env$start_year)
  
  # After trimming data
  n_years <- max(meta_years$yr) - min(meta_years$yr) + 1
  
  # Backcast
  if(is.null(max_backcast)) max_backcast <- n_years
  
  # Posterior draws
  n <- samples_to_array(model_fit, alternate_n,
                        meta_strata = meta_strata,
                        meta_years = meta_years,
                        years_to_keep = start_year:end_year)
  
  weeks_sample <- raw_data %>% 
    dplyr::select(week) %>% 
    drop_na() %>% 
    dplyr::slice_sample(n = dim(n)[1],
                        replace = TRUE)
  


  # Meta strata data
  meta_strata <- meta_strata %>%
    dplyr::mutate(strata_level = .data$strata_name)
  
  
  # Adding extra regions
  # if(!is.null(regions_index)) {
  #   
  #   # Check if strata_names don't match
  #   if(!all(meta_strata$strata_name %in% regions_index$strata_name)){
  #     stop("'strata_name's in the `regions_index` don't match 'strata_name's ",
  #          "in the data. ",
  #          "See `model_output$meta_strata` for the strata to match",
  #          call. = FALSE)
  #   }
  #   
  #   # Keep only relevant regions
  #   r <- regions[!regions %in% c("survey_wide", "strata_level")]
  #   regions_index <- regions_index %>%
  #     dplyr::select("strata_name", dplyr::all_of(r)) %>%
  #     dplyr::arrange(.data$strata_name) %>%
  #     dplyr::distinct() %>%
  #     dplyr::mutate(strata_name = as.character(.data$strata_name))
  #   
  #   # Add new regional definitions to existing meta_strata
  #   meta_strata <- meta_strata %>%
  #     dplyr::select(-dplyr::any_of(r)) %>% # Remove any existing regions
  #     dplyr::left_join(regions_index, by = "strata_name") # Join in new
  # }
  
  
  # Calculate strata/year-level observation statistics
  obs_strata <- raw_data_sel %>%
    dplyr::select("stratum", "year", "first_year", "count") %>%
    dplyr::group_by(.data$stratum) %>%
    tidyr::complete(year = seq(.env$start_year, .env$end_year),
                    .data$first_year) %>%
    dplyr::arrange(.data$stratum, .data$year, .data$count) %>%
    dplyr::group_by(.data$stratum, .data$year, .data$first_year) %>%
    dplyr::summarize(obs_mean = mean(.data$count, na.rm = TRUE),
                     n_routes = sum(!is.na(.data$count)),
                     n_non_zero = sum(.data$count > 0, na.rm = TRUE),
                     strata_remove_flag = 0, .groups = "drop")
  
  n_routes_total <- raw_data_sel %>%
    dplyr::select("stratum", "route_id", "count") %>%
    dplyr::group_by(.data$stratum) %>% 
    dplyr::summarize(n_routes_total = length(unique(route_id)), .groups = "drop")
  
  obs_strata <- obs_strata %>% 
    dplyr::inner_join(n_routes_total,
                       by = "stratum") #%>% 
    # dplyr::mutate(obs_mean = obs_mean) ## obs_mean is a proportion value, proportion of this year's mean count relative to the mean count across all years in the stratum
    # 
  indices <- dplyr::tibble()
  N_all <- list()
  
  

# strata level indices ----------------------------------------------------

  
  rr <- "strata_level"
  if(!quiet) message("Processing region ", rr)
  rel_abund <- rel_abundance_weights[[rr]]
  weight_array <- as.matrix(t(rel_abund[,weeks_sample$week]))
  dimnames(weight_array) <- dimnames(n)[1:2]
  
  
    meta_strata_sub <- meta_strata %>%
      dplyr::mutate(area_weight_non_zero = 1) 
 
  
  # Calculate observation statistics for this composite region
  obs_region <- obs_strata %>%
    dplyr::inner_join(meta_strata_sub, by = "stratum") %>%
    dplyr::mutate(obs_mean = .data$obs_mean * .data$area_weight_non_zero) %>%
    dplyr::group_by(.data[[rr]], .data$stratum)
  
  # Flag strata to remove due to max_backcast
  # - Flag first max_backcast no. years IF:
  #    - If no obs in those years, AND
  #    - first_year is AFTER the current start of the data range
  #      (i.e. flag data that has no true counts in it.)
  
  obs_region <- obs_region %>%
    dplyr::mutate(
      flag_remove = sum(.data$n_non_zero[seq_len(.env$max_backcast)]) < 1 &
        .data$first_year > .env$start_year,
      flag_year = dplyr::if_else(.data$flag_remove &
                                   .data$year < .data$first_year,
                                 .data$area_weight_non_zero, 0))
  
  # Mark strata included/excluded
  obs_region <- obs_region %>%
    dplyr::ungroup() %>% 
    dplyr::group_by(.data[[rr]], .data$year) %>%
    dplyr::mutate(
      strata_included = paste0(.data$strata_name[!.data$flag_remove],
                               collapse = " ; "),
      strata_excluded = paste0(.data$strata_name[.data$flag_remove],
                               collapse = " ; "))
  
  # Exclude if requested
  if(drop_exclude) {
    rm <- unique(obs_region$strata_name[obs_region$flag_remove])
    
    obs_region <- dplyr::filter(obs_region, !.data$flag_remove)
    meta_strata_sub <- dplyr::filter(meta_strata_sub,
                                     !.data$strata_name %in% rm)
    n_sub <- n[, unique(obs_region$strata_name), ] # Keep only good
  } else n_sub <- n
  
  
  # Missing data (a missing year identified by n_not_missing = 0)
  missing_yrs <- obs_region %>%
    dplyr::ungroup() %>%
    dplyr::select("stratum", "year", "obs_mean") %>%
    dplyr::distinct() %>%
    dplyr::group_by(.data$year) %>%
    dplyr::summarize(n_not_missing = sum(.data$obs_mean, na.rm = TRUE)) %>%
    dplyr::filter(.data$n_not_missing == 0) %>%
    dplyr::pull(.data$year)
  
  obs_region <- obs_region %>%
    dplyr::group_by(.data$strata_included, .data$strata_excluded,
                    .add = TRUE) %>%
    dplyr::summarize(
      dplyr::across(.cols = c("n_routes",
                              "n_routes_total",
                              "n_non_zero", "flag_year"),
                    ~ sum(.x, na.rm = TRUE)),
      dplyr::across(.cols = c("obs_mean"),
                    ~ mean(.x, na.rm = TRUE)),
      .groups = "drop")
  
  
  if(hpdi){
    calc_quantiles <- calc_quantiles_hpdi
  }else{
    calc_quantiles <- calc_quantiles_original
  }
  # Calculate sample statistics for this composite region
  samples <- meta_strata_sub %>%
    # Create back up col for use in calculations
    tidyr::nest(data = -dplyr::any_of(rr)) %>%
    dplyr::group_by(.data[[rr]]) %>%
    dplyr::summarize(N = purrr::map(.data$data, calc_weights, .env$n_sub, .env$weight_array),
                     N_names = paste0(rr, "_", .data[[rr]]),
                     Q = purrr::map(.data$N, calc_quantiles,
                                    .env$quantiles)) %>%
    dplyr::mutate(r = .env$rr)
  
  # Save sample stats for output
  N_all <- append(N_all, stats::setNames(samples$N, samples$N_names))
  
  # Calculate data summaries for output
  indices <- obs_region %>%
    dplyr::mutate(
      backcast_flag = 1 - .data$flag_year,
      region_type = .env$rr,
      # Replace with NA, if entire year missing
      obs_mean = dplyr::if_else((.data$year %in% .env$missing_yrs) |
                                  (.data$n_routes == 0),
                                NA_real_,
                                .data$obs_mean))  %>% 
    # Add in quantiles
    dplyr::left_join(tidyr::unnest(samples, "Q"), by = c(rr, "year")) %>%
    # Clean up
    dplyr::rename(region = dplyr::any_of(rr)) %>%
    dplyr::select("year", "region", "region_type",
                  "strata_included", "strata_excluded",
                  "index", dplyr::contains("index_q"),
                  "obs_mean", "n_routes", "n_routes_total", "n_non_zero",
                  "backcast_flag") %>%
    dplyr::bind_rows(indices, .)
  
  
  

# End of strata calculations ----------------------------------------------

# calculate indices for bbs_strata ----------------------------------------------
  
  rr <- "bbs_strata"
  
  
  if(!quiet) message("Processing region ", rr)
  
  rel_abund <- rel_abundance_weights[[rr]]
  weight_array <- as.matrix(t(rel_abund[,weeks_sample$week]))
  dimnames(weight_array) <- dimnames(n)[1:2]
  
 
  meta_strata_sub <- fractional_weights %>%
    dplyr::mutate(area_weight_non_zero = weights) 
  
  
  
  # Calculate observation statistics for this composite region
  obs_region <- obs_strata %>%
    dplyr::inner_join(meta_strata_sub, by = "stratum",
                      relationship = "many-to-many") %>%
    dplyr::mutate(obs_mean = .data$obs_mean) |> # * .data$area_weight_non_zero)  %>%
    dplyr::group_by(.data[[rr]]) 
  
  # Flag strata to remove due to max_backcast
  # - Flag first max_backcast no. years IF:
  #    - If no obs in those years, AND
  #    - first_year is AFTER the current start of the data range
  #      (i.e. flag data that has no true counts in it.)
  
  obs_region <- obs_region %>%
    dplyr::mutate(
      flag_remove = sum(.data$n_non_zero[seq_len(.env$max_backcast)]) < 1 &
        .data$first_year > .env$start_year,
      flag_year = dplyr::if_else(.data$flag_remove &
                                   .data$year < .data$first_year,
                                 .data$area_weight_non_zero, 0))
  
  # Mark strata included/excluded
  obs_region <- obs_region %>%
    dplyr::ungroup() %>% 
    dplyr::group_by(.data[[rr]], .data$year) %>%
    dplyr::mutate(
      strata_included = paste0(.data$strata_name[!.data$flag_remove],
                               collapse = " ; "),
      strata_excluded = paste0(.data$strata_name[.data$flag_remove],
                               collapse = " ; "))
  
  # Exclude if requested
  if(drop_exclude) {
    rm <- unique(obs_region$strata_name[obs_region$flag_remove])
    
    obs_region <- dplyr::filter(obs_region, !.data$flag_remove)
    meta_strata_sub <- dplyr::filter(meta_strata_sub,
                                     !.data$strata_name %in% rm)
    n_sub <- n[, unique(obs_region$strata_name), ] # Keep only good
  } else n_sub <- n
  
  
  # Missing data (a missing year identified by n_not_missing = 0)
  missing_yrs <- obs_region %>%
    dplyr::ungroup() %>%
    dplyr::select(.env$rr, "year", "obs_mean") %>%
    dplyr::distinct() %>%
    dplyr::group_by(.data$year) %>%
    dplyr::summarize(n_not_missing = sum(.data$obs_mean, na.rm = TRUE)) %>%
    dplyr::filter(.data$n_not_missing == 0) %>%
    dplyr::pull(.data$year)
  

  obs_region <- obs_region %>%
    dplyr::group_by(.data$strata_included, .data$strata_excluded,
                    .add = TRUE) %>%
    dplyr::summarize(
      dplyr::across(.cols = c("n_routes",
                              "n_routes_total",
                              "n_non_zero", "flag_year"),
                    ~ sum(.x, na.rm = TRUE)),
      dplyr::across(.cols = c("obs_mean"),
                    ~ weighted_mean(.x, area_weight_non_zero, na.rm = TRUE)),
      .groups = "drop")
  
  

  # Calculate sample statistics for this composite region
  samples <- meta_strata_sub %>%
    # Create back up col for use in calculations
    tidyr::nest(data = -dplyr::any_of(rr)) %>%
    dplyr::group_by(.data[[rr]]) %>%
    dplyr::summarize(N = purrr::map(.data$data, calc_weights, .env$n_sub, .env$weight_array),
                     N_names = paste0(rr, "_", .data[[rr]]),
                     Q = purrr::map(.data$N, calc_quantiles,
                                    .env$quantiles)) %>%
    dplyr::mutate(r = .env$rr)
  
  # Save sample stats for output
  N_all <- append(N_all, stats::setNames(samples$N, samples$N_names))
  
  # Calculate data summaries for output
  indices <- obs_region %>%
    dplyr::mutate(
      backcast_flag = 1 - .data$flag_year,
      region_type = .env$rr,
      # Replace with NA, if entire year missing
      obs_mean = dplyr::if_else((.data$year %in% .env$missing_yrs) |
                                  (.data$n_routes == 0),
                                NA_real_,
                                .data$obs_mean))  %>% 
    # Add in quantiles
    dplyr::left_join(tidyr::unnest(samples, "Q"), by = c(rr, "year")) %>%
    # Clean up
    dplyr::rename(region = dplyr::any_of(rr)) %>%
    dplyr::select("year", "region", "region_type",
                  "strata_included", "strata_excluded",
                  "index", dplyr::contains("index_q"),
                  "obs_mean", "n_routes", "n_routes_total", "n_non_zero",
                  "backcast_flag") %>%
    dplyr::bind_rows(indices, .)
  
  
  
  
  
  
  
  
  ### as of here it seems to be working! now to calculate the weighted summaries of the
  ### bbs_strata indices to create the broader regional estimates
  
  
  
  
  
    
  for(rr in regions[-c(1,2)]) { #selecting the type of composite region
    
    if(!quiet) message("Processing region ", rr)
    
    sub_regions <- regional_sub[[rr]]
    
    rel_abund <- rel_abundance_weights[[rr]] 
    


    # Calculate strata-level info for sub-regions in this composite region
    sub_regions_incl <- stringr::str_replace(names(N_all)[which(stringr::str_detect(names(N_all),sub_regions))],
                                             paste0(sub_regions,"_"),"")
    meta_strata_sub <- data.frame(x = sub_regions_incl)
    names(meta_strata_sub) <- sub_regions 
    meta_strata_sub <- meta_strata_sub %>%
      dplyr::mutate(area_weight_non_zero = 1) |> 
      dplyr::left_join(regions_index[[rr]], # left join drops sub_regions with no estimates
                       by = sub_regions)
    
    # Calculate observation statistics for this composite region
    obs_region <- indices |> 
      dplyr::filter(region %in% sub_regions_incl) |> 
      dplyr::select(year, region, strata_included, strata_excluded,
                    obs_mean,n_routes,n_routes_total,
                    n_non_zero,backcast_flag) %>%
      dplyr::mutate(
        flag_remove = sum(.data$n_non_zero[seq_len(.env$max_backcast)]) < 1)
    
    dnam <- dimnames(n)
    names(dnam)[2]<- sub_regions
    dnam[[sub_regions]] <- sub_regions_incl
    
    n_sub <- array(data = NA,dim = c(dim(n)[1],length(sub_regions_incl),dim(n)[3]),
                   dimnames = dnam)
                   
    for(jj in sub_regions_incl){
      n_sub[,jj,] <- N_all[[paste0(sub_regions,"_",jj)]]
    }
    
    
    # Exclude if requested
    if(drop_exclude) {
      rm <- unique(obs_region$region[obs_region$flag_remove])
      
      obs_region <- dplyr::filter(obs_region, !.data$flag_remove)
      meta_strata_sub <- dplyr::filter(meta_strata_sub,
                                       !.env$sub_regions %in% rm)
        
        n_sub <- n_sub[, unique(obs_region$region), ] # Keep only good
    } 
    
    weeks_sample <- raw_data %>% 
      dplyr::select(week) %>% 
      drop_na() %>% 
      dplyr::slice_sample(n = dim(n)[1],
                          replace = TRUE)
    
    weeks <- unique(raw_data$week)
    
    rel_abund <- rel_abund |> 
      dplyr::select(dplyr::starts_with(weeks),
                    dplyr::starts_with(sub_regions)) |> 
      dplyr::filter(.data[[sub_regions]] %in% sub_regions_incl)
      
      
    obs_region <- obs_region |> 
      dplyr::left_join(regions_index[[rr]], # left join drops sub_regions with no estimates
                       by = c("region" = sub_regions)) |> 
      dplyr::group_by(.data[[rr]],year) |> 
      dplyr::summarise(strata_included = paste(unique(strata_included),collapse = "; "),
                       strata_excluded = paste(unique(strata_excluded),collapse = "; "),
                       obs_mean = sum(obs_mean * n_routes,na.rm = TRUE)/sum(n_routes,na.rm = TRUE),
                       n_routes = sum(n_routes,na.rm = TRUE),
                       n_routes_total = sum(n_routes_total,na.rm = TRUE),
                       n_non_zero = sum(n_non_zero,na.rm = TRUE),
                       backcast_flag = mean(backcast_flag),
                       .groups = "drop") |> 
      dplyr::mutate(region_type = rr)
    

      
    weight_array <- as.matrix(t(rel_abund[,weeks_sample$week]))
    
    dimnames(weight_array) <- dimnames(n_sub)[1:2]
    


    # Calculate sample statistics for this composite region
    samples <- meta_strata_sub %>%
      # Create back up col for use in calculations
      tidyr::nest(data = -dplyr::any_of(rr)) %>%
      dplyr::group_by(.data[[rr]]) %>%
      dplyr::summarize(N = purrr::map(.data$data, calc_weights_shared, .env$n_sub, .env$weight_array, .env$sub_regions),
                       N_names = paste0(rr, "_", .data[[rr]]),
                       Q = purrr::map(.data$N, calc_quantiles,
                                      .env$quantiles)) %>%
      dplyr::mutate(r = .env$rr)
    
    # Save sample stats for output
    N_all <- append(N_all, stats::setNames(samples$N, samples$N_names))
    
    # Calculate data summaries for output
    indices <- obs_region %>%
      # Add in quantiles
      dplyr::left_join(tidyr::unnest(samples, "Q"), by = c(rr, "year")) %>%
      # Clean up
      dplyr::rename(region = dplyr::any_of(rr)) %>%
      dplyr::select("year", "region", "region_type",
                    "strata_included", "strata_excluded",
                    "index", dplyr::contains("index_q"),
                    "obs_mean", "n_routes", "n_routes_total", "n_non_zero",
                    "backcast_flag") %>%
      dplyr::bind_rows(indices, .)
  }
  

  
  if(gam_smooths){
    
    N_all_smooth <- N_all
    for(j in names(N_all)){
      if(!quiet) message("Creating gam smooths for region ",j )
      N_all_smooth[[j]] <- t(apply(N_all[[j]],1,gam_sm))
      dimnames(N_all_smooth[[j]]) <- dimnames(N_all[[j]])
    }
  }else{
    N_all_smooth <- NA
  }
  
  
  list("indices" = indices,
       "samples" = N_all,
       "meta_data" = (list("regions" = regions,
                                 "start_year" = start_year,
                                 "n_years" = n_years,
                                 "hpdi_indices" = hpdi)),
       "meta_strata" = meta_strata,
       "raw_data" = raw_data, # Original data before trimming
       "gam_smooth_samples" = N_all_smooth
  )
}









calc_quantiles_hpdi <- function(N, quantiles) {
  apply(N, 2, interval_function_hpdi, probs = c(quantiles, 0.5)) %>%
    t() %>%
    as.data.frame() %>%
    stats::setNames(c(paste0("index_q_", quantiles), "index")) %>%
    dplyr::bind_cols(year = as.numeric(dimnames(N)$year))
}

calc_quantiles_original <- function(N, quantiles) {
  apply(N, 2, stats::quantile, probs = c(quantiles, 0.5)) %>%
    t() %>%
    as.data.frame() %>%
    stats::setNames(c(paste0("index_q_", quantiles), "index")) %>%
    dplyr::bind_cols(year = as.numeric(dimnames(N)$year))
}


calc_alt_names <- function(r, region_names) {
  col_region_name <- dplyr::case_when(r == "prov_state" ~ "province_state",
                                      TRUE ~ r)
  
  region_alt_name <- dplyr::bind_cols(
    {{r}} := region_names[[r]],
    region_alt = region_names[[col_region_name]]) %>%
    dplyr::distinct()
  
  if(r == "bcr") {
    region_alt_name <- dplyr::mutate(
      region_alt_name, region_alt = paste0("BCR_", .data$region_alt))
  }
  
  region_alt_name
}


gam_sm <- function(i) {
  df <- data.frame(i = log(i),
                   y = as.integer(names(i)))
  sm <- mgcv::gam(data = df,
                  formula = i~s(y))
  smf <- exp(sm$fitted.values)
}




samples_to_array <- function(model_fit, alternate_n, years_to_keep,
                             meta_strata, meta_years) {
  
  # Extract samples
  n <- model_fit$draws(variables = alternate_n,
                                    format = "draws_matrix")
  # Determine dim names
  strata_name <- meta_strata$strata_name
  year <- meta_years$year
  
  
  # Transform samples to array with appropriate dimnames
  n <- array(as.vector(n),
             dim = c(posterior::ndraws(n), length(strata_name), length(year)),
             dimnames = list("iter" = 1:posterior::ndraws(n),
                             "strata_name" = strata_name,
                             "year" = year))
  
  # Filter to years selected
  years_to_keep <- years_to_keep[years_to_keep %in% year]
  n[ , , as.character(years_to_keep)]
}



calc_weights_shared <- function(data, n, weight_array, sub_regions) {
  
  # Weight each sampled n
  n_weight <- n[, data[[sub_regions]], , drop = FALSE]
  weight_sub <- weight_array[,data[[sub_regions]], drop = FALSE]
  area_weights <- data$area_weight_non_zero
  
    # Use numbers for indexing as is slightly faster
  for (i in seq_len(dim(n_weight)[1])) {       # iter
    wts <- weight_sub[i,] *  area_weights
    wts <- wts/sum(wts)
    for (j in seq_len(dim(n_weight)[2])) {     # sub_regions
      n_weight[i, j, ] <- n_weight[i, j, ] * wts[j]
      
    }
  }
  
  # Sum over strata
  apply(n_weight, c(1, 3), sum)
}

 # purrr::map(.data$data, calc_weights, .env$n_sub, .env$meta_strata_sub, .env$weight_array)

calc_weights <- function(data, n, weight_array) {

  # Weight each sampled n
  n_weight <- n[, data$strata_name, , drop = FALSE]
  weight_sub <- weight_array[,data$strata_name, drop = FALSE]
  weight_sub <- weight_sub/rowSums(weight_sub)
  area_weights <- data$area_weight_non_zero
  # Use numbers for indexing as is slightly faster
  for (i in seq_len(dim(n_weight)[1])) {       # iter
      wts <- weight_sub[i,] *  area_weights
      wts <- wts/sum(wts)
      for (j in seq_len(dim(n_weight)[2])) {     # sub_regions
        n_weight[i, j, ] <- n_weight[i, j, ] * wts[j]
        
    }
  }

  # Sum over strata
  apply(n_weight, c(1, 3), sum)
}


weighted_sum <- function(x,w,na.rm = TRUE){
  ss <- sum(x*w,na.rm = na.rm)
}

weighted_mean <- function(x,w,na.rm = TRUE){
  wt <- w/sum(w,na.rm = na.rm)
  ss <- sum(x*wt,na.rm = na.rm)
  return(ss)
}


