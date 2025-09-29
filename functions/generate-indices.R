
generate_indices <- function(
    model_fit = NULL,
    meta_strata = NULL, # df with columns strata, strata_name, and optional weights
    meta_years = NULL,
    quantiles = c(0.025, 0.05, 0.25, 0.75, 0.95, 0.975),
    regions = c("survey_wide","strata_level"),
    regions_index = NULL, # alternate post-hoc combinations of strata to form new regions - data frame with strata_name and region
    alternate_n = "n",
    gam_smooths = FALSE,
    start_year = NULL,
    max_backcast = NULL,
    drop_exclude = FALSE,
    hpdi = FALSE,
    quiet = FALSE) {
  
  # Checks
  check_numeric(quantiles)
  check_numeric(start_year, max_backcast, allow_null = TRUE)
  check_logical(drop_exclude, quiet, hpdi)
  
  # Get data

  if(! "weights" %in% names(meta_strata)){
    meta_strata <- meta_strata %>% 
      mutate(weights = 1)
  }
meta_strata <- meta_strata %>% 
  arrange(stratum)
 
  
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
  
  # Meta strata data
  meta_strata <- meta_strata %>%
    dplyr::mutate(survey_wide = "survey_wide",
                  strata_level = .data$strata_name)
  
  # Adding extra regions
  if(!is.null(regions_index)) {
    
    # Check if strata_names don't match
    if(!all(meta_strata$strata_name %in% regions_index$strata_name)){
      stop("'strata_name's in the `regions_index` don't match 'strata_name's ",
           "in the data. ",
           "See `model_output$meta_strata` for the strata to match",
           call. = FALSE)
    }
    
    # Keep only relevant regions
    r <- regions[!regions %in% c("survey_wide", "strata_level")]
    regions_index <- regions_index %>%
      dplyr::select("strata_name", dplyr::all_of(r)) %>%
      dplyr::arrange(.data$strata_name) %>%
      dplyr::distinct() %>%
      dplyr::mutate(strata_name = as.character(.data$strata_name))
    
    # Add new regional definitions to existing meta_strata
    meta_strata <- meta_strata %>%
      dplyr::select(-dplyr::any_of(r)) %>% # Remove any existing regions
      dplyr::left_join(regions_index, by = "strata_name") # Join in new
  }
  
  
  # Calculate strata/year-level observation statistics
  obs_strata <- raw_data_sel %>%
    dplyr::select("stratum", "year", "first_year", "count") %>%
    dplyr::arrange(.data$stratum, .data$year, .data$count) %>%
    dplyr::group_by(.data$stratum, .data$year, .data$first_year) %>%
    dplyr::summarize(obs_mean = mean(.data$count, na.rm = TRUE),
                     n_routes = sum(!is.na(.data$count)),
                     n_non_zero = sum(.data$count > 0, na.rm = TRUE),
                     strata_remove_flag = 0, .groups = "drop")
  
  n_routes_total <- raw_data_sel %>%
    dplyr::select("stratum", "route_id") %>%
    dplyr::group_by(.data$stratum) %>% 
    dplyr::summarize(n_routes_total = length(unique(route_id)), .groups = "drop")
  
  obs_strata <- obs_strata %>% 
    dplyr::inner_join(n_routes_total,
                      by = "stratum")
  
  indices <- dplyr::tibble()
  N_all <- list()
  
  for(rr in regions) { #selecting the type of composite region
    
    if(!quiet) message("Processing region ", rr)
    
    # Calculate strata-level info for sub-regions in this composite region
    meta_strata_sub <- meta_strata %>%
      # Ensure region columns are character
      dplyr::mutate("{rr}" := as.character(.data[[rr]])) %>%
      dplyr::group_by(.data[[rr]]) %>%
      dplyr::mutate(area_weight_non_zero = .data$weights/sum(.data$weights)) %>%
      dplyr::ungroup()
    
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
        dplyr::across(.cols = c("obs_mean", "n_routes",
                                "n_routes_total",
                                "n_non_zero", "flag_year"),
                      ~ sum(.x, na.rm = TRUE)),
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
      dplyr::summarize(N = purrr::map(.data$data, calc_weights, .env$n_sub),
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
                                  .data$obs_mean)) %>%
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
  
  meta_strata <- dplyr::select(meta_strata,
                               "strata_name", "stratum", "weights",
                               dplyr::all_of(.env$regions))
  
  
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



calc_weights <- function(data, n) {
  
  # Weight each sampled n
  n_weight <- n[, data$strata_name, , drop = FALSE]
  
  # Use numbers for indexing as is slightly faster
  for (i in seq_len(dim(n_weight)[1])) {       # iter
    for (j in seq_len(dim(n_weight)[2])) {     # strata_name
      n_weight[i, j, ] <- n_weight[i, j, ] * meta_strata_sub$area_weight_non_zero[j]
      
    }
  }
  
  # Sum over strata
  apply(n_weight, c(1, 3), sum)
}
