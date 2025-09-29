## generic trend function
## input = posterior distribution from model with annual indices




generate_trends <- function(indices,
                            min_year = NULL,
                            max_year = NULL,
                            quantiles = c(0.025, 0.05, 0.25, 0.75, 0.95, 0.975),
                            slope = FALSE,
                            gam = FALSE,
                            prob_decrease = NULL,
                            prob_increase = NULL,
                            hpdi = FALSE) {
  
  # Checks
  check_logical(slope, hpdi, gam)
  check_numeric(quantiles)
  check_numeric(min_year, max_year, quantiles, prob_decrease, prob_increase,
                allow_null = TRUE)
  check_range(quantiles, c(0, 1))
  check_range(prob_decrease, c(0, 100))
  check_range(prob_increase, c(0, Inf))
  
  if(hpdi){
    calc_quantiles <- interval_function_hpdi
  }else{
    calc_quantiles <- stats::quantile
  }
  
  start_year <- indices[["meta_data"]]$start_year
  n_years <- indices[["meta_data"]]$n_years
  indx <- indices[["indices"]]
  
  if(!gam){
    ind_samples <- indices[["samples"]]
  }
  
  if(gam){
    ind_samples <- indices[["gam_smooth_samples"]]
    if(all(is.na(ind_samples))){
      stop("indices object has no gam_smooth_samples. ",
           "To calculate gam-based trends, generate_trends requires the output ",
           "from generate_indices with gam_smooths = TRUE")
    }
  }
  if(is.null(min_year)) {
    min_year <- start_year
  } else {
    
    if(min_year < start_year) {
      message("`min_year` is before the date range, using minimum year of ",
              "the data (", min_year <- start_year, ") instead.")
    }
  }
  
  if (is.null(max_year)) {
    max_year <- max(indx$year)
  } else if(max_year > max(indx$year)) {
    message("`max_year` is beyond the date range, using maximum year of ",
            "the data (", max_year <- max(indx$year), ") instead.")
  }
  
  # For indexing
  min_year_num <- min_year - start_year + 1
  max_year_num <- max_year - start_year + 1
  
  
  trends <- indx %>%
    dplyr::filter(.data$year %in% min_year:max_year) %>%
    dplyr::group_by(.data$region, .data$region_type,
                    .data$strata_included, .data$strata_excluded) %>%
    dplyr::summarize(
      # Basic statistics
      rel_abundance = mean(.data$index),
      obs_rel_abundance = mean(.data$obs_mean,na.rm = TRUE),
      mean_n_routes = mean(.data$n_routes),
      n_routes = mean(.data$n_routes_total),
      backcast_flag = mean(.data$backcast_flag),
      
      # Metadata
      start_year = .env$min_year,
      end_year = .env$max_year, .groups = "keep") %>%
    
    dplyr::mutate(
      n_strata_included = purrr::map_dbl(
        .data$strata_included, ~length(unlist(stringr::str_split(.x, " ; ")))),
      
      # Add in samples
      
      n = purrr::map2(.data$region_type, .data$region,
                      ~ind_samples[[paste0(.x, "_", .y)]]),
      # Calculate change start to end for each iteration
      ch = purrr::map(.data$n,
                      ~.x[, .env$max_year_num] / .x[, .env$min_year_num]),
      # Calculate change as trend for each iteration
      tr = purrr::map(
        .data$ch,
        ~100 * ((.x^(1/(.env$max_year_num - .env$min_year_num))) - 1)),
      
      # Median and percentiles of trend per region
      trend = purrr::map_dbl(.data$tr, stats::median),
      trend_q = purrr::map_df(
        .data$tr,
        ~stats::setNames(calc_quantiles(.x, quantiles, names = FALSE),
                         paste0("trend_q_", quantiles))),
      
      # Percent change and quantiles thereof per region
      percent_change = purrr::map_dbl(.data$ch, ~100 * (stats::median(.x) - 1)),
      pc_q = purrr::map_df(
        .data$ch, ~stats::setNames(
          100 * (calc_quantiles(.x, quantiles, names = FALSE) - 1),
          paste0("percent_change_q_", quantiles)))) %>%
    dplyr::ungroup() %>%
    tidyr::unnest(cols = c("trend_q", "pc_q")) %>%
    dplyr::arrange(.data$region_type, .data$region)
  
  
  
  # Reliability Criteria
  q1 <- quantiles[1]
  q2 <- quantiles[length(quantiles)]
  q <- (q2 - q1) * 100
  
  trends <- trends %>%
    dplyr::mutate(
      "width_of_{{q}}_percent_credible_interval" :=
        .data[[paste0("trend_q_", q2)]] - .data[[paste0("trend_q_", q1)]])
  
  # Optional slope based trends
  if(slope) {
    trends <- trends %>%
      dplyr::mutate(
        sl_t = purrr::map(.data$n, calc_slope,
                          .env$min_year_num, .env$max_year_num),
        slope_trend = purrr::map_dbl(.data$sl_t, stats::median),
        slope_trend_q = purrr::map_df(
          .data$sl_t, ~stats::setNames(
            calc_quantiles(.x, quantiles, names = FALSE),
            paste0("slope_trend_q_", quantiles)))) %>%
      tidyr::unnest("slope_trend_q") %>%
      dplyr::mutate(
        "width_of_{q}_percent_credible_interval_slope" :=
          .data[[paste0("slope_trend_q_", q2)]] -
          .data[[paste0("slope_trend_q_", q1)]])
  }
  
  
  # # Optional gam based trends
  # if(gam) {
  #   trends <- trends %>%
  #     dplyr::mutate(
  #       sl_t = purrr::map(.data$n, calc_gam,
  #                         .env$min_year_num, .env$max_year_num),
  #       gam_trend = purrr::map_dbl(.data$sl_t, stats::median),
  #       gam_trend_q = purrr::map_df(
  #         .data$sl_t, ~stats::setNames(
  #           calc_quantiles(.x, quantiles, names = FALSE),
  #           paste0("gam_trend_q_", quantiles)))) %>%
  #     tidyr::unnest("gam_trend_q") %>%
  #     dplyr::mutate(
  #       "width_of_{q}_percent_credible_interval_gam" :=
  #         .data[[paste0("gam_trend_q_", q2)]] -
  #         .data[[paste0("gam_trend_q_", q1)]])
  # }
  
  
  # Model conditional probabilities of population change during trends period
  if(!is.null(prob_decrease)) {
    trends <- trends %>%
      dplyr::mutate(
        pch = purrr::map(.data$ch, ~100 * (.x - 1)),
        pch_pp = purrr::map_df(.data$pch, calc_prob_crease,
                               .env$prob_decrease, type = "decrease")) %>%
      tidyr::unnest("pch_pp") %>%
      dplyr::select(-"pch")
  }
  
  if(!is.null(prob_increase)){
    trends <- trends %>%
      dplyr::mutate(
        pch = purrr::map(.data$ch, ~100 * (.x - 1)),
        pch_pp = purrr::map_df(.data$pch, calc_prob_crease,
                               .env$prob_increase, type = "increase")) %>%
      tidyr::unnest("pch_pp") %>%
      dplyr::select(-"pch")
  }
  
  trends <- trends %>%
    dplyr::select(-"n", -"ch", -"tr") %>%
    dplyr::select(
      "start_year", "end_year", "region", "region_type",
      "strata_included", "strata_excluded",
      dplyr::starts_with("trend"),
      dplyr::starts_with("percent_change"),
      dplyr::starts_with("slope"),
      dplyr::starts_with("gam"),
      dplyr::starts_with("width"),
      dplyr::starts_with("prob"),
      "rel_abundance", "obs_rel_abundance",
      "n_routes", "mean_n_routes",
      "n_strata_included", "backcast_flag")
  
  
  list("trends" = trends,
       "meta_data" = append(indices[["meta_data"]],
                            list("hpdi_trends" = hpdi,
                                 "gam_smooth_trends" = gam)),
       "meta_strata" = indices[["meta_strata"]],
       "raw_data" = indices[["raw_data"]])
}

bsl <- function(i, wy) {
  n <- length(wy)
  sy <- sum(i)
  sx <- sum(wy)
  ssx <- sum(wy^2)
  sxy <- sum(i*wy)
  
  (n * sxy - sx * sy) / (n * ssx - sx^2)
}

calc_slope <- function(n, min_year_num, max_year_num) {
  wy <- min_year_num:max_year_num
  
  ne <- log(n[, wy])
  m <-  t(apply(ne, 1, FUN = bsl, wy))
  
  as.vector((exp(m) - 1) * 100)
}
#
# gam_sl <- function(i, wy) {
#   df <- data.frame(i = i,
#                    y = 1:length(i))
#   sm <- mgcv::gam(data = df,
#                   formula = i~s(y))
#   smf <- sm$fitted.values[wy]
#   ny <- wy[2]-wy[1]
#   smt <- (smf[2]-smf[1])*(1/ny)
# }
#
# calc_gam <- function(n, min_year_num, max_year_num) {
#   wy <- c(min_year_num,max_year_num)
#
#   ne <- log(n)
#   m <-  t(apply(ne, 1, FUN = gam_sl, wy))
#
#   as.vector((exp(m) - 1) * 100)
# }

calc_prob_crease <- function(x, p, type = "decrease") {
  if(type == "decrease") f <- function(p) length(x[x < (-1 * p)]) / length(x)
  if(type == "increase") f <- function(p) length(x[x > p]) / length(x)
  
  vapply(p, FUN = f, FUN.VALUE = 1.1) %>%
    stats::setNames(paste0("prob_", type, "_", p, "_percent"))
}
