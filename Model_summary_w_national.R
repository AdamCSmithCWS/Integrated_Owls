## model exploration
## 

library(bbsBayes2)
library(tidyverse)
library(cmdstanr)
library(sf)
library(ebirdst)
library(patchwork)


source("functions/checks.R")
source("functions/generate_trends.R")
source("functions/generate_indices.R")

source("functions/generate_indices_shared.R")
source("functions/utils.R")


output_dir <- "output/"
year_end <- 2025
year_start <- 2000

# owl observations --------------------------------------------------------

obs_owl <- read_csv("data/owls2/NationalOwlData.csv") %>% 
  filter(RouteIdentifier != "NS050") %>%  # temporary removal of one problematic route that is replicated in two provinces%>% 
  mutate(unique_survey = paste(RouteIdentifier,survey_year,survey_month,survey_day, sep = "-"))

owl_species <- obs_owl %>% 
  group_by(CommonName) %>% 
  summarise(n_years = length(unique(survey_year)),
            n_routes = length(unique(RouteIdentifier)),
            n_obs = length(which(Count > 0))) %>% 
  arrange(n_obs) %>% 
  filter(n_years > 24) # species must have observations for 25 years (2001 through 2025)


# explore results ---------------------------------------------------------

trends_out <- NULL
converge <- NULL

#re_fit <- FALSE
for(sp in owl_species$CommonName){
  
  
  # sp <- "Barred Owl"
  
  #sp_id <- unique(sp_obs_owl$species_id)
  sp_ebird <- ebirdst::get_species(sp)
  if(!file.exists(paste0("data/pre_fit_data_",sp_ebird,".RData"))){
    next
  }
  load(paste0("data/pre_fit_data_",sp_ebird,".RData"))
  
  
  fit2 <- readRDS(paste0(output_dir,"fit_","new_scale","_",sp_ebird,".rds"))
  summ2 <- readRDS(paste0(output_dir,"fit_summary_","new_scale","_",sp_ebird,".rds"))
  
  
  summa <- summ2 %>% 
    mutate(species = sp,
           rhat_fail = ifelse(rhat > 1.05,TRUE,FALSE),
           ess_fail = ifelse(ess_bulk < 100, TRUE, FALSE)) 
  
  converge <- bind_rows(converge,summa)
  
  
  
  
  
  
  provs <- bbsBayes2::load_map("prov_state") %>% 
    filter(country_code == "CA") %>% 
    mutate(survey_region = ifelse(province_state %in% c("New Brunswick",
                                                        "Nova Scotia",
                                                        "Prince Edward Island"),
                                  "Maritimes",
                                  province_state)) %>% 
    select(province_state,survey_region,prov_state)
  
  provs_join <- provs |> 
    st_drop_geometry()
  
  strata_bbs <- bbsBayes2::load_map("bbs") |> 
    filter(country_code == "CA") |> 
    select(prov_state,strata_name) |> 
    rename(bbs_strata_name = strata_name) 
  
  
  
  
  
  strata_bbs_regions <- strata_bbs |> 
    inner_join(provs_join,
               by = "prov_state")

  
  n_owls <- df_full %>% 
    mutate(Survey = ifelse(dataset == "bbs",
                           "BBS",
                           "Nocturnal Owl Survey"),
           Survey = factor(Survey,
                           levels = rev(c("BBS","Nocturnal Owl Survey")),
                           ordered = TRUE)) %>% 
    ungroup() %>% 
    group_by(Survey) %>% 
    summarise(n_surveys = n(),
              n_sites = length(unique(route_id)),
              n_non_zeros = length(which(count > 0))) %>% 
    mutate(proportion_non_zero = round(n_non_zeros/n_surveys,2))
  
  ns_o <- n_owls %>% 
    filter(Survey == "Nocturnal Owl Survey") 
  ns_b <- n_owls %>% 
    filter(Survey == "BBS") 
  if(nrow(ns_b) == 0){
    ns_b <- data.frame(n_surveys = 0,
                       n_sites = 0,
                       proportion_non_zero = NA)
  }
  
  
  # 
  # obs_owls1 <- df_full %>% 
  #   mutate(Survey = ifelse(dataset == "bbs",
  #                          "BBS",
  #                          "Nocturnal Owl Survey"),
  #          Survey = factor(Survey,
  #                          levels = rev(c("BBS","Nocturnal Owl Survey")),
  #                          ordered = TRUE)) %>% 
  #   ungroup() %>% 
  #   group_by(Survey,strata_name,stratum) %>%
  #   summarise(max_obs = max(count))
  # 
  # 
  # obs_owls <- df_full %>% 
  #   mutate(Survey = ifelse(dataset == "bbs",
  #                          "BBS",
  #                          "Nocturnal Owl Survey"),
  #          Survey = factor(Survey,
  #                          levels = rev(c("BBS","Nocturnal Owl Survey")),
  #                          ordered = TRUE)) %>% 
  #   ungroup() %>% 
  #   group_by(Survey,strata_name,stratum,year) %>%
  #   summarise(obs_y = mean(count),
  #             n_surveys = n()) %>% 
  #   left_join(obs_owls1,
  #             by = c("Survey","strata_name","stratum")) %>% 
  #   mutate(p_obsy = obs_y/max_obs) %>% 
  #   ungroup() %>% 
  #   group_by(strata_name,year) %>% 
  #   summarise(mean_prop_obs = mean(p_obsy),
  #             n_datasources = n(),
  #             datasources = paste(unique(Survey),
  #                                 collapse = "-"))
  # 
  # 
  # 
  # 
  
  # obs_owls1 <- df_full %>% 
  #   left_join(alt_regs) %>% 
  #   mutate(Survey = ifelse(dataset == "bbs",
  #                          "BBS",
  #                          "Nocturnal Owl Survey"),
  #          Survey = factor(Survey,
  #                          levels = rev(c("BBS","Nocturnal Owl Survey")),
  #                          ordered = TRUE)) %>% 
  #   ungroup() %>% 
  #   group_by(Survey,survey_region) %>%
  #   summarise(max_obs = max(count))
  # 
  # obs_owls_region <- df_full %>% 
  #   left_join(alt_regs) %>% 
  #   mutate(Survey = ifelse(dataset == "bbs",
  #                          "BBS",
  #                          "Nocturnal Owl Survey"),
  #          Survey = factor(Survey,
  #                          levels = rev(c("BBS","Nocturnal Owl Survey")),
  #                          ordered = TRUE)) %>% 
  #   ungroup() %>% 
  #   group_by(Survey,survey_region,year) %>%
  #   summarise(obs_y = mean(count),
  #             n_surveys = n()) %>% 
  #   left_join(obs_owls1,
  #             by = c("Survey","survey_region")) %>% 
  #   mutate(p_obsy = obs_y/max_obs) %>% 
  #   ungroup() %>% 
  #   group_by(survey_region,year) %>% 
  #   summarise(mean_prop_obs = mean(p_obsy),
  #             n_datasources = n(),
  #             datasources = paste(unique(Survey),
  #                                 collapse = "-"))
  # 
  # 
  pdf(paste0("figures/new_fit_summary_",sp_ebird,".pdf"),
      width = 11,
      height = 8.5) 
  
  
  
  
  breaks <- c(-7, -4, -2, -1, -0.5, 0.5, 1, 2, 4, 7)
  labls <- c(paste0("< ", breaks[1]), paste0(breaks[-c(length(breaks))], 
                                             ":", breaks[-c(1)]), paste0("> ", breaks[length(breaks)]))
  labls <- paste0(labls, " %")
  
  
  pal <- stats::setNames(c("#a50026", "#d73027", "#f46d43", 
                           "#fdae61", "#fee090", "#ffffbf", "#e0f3f8", "#abd9e9", 
                           "#74add1", "#4575b4", "#313695"), labls)
  
  
  # pal <- stats::setNames(scico(11, palette = "roma"),
  #                        labls)
  # 
  
  
  # join strata with regions based on %overlap ------------------------------
  
  
  
  alt_regs_bbs_strata <- strata_used %>% 
    st_intersection(strata_bbs_regions) |> 
    st_make_valid() |> 
    mutate(area_km2 = as.numeric(st_area(geom)/1e6)) |> 
    filter(area_km2 > 10) |> 
    group_by(strata_name,stratum) |> 
    mutate(weights = area_km2/sum(area_km2))|> 
    st_drop_geometry() |> 
    select(strata_name,stratum,bbs_strata_name,weights) |> 
    rename(bbs_strata = bbs_strata_name)
  
  
  alt_regs_survey_regions <- strata_bbs_regions |> 
    select(bbs_strata_name,survey_region) |> 
    st_drop_geometry() |> 
    distinct()|> 
    rename(bbs_strata = bbs_strata_name)
  
  alt_regs_survey_wide <- strata_bbs_regions |> 
    select(bbs_strata_name) |> 
    st_drop_geometry() |> 
    mutate(survey_wide = "survey_wide") |> 
    distinct()|> 
    rename(bbs_strata = bbs_strata_name)
  
  # 
  # alt_regs <- strata_used %>%
  #   #st_intersection(provs) %>%
  #   st_join(provs,
  #           largest = TRUE,
  #           left = TRUE) |>
  #   st_make_valid() |>
  #   mutate(area_km2 = as.numeric(st_area(geom)/1e6)) |>
  #   filter(area_km2 > 10) |>
  #   group_by(strata_name,stratum, survey_region) |>
  #   mutate(weights = area_km2/sum(area_km2)) |>
  #   st_drop_geometry()
  # 
  reg_weights_all <- vector("list",
                                  length = 4)
  names(reg_weights_all) <- c("strata_level","bbs_strata","survey_wide","survey_region")
  
  reg_weights_all[["strata_level"]] <- prop_in_strata_used
  reg_weights_all[["bbs_strata"]] <- prop_in_strata_used
  reg_weights_all[["survey_wide"]] <- prop_in_bcr_prov
  reg_weights_all[["survey_region"]] <- prop_in_bcr_prov
  
  regions_index_all <- vector("list",
                              length = 4)
  names(regions_index_all) <- c("strata_level","bbs_strata","survey_wide","survey_region")
  regions_index_all[["strata_level"]] <- NA
  regions_index_all[["bbs_strata"]] <- alt_regs_bbs_strata
  regions_index_all[["survey_wide"]] <- alt_regs_survey_wide
  regions_index_all[["survey_region"]] <- alt_regs_survey_regions
  
  weeks_keep <- names(prop_in_strata_used)
  
  
 
  
  indices <- generate_indices_shared(model_fit = fit2,
                              meta_strata = meta_strata, # df with columns strata, strata_name, and optional weights
                              meta_years = meta_years,
                              raw_data = df_full,
                              quantiles = c(0.025, 0.05, 0.25, 0.75, 0.95, 0.975),
                              regions = c("strata_level","bbs_strata","survey_wide","survey_region"),
                              regional_sub = c("strata_level","strata_level","bbs_strata","bbs_strata"),# same length as regions indicates which regions form the subunits of each level
                              regions_index = regions_index_all, # alternate post-hoc combinations of strata to form new regions - list of data frames that link each value of the lower level in regional_sub to each value of the upper level of regions
                              alternate_n = "n",
                              gam_smooths = FALSE,
                              start_year = NULL,
                              max_backcast = NULL,
                              drop_exclude = FALSE,
                              hpdi = TRUE,
                              quiet = FALSE,
                              fractional_weights = alt_regs_bbs_strata,
                              rel_abundance_weights = reg_weights_all)
  
  
  indices_bbs <- generate_indices(model_fit = fit2,
                              meta_strata = meta_strata, # df with columns strata, strata_name, and optional weights
                              meta_years = meta_years,
                              raw_data = df_full,
                              quantiles = c(0.025, 0.05, 0.25, 0.75, 0.95, 0.975),
                              regions = c("bbs_strata_name"),
                              regions_index = alt_regs, # alternate post-hoc combinations of strata to form new regions - data frame with strata_name and region
                              alternate_n = "n",
                              gam_smooths = FALSE,
                              start_year = NULL,
                              max_backcast = NULL,
                              drop_exclude = FALSE,
                              hpdi = TRUE,
                              quiet = FALSE,
                              weighted = TRUE,
                              rel_abundance_weights = prop_in_strata_used)
  
  
  tt_all <- NULL
  
  earliest_trend <- min(df_full$year, na.rm = TRUE)
  
  
  
  yrs <- data.frame(st_year = c(earliest_trend,
                                earliest_trend,
                                2012,
                                2014),
                    en_year = c(2024,
                                2012,
                                2024,
                                2024))
  for(i in 1:nrow(yrs)){
    
    st_year <- yrs[i,"st_year"]
    en_year <- yrs[i,"en_year"]
    
    
    trends <- generate_trends(indices,
                              min_year = st_year,
                              max_year = en_year,
                              quantiles = c(0.025, 0.05, 0.1,0.25, 0.75,0.9, 0.95, 0.975),
                              slope = FALSE,
                              gam = FALSE,
                              prob_decrease = NULL,
                              prob_increase = NULL,
                              hpdi = TRUE)
    
    
    tt <- trends$trends %>% 
      mutate(trend_period = paste(start_year,end_year,sep = "-"))
    
    tt_all <- bind_rows(tt_all,tt)
  }
  
  tt_all <- tt_all %>% 
    mutate(trends = cut(trend, breaks = c(-Inf, 
                                          breaks, Inf), labels = labls),
           trends_lci = cut(trend_q_0.1, breaks = c(-Inf, 
                                                    breaks, Inf), labels = labls),
           trends_uci = cut(trend_q_0.9, breaks = c(-Inf, 
                                                    breaks, Inf), labels = labls))
  
  map_strat_trends <- strata_used %>% 
    inner_join(tt_all,
               by = c("strata_name" = "region"))
  
  bb <- st_bbox(map_strat_trends)
  
  t_map <- ggplot()+
    geom_sf(data = map_strat_trends,
            aes(fill = trends))+
    #scale_fill_viridis_c()+
    # scico::scale_fill_scico_d(direction = -1,
    #                           palette = "roma",
    #                           name = "Trend %/year")+
    scale_fill_manual(values = pal, 
                      na.value = "white",
                      name = "Trend %/year")+
    labs(title = paste(sp,"selected trends by strata with data"),
         subtitle = "Point estimates of the trends for each year by strata,
       showing the spatial pattern of annual changes.")+
    theme_bw()+
    geom_sf(data = provs, fill = NA)+
    coord_sf(xlim = bb[c("xmin","xmax")],
             ylim = bb[c("ymin","ymax")])+
    facet_wrap(vars(trend_period))
  
  
  print(t_map)
  
  
  
  
  t_map_l <- ggplot()+
    geom_sf(data = map_strat_trends,
            aes(fill = trends_lci))+
    #scale_fill_viridis_c()+
    # scico::scale_fill_scico_d(direction = -1,
    #                           palette = "roma",
    #                           name = "Trend %/year")+
    scale_fill_manual(values = pal, 
                      na.value = "white",
                      name = "Trend %/year")+
    labs(title = paste(sp,"selected trends by strata with data"),
         subtitle = "Lower 80% credible limit of the trends for each year by strata,
       showing the spatial pattern of annual changes.",
         caption = "The model implies an 90% probability that the trend is higher (more positive) than this value.")+
    theme_bw()+
    geom_sf(data = provs, fill = NA)+
    coord_sf(xlim = bb[c("xmin","xmax")],
             ylim = bb[c("ymin","ymax")])+
    facet_wrap(vars(trend_period))
  
  
  
  
  
  t_map_u <- ggplot()+
    geom_sf(data = map_strat_trends,
            aes(fill = trends_uci))+
    #scale_fill_viridis_c()+
    # scico::scale_fill_scico_d(direction = -1,
    #                           palette = "roma",
    #                           name = "Trend %/year")+
    scale_fill_manual(values = pal, 
                      na.value = "white",
                      name = "Trend %/year")+
    labs(title = paste(sp,"selected trends by strata with data"),
         subtitle = "Upper 80% credible limit of the trends for each year by strata,
       showing the spatial pattern of annual changes.",
         caption = "The model implies an 90% probability that the trend is lower (more negative) than this value.")+
    theme_bw()+
    geom_sf(data = provs, fill = NA)+
    coord_sf(xlim = bb[c("xmin","xmax")],
             ylim = bb[c("ymin","ymax")])+
    facet_wrap(vars(trend_period))
  
  map_uncert <- (t_map_u/t_map_l  + plot_layout(guides = "collect") )#& theme(legend.position = "bottom"))
  print(map_uncert)
  
  
  tt_all <- NULL
  
  yrs <- data.frame(st_year = c(seq(2000,2023)),
                    en_year = c(seq(2001,2024)))
  for(i in 1:nrow(yrs)){
    
    st_year <- yrs[i,"st_year"]
    en_year <- yrs[i,"en_year"]
    
    
    trends <- generate_trends(indices,
                              min_year = st_year,
                              max_year = en_year,
                              quantiles = c(0.025, 0.05, 0.1,0.25, 0.75,0.9, 0.95, 0.975),
                              slope = FALSE,
                              gam = FALSE,
                              prob_decrease = NULL,
                              prob_increase = NULL,
                              hpdi = TRUE)
    
    
    tt <- trends$trends %>% 
      mutate(trend_period = paste(start_year,end_year,sep = "-"))
    
    tt_all <- bind_rows(tt_all,tt)
  }
  
  tt_all <- tt_all %>% 
    mutate(trends = cut(trend, breaks = c(-Inf, 
                                          breaks, Inf), labels = labls))
  
  map_strat_trends <- strata_used %>% 
    inner_join(tt_all,
               by = c("strata_name" = "region"))
  
  t_map <- ggplot()+
    geom_sf(data = map_strat_trends,
            aes(fill = trends))+
    #scale_fill_viridis_c()+
    # scico::scale_fill_scico_d(direction = -1,
    #                           palette = "roma",
    #                           name = "Trend %/year")+
    scale_fill_manual(values = pal, 
                      na.value = "white",
                      name = "Trend %/year")+
    labs(title = paste(sp,"annual trends by strata with data"),
         subtitle = "Point estimates of the trends for each year by strata,
       showing the spatial pattern of annual changes.")+
    theme_bw()+
    geom_sf(data = provs, fill = NA)+
    coord_sf(xlim = bb[c("xmin","xmax")],
             ylim = bb[c("ymin","ymax")])+
    facet_wrap(vars(trend_period))
  
  
  print(t_map)
  
  
  
  trends <- generate_trends(indices,
                            # min_year = 2011,
                            # max_year = 2022,
                            quantiles = c(0.025, 0.05, 0.1,0.25, 0.75,0.9, 0.95, 0.975),
                            slope = FALSE,
                            gam = FALSE,
                            prob_decrease = NULL,
                            prob_increase = NULL,
                            hpdi = TRUE)
  
  
  tt <- trends$trends
  
  tt <- tt %>% 
    mutate(trends = cut(trend, breaks = c(-Inf, 
                                          breaks, Inf), labels = labls),
           species = sp)
  
  
  trends_out <- bind_rows(tt,trends_out)
  
  
  country <- rnaturalearth::ne_countries(continent = "North America") %>%
    filter(admin %in% c("Canada")) %>%
    sf::st_transform(crs = sf::st_crs(strata_bbs)) %>%
    rename(country = admin) %>%
    select(country)
  
  map_strat_trends <- strata_used %>% 
    sf::st_intersection(country) %>% 
    inner_join(tt,
               by = c("strata_name" = "region"))
  
  
  bb <- st_bbox(map_strat_trends)
  
  
  survey_sites <- df_full %>% 
    mutate(Survey = ifelse(dataset == "bbs",
                           "BBS",
                           "Nocturnal Owl Survey"),
           Survey = factor(Survey,
                           levels = rev(c("BBS","Nocturnal Owl Survey")),
                           ordered = TRUE)) %>% 
    select(route_id,Survey,longitude,latitude) %>% 
    distinct() %>% 
    st_as_sf(coords = c("longitude","latitude"),
             crs = 4326)
  
  
  t_map_long <- ggplot()+
    geom_sf(data = provs,
            fill = NA)+
    geom_sf(data = map_strat_trends,
            aes(fill = trends))+
    geom_sf(data = survey_sites,
            aes(colour = Survey),
            inherit.aes = FALSE,
            size = 0.05)+
    labs(title = paste(sp,"Trends across Canada, 2000-2024"),
         subtitle = paste0(ns_o$n_sites," Nocturnal Owl Survey routes with ", ns_o$n_surveys, " surveys and species observed during ",ns_o$proportion_non_zero*100,"% \n",
                           ns_b$n_sites, " BBS routes with ",ns_b$n_surveys," surveys and species observed during ",ns_b$proportion_non_zero*100,"%"),
         caption = paste("Population trends from an integrated analysis of data from Nocturnal Owl Surveys and BBS\n",
                         "using a spatially explicit model to estimate annual changes in counts within grid cells and\n",
                         "eBird seasonal relative abundance in 2023 to weight local changes based on relative population size."))+
    coord_sf(xlim = bb[c("xmin","xmax")],
             ylim = bb[c("ymin","ymax")])+
    #scale_fill_viridis_c()+
    scale_colour_viridis_d(direction = 1,
                           name = "Survey",
                           end = 0.5)+
    # scico::scale_fill_scico_d(direction = -1,
    #                           palette = "roma",
    #                           name = "Trend %/year")+
    scale_fill_manual(values = pal, 
                      na.value = "white",
                      name = "Trend %/year")+
    theme_bw()
  
  
  
  # pdf("temp.pdf")
  print(t_map_long)
  # dev.off()
  
  
  
  
  trends <- generate_trends(indices,
                            min_year = 2010,
                            # max_year = 2022,
                            quantiles = c(0.025, 0.05, 0.1,0.25, 0.75,0.9, 0.95, 0.975),
                            slope = FALSE,
                            gam = FALSE,
                            prob_decrease = NULL,
                            prob_increase = NULL,
                            hpdi = TRUE)
  
  
  tt <- trends$trends
  
  tt <- tt %>% 
    mutate(trends = cut(trend, breaks = c(-Inf, 
                                          breaks, Inf), labels = labls),
           species = sp)
  
  trends_out <- bind_rows(tt,trends_out)
  
  map_strat_trends <- strata_used %>% 
    sf::st_intersection(country) %>% 
    inner_join(tt,
               by = c("strata_name" = "region"))
  
  
  bb <- st_bbox(map_strat_trends)
  
  
  survey_sites <- df_full %>% 
    mutate(Survey = ifelse(dataset == "bbs",
                           "BBS",
                           "Nocturnal Owl Survey"),
           Survey = factor(Survey,
                           levels = rev(c("BBS","Nocturnal Owl Survey")),
                           ordered = TRUE)) %>% 
    select(route_id,Survey,longitude,latitude) %>% 
    distinct() %>% 
    st_as_sf(coords = c("longitude","latitude"),
             crs = 4326)
  
  
  t_map_short <- ggplot()+
    geom_sf(data = provs,
            fill = NA)+
    geom_sf(data = map_strat_trends,
            aes(fill = trends))+
    geom_sf(data = survey_sites,
            aes(colour = Survey),
            inherit.aes = FALSE,
            size = 0.05)+
    labs(title = paste(sp,"Trends across Canada, 2010-2024"),
         subtitle = paste0(ns_o$n_sites," Nocturnal Owl Survey routes with ", ns_o$n_surveys, " surveys and species observed during ",ns_o$proportion_non_zero*100,"% \n",
                           ns_b$n_sites, " BBS routes with ",ns_b$n_surveys," surveys and species observed during ",ns_b$proportion_non_zero*100,"%"),
         caption = paste("Population trends from an integrated analysis of data from Nocturnal Owl Surveys and BBS\n",
                         "using a spatially explicit model to estimate annual changes in counts within grid cells and\n",
                         "eBird seasonal relative abundance in 2023 to weight local changes based on relative population size."))+
    coord_sf(xlim = bb[c("xmin","xmax")],
             ylim = bb[c("ymin","ymax")])+
    #scale_fill_viridis_c()+
    scale_colour_viridis_d(direction = 1,
                           name = "Survey",
                           end = 0.7)+
    # scico::scale_fill_scico_d(direction = -1,
    #                           palette = "roma",
    #                           name = "Trend %/year")+
    scale_fill_manual(values = pal, 
                      na.value = "white",
                      name = "Trend %/year")+
    theme_bw()
  
  
  
  # pdf("temp.pdf")
  print(t_map_short)
  
  
  t_map <- ggplot()+
    geom_sf(data = provs,
            fill = NA)+
    geom_sf(data = map_strat_trends,
            fill = NA)+
    geom_sf(data = survey_sites,
            aes(colour = Survey),
            inherit.aes = FALSE,
            size = 0.2)+
    labs(title = paste0(sp," Survey routes contributing to trends"),
         subtitle = paste0(ns_o$n_sites," Nocturnal Owl Survey routes with ", ns_o$n_surveys, " surveys and species observed during ",ns_o$proportion_non_zero*100,"% \n",
                           ns_b$n_sites, " BBS routes with ",ns_b$n_surveys," surveys and species observed during ",ns_b$proportion_non_zero*100,"%"),
         caption = paste("Population trends from an integrated analysis of data from Nocturnal Owl Surveys and BBS\n",
                         "using a spatially explicit model to estimate annual changes in counts within grid cells and\n",
                         "eBird seasonal relative abundance in 2023 to weight local changes based on relative population size."))+
    coord_sf(xlim = bb[c("xmin","xmax")],
             ylim = bb[c("ymin","ymax")])+
    #scale_fill_viridis_c()+
    scale_colour_viridis_d(direction = 1,
                           name = "Survey",
                           end = 0.7)+
    theme_bw()
  
  
  
  # pdf("temp.pdf")
  print(t_map)
  
  # ch_map <- ggplot()+
  #   geom_sf(data = map_strat_trends,
  #           aes(fill = percent_change))+
  #   #scale_fill_viridis_c()+
  #   colorspace::scale_fill_binned_diverging(rev = TRUE,
  #                                           guide = "colorsteps",
  #                                               palette = "Blue-Red 3",
  #                                               breaks = c(-Inf,-50,-25,-10,10,33,100,Inf))
  # 
  # 
  # print(ch_map)
  
  # ch_map_up <- ggplot()+
  #   geom_sf(data = map_strat_trends,
  #           aes(fill = percent_change_q_0.75))+
  #   #scale_fill_viridis_c()+
  #   colorspace::scale_fill_continuous_diverging(rev = TRUE,
  #                                               palette = "Blue-Red 3")
  # 
  # 
  # ch_map_up
  # ch_map_do <- ggplot()+
  #   geom_sf(data = map_strat_trends,
  #           aes(fill = percent_change_q_0.25))+
  #   #scale_fill_viridis_c()+
  #   colorspace::scale_fill_continuous_diverging(rev = TRUE,
  #                                               palette = "Blue-Red 3")
  
  # 
  # print(ch_map_do)
  
  
  # mapping abundance -------------------------------------------------------
  
  
  # 
  # a_map <- ggplot()+
  #    geom_sf(data = provs,
  #            fill = NA)++
  #   geom_sf(data = map_strat_trends,
  #           aes(fill = rel_abundance))+
  #   scale_fill_viridis_c()+
  #    coord_sf(xlim = bb[c("xmin","xmax")],
  #             ylim = bb[c("ymin","ymax")])+
  #    labs(title = paste(sp,"Estimated relative abundance averaged "))
  # 
  # 
  # print(a_map)
  
  
  # trajectories ------------------------------------------------------------
  
  
  ii <- indices$indices 
  yups <- ii %>% 
    group_by(region) %>% 
    summarise(yup = max(index_q_0.75))
  
  ii_strat <- indices$indices #%>% 
    # filter(region_type == "strata_level") %>% 
    # left_join(obs_owls,
    #           by = c("region" = "strata_name",
    #                  "year")) %>% 
    # left_join(yups,
    #           by = "region") %>% 
    # mutate(mean_prop_obs_plot = mean_prop_obs*yup) %>% 
    # inner_join(alt_regs,by = c("region" = "strata_name"))
    # 
  nstrat <- nrow(yups)
  if(nstrat > 70){
    
    
    ii_test <- ggplot(data = filter(ii_strat,
                                    survey_region %in% c("Maritimes",
                                                         "Quebec",
                                                         "Ontario")),
                      aes(x = year, y = index))+
      geom_line()+
      geom_ribbon(aes(ymin = index_q_0.025,
                      ymax = index_q_0.975),
                  alpha = 0.2)+
      #geom_point(aes(y = mean_prop_obs_plot))+
      geom_point(aes(y = obs_mean, 
                     colour = n_routes, 
                     shape = datasources))+
      scale_colour_viridis_b(breaks = c(1,2,3,5,10,16))+
      facet_wrap(vars(region),
                 scales = "free_y")+
      theme(legend.position = "bottom",
            axis.text.x = element_text(angle = 90))+
      labs(subtitle = paste(sp,"Strata-level trajectories for strata in",
                            "Maritimes, ",
                            "Quebec, ",
                            "Ontario, "))
    
    
    print(ii_test)
    
    
    ii_test <- ggplot(data = filter(ii_strat,
                                    survey_region %in% c("Manitoba","Northwest Territories",
                                                         "Yukon",
                                                         "Saskatchewan",
                                                         "Alberta",
                                                         "British Columbia")),
                      aes(x = year, y = index))+
      geom_line()+
      geom_ribbon(aes(ymin = index_q_0.025,
                      ymax = index_q_0.975),
                  alpha = 0.2)+
      #geom_point(aes(y = mean_prop_obs_plot))+
      geom_point(aes(y = obs_mean, 
                     colour = n_routes, 
                     shape = datasources))+
      scale_colour_viridis_b(breaks = c(1,2,3,5,10,16))+
      facet_wrap(vars(region),
                 scales = "free_y")+
      theme(legend.position = "bottom",
            axis.text.x = element_text(angle = 90))+
      labs(subtitle = paste(sp,"Strata-level trajectories for strata in",
                            "Manitoba","Northwest Territories",
                            "Yukon",
                            "Saskatchewan",
                            "Alberta",
                            "British Columbia"))
    
    
    print(ii_test)
    
    
    
  }else{
    ii_test <- ggplot(data = ii_strat,
                      aes(x = year, y = index))+
      geom_line()+
      geom_ribbon(aes(ymin = index_q_0.025,
                      ymax = index_q_0.975),
                  alpha = 0.2)+
      #geom_point(aes(y = mean_prop_obs_plot))+
      geom_point(aes(y = obs_mean, 
                     colour = n_routes, 
                     shape = datasources))+
      scale_colour_viridis_b(breaks = c(1,2,3,5,10,16))+
      facet_wrap(vars(region),
                 scales = "free_y")+
      theme(legend.position = "bottom",
            axis.text.x = element_text(angle = 90))+
      labs(subtitle = paste(sp,"Strata-level trajectories"))
    
    
    print(ii_test)
  }
  
  indices_map <- strata_used %>% 
    inner_join(ii,
               by = c("strata_name" = "region"))
  
  
  
  indices_map_sel <- indices_map %>% 
    filter(year %in% seq(1995,2023,by = 5))
  
  
  a_map <- ggplot()+
    geom_sf(data = indices_map_sel,
            aes(fill = index))+
    scale_fill_viridis_c()+
    facet_wrap(vars(year))
  
  
  print(a_map)
  
  
  
  ii_surv <- ii %>% 
    filter(region == "survey_wide") 
  
  
  traj_nat <- ggplot(data = ii_surv,
                     aes(x = year, 
                         y = index))+
    geom_ribbon(aes(ymin = index_q_0.05,
                    ymax = index_q_0.95),
                alpha = 0.3)+
    coord_cartesian(ylim = c(0,NA))+
    scale_x_continuous(breaks = c(1995,
                                  2000,
                                  2005,
                                  2010,
                                  2015,
                                  2020,
                                  2025))+
    geom_line()+
    theme_bw()+
    ylab("Index of relative abundance")+
    xlab("")+
    labs(title = paste(sp,"Population trajectory for Canada"),
         caption = paste("Population trajectory from an integrated analysis of data from Nocturnal Owl Surveys and BBS\n",
                         "using a spatially explicit model to estimate annual changes in counts within grid cells and\n",
                         "eBird seasonal relative abundance in 2023 to weight local changes based on relative population size."))+
    theme(plot.caption = element_text(hjust = 0))
  
  
  print(traj_nat)
  
  
  
  ii_regs <- ii %>% 
    filter(region_type == "survey_region") #%>% 
    # left_join(obs_owls_region,
    #           by = c("region" = "survey_region",
    #                  "year")) %>% 
    # left_join(yups,
    #           by = "region") %>% 
    # mutate(mean_prop_obs_plot = mean_prop_obs*yup)
  
  traj_test <- ggplot(data = ii_regs,
                      aes(x = year, y = index))+
    geom_line()+
    geom_ribbon(aes(ymin = index_q_0.025,
                    ymax = index_q_0.975),
                alpha = 0.2)+
    #geom_point(aes(y = mean_prop_obs_plot))+
    geom_point(aes(y = obs_mean, colour = n_routes, shape = datasources))+
    scale_colour_viridis_c()+#b(breaks = c(1,2,3,5,10,16))+
    facet_wrap(vars(region),
               scales = "free_y")
  
  
  print(traj_test)
  
  
  traj <- ggplot(data = ii_regs,
                 aes(x = year, 
                     y = index))+
    geom_ribbon(aes(ymin = index_q_0.05,
                    ymax = index_q_0.95),
                alpha = 0.3)+
    geom_line()+
    scale_y_continuous(limits = c(0,NA))+
    facet_wrap(vars(region),
               scales = "free_y")
  
  
  print(traj)
  
  
  
  
  
  
  survey_explore <- df_full %>% 
    group_by(strata_name,dataset) %>% 
    summarise(n_surveys = n()) %>% 
    arrange(n_surveys) %>% 
    group_by(strata_name) %>% 
    sample_n(1)
  
  
  survey_expl_map <- strata_used %>% 
    left_join(survey_explore) 
  
  
  survey_map <- ggplot()+
    geom_sf(data = survey_expl_map,
            aes(fill = dataset))+
    scale_fill_viridis_d()
  print(survey_map)
  
  
  dev.off()
  
  
  png(paste0("Figures/",sp_ebird,"_Canada_trajectory.png"),
      res = 300, height = 5, width = 8,
      units = "in")
  print(traj_nat)
  dev.off()
  
  png(paste0("Figures/",sp_ebird,"_Canada_trend_long.png"),
      res = 300, height = 5, width = 8,
      units = "in")
  print(t_map_long)
  dev.off()
  
  png(paste0("Figures/",sp_ebird,"_Canada_trend_short.png"),
      res = 300, height = 5, width = 8,
      units = "in")
  print(t_map_short)
  dev.off()
  
} #end species loop


write_csv(trends_out,"All_owl_trends.csv")
write_csv(converge,"convergence_all_sp.csv")

fail <- converge %>% 
  filter(ess_fail | rhat_fail) %>% 
  group_by(species) %>% 
  summarise(n_fail = n())

fail


trends_out_broad <- trends_out %>% 
  filter(region_type != "strata_level")

write_csv(trends_out_broad,"All_broad_scale_owl_trends.csv")











