
library(bbsBayes2)
library(tidyverse)
library(cmdstanr)
library(sf)
library(ebirdst)
library(SurveyCoverage)

use_uncertainty <- FALSE
re_ebird <- FALSE
re_fit <- TRUE


source("functions/checks.R")
source("functions/generate_trends.R")
source("functions/generate_indices.R")
source("functions/utils.R")
source("functions/neighbours_define.R")

output_dir <- "output/"
# The code is here: https://github.com/BirdsCanada/National_NOS_Clean
# 
# The counts are summed per route for the target species list. 
# sp.list<-c("Great Horned Owl", "Boreal Owl", "Northern Saw-whet Owl", "Barred Owl", "Great Grey Owl")
# 
# Routes that had less than 5 stops were considered incomplete and removed. The data were also filtered for routes that were run under appropriate survey conditions (temp, wind, precip), as specified in the survey protocol. If a route was run more than one per year, I have retained it. 

# Load Owl Data -----------------------------------------------------------

events_owl <- read_csv("data/owls/samplingEvents.csv")%>% 
  distinct() %>% 
  mutate(unique_survey = paste(RouteIdentifier,survey_year,survey_month,survey_day, sep = "-")) %>% 
  arrange(unique_survey) %>% 
  filter(RouteIdentifier != "NS050") # temporary removal of one problematic route that is replicated in two provinces

# sel <- c(which(duplicated(events_owl$unique_survey))-1,
#          which(duplicated(events_owl$unique_survey)),
#          which(duplicated(events_owl$unique_survey))+1)
# tmp <- events_owl[sel,] %>% 
#   arrange(unique_survey)
# 
# tmp2 <- events_owl %>% 
#   filter(RouteIdentifier == "NS050")

protocols_owl <- events_owl %>% 
  group_by(collection,protocol_id) %>% 
  summarise(n_events = n())

protocols_owl

observers_owl <- events_owl %>% 
  group_by(RouteIdentifier,CollectorNumber) %>% 
  summarise(n_events = n())

hist(observers_owl$n_events)
length(which(observers_owl$n_events > 1))/nrow(observers_owl)

## 35% of surveys are single route-observer combinations
## may be very difficult to fit observer and route effects

routes_owl <- events_owl %>% 
  group_by(RouteIdentifier,protocol_id) %>% 
  summarise(n_events = n())

hist(routes_owl$n_events)
length(which(routes_owl$n_events > 1))/nrow(routes_owl)


# Mapping owl surveys --------------------------------------------------------

event_map <- events_owl %>% 
  group_by(RouteIdentifier,StateProvince,
           latitude,longitude,protocol_id) %>% 
  summarise(n_surveys = n(),
            min_year = min(survey_year,na.rm = TRUE),
            max_year = max(survey_year,na.rm = TRUE)) %>% 
  mutate(span_years = max_year - min_year) %>% 
  sf::st_as_sf(., coords = c("longitude","latitude")) #%>% 
sf::st_crs(event_map) <- 4269

tst_map <- ggplot()+
  geom_sf(data = event_map,
          aes(colour = span_years))

tst_map

event_map10 <- event_map %>% 
  filter(span_years > 9)

tst_map <- ggplot()+
  geom_sf(data = event_map10,
          aes(colour = span_years))+
  scale_colour_viridis_c()

tst_map

tst_map <- ggplot()+
  geom_sf(data = event_map,
          aes(colour = min_year))+
  scale_colour_viridis_c()

tst_map

progr_summary <- events_owl %>% 
  group_by(collection,protocol_id) %>% 
  summarise(first_year = min(survey_year),
            latest_year = max(survey_year),
            n_total = n())

# owl observations --------------------------------------------------------

obs_owl <- read_csv("data/owls/NationalOwlData.csv") %>% 
  filter(RouteIdentifier != "NS050") %>%  # temporary removal of one problematic route that is replicated in two provinces%>% 
  mutate(unique_survey = paste(RouteIdentifier,survey_year,survey_month,survey_day, sep = "-"))
  

# # identifying replicated survey events ---------------------------------
# 
# 
# tmp3 <- obs_owl %>% filter(RouteIdentifier == "NS050")
# 
#  write_csv(tmp2, "replicated_ATOWLS_NS050_two_province.csv")
#  write_csv(tmp3, "non-matching_obs_on_replicated_ATOWLS_NS050_two_province.csv")
# 
#  routes_GT1_years <- events_owl %>% 
#    group_by(RouteIdentifier) %>% 
#    summarise(n_years = length(unique(survey_year))) %>% 
#    filter(n_years > 1)
#  
#  events_owl <- events_owl %>% 
#    filter(RouteIdentifier != "NS050",
#           RouteIdentifier %in% routes_GT1_years$RouteIdentifier)

# step through species to build full dataset with zero-fill

# load all of the bbs data
full_bbs <- bbsBayes2::load_bbs_data()



# regional maps for coverage calculations ---------------------------------



# load maps of regions ----------------------------------------------------


stratum <- load_map("bbs_cws")

prov_state <- load_map("prov_state")  %>% 
  filter(country_code == "CA")%>%
  group_by(prov_state) %>%
  summarise() %>%
  mutate(strata_name = prov_state)


bcr_by_country <- load_map("bbs_cws") %>%
  filter(country_code == "CA")%>%
  group_by(bcr_by_country) %>%
  summarise() %>%
  rename(strata_name = bcr_by_country)

country <- rnaturalearth::ne_countries(continent = "North America") %>%
  filter(admin %in% c("Canada")) %>%
  sf::st_transform(crs = sf::st_crs(stratum)) %>%
  rename(strata_name = admin) %>%
  select(strata_name)

regs_maps <- list(country = country,
                  prov_state = prov_state,
                  stratum = stratum,
                  bcr_by_country = bcr_by_country)


# Species loop ------------------------------------------------------------


for(sp in c("Great Horned Owl","Barred Owl","Northern Saw-whet Owl","Boreal Owl")){
  
 #sp <- "Great Horned Owl"
  
  sp_obs_owl <- obs_owl %>% 
    filter(CommonName == sp) 
  
  sp_id <- unique(sp_obs_owl$species_id)
  
  sp_obs_owl <- sp_obs_owl %>% 
    select(-c(species_id,CommonName,doy))
  
  routes_w_species_owl <- sp_obs_owl %>% 
    group_by(RouteIdentifier) %>% 
    summarise(n_years = length(unique(survey_year))) %>% 
    filter(n_years > 1) ## ignore routes where species has been observed only once
  
  
  events_on_routes_w_sp_owl <- events_owl %>% 
    filter(RouteIdentifier %in% routes_w_species_owl$RouteIdentifier)
  
  obs_w_zeros_owl <- events_on_routes_w_sp_owl %>% 
    left_join(sp_obs_owl,
              by = c("collection", 
                     "RouteIdentifier", 
                     "StateProvince", 
                     "survey_year",
                     "survey_month", 
                     "survey_day", 
                     "protocol_id")) %>% 
    mutate(common_name = sp,
           species_id = sp_id,
           count = ifelse(is.na(Count),
                          0,
                          Count),
           dataset = "owl")
  
  check_owl_n_obs <- obs_w_zeros_owl %>% 
    filter(count > 0) %>% 
    group_by(RouteIdentifier) %>% 
    summarise(n_years = length(unique(survey_year))) %>% 
    filter(n_years > 1) ## ignore routes where species has been observed only once
  
  obs_w_zeros_owl <- obs_w_zeros_owl %>% 
    filter(RouteIdentifier %in% check_owl_n_obs$RouteIdentifier)
  
  
  
  







# BBS data ----------------------------------------------------------------


sp_aou <- as.integer(bbsBayes2::search_species(sp)["aou"])

obs_bbs <- full_bbs$birds %>% 
  filter(year >= min(obs_w_zeros_owl$survey_year),
         country_num == 124,
         aou == sp_aou) %>% # Canada only  
  mutate(route_id = paste(state_num,route,sep = "-")) %>% 
  select(route_id,
         year,
         route_data_id,
         species_total)

routes_w_species_bbs <- obs_bbs %>% 
  group_by(route_id) %>% 
  summarise(n_years = length(unique(year))) %>% 
  filter(n_years > 1) ## ignore routes where species has been observed only once

  
  
events_bbs <- full_bbs$routes %>% 
  mutate(route_id = paste(state_num,route,sep = "-")) %>% 
  filter(year >= min(obs_w_zeros_owl$survey_year),
         country_num == 124,
         route_id %in% routes_w_species_bbs$route_id) #Canada only  


obs_w_zeros_bbs <- events_bbs %>% 
  left_join(obs_bbs,
            by = c("route_id","route_data_id","year")) %>% 
  mutate(common_name = sp,
         species_id = sp_id,
         count = ifelse(is.na(species_total),
                        0,
                        species_total),
         dataset = "bbs")
  

check_bbs_n_obs <- obs_w_zeros_bbs %>% 
  filter(count > 0) %>% 
  group_by(route_id) %>% 
  summarise(n_years = length(unique(year))) %>% 
  filter(n_years > 1) ## ignore routes where species has been observed only once

obs_w_zeros_bbs <- obs_w_zeros_bbs %>% 
  filter(route_id %in% check_bbs_n_obs$route_id)

# Reconcile and combine two datasets --------------------------------------


# owl data
df_owl <- obs_w_zeros_owl %>% 
  select(RouteIdentifier,survey_year,survey_month,survey_day,
         dataset,protocol_id,
         nstop,latitude,longitude,
         common_name,species_id,
         count) %>% 
  mutate(route_id = paste(dataset,RouteIdentifier,sep = "-")) %>% 
  rename(year = survey_year,
         day = survey_day,
         month = survey_month) %>% 
  select(-RouteIdentifier)


#bbs data
df_bbs <- obs_w_zeros_bbs %>% 
  select(route_id,year,
         dataset,
         latitude,longitude,
         common_name,species_id,
         count) 

#combine into single dataframe
df_full <- df_owl %>% 
  bind_rows(df_bbs)




# stratification ----------------------------------------------------------

# benefit = can combine multiple surveys to share on local trends
# benefit = can drop strata-level intercept and only use route-level intercepts

strata_map <- bbsBayes2::load_map("latlong") %>% 
  select(-area_sq_km)


## Can treat NAD83 coordinates and WGS84 coordinates as equivalent, given the 
## differences between the two datum are ~1-2m, and therefore far smaller than
## the precision of the original coordinates and even less relevant given the 
## length of each route.
## 

route_coords <- df_full %>% 
  select(route_id,longitude,latitude) %>% 
  distinct() 

if(length(unique(route_coords$route_id)) != nrow(route_coords)){
  warning("Some routes have more than 1 set of coordinates")
  route_coords <- route_coords %>% 
    arrange(route_id) %>% 
    group_by(route_id) %>% 
    sample_n(1)
}
 route_coords <- route_coords %>% 
  sf::st_as_sf(., coords = c("longitude","latitude"),
               crs = 4269) %>% # NAD83
  sf::st_transform(crs = sf::st_crs(strata_map))

route_bb <- sf::st_bbox(route_coords)

tst <- ggplot()+
  geom_sf(data = strata_map)+
  geom_sf(data = route_coords)+
  coord_sf(xlim = route_bb[c("xmin","xmax")],
           ylim = route_bb[c("ymin","ymax")])
tst


route_by_strata <- route_coords %>% 
  sf::st_join(strata_map,
              left = TRUE,
              largest = TRUE) %>% 
  filter(!is.na(strata_name))

### this currently drops routes that don't fall within a stratum
## we could reasonably link all routes to the nearest stratum if we wanted to
## retain all of the observations


tst <- ggplot()+
  geom_sf(data = strata_map)+
  geom_sf(data = route_by_strata,
          aes(colour = strata_name))+
  coord_sf(xlim = route_bb[c("xmin","xmax")],
           ylim = route_bb[c("ymin","ymax")])+
  theme(legend.position = "none")
tst


route_strata_join <- route_by_strata %>% 
  sf::st_drop_geometry() 

df_full <- df_full %>% 
  inner_join(route_strata_join,
            by = "route_id") %>% 
  mutate(yr = year - (min(year)-1))
  


# filter out strata that don't meet minimum-span criterion ----------------

min_span <- 15 #minimum number of years with surveys

strata_span <- df_full %>% 
  group_by(strata_name) %>% 
  summarise(min_year = min(year),
            max_year = max(year)) %>% 
  mutate(span = max_year - min_year) %>% 
  filter(span >= min_span)

df_full <- df_full %>% 
  filter(strata_name %in% strata_span$strata_name) %>% 
  mutate(stratum = as.integer(factor(strata_name)))
 

route_strata_join <- df_full %>% 
  select(stratum,strata_name) %>% 
  distinct() %>% 
  arrange(stratum)




meta_years <- df_full %>% 
  select(year,yr) %>% 
  distinct() %>% 
  arrange(yr)

if(nrow(meta_years) != max(meta_years$yr)){
  stop("There is at least one year with no surveys.")
}

df_bbs_final <- df_full %>% 
  filter(dataset == "bbs") %>% 
  mutate(route = as.integer(factor(route_id))) 
  
  

df_owl_final <- df_full %>% 
  filter(dataset == "owl") %>% 
  mutate(route = as.integer(factor(route_id)),
         protocol = as.integer(factor(protocol_id)))

# prepare spatial components ----------------------------------------------
meta_strata <- route_strata_join %>% 
  select(strata_name,stratum) %>% 
  distinct()
  

strata_used <- strata_map %>% 
  inner_join(meta_strata,
             by = c("strata_name"))
  
route_by_strata <- route_by_strata %>% 
  filter(strata_name %in% df_full$strata_name)

tst <- ggplot()+
  geom_sf(data = strata_used)+
  geom_sf(data = route_by_strata,
          aes(colour = strata_name))+
  theme(legend.position = "none")
tst
 
neighbours <- neighbours_define(strata_used,
                                strat_link_fill = 1000, #distance to fill if strata are not connected
                                buffer = TRUE,
                                convex_hull = TRUE,
                                plot_neighbours = TRUE,
                                species = sp,
                                plot_dir = "neighbour_maps/",
                                plot_file = "_neighbour_map_vor_",
                                save_plot_data = TRUE,
                                voronoi = TRUE,
                                nn_fill = FALSE,
                                add_map = NULL,
                                strat_indicator = "stratum",
                                island_link_dist_factor = 1.2) #consider nearest strata neighbours if distances are within this factor of each other, when linking otherwise isolated islands of strata



# eBird relative abundance ------------------------------------------------

sp_ebird <- ebirdst::get_species(sp)
metric_used <- "mean"

qual_sel <- ebirdst_runs[which(ebirdst_runs$species_code == sp_ebird),]
breed_qual <- unname(unlist(qual_sel[,"breeding_quality"]))
resident_qual <- unname(unlist(qual_sel[,"resident_quality"]))
resident <- unname(unlist(qual_sel[,"is_resident"]))

breeding_start <- unname(unlist(qual_sel[,"breeding_start"]))
breeding_end <- unname(unlist(qual_sel[,"breeding_end"]))

season <- ifelse(resident,"resident","breeding")
yr_ebird <- as.integer(unname(unlist(qual_sel[,"status_version_year"])))

if(re_ebird){
  

down <- try(ebirdst::ebirdst_download_status(sp_ebird,
                                             download_ranges = TRUE,
                                             download_abundance = TRUE,
                                             download_occurrence = FALSE,
                                             force = FALSE,
                                             pattern = paste0("abundance_seasonal_",metric_used,"_3km_")),
            silent = TRUE)



if(use_uncertainty){
  down_lower <- try(ebirdst::ebirdst_download_status(sp_ebird,
                                                     download_ranges = FALSE,
                                                     download_abundance = TRUE,
                                                     download_occurrence = FALSE,
                                                     force = FALSE,
                                                     pattern = "abundance_lower_3km_"),
                    silent = TRUE)
  
  
  down_upper <- try(ebirdst::ebirdst_download_status(sp_ebird,
                                                     download_ranges = FALSE,
                                                     download_abundance = TRUE,
                                                     download_occurrence = FALSE,
                                                     force = FALSE,
                                                     pattern = "abundance_upper_3km_"),
                    silent = TRUE)
  
  
}
if(class(down) == "try-error"){
  
  warning(paste("Download fail", sp))

  next}


abd_seasonal_abundance <- ebirdst::load_raster(species = sp_ebird,
                                               resolution = "3km",
                                               period = "seasonal",
                                               product = "abundance",
                                               metric = metric_used)  #3km high resolution

breed_abundance <- abd_seasonal_abundance[[season]]

if(use_uncertainty){


abd_seasonal_abundance_lower <- ebirdst::load_raster(species = sp_ebird,
                                               resolution = "3km",
                                               period = "seasonal",
                                               product = "abundance",
                                               metric = "lower")  #3km high resolution

breed_abundance_lower <- abd_seasonal_abundance_lower[[season]]


abd_seasonal_abundance_upper <- ebirdst::load_raster(species = sp_ebird,
                                                   resolution = "3km",
                                                   period = "seasonal",
                                                   product = "abundance",
                                                   metric = metric_used)  #3km high resolution

breed_abundance_upper <- abd_seasonal_abundance_upper[[season]]

}

saveRDS(breed_abundance,paste0("data/species_relative_abundance_2023/",sp_ebird,"_derived_breeding_relative_abundance.rds"))

if(use_uncertainty){
  saveRDS(breed_abundance_upper,paste0("data/species_relative_abundance_2023/",sp_ebird,"_derived_breeding_relative_abundance_upper.rds"))
  saveRDS(breed_abundance_lower,paste0("data/species_relative_abundance_2023/",sp_ebird,"_derived_breeding_relative_abundance_lower.rds"))
}

# Calculate mean relative abundance in each stratum -----------------------

# both mean relative abundance in each of the full continental grid cell
# mean relative abundance of surveyed population in each grid-cell included in the analysis (grid cells with data)
# 
 
## crop relative abundance to US-Canada
bbs_strata_buffer <- bbsBayes2::load_map("bbs_usgs") %>%
  st_buffer(., dist = 10000) #add 10km buffer for clipping eBird data

# project boundary to match raster data
region_boundary_proj <- st_transform(bbs_strata_buffer, st_crs(breed_abundance))

## projected bbs strata
bbs_strata <- bbsBayes2::load_map("bbs_usgs")
bbs_strata_proj <- st_transform(bbs_strata, st_crs(breed_abundance))


breed_abundance <- terra::crop(breed_abundance, region_boundary_proj) |>
  terra::mask(region_boundary_proj)

if(use_uncertainty){
  
  breed_abundance_sd <- terra::crop(breed_abundance_sd, region_boundary_proj) |>
    terra::mask(region_boundary_proj)
  
  breed_abundance_ci <- terra::crop(breed_abundance_ci, region_boundary_proj) |>
    terra::mask(region_boundary_proj)
}


strata_used_proj <-  st_transform(strata_used, st_crs(breed_abundance))
  
strata_map_proj <-  st_transform(strata_map, st_crs(breed_abundance))

abundance_in_strata <- terra::extract(breed_abundance,
                                      strata_map_proj,
                                       fun = sum,
                                       na.rm = TRUE,
                                       ID = FALSE,
                                       exact = TRUE)

abundance_in_strata_used <- terra::extract(breed_abundance,
                                           strata_used_proj,
                                      fun = sum,
                                      na.rm = TRUE,
                                      ID = FALSE,
                                      exact = TRUE)

abundance_in_strata_mean <- terra::extract(breed_abundance,
                                      strata_map_proj,
                                      fun = mean,
                                      na.rm = TRUE,
                                      ID = FALSE,
                                      exact = TRUE)

abundance_in_strata_used_mean <- terra::extract(breed_abundance,
                                           strata_used_proj,
                                           fun = mean,
                                           na.rm = TRUE,
                                           ID = FALSE,
                                           exact = TRUE)



abundance_df <- data.frame(strata_name = strata_used_proj$strata_name,
                           stratum = strata_used_proj$stratum,
                           ebird_abund = abundance_in_strata_used[[1]],
                           ebird_abund_mean = abundance_in_strata_used_mean[[1]],
                            area_100km2 = units::drop_units(sf::st_area(strata_used_proj))/1e8) %>% 
  arrange(stratum) #

min_abund <- abundance_df %>%
  filter(ebird_abund > 0) %>% 
  summarise(min = min(ebird_abund),
            min_m = min(ebird_abund_mean))


abundance_df <- abundance_df %>% 
  mutate(ebird_abund = ifelse(ebird_abund == 0,
                              as.numeric(min_abund$min),
                              ebird_abund),
         ebird_abund_mean = ifelse(ebird_abund_mean == 0,
                              as.numeric(min_abund$min_m),
                              ebird_abund_mean),
         log_mean_rel_abund = log(ebird_abund/area_100km2))


abundance_df_all_strata <- data.frame(strata_name = strata_map_proj$strata_name,
                           ebird_abund = abundance_in_strata[[1]],
                           ebird_abund_mean = abundance_in_strata_mean[[1]],
                           area_100km2 = units::drop_units(sf::st_area(strata_map_proj))/1e8) 

min_abund <- abundance_df_all_strata %>%
  filter(ebird_abund > 0) %>% 
  summarise(min = min(ebird_abund),
            min_m = min(ebird_abund_mean))

abundance_df_all_strata <- abundance_df_all_strata %>% 
  mutate(ebird_abund = ifelse(ebird_abund == 0,
                              as.numeric(min_abund$min),
                              ebird_abund),
         ebird_abund_mean = ifelse(ebird_abund_mean == 0,
                                   as.numeric(min_abund$min_m),
                                   ebird_abund_mean),
         log_mean_rel_abund = log(ebird_abund/area_100km2))


saveRDS(abundance_df,paste0("data/",sp_ebird,"_strata_w_data_level_rel_abund.rds"))
saveRDS(abundance_df_all_strata,paste0("data/",sp_ebird,"_strata_level_rel_abund.rds"))
}else{
  abundance_df <- readRDS(paste0("data/",sp_ebird,"_strata_w_data_level_rel_abund.rds"))
abundance_df_all_strata <- readRDS(paste0("data/",sp_ebird,"_strata_level_rel_abund.rds"))
}

# Build Stan data object --------------------------------------------------

# df_bbs_final
# df_owl_final
# df_full

if(nrow(df_bbs_final) > 0){
train_bbs <- c(1:nrow(df_bbs_final))
count_bbs <- df_bbs_final$count
site_bbs <- df_bbs_final$route

year_bbs <- df_bbs_final$yr
strat_bbs <- df_bbs_final$stratum
}else{
  train_bbs <- 1
  count_bbs <- 1
  site_bbs <- 1
  year_bbs <- 1
  strat_bbs <- 1
}
stan_data <- list(n_sites_owl = max(df_owl_final$route),
                  n_sites_bbs = max(1,max(df_bbs_final$route)),
                  n_counts_bbs = max(1,nrow(df_bbs_final)),
                  n_counts_owl = nrow(df_owl_final),
                  n_strata = max(df_full$stratum),
                  n_years = max(df_full$yr),
                  n_protocols = max(df_owl_final$protocol),
                  ebird_year = yr_ebird-(min(df_full$year)-1),
                  yrev = seq(from = (yr_ebird-(min(df_full$year))),to = 1, by = -1),
                  
                  
                  count_bbs = count_bbs,
                  site_bbs = site_bbs,
                  year_bbs = year_bbs,
                  strat_bbs = strat_bbs,
                  
                  count_owl = df_owl_final$count,
                  site_owl = df_owl_final$route,
                  year_owl = df_owl_final$yr,
                  strat_owl = df_owl_final$stratum,
                  proto = df_owl_final$protocol,
                  off_set = log(df_owl_final$nstop), #effort offset for owl surveys
                  
                  n_train_bbs = ifelse(nrow(df_bbs_final) > 0,
                                       nrow(df_bbs_final),
                                       1),
                  n_train_owl = nrow(df_owl_final),
                  n_test_bbs = 1,
                  n_test_owl = 1,
                  train_bbs = train_bbs,
                  test_bbs = 1,
                  train_owl = 1:nrow(df_owl_final),
                  test_owl = 1,
                  
                  log_mean_rel_abund = abundance_df$log_mean_rel_abund,
                  
                  zero_betas = rep(0,max(abundance_df$stratum)),
                  
                  n_edges = neighbours$N_edges,
                  node1 = neighbours$node1,
                  node2 = neighbours$node2,
                  
                  # Conditional statments
                  use_pois = 0,
                  use_t = 1,
                  calc_cv = 0,
                  calc_nu = 0,
                  calc_log_lik = 0,
                  heavy_tailed = 0
)

for(j in names(stan_data)){
  if(!j %in% c("log_mean_rel_abund",
              "off_set")){
    stan_data[[j]] <- as.integer(stan_data[[j]])
  }
}


# coverage comparison -----------------------------------------------------

strata_cov <- strata_map %>%
  rename(grid_cell_name = strata_name) %>%  
  sf::st_intersection(sf::st_buffer(country,10000))%>% 
  filter(!is.na(strata_name))#

strat_area <- as.numeric(sf::st_area(strata_cov)/1e6)

strata_cov <- strata_cov %>% 
  mutate(area_km2 = strat_area) 

range_info <- grid_range(sp,
                             coverage_grid_custom = strata_cov,
                             seasonal_range = season)


survey_data <- df_full %>% 
  select(year,route_id,longitude,latitude)
  
survey_data_sf <- df_full %>% 
  select(year,route_id,longitude,latitude) %>% 
  st_as_sf(coords = c("longitude","latitude"),crs = 4326)

sp_coverage <- overlay_range_data(range = range_info,
                                  survey_sites = survey_data,
                                  sites = "route_id",
                                  years = "year",
                                  x_coord = "longitude",
                                  y_coord = "latitude",
                                  crs_site_coordinates = 4326,
                                  add_survey_sites_to_range = TRUE,
                                  proportion_area_included = 0.5)



cumulative_coverage_map <- sp_coverage$cumulative_coverage_map
overall_coverage_estimate <- sp_coverage$cumulative_coverage_estimate

coverage_overall <- ggplot()+
  geom_sf(data = cumulative_coverage_map,
          aes(fill = coverage))+
  geom_sf(data = survey_data_sf,
          alpha = 0.3,
          size = 0.3)+
  geom_sf(data = range_info$range_map,
          fill = NA)+
  scale_fill_viridis_d()+
  labs(title = paste(sp,"proportion covered = ",round(overall_coverage_estimate$coverage_proportion,2)))

pdf(paste0("figures/coverage_",sp_ebird,".pdf"))
print(coverage_overall)
dev.off()


ann_coverage <- NULL
cumulative_coverage <- NULL

for(reg in names(regs_maps)){
  
  mp_tmp <- regs_maps[[reg]]
  tmp_coverage <- regional_summary(sp_coverage,
                                   regions = mp_tmp,
                                   region_name = "strata_name")
  
  ann_tmp <- tmp_coverage$regional_annual_coverage_estimate %>%
    filter(coverage) %>%
    mutate(region_type = reg,
           species = sp)
  
  cumulative_tmp <- tmp_coverage$regional_cumulative_coverage_estimate %>%
    filter(coverage) %>%
    mutate(region_type = reg,
           species = sp)
  
  ann_coverage <- bind_rows(ann_coverage,ann_tmp)
  cumulative_coverage <- bind_rows(cumulative_coverage,cumulative_tmp)
  
} #end of regions loop





save(list = c("stan_data",
              "df_full",
              "meta_strata",
              "meta_years","strata_used",
              "df_owl_final",
              "df_bbs_final",
              "ann_coverage",
              "cumulative_coverage",
              "sp_coverage"),
     file = paste0("data/pre_fit_data_",sp_ebird,".RData"))
}# end data prep species loop
# Model fit ---------------------------------------------------------------


re_fit <- FALSE
for(sp in c("Great Horned Owl","Barred Owl","Northern Saw-whet Owl","Boreal Owl")[-1]){
 
  

  #sp_id <- unique(sp_obs_owl$species_id)
  sp_ebird <- ebirdst::get_species(sp)
 
model <- cmdstanr::cmdstan_model("models/first_difference_spatial_owls_integrated.stan")


load(paste0("data/pre_fit_data_",sp_ebird,".RData"))

if(re_fit){

fit2 <- model$sample(data = stan_data,
                    parallel_chains = 4,
                    refresh = 500,
                    iter_warmup = 1000,
                    iter_sampling = 1000,
                    adapt_delta = 0.8,
                    max_treedepth = 11,
                    show_exceptions = TRUE)

fit2$save_object(paste0(output_dir,"fit_","temp","_",sp_ebird,".rds"))

summ2 <- fit2$summary()

saveRDS(summ2, paste0(output_dir,"fit_summary_","temp","_",sp_ebird,".rds"))

}else{
  fit2 <- readRDS(paste0(output_dir,"fit_","temp","_",sp_ebird,".rds"))
  summ2 <- readRDS(paste0(output_dir,"fit_summary_","temp","_",sp_ebird,".rds"))
}



# 
# betas <- summ2 %>% filter(grepl("beta[",variable, fixed = TRUE)) %>% 
#   mutate(year = as.integer(str_extract(variable,"([[:digit:]]{1,})(?=]$)")),
#          stratum = as.integer(str_extract(variable,"(?<=\\[)[[:digit:]]{1,}")))
# 
# beta_by_year <- betas %>% 
#   group_by(year) %>% 
#   summarise(min_ess = min(ess_bulk),
#             max_rhat = max(rhat),
#             mean_ess = mean(ess_bulk),
#             mean_rhat = mean(rhat))
# 
# 
# beta_by_stratum <- betas %>% 
#   filter(!is.na(rhat)) %>% 
#   group_by(stratum) %>% 
#   summarise(min_ess = min(ess_bulk),
#             max_rhat = max(rhat),
#             mean_ess = mean(ess_bulk),
#             mean_rhat = mean(rhat))
# 
# 
# ns <- summ2 %>% filter(grepl("n[",variable,fixed = TRUE))%>% 
#   mutate(year = as.integer(str_extract(variable,"([[:digit:]]{1,})(?=]$)")),
#          stratum = as.integer(str_extract(variable,"(?<=\\[)[[:digit:]]{1,}")))
# 

protocols <- summ2 %>% filter(grepl("protocol[",variable, fixed = TRUE))


pdf(paste0("figures/temp_fit_summary_",sp_ebird,".pdf"),
    width = 11,
    height = 8.5) 




# explore results ---------------------------------------------------------
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


provs <- bbsBayes2::load_map("prov_state") %>% 
  filter(country_code == "CA") %>% 
  mutate(survey_region = ifelse(province_state %in% c("New Brunswick",
                                                      "Nova Scotia",
                                                      "Prince Edward Island"),
                                "Maritimes",
                                province_state)) %>% 
  select(province_state,survey_region)

alt_regs <- strata_used %>% 
  st_join(provs,
          left = TRUE,
          largest = TRUE) %>% 
  st_drop_geometry() %>% 
  select(strata_name,survey_region)
  

indices <- generate_indices(model_fit = fit2,
                           meta_strata = meta_strata, # df with columns strata, strata_name, and optional weights
                            meta_years = meta_years,
                           raw_data = df_full,
                            quantiles = c(0.025, 0.05, 0.25, 0.75, 0.95, 0.975),
                            regions = c("survey_wide","strata_level","survey_region"),
                            regions_index = alt_regs, # alternate post-hoc combinations of strata to form new regions - data frame with strata_name and region
                            alternate_n = "n",
                            gam_smooths = FALSE,
                            start_year = NULL,
                            max_backcast = NULL,
                            drop_exclude = FALSE,
                            hpdi = TRUE,
                            quiet = FALSE,
                           weighted = TRUE)


tt_all <- NULL

yrs <- data.frame(st_year = c(1995,
                                1995,
                                2005,
                                2015,
                                2011),
                  en_year = c(2024,
                              2005,
                              2015,
                              2024,
                              2024))
for(i in 1:nrow(yrs)){
  
  st_year <- yrs[i,"st_year"]
  en_year <- yrs[i,"en_year"]
  

trends <- generate_trends(indices,
                          min_year = st_year,
                          max_year = en_year,
                          quantiles = c(0.025, 0.05, 0.25, 0.75, 0.95, 0.975),
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
  theme_bw()+
  facet_wrap(vars(trend_period))


print(t_map)



tt_all <- NULL

yrs <- data.frame(st_year = c(seq(1995,2023)),
                  en_year = c(seq(1996,2024)))
for(i in 1:nrow(yrs)){
  
  st_year <- yrs[i,"st_year"]
  en_year <- yrs[i,"en_year"]
  
  
  trends <- generate_trends(indices,
                            min_year = st_year,
                            max_year = en_year,
                            quantiles = c(0.025, 0.05, 0.25, 0.75, 0.95, 0.975),
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
  theme_bw()+
  facet_wrap(vars(trend_period))


print(t_map)



trends <- generate_trends(indices,
                          # min_year = 2011,
                          # max_year = 2022,
                          quantiles = c(0.025, 0.05, 0.25, 0.75, 0.95, 0.975),
                          slope = FALSE,
                          gam = FALSE,
                          prob_decrease = NULL,
                          prob_increase = NULL,
                          hpdi = FALSE)


tt <- trends$trends

tt <- tt %>% 
  mutate(trends = cut(trend, breaks = c(-Inf, 
                                        breaks, Inf), labels = labls))

country <- rnaturalearth::ne_countries(continent = "North America") %>%
  filter(admin %in% c("Canada")) %>%
  sf::st_transform(crs = sf::st_crs(stratum)) %>%
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


t_map <- ggplot()+
  geom_sf(data = provs,
          fill = NA)+
  geom_sf(data = map_strat_trends,
          aes(fill = trends))+
  geom_sf(data = survey_sites,
          aes(colour = Survey),
          inherit.aes = FALSE,
          size = 0.05)+
  labs(title = paste(sp,"Trends across Canada, 1995-2024"),
       caption = paste("Population trends from an integrated analysis of Nocturnal Owl Monitoring data and BBS"))+
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
 print(t_map)
# dev.off()

 
 
 
 trends <- generate_trends(indices,
                            min_year = 2010,
                           # max_year = 2022,
                           quantiles = c(0.025, 0.05, 0.25, 0.75, 0.95, 0.975),
                           slope = FALSE,
                           gam = FALSE,
                           prob_decrease = NULL,
                           prob_increase = NULL,
                           hpdi = FALSE)
 
 
 tt <- trends$trends
 
 tt <- tt %>% 
   mutate(trends = cut(trend, breaks = c(-Inf, 
                                         breaks, Inf), labels = labls))
 
 
 
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
 
 
 t_map <- ggplot()+
   geom_sf(data = provs,
           fill = NA)+
   geom_sf(data = map_strat_trends,
           aes(fill = trends))+
   geom_sf(data = survey_sites,
           aes(colour = Survey),
           inherit.aes = FALSE,
           size = 0.05)+
   labs(title = paste(sp,"Trends across Canada, 2010-2024"),
        caption = paste("Population trends from an integrated analysis of Nocturnal Owl Monitoring data and BBS"))+
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



a_map <- ggplot()+
  geom_sf(data = map_strat_trends,
          aes(fill = rel_abundance))+
  scale_fill_viridis_c()


print(a_map)


ii <- indices$indices 

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


traj <- ggplot(data = ii_surv,
               aes(x = year, 
                   y = index))+
  #geom_point(aes(x = year, y = obs_mean))+
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
             caption = paste("Population trajectory from an integrated analysis of Nocturnal Owl Monitoring data and BBS"))


print(traj)



ii_regs <- ii %>% 
  filter(region_type == "survey_region")



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





} #end species loop

