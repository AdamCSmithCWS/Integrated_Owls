
library(bbsBayes2)
library(tidyverse)
library(cmdstanr)
library(sf)
library(ebirdst)
library(SurveyCoverage)

use_uncertainty <- TRUE
re_ebird <- TRUE #should it re-download eBird
# re_fit <- TRUE


source("functions/neighbours_define.R")
source("functions/which_min_LT.R")

output_dir <- "output/"
# The code is here: https://github.com/BirdsCanada/National_NOS_Clean
# 
# The counts are summed per route for the target species list. 
# sp.list<-c("Great Horned Owl", "Boreal Owl", "Northern Saw-whet Owl", "Barred Owl", "Great Grey Owl")
# 
# Routes that had less than 5 stops were considered incomplete and removed. 
# The data were also filtered for routes that were run under appropriate survey conditions
#  (temp, wind, precip), as specified in the survey protocol. 
#  If a route was run more than once per year, I have retained it. 

 year_end <- 2025
 year_start <- 2000
 # Load Owl Data -----------------------------------------------------------

events_owl <- read_csv("data/owls2/samplingEvents.csv")%>% 
  distinct() %>% 
  mutate(unique_survey = paste(RouteIdentifier,survey_year,survey_month,survey_day, sep = "-")) %>% 
  arrange(unique_survey) %>% 
  filter(RouteIdentifier != "NS050", # temporary removal of one problematic route that is replicated in two provinces
         survey_year <= year_end, # exclude 2026
         survey_year >= year_start) #only since 2000

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

obs_owl <- read_csv("data/owls2/NationalOwlData.csv") %>% 
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
full_bbs <- bbsBayes2::load_bbs_data(release = 2025)



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


regen_map <- FALSE
if(regen_map){
event_map2 <- event_map %>% 
  #filter(min_year < 2016) %>% 
  st_transform(crs = st_crs(prov_state))

bb <- st_bbox(event_map2)

survey_map <- ggplot()+
  geom_sf(data = prov_state,
          fill = NA)+
  geom_sf(data = event_map2,
          aes(colour = min_year),
          alpha = 0.3)+
  scale_colour_viridis_b(direction = -1,
                         name = "First year with\ndata in analyses",
                         breaks = c(2000,2005,2010,2015,2020,2025),
                         right = FALSE)+
  labs(title = paste0("Nocturnal Owl Survey routes"),
       subtitle = paste(nrow(events_owl), "surveys, at",length(unique(event_map$RouteIdentifier)),"routes\n",
                         "by",length(unique(events_owl$CollectorNumber)),"volunteer observers"))+
  coord_sf(xlim = bb[c("xmin","xmax")],
           ylim = bb[c("ymin","ymax")])+
  theme_bw()
  

png(paste0("Figures/Nocturnal_owl_survey_routes.png"),
    res = 300, height = 5, width = 8,
    units = "in")
print(survey_map)
dev.off()
}

owl_species <- obs_owl %>% 
  group_by(CommonName) %>% 
  summarise(n_years = length(unique(survey_year)),
            n_routes = length(unique(RouteIdentifier)),
            n_obs = length(which(Count > 0))) %>% 
  arrange(n_obs) %>% 
  filter(n_years > 24) # species must have observations for 25 years (2001 through 2025)


owls_ebird <- ebirdst::get_species(owl_species$CommonName)


ebird_owls_status <- ebirdst_runs[which(ebirdst_runs$species_code %in% owls_ebird),]






# Species data prep loop ------------------------------------------------------------


for(sp in owl_species$CommonName){
  
 #sp <- "Barred Owl"
  
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
  
  
  
  
  obs_check <- obs_w_zeros_owl %>% 
    group_by(survey_year,StateProvince) %>% 
    summarise(n_pos = length(which(count > 0)),
              n_surveys = n())



  if(sp %in% c("American Woodcock",
               "Wilson's Snipe",
               "Ruffed Grouse")){
    obs_w_zeros_owl <- obs_w_zeros_owl %>% 
      filter(survey_year > 2000,
             !(survey_year < 2022 & StateProvince == "BCY"))
    
  }
  



# BBS data ----------------------------------------------------------------


sp_aou <- as.integer(bbsBayes2::search_species(sp)["aou"])

obs_bbs <- full_bbs$birds %>% 
  filter(year >= min(obs_w_zeros_owl$survey_year), # only include BBS data from years with owl data
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
         dataset = "bbs",
         prov_route = st_abrev)
  

check_bbs_n_obs <- obs_w_zeros_bbs %>% 
  filter(count > 0) %>% 
  group_by(route_id) %>% 
  summarise(n_years = length(unique(year)))  %>% 
  filter(n_years > 1) ## ignore routes where species has been observed only once

obs_w_zeros_bbs <- obs_w_zeros_bbs %>% 
  filter(route_id %in% check_bbs_n_obs$route_id)

# Reconcile and combine two datasets --------------------------------------


# owl data
df_owl <- obs_w_zeros_owl %>% 
  select(RouteIdentifier,survey_year,survey_month,survey_day,
         dataset,protocol_id,
         nstop,latitude,longitude,StateProvince,
         common_name,species_id,
         count) %>% 
  mutate(route_id = paste(dataset,RouteIdentifier,sep = "-"),
         prov_route = ifelse(StateProvince == "BCY",
                             "BC",StateProvince)) %>% 
  rename(year = survey_year,
         day = survey_day,
         month = survey_month) %>% 
  select(-RouteIdentifier)


#bbs data
df_bbs <- obs_w_zeros_bbs %>% 
  select(route_id,year,
         day,
         month,
         dataset,
         latitude,longitude,prov_route,
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
## the precision of the original coordinates, or the area surveyed during a 
## stop, and even less relevant given the length of each route.
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
              largest = TRUE) 

missed_routes <- route_by_strata |> 
  filter(is.na(strata_name)) |> 
  select(-strata_name)

if(nrow(missed_routes) > 0){
w_miss_join <- which_min_LT(outside = missed_routes,
                            strata_map,
                            distance_to_strata = 5000)
route_by_strata <- route_by_strata %>%
  dplyr::filter(!is.na(.data$strata_name)) %>%
  dplyr::bind_rows(w_miss_join)

}
### If a route does not fall directly within a stratum, the if loop above 
### currently links routes to the nearest strata based on their start locations being
### within 5km of the strata, 
bbs_strata <- bbsBayes2::load_map("bbs") |> 
  filter(country_code == "CA",
         area_sq_km > 50) # removes the sliver of BCR 23 that seems to have been allocated to Ontario


route_by_bbs_strata <- route_coords %>% 
  sf::st_join(bbs_strata,
              left = TRUE,
              largest = TRUE) 

missed_routes2 <- route_by_bbs_strata |> 
  filter(is.na(strata_name)) |> 
  select(route_id)
if(nrow(missed_routes2) > 0){
  w_miss_join2 <- which_min_LT(outside = missed_routes2,
                              bbs_strata,
                              distance_to_strata = 5000)
  route_by_bbs_strata <- route_by_bbs_strata %>%
    dplyr::filter(!is.na(.data$strata_name)) %>%
    dplyr::bind_rows(w_miss_join2)
  
}

route_by_bbs_strata <- route_by_bbs_strata |> 
  rename(bbs_strata = strata_name) |> 
  select(route_id,bbs_strata) |> 
  st_drop_geometry()


  


route_strata_join <- route_by_strata %>% 
  sf::st_drop_geometry()  |> 
  left_join(route_by_bbs_strata,
            by = "route_id")

df_full <- df_full %>% 
  inner_join(route_strata_join,
            by = "route_id") %>% 
  mutate(yr = year - (min(year)-1))
  
tst_strata_allocation <- df_full |> 
  group_by(strata_name,bbs_strata,prov_route,dataset) |> 
  summarise(n = n(),
            n_route = length(unique(route_id))) |> 
  mutate(prov_strat = str_extract(bbs_strata,pattern = "(?<=-)[[:alpha:]]{2}"),
         prov_check = ifelse(prov_strat == prov_route,TRUE,FALSE))

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
  warning("There is at least one year with no surveys.")
  next
}

df_bbs_final <- df_full %>% 
  filter(dataset == "bbs") %>% 
  mutate(route = as.integer(factor(route_id))) 
  
  

df_owl_final <- df_full %>% 
  filter(dataset == "owl") %>% 
  mutate(route = as.integer(factor(route_id)),
         protocol = as.integer(factor(protocol_id)))
# 
# 
# tst <- ggplot(data = df_owl_final,
#               aes(group = factor(protocol),
#                   colour = factor(protocol),
#                   x = nstop, y = log(count+1,10)))+
#   geom_point(alpha = 0.2,position = position_jitter(width = 0.1, height = 0.1))+
#   geom_smooth(method = "lm",se = FALSE)+
#   scale_x_continuous(transform = "log10")+
#   scale_colour_viridis_d()+
#   facet_wrap(vars(protocol_id))
#         
# 
# tst        
#   
# 
# rt_var <- df_owl_final %>% 
#   group_by(route_id) %>% 
#   summarise(min_nstop = min(nstop),
#             max_nstop = max(nstop),
#             mean_nstop = mean(nstop),
#             sd_nstop = sd(nstop),
#             nsurveys = n(),
#             nyears = length(unique(year)),
#             mean_count = mean(count),
#             sd_count = sd(count))


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
                                             download_ranges = FALSE,
                                             download_abundance = TRUE,
                                             download_occurrence = FALSE,
                                             force = FALSE,
                                             pattern = "abundance_(median|upper|lower)_3km"),
            silent = TRUE)

}

# if(use_uncertainty){
#   down_lower <- try(ebirdst::ebirdst_download_status(sp_ebird,
#                                                      download_ranges = FALSE,
#                                                      download_abundance = TRUE,
#                                                      download_occurrence = FALSE,
#                                                      force = FALSE,
#                                                      pattern = "abundance_lower_3km_"),
#                     silent = TRUE)
#   
#   
#   down_upper <- try(ebirdst::ebirdst_download_status(sp_ebird,
#                                                      download_ranges = FALSE,
#                                                      download_abundance = TRUE,
#                                                      download_occurrence = FALSE,
#                                                      force = FALSE,
#                                                      pattern = "abundance_upper_3km_"),
#                     silent = TRUE)
#   
#   
# }
if(class(down) == "try-error"){
  
  warning(paste("Download fail", sp))

  next}

### insert code to pull weekly abundance for the relevant season
### 
df_full <- df_full %>% 
  mutate(doy = yday(as_date(paste(year,month,day,sep = "-"))))
  
  
survey_days <- df_full %>% 
  select(month,day) %>% 
  distinct() %>% 
  mutate(year = yr_ebird,
         doy = as_date(paste(year,month,day,sep = "-"))) %>% 
  drop_na() ## drops routes without day of survey information (or with incorrect information, e.g., April 31)

survey_day_range <- as_date(min(survey_days$doy):max(survey_days$doy))


abd_seasonal_abundance <- ebirdst::load_raster(species = sp_ebird,
                                               resolution = "3km",
                                               #period = "seasonal",
                                               product = "abundance",
                                               metric = "median")  #3km high resolution

weeks_keep <- names(abd_seasonal_abundance)
weeks_keep <- weeks_keep[weeks_keep %in% survey_day_range]

## keep only estimates for weeks with survey data
breed_abundance <- abd_seasonal_abundance[[weeks_keep]]

## translate the weeks into breaks for the doy values
doy_weeks_keep <- c((yday(weeks_keep)-3),(yday(weeks_keep)[length(weeks_keep)]+4))




## crop relative abundance to US-Canada
bbs_strata_buffer <- bbsBayes2::load_map("bbs") %>%
  st_buffer(., dist = 10000) #add 10km buffer for clipping eBird data

# project boundary to match raster data
region_boundary_proj <- st_transform(bbs_strata_buffer, st_crs(breed_abundance))

# bbs_strata_o <- bbs_strata

bbs_strata <- bbs_strata |>
#   select(-area_sq_km) |>
  rename(bbs_strata = strata_name) 
#   st_intersection(strata_map) |> 
#   mutate(area_km2 = as.numeric(st_area(geom))/1e6) |> 
#   filter(area_km2 > 10) |> 
#   sf::st_make_valid()


## projected bbs strata
bbs_strata_proj <- st_transform(bbs_strata, st_crs(breed_abundance))%>% 
  st_make_valid()


strata_used_proj <-  st_transform(strata_used, st_crs(breed_abundance)) %>% 
  st_make_valid()

strata_map_proj <-  st_transform(strata_map, st_crs(breed_abundance)) %>% 
  st_make_valid()


replace_zeros <- function(x){
  min_x <- min(x[which(x>0)],na.rm = TRUE)
  y <- ifelse(x == 0,min_x*0.5, # currently forces a non-zero abundance value that is half of the lowest positive abundance
              x)
  y <- ifelse(is.na(y),
              min_x*0.1, # also forces a non-zero abundance value that is 1/10th of the lowest positive abundance if missing, but also has monitoring observations
              y)
  return(y)
} 


abundance_in_strata_used_wide <- terra::extract(breed_abundance,
                                                strata_used_proj,
                                                fun = sum,
                                                na.rm = TRUE,
                                                ID = FALSE,
                                                exact = TRUE) %>% 
  mutate(strata_name = strata_used_proj$strata_name,
         stratum = strata_used_proj$stratum,
         across(.cols = starts_with(as.character(yr_ebird)),
                .fns = ~replace_zeros(.x))) # forces a non-zero abundance value for regions with monitoring data

if(any(!is.finite(colSums(abundance_in_strata_used_wide[,weeks_keep])))){
  weeks_drop <- which(!is.finite(colSums(abundance_in_strata_used_wide[,weeks_keep])))
  
  weeks_keep <- weeks_keep[-weeks_drop]
  
  # doy_keep <- min(yday(as_date(weeks_keep))):max(yday(as_date(weeks_keep)))
  # df_full <- df_full |> 
  #   filter(doy %in% doy_keep)
  # 
  
  ## keep only estimates for weeks with survey data
    breed_abundance <- abd_seasonal_abundance[[weeks_keep]]
  
  ## translate the weeks into breaks for the doy values
  doy_weeks_keep <- c((yday(weeks_keep)-3),(yday(weeks_keep)[length(weeks_keep)]+4))
  
  
  abundance_in_strata_used_wide <- terra::extract(breed_abundance,
                                                  strata_used_proj,
                                                  fun = sum,
                                                  na.rm = TRUE,
                                                  ID = FALSE,
                                                  exact = TRUE) %>% 
    mutate(strata_name = strata_used_proj$strata_name,
           stratum = strata_used_proj$stratum,
           across(.cols = starts_with(as.character(yr_ebird)),
                  .fns = ~replace_zeros(.x))) # forces a non-zero abundance value for regions with monitoring data
  
  
}


saveRDS(breed_abundance,paste0("data/species_relative_abundance_",yr_ebird,"/",sp_ebird,"_derived_breeding_weekly_relative_abundance.rds"))


 
  # abundance_in_strata_used_wide <- terra::extract(breed_abundance,
  #                                            strata_used_proj,
  #                                            fun = sum,
  #                                            na.rm = TRUE,
  #                                            ID = FALSE,
  #                                            exact = TRUE) %>% 
  #   mutate(strata_name = strata_used_proj$strata_name,
  #          stratum = strata_used_proj$stratum,
  #          across(.cols = starts_with(as.character(yr_ebird)),
  #                 .fns = ~replace_zeros(.x))) # forces a non-zero abundance value for regions with monitoring data
  # 

  
  prop_in_strata_used <- abundance_in_strata_used_wide %>% 
    mutate(across(.cols = starts_with(as.character(yr_ebird)), # calculates the proportion of the total abundance in each polygon
                  .fns = ~.x/sum(.x,na.rm = TRUE))) %>% 
    arrange(stratum)
  
  abundance_in_strata_used_summary <- abundance_in_strata_used_wide %>% 
    pivot_longer(cols = starts_with(as.character(yr_ebird)),
                 values_to = "rel_abundance",
                 names_to = "week") %>% 
    group_by(week) %>% 
    mutate(prop_population = rel_abundance/sum(rel_abundance,na.rm = TRUE)) %>% 
    ungroup() %>% 
    group_by(strata_name,stratum) %>% 
    summarise(mean_rel_abundance = mean(rel_abundance),
              median_rel_abundance = median(rel_abundance),
              min_rel_abundance = min(rel_abundance),
              max_rel_abundance = max(rel_abundance),
              sd_rel_abundance = sd(rel_abundance),
              mean_prop_population = mean(prop_population),
              median_prop_population = median(prop_population),
              min_prop_population = min(prop_population),
              max_prop_population = max(prop_population),
              sd_prop_population = sd(prop_population)) %>% 
    mutate(cv_rel_abundance = sd_rel_abundance/median_rel_abundance,
           cv_prop_population = sd_prop_population/mean_prop_population)
  
  

  

# Abundance in BCR-prov ---------------------------------------------------


  abundance_in_bcr_prov_wide <- terra::extract(breed_abundance,
                                                  bbs_strata_proj,
                                                  fun = sum,
                                                  na.rm = TRUE,
                                                  ID = FALSE,
                                                  exact = TRUE) %>% 
    mutate(bbs_strata = bbs_strata$bbs_strata,
           country = bbs_strata_proj$country,
           prov_state = bbs_strata_proj$prov_state,
           bcr = bbs_strata_proj$bcr,
           across(.cols = starts_with(as.character(yr_ebird)),
                  .fns = ~replace_zeros(.x))) # forces a non-zero abundance value for regions with monitoring data
  
  
  prop_in_bcr_prov <- abundance_in_bcr_prov_wide %>% 
    mutate(across(.cols = starts_with(as.character(yr_ebird)),
                  .fns = ~.x/sum(.x,na.rm = TRUE))) # calculates the proportion of the weekly abundance in each cell
  
  abundance_in_bcr_prov_summary <- abundance_in_bcr_prov_wide %>% 
    pivot_longer(cols = starts_with(as.character(yr_ebird)),
                 values_to = "rel_abundance",
                 names_to = "week") %>% 
    group_by(week) %>% 
    mutate(prop_population = rel_abundance/sum(rel_abundance,na.rm = TRUE)) %>% 
    ungroup() %>% 
    group_by(bbs_strata, 
             country,
             prov_state,
             bcr) %>% 
    summarise(mean_rel_abundance = mean(rel_abundance),
              median_rel_abundance = median(rel_abundance),
              min_rel_abundance = min(rel_abundance),
              max_rel_abundance = max(rel_abundance),
              sd_rel_abundance = sd(rel_abundance),
              mean_prop_population = mean(prop_population),
              median_prop_population = median(prop_population),
              min_prop_population = min(prop_population),
              max_prop_population = max(prop_population),
              sd_prop_population = sd(prop_population)) %>% 
    mutate(cv_rel_abundance = sd_rel_abundance/median_rel_abundance,
           cv_prop_population = sd_prop_population/mean_prop_population)
  
  
  
    
# Abundance in all strata -------------------------------------------------

  
  
  abundance_in_strata_wide <- terra::extract(breed_abundance,
                                             strata_map_proj,
                                                  fun = sum,
                                                  na.rm = TRUE,
                                                  ID = FALSE,
                                                  exact = TRUE) %>% 
    mutate(strata_name = strata_map_proj$strata_name,
           across(.cols = starts_with(as.character(yr_ebird)),
                  .fns = ~replace_zeros(.x))) # forces a non-zero abundance value for regions with monitoring data
  
  
  prop_in_strata <- abundance_in_strata_wide %>% 
    mutate(across(.cols = starts_with(as.character(yr_ebird)),
                  .fns = ~.x/sum(.x,na.rm = TRUE))) 
  
  abundance_in_strata_summary <- abundance_in_strata_wide %>% 
    pivot_longer(cols = starts_with(as.character(yr_ebird)),
                 values_to = "rel_abundance",
                 names_to = "week") %>% 
    group_by(week) %>% 
    mutate(prop_population = rel_abundance/sum(rel_abundance,na.rm = TRUE)) %>% 
    ungroup() %>% 
    group_by(strata_name) %>% 
    summarise(mean_rel_abundance = mean(rel_abundance),
              median_rel_abundance = median(rel_abundance),
              min_rel_abundance = min(rel_abundance),
              max_rel_abundance = max(rel_abundance),
              sd_rel_abundance = sd(rel_abundance),
              mean_prop_population = mean(prop_population),
              median_prop_population = median(prop_population),
              min_prop_population = min(prop_population),
              max_prop_population = max(prop_population),
              sd_prop_population = sd(prop_population)) %>% 
    mutate(cv_rel_abundance = sd_rel_abundance/median_rel_abundance,
           cv_prop_population = sd_prop_population/mean_prop_population)
  
  
 
  
  df_full <- df_full %>% 
    mutate(week = as.integer(cut(doy,breaks = doy_weeks_keep))) |> 
    filter(!is.na(week)) |> 
    mutate(week = weeks_keep[week])
  
  # 


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
                  
                  zero_betas = rep(0,max(df_full$stratum)),
                  
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
  if(!j %in% c("off_set")){
    stan_data[[j]] <- as.integer(stan_data[[j]])
  }
}


do_coverage <- TRUE

if(do_coverage){
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


}# end coverage loop


save(list = c("stan_data",
              "df_full",
              "meta_strata",
              "meta_years","strata_used",
              "df_owl_final",
              "df_bbs_final",
              "ann_coverage",
              "cumulative_coverage",
              "sp_coverage",
              "prop_in_strata_used",
              "prop_in_strata",
              "prop_in_bcr_prov",
              "abundance_in_strata_used_summary",
              "abundance_in_strata_summary",
              "abundance_in_bcr_prov_summary"),
     file = paste0("data/pre_fit_data_",sp_ebird,".RData"))

print(sp_ebird)

}# end data prep species loop
# Model fit ---------------------------------------------------------------


#re_fit <- FALSE
for(sp in owl_species$CommonName){

# for(sp in fail$species){  

  #sp_id <- unique(sp_obs_owl$species_id)
  sp_ebird <- ebirdst::get_species(sp)
 
model <- cmdstanr::cmdstan_model("models/first_difference_spatial_owls_integrated_no_scale.stan")

if(!file.exists(paste0("data/pre_fit_data_",sp_ebird,".RData"))){
  next
}
load(paste0("data/pre_fit_data_",sp_ebird,".RData"))

# if(any(is.na(stan_data$log_mean_rel_abund))){
#   stan_data$log_mean_rel_abund[which(is.na(stan_data$log_mean_rel_abund))] <-
#     mean(stan_data$log_mean_rel_abund,na.rm = TRUE)
# }
re_fit <- TRUE
if(re_fit){

  # if(sp == "Great Horned Owl"){
    ni <- 6000
    th <- ni/1000
  # }else{
  #   ni <- 1000
  #   th <- 1
  # }
fit2 <- model$sample(data = stan_data,
                    parallel_chains = 4,
                    refresh = 500,
                    iter_warmup = 3000,
                    iter_sampling = ni,
                    thin = th,
                    adapt_delta = 0.95,
                    max_treedepth = 11,
                    show_exceptions = FALSE)

fit2$save_object(paste0(output_dir,"fit_","new_scale","_",sp_ebird,".rds"))

summ2 <- fit2$summary()

saveRDS(summ2, paste0(output_dir,"fit_summary_","new_scale","_",sp_ebird,".rds"))

}else{
  fit2 <- readRDS(paste0(output_dir,"fit_","new_scale","_",sp_ebird,".rds"))
  summ2 <- readRDS(paste0(output_dir,"fit_summary_","new_scale","_",sp_ebird,".rds"))
}


}


