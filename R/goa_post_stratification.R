##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  Recalculation of total abundance and biomass of major FMP Gulf of Alaska
##  species after post-stratifying historical stations into new 2025 strata 
##
##  Use the "survey" package to calculate post-stratification weights and
##  post-stratified total abundances and biomass at both regional and Western/
##  Central/Eastern GOA management levels.
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Install/Import libraries
library(gapindex) #devtools::install_github("afsc-gap-products/gapindex")
library(akgfmaps) #devtools::install_github("afsc-gap-products/akgfmaps", build_vignettes = TRUE)
library(survey) # install.packages("survey)
library(cowplot)

## Survey change from x strata to y strata. From x number of stations to y number of stations. Citation for unbiasness. 

## Connect to Oracle either using AFSC credentials
chl <- gapindex::get_connected(db = "AFSC")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  Import Data
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Import historical (DESIGN_YEAR 2024) and new (DESIGN_YEAR 2025) GOA strata
goa_strata_2025 <- 
  akgfmaps::get_base_layers(select.region = "goa", 
                            design.year = 2025, 
                            set.crs = "EPSG:4326")$survey.strata
goa_strata_2024 <- 
  akgfmaps::get_base_layers(select.region = "goa", 
                            design.year = 2024, 
                            set.crs = "EPSG:4326")$survey.strata

## Import historical stations pre-2025 and reassign them to the new 2025 strata
goa_stations_hist <- RODBC::sqlQuery(
  channel = chl,
  query = "
  select hauljoin, year, 
  (LATITUDE_DD_START + LATITUDE_DD_END)/2 AS LAT, 
  (LONGITUDE_DD_START + LONGITUDE_DD_END) / 2 AS LON
  from gap_products.akfin_haul
  join (gap_products.akfin_cruise) using (cruisejoin)
  where survey_definition_id = 47
  and year <= 2023
  "
) |> sf::st_as_sf(coords = c("LON", "LAT"),
                  crs = "EPSG:4326") |>
  sf::st_intersection(y = goa_strata_2025[, c("STRATUM")])

## Import CPUE station data
spp_codes <- #c(21720, 21740, 30060, 10110, 30420)
  c(10110, 10130, 10180, 20510, 21720, 21740, 30060, 30420, 30050, 
    30051, 30052, 30150, 30152, 10261, 10262, 10200, 310, 406, 410, 
    420, 425, 435, 440, 445, 450, 455, 460, 471, 472, 475, 477, 480, 
    483, 485, 490, 495, 21230, 10170, 10210, 10220, 10250, 10270, 10285, 
    30020, 30100, 30430, 30475, 30535, 30560, 30576)

years <- c(1990, 1993, 1996, 1999, 2001, 2003, 2005, 2007, 2009, 
           2011, 2013, 2015, 2017, 2019, 2021, 2023)

gp_data <- gapindex::get_data(survey_set = "GOA", 
                              year_set = years, 
                              spp_codes = spp_codes,
                              channel = chl)
gp_cpue <- gapindex::calc_cpue(gapdata = gp_data)

gp_data_2025 <- gapindex::get_data(survey_set = "GOA", 
                                   year_set = 2025, 
                                   spp_codes = 21720,
                                   channel = chl)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  Import the original estimates from GAP_PRODUCTS
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
original_est <- RODBC::sqlQuery(
  channel = chl,
  query = paste("
  select YEAR, SPECIES_CODE, 
  'ORIG' as EST_TYPE,
  AREA_ID,
  BIOMASS_MT, 
  BIOMASS_MT - 1.96 * sqrt(BIOMASS_VAR) as BIOMASS_LCI,
  BIOMASS_MT + 1.96 * sqrt(BIOMASS_VAR) as BIOMASS_HCI,
  POPULATION_COUNT, 
  POPULATION_COUNT - 1.96 * sqrt(POPULATION_VAR) as POPULATION_LCI,
  POPULATION_COUNT + 1.96 * sqrt(POPULATION_VAR) POPULATION_HCI
  
  FROM GAP_PRODUCTS.BIOMASS
  where species_code in ", gapindex::stitch_entries(spp_codes), "
  and area_id in (803, 804, 805, 99903)
  order by YEAR
  ")
) 

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  Naive estimator: Reassign the 2025 survey onto the historical stations 
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# naive_data <- gp_data
# naive_data$haul <- 
#   merge(x = subset(x = naive_data$haul, select = -STRATUM),
#         y = as.data.frame(x = goa_stations_hist)[, c("HAULJOIN", "STRATUM")],
#         by = "HAULJOIN")
# naive_data$survey$DESIGN_YEAR <- 
#   naive_data$survey_design$DESIGN_YEAR <- 
#   naive_data$cruise$DESIGN_YEAR <- 2025
# naive_data$strata <- gp_data_2025$strata
# naive_data$stratum_groups <- gp_data_2025$stratum_groups
# naive_data$subarea <- gp_data_2025$subarea
# 
# naive_cpue <- 
#   merge(x = subset(x = gp_cpue, select = -STRATUM),
#         y = as.data.frame(x = goa_stations_hist)[, c("HAULJOIN", "STRATUM")],
#         by = "HAULJOIN")
# naive_cpue$DESIGN_YEAR <- 2025
# 
# naive_biomass_stratum <- 
#   gapindex::calc_biomass_stratum(gapdata = naive_data, 
#                                  cpue = naive_cpue)
# naive_biomass_region <- 
#   gapindex::calc_biomass_subarea(gapdata = naive_data, 
#                                  biomass_stratum = naive_biomass_stratum) |>
#   subset(subset = AREA_ID == 99903) |> 
#   transform(
#     EST_TYPE = 'NAIVE',
#     BIOMASS_LCI = BIOMASS_MT - 1.96 * sqrt(x = BIOMASS_VAR),
#     BIOMASS_HCI = BIOMASS_MT + 1.96 * sqrt(x = BIOMASS_VAR),
#     POPULATION_LCI = POPULATION_COUNT - 1.96 * sqrt(x = POPULATION_VAR),
#     POPULATION_HCI = POPULATION_COUNT + 1.96 * sqrt(x = POPULATION_VAR)
#   ) |>
#   transform(
#     BIOMASS_LCI = ifelse(test = BIOMASS_LCI < 0, 0, BIOMASS_LCI),
#     POPULATION_LCI = ifelse(test = POPULATION_LCI < 0, 0, POPULATION_LCI)
#   ) |>
#   subset( select = c(YEAR, SPECIES_CODE, EST_TYPE,
#                      BIOMASS_MT, BIOMASS_LCI, BIOMASS_HCI,
#                      POPULATION_COUNT, POPULATION_LCI, POPULATION_HCI) )

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  Assign original and reclassified strata along with the stratum weights
##  to 
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

new_strata_areas <- data.frame(STRATUM_NEW = gp_data_2025$strata$STRATUM,
                               Freq = gp_data_2025$strata$AREA_KM2)

years <- sort(unique(gp_cpue$YEAR))
restrat_df <- data.frame()

for (ispp in spp_codes) {
  for (iyear in years) {
    cod <- subset(
      x = gp_cpue,
      subset = SPECIES_CODE == ispp & YEAR == iyear & !is.na(x = CPUE_NOKM2),
      select = c(HAULJOIN, YEAR, SPECIES_CODE, STRATUM,
                 CPUE_KGKM2, CPUE_NOKM2)
    ) 
    if (nrow(x = cod) == 0) next
    cod2 <- 
      merge(x = cod, 
            y = as.data.frame(goa_stations_hist)[, c("HAULJOIN", "STRATUM")],
            by = "HAULJOIN", suffixes = c("", "_NEW")) |>
      merge(y = gp_data$strata[, c("STRATUM", "AREA_KM2")],
            by = "STRATUM") |>
      merge(y = new_strata_areas,
            by = "STRATUM_NEW")
    
    ## Impute stations where stratum effort is 1
    singleton_new_strata <-
      which(table(cod2$STRATUM_NEW) == 1) |> names() |> as.numeric()
    singleton_old_strata <-
      which(table(cod2$STRATUM) == 1) |> names() |> as.numeric()
    
    cod2 <- 
      rbind(cod2,
            subset(x = cod2, subset = STRATUM_NEW %in% singleton_new_strata),
            subset(x = cod2, subset = STRATUM %in% singleton_old_strata))
    # cod2$Freq <- cod2$
    
    orig_design <- survey::svydesign(
      id = ~1,
      strata = ~STRATUM,
      data = cod2,
      fpc = ~AREA_KM2
    )
    
    post_design <- survey::postStratify(design = orig_design,
                                        strata = ~STRATUM_NEW,
                                        population = new_strata_areas,
                                        partial = TRUE)
    
    stratum_estimates <- merge(
      ## Calculate stratum-level mean/var weight CPUE
      x = svyby(
        formula = ~CPUE_KGKM2,          
        by = ~STRATUM_NEW,        
        design = post_design,     
        FUN = svymean             
      ) |> 
        transform(
          STRATUM = STRATUM_NEW,
          CPUE_KGKM2_MEAN = CPUE_KGKM2 ,
          CPUE_KGKM2_VAR = se^2 
        ) |>
        subset(select = c(STRATUM, CPUE_KGKM2_MEAN, CPUE_KGKM2_VAR)),
      
      ## Calculate stratum-level biomass with variances
      y = svyby(
        formula = ~CPUE_KGKM2,          
        by = ~STRATUM_NEW,        
        design = post_design,     
        FUN = svytotal            
      ) |> 
        transform(
          STRATUM = STRATUM_NEW,
          BIOMASS_MT = CPUE_KGKM2 * 0.001,
          BIOMASS_VAR = se^2 * 1e-6
        ) |>
        subset(select = c(STRATUM, BIOMASS_MT, BIOMASS_VAR)),
      by = "STRATUM") |>
      
      ## Calculate stratum-level mean/var numbers CPUE
      merge(y =  svyby(
        formula = ~CPUE_NOKM2,          
        by = ~STRATUM_NEW,        
        design = post_design,     
        FUN = svymean             
      ) |> 
        transform(
          STRATUM = STRATUM_NEW,
          CPUE_NOKM2_MEAN = CPUE_NOKM2,
          CPUE_NOKM2_VAR = se^2
        ) |>
        subset(select = c(STRATUM, CPUE_NOKM2_MEAN, CPUE_NOKM2_VAR)),
      by = "STRATUM") |>
      
      ## Calculate stratum-level abundance with variances
      merge(y = svyby(
        formula = ~CPUE_NOKM2,          # The variable you are measuring
        by = ~STRATUM_NEW,        # The variable to group by
        design = post_design,     # Your post-stratified design object
        FUN = svytotal             # The estimator function
      ) |> 
        transform(
          STRATUM = STRATUM_NEW,
          POPULATION_COUNT = CPUE_NOKM2,
          POPULATION_VAR = se^2
        ) |>
        subset(select = c(STRATUM, POPULATION_COUNT, POPULATION_VAR))) |>
      transform (SURVEY_DEFINITION_ID = 47,
                 SURVEY = "GOA",
                 YEAR = iyear,
                 SPECIES_CODE = ispp)
    
    ## Aggregate stratum-level estimates up to the Central/Western/Eastern
    ## management areas as well as GOA-wide
    gp_data_2025$survey$YEAR = iyear
    subarea_estimates <- 
      gapindex::calc_biomass_subarea(gapdata = gp_data_2025,  
                                     biomass_stratum = stratum_estimates) 
    
    ## Add to restrat_df, express variance as 95% CIs. Truncate lower CIs to 
    ## zero if they are negative. 
    restrat_df <- 
      rbind(
        restrat_df, 
        with(subarea_estimates, 
             data.frame(
               YEAR = YEAR,
               SPECIES_CODE = SPECIES_CODE ,
               EST_TYPE = "PS",
               AREA_ID = AREA_ID,
               BIOMASS_MT = BIOMASS_MT,
               BIOMASS_LCI = BIOMASS_MT - sqrt(BIOMASS_VAR) * 1.96,
               BIOMASS_HCI = BIOMASS_MT + sqrt(BIOMASS_VAR) * 1.96,
               POPULATION_COUNT = POPULATION_COUNT,
               POPULATION_LCI = POPULATION_COUNT - sqrt(POPULATION_VAR) * 1.96,
               POPULATION_HCI = POPULATION_COUNT + sqrt(POPULATION_VAR) * 1.96
             )
        ) |>
          transform(
            BIOMASS_LCI = ifelse(test = BIOMASS_LCI < 0, 
                                 yes = 0, 
                                 no = BIOMASS_LCI),
            POPULATION_LCI = ifelse(test = POPULATION_LCI < 0, 
                                    yes = 0, 
                                    no = POPULATION_LCI)
          )
      )
    
  }
  cat("Finished with", ispp, "\n")
}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  Merge all estimators together
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
all_ests <- 
  rbind(original_est,
        restrat_df) |>
  transform(REGION = dplyr::case_when(
    AREA_ID == 803 ~ "Central GOA",
    AREA_ID == 804 ~ "Eastern GOA",
    AREA_ID == 805 ~ "Western GOA",
    AREA_ID == 99903 ~ "GOA"
  )
  ) |>
  transform(YEAR = dplyr::case_when(
    EST_TYPE == "PS" ~ (YEAR + 0.25),
    .default = YEAR
  )
  ) |>
  subset(subset = !(REGION == "Eastern GOA" & YEAR %in% c(2001, 2001.25)))
all_ests$BIOMASS_LCI <- ifelse(test = all_ests$BIOMASS_LCI < 0, 
                               yes = 0,
                               no = all_ests$BIOMASS_LCI)
all_ests$POPULATION_LCI <- ifelse(test = all_ests$POPULATION_LCI < 0, 
                                  yes = 0,
                                  no = all_ests$POPULATION_LCI)

## Order regions so that it goes from W -> C -> E GOA
all_ests$REGION <- factor(x = all_ests$REGION, 
                          levels = c("GOA", "Western GOA", "Central GOA", "Eastern GOA"))

## Save output
write.csv(x = all_ests,
          file = "plots/all_ests.csv", row.names = F)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  Output plots
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for (ispp in spp_codes) {
  
  stock_name <- 
    gp_data$species$REPORT_NAME_SCIENTIFIC[gp_data$species$SPECIES_CODE == ispp]
  
  if (all_ests %>% 
      filter(REGION != 'GOA' & 
             SPECIES_CODE == ispp) %>% nrow == 0) next 
  
  biomass <- ggplot(
    data = all_ests %>% filter(REGION != 'GOA' & SPECIES_CODE == ispp), 
    aes(x = YEAR, 
        y = BIOMASS_MT / 1000, 
        colour = EST_TYPE, 
        shape = EST_TYPE)) +
    geom_line(linewidth = 0.75, linetype = "dashed") +
    geom_pointrange(aes(ymin = BIOMASS_LCI / 1000, ymax = BIOMASS_HCI  / 1000), 
                    linewidth = 0.75, 
                    size = 0.5) +
    facet_wrap(facets = ~REGION, 
               ncol = 1, 
               scales = 'free_y') +
    theme_bw(base_size = 14) + 
    labs(x = "Year", y = "Survey biomass (1000s mt)", 
         colour = "Index type", 
         shape = "Index type", 
         title = stock_name) +
    scico::scale_color_scico_d(palette = 'roma')
  
  numbers <- ggplot(
    data = all_ests %>% filter(REGION != 'GOA' & SPECIES_CODE == ispp), 
    aes(x = YEAR, y = POPULATION_COUNT / 1000000, 
        colour = EST_TYPE, 
        shape = EST_TYPE)) +
    geom_line(linewidth = 0.75, 
              linetype = "dashed") +
    geom_pointrange(aes(ymin = POPULATION_LCI / 1000000, 
                        ymax = POPULATION_HCI / 1000000), 
                    linewidth = 0.75, 
                    size = 0.5) +
    facet_wrap(~REGION, 
               ncol = 1, 
               scales = 'free_y') +
    theme_bw(base_size = 14) + 
    labs(x = "Year", y = "Survey numbers (millions)", 
         colour = "Index type", 
         shape = "Index type", 
         title = stock_name) +
    scico::scale_color_scico_d(palette = 'roma')
  
  ggsave(filename = here::here('plots', 
                               paste0(gsub(x = stock_name, 
                                           pattern = " ", 
                                           replacement = "_"), 
                                      '_subreg_num.png')),
         plot = numbers,
         width = 11, height = 7, units = "in")
  
  
  ggsave(filename = here::here('plots', paste0(gsub(x = stock_name, 
                                                    pattern = " ", 
                                                    replacement = "_"), 
                                               '_subreg_biom.png')),
         plot = biomass,
         width = 11, height = 7, units = "in")
  
  biomass_region <- ggplot(
    data = all_ests %>% filter(REGION == 'GOA' & SPECIES_CODE == ispp), 
    aes(x = YEAR, 
        y = BIOMASS_MT * 1e-3, 
        colour = EST_TYPE, shape = EST_TYPE)) +
    geom_line(linewidth = 0.75, linetype = "dashed") +
    geom_pointrange(aes(ymin = BIOMASS_LCI * 1e-3, 
                        ymax = BIOMASS_HCI * 1e-3), 
                    linewidth = 0.75, size = 0.5) +
    theme_bw(base_size = 14) + 
    labs(x = "", y = "Survey biomass (thousand mt)", 
         colour = "Index type", shape = "Index type", title = stock_name) +
    scico::scale_color_scico_d(palette = 'roma')
  
  num_region <- ggplot(
    data = all_ests %>% filter(REGION == 'GOA' & 
                                 SPECIES_CODE == ispp & 
                                 POPULATION_COUNT > 0), 
    aes(x = YEAR, 
        y = POPULATION_COUNT * 1e-6, 
        colour = EST_TYPE, shape = EST_TYPE)) +
    geom_line(linewidth = 0.75, linetype = "dashed") +
    geom_pointrange(aes(ymin = POPULATION_LCI * 1e-6, 
                        ymax = POPULATION_HCI * 1e-6), 
                    linewidth = 0.75, size = 0.5) +
    theme_bw(base_size = 14) + 
    labs(x = "Year", y = "Survey numbers (millions)", 
         colour = "Index type", 
         shape = "Index type") +
    scico::scale_color_scico_d(palette = 'roma')
  
  ggsave(filename = here::here('plots', paste0(gsub(x = stock_name, 
                                                    pattern = " ", 
                                                    replacement = "_"), 
                                               '_region.png')),
         plot =   plot_grid(biomass_region, num_region, ncol = 1),
         width = 11, height = 7, units = "in")
}
