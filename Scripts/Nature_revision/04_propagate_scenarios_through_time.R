#propagate the uncertainty and estimates of scenario ACD 
#NOTE - need to go back and correct belowground process.
library(dplyr)
library(purrr)
library(tidyverse)
library(data.table)

#___________________________
#hardcoded decision
#___________________________
#SELECT TWICE-LOGGED RECOVERY TRAJECTORY- HARD-CODED
#1 = same slope as 1L, 0.8 is 20% slower recovery, 1.2 = 20% faster recovery

twice_logged_slope_trajectory = 1

#___________________________
#Inputs 
#___________________________
#Inputs params including the makeup of starting landscapes
source('Inputs/FixedScenarioParmams.R')
all_start_landscape <- all_start_landscape %>%  
  mutate(functional_habitat = habitat) %>%  
  mutate(
    functional_habitat = case_when(
      #deals with our rule that scenarios starting once-logged were logged at t-15 and twice-logged were at t-1 (see rationale below)
      habitat == "once-logged"  ~ "once-logged_start",
      habitat == "twice-logged" ~ "twice-logged_start",
      TRUE                       ~ habitat  # keep all other values unchanged
    )
  )

#social discount rate caluclate for 2,4,6% 
#built in the CalculateSocialDiscountRates.R script 
socialDR <- read.csv("Outputs/SocialDiscountRates_2_4_6pc_DR.csv")

 
#----------------read in scenarios -------------------------------

#yield matched scenarios where 1/30th of plantation conversion happens annually - with no time delay
#temporarily read in to get composition
scenarios_rm <- readRDS("Inputs/MasterAllScenarios.rds")
scenario_composition <- rbindlist(scenarios_rm, use.names=TRUE)

#yield matched scenarios where 1/30th of plantation conversion happens annually - WITH TIME DELAY 
scenarios <- readRDS("Inputs/MasterAllScenarios_withHarvestDelays.rds")

#function for making sure we can visualise the order of scenarios correctly
scenario_row_order_visualisation_fun <- function(x){
  x %>%  
    group_by(index, production_target, habitat, original_habitat, harvest_delay) %>%  
    arrange(true_year, .by_group = TRUE) %>%
    ungroup()
}

#___________________________
#read in posterior draws
#___________________________

#Build final posterior draws across habs ####
p_l_r <- readRDS("Outputs/primary_restored_once_logged_ACD_draws.rds") 
#currently has different twice-logged trajectories - from either clearing 15 yr or 30 once-logged.
twice <- readRDS("Outputs/twice_logged_draws_diff_assumptions.rds")
plant <- readRDS("Outputs/plantation_carbon_draws.rds") %>%  rename(habitat = species, functionalhabAge = plantationAge)

#_______________________________________________________________________________
#transition rules 
#_______________________________________________________________________________
#primary - no rules 

#once logged
#a. P -> 1L  If parcel goes to from primary to once-logged, we start at 1L yr 0
#b. 1L -> 1L if parcel starts as once-logged and stays 1L, assume logging happened 15 yrs befor scnanario start
#c. 1L -> 2L if parcel starts as once logged before re-harvest, we assume 1L recovery (already explicit in scenarios) that resets to yr 0 once-logged values.

#for b. we use the function update_1L_functional_habitat to ensure that scenarios beginning the scenario as once-logged were harvest 15yrs before

#if parcel starts as twice-logged, use the 15yr after 1L curves - with different slope factors, tracking the different slope assumption rates for twice-logged forest 
#if parcel starts as primary and goes to twice-logged, use the 30yrafter 1L, with different slope factors 
#_________________________________________________________________
# modify scenarios to accommodate different harvest trajectories
#_____________________________________________________________

#if original habitat remains as habitat (ie no transition), select one delay year 
filter_original_eq_habitat <- function(df_list) {
  df_list %>%
    map(~ .x %>%
          filter(!(original_habitat == habitat & harvest_delay != "delay 0"))
    )
}


#___________________________________________________________________
#A) 1L- 1L modification. Forest that starts a scenario as once-logged was first harvested at yr t-15

# changes `functional_habitat` to "once_logged_start" when both
# `original_habitat` and `functional_habitat` are "once_logged".
# Otherwise, leaves `functional_habitat` unchanged.

# Function: update_functional_habitat()
#_______________________________________
# Takes a list of data frames and applies two operations to each:
# (1) If both `original_habitat` and `functional_habitat` are "once_logged",
#     change `functional_habitat` to "once_logged_start".
# (2) If `original_habitat` is "once_logged", remove rows where `harvest_delay > 15`.
update_1L_functional_habitat <- function(df_list) {
  df_list %>%
    map(~ .x %>%
          mutate(
            functional_habitat = if_else(
              original_habitat == "once-logged" & functional_habitat == "once-logged",
              "once-logged_start",
              functional_habitat
            )
          ) %>%
          filter(!(original_habitat == "once-logged" & harvest_delay < 15))
    )
}

#___________________________________________________________________
# B)2L - 2L modification. Forest that starts with twice-logged was re-harvest just before scenario start
# Takes a list of data frames and applies two operations to each:
# (1) If both `original_habitat` and `functional_habitat` are "twice-logged",
#     change `functional_habitat` to "twice_logged_start".
# (2) If `original_habitat` is "once_logged", remove rows where `harvest_delay < 15`.- because you have to wait 30 yrs before can reharvest once-logged forest 
update_2L_functional_habitat <- function(df_list) {
  df_list %>%
    map(~ .x %>%
          mutate(
            functional_habitat = if_else(
              original_habitat == "twice-logged" & functional_habitat == "twice-logged",
              "twice-logged_start",
              functional_habitat
            )
          ) %>%
          filter(!(original_habitat == "twice-logged" & harvest_delay < 15))
    )
}


#remove scenarios that produce zero timber 

remove_no_timber_scenarios <- function(df_list) {
  df_list %>%
    map(~ .x %>%
          filter(production_target > 0)  # keep only scenarios producing timber
    )
}

scenarios <- scenarios %>%  update_1L_functional_habitat() %>% 
  update_2L_functional_habitat() %>%  
  filter_original_eq_habitat() %>% 
  remove_no_timber_scenarios()

#convert to datatable format
scenarios <- lapply(scenarios, as.data.table)

#________________________________________________________________
#make correct modifications to each draw
#__________________________________________________________________
primary <- p_l_r %>% filter((habitat == "primary"))
restored <- p_l_r %>% filter((habitat == "restored"))

once_logged <-  p_l_r %>% filter((habitat == "once_logged")) %>%  mutate(habitat = if_else(habitat == "once_logged", "once-logged", habitat))
#for scenarios starting once-logged assume the forest was logged 15yrs before scneario data
once_logged_start <- p_l_r %>% filter((habitat == "once_logged")) %>% 
  mutate(functionalhabAge = functionalhabAge -15) %>%  
  filter(functionalhabAge >-1) %>%  
  mutate(habitat = "once-logged_start")

#2L ->2L forest that starts as twice-logged was first logged fifteen years previously then reharvested just before scenario start
twice_logged_start <-twice %>% filter(start_age == "15yrAfter1L - e.g if parcel starts scenario 2L") %>%  
  mutate(functionalhabAge = functionalhabAge -15) %>%  
  filter(functionalhabAge >-1) %>% 
  mutate(habitat = "twice-logged_start")

#1L-2L forest that goes from once-logged to twice-logged during scenarios
#not allowed to be harvest for 15yrs)

#second harvest is not allowed for first 15 yrs, then it goes to twice-logged(30 yr after once-logging version)
twice_logged <-twice %>% filter(start_age == "30yrAfter1L - e.g if primary goes to 2L during scenario") %>% 
  mutate(habitat = "twice-logged")

#____________________________________________
#combine all draws
#____________________________________________
#Note primary forest has NA for functionalhabAge
all_draws <- primary %>% rbind(restored) %>% rbind(once_logged) %>% rbind(once_logged_start) %>%  
   rbind(plant)

#SELECT TWICE-LOGGED RECOVERY TRAJECTORY- HARD-CODE
#!!!!!select which twice-logged recovery trajectory to consider using 'slope factor' - 
#1 = same slope as 1L, 0.8 is 20% slower recovery, 1.2 = 20% faster recovery
diff_2Ls <- twice_logged %>% rbind(twice_logged_start) %>% 
  filter(slope_factor == twice_logged_slope_trajectory) %>% 
  select( draw, functionalhabAge, ACD, habitat) %>% unique()

all_draws<- all_draws %>% rbind(diff_2Ls)

#turn draws into a list where each draw is it's own df
all_draws <- all_draws %>%
  rename(functional_habitat = habitat) %>%  # rename first
  group_split(draw)     

#convert to datatable
all_draws <- lapply(all_draws, as.data.table)

#____________________________________________________________________________________________________________________________________
#We are now ready to run the full workflow below
# For each scenario set i, we extract the scenario schedule and sequentially apply the full workflow:
#   1) Add carbon estimates from all posterior draws to the scenario landscape.
#   2) Calculate ACD (aboveground carbon density) at the parcel, habitat-transition, and scenario-year scales.
#   3) Compute carbon stock-years for the scenario landscape.
#   4) Add starting landscape carbon trajectories and calculate annual ACD changes for both scenario and starting landscapes.
#   5) Convert annual ACD changes into CO2-equivalent fluxes and apply Social Cost of Carbon (SCC) discounting for multiple discount rates (2%, 4%, 6%).
#   6) Summarize posterior draws to calculate mean and 95% / 80% credible intervals of carbon fluxes.
# At the end of each iteration, the processed results for scenario i are stored in the list `posterior_summary_all[[i]]`,
# containing SCC-discounted carbon impacts with uncertainty and cumulative carbon stock-years.
#___________________________________________________________________________________________________________________________________

#to test for one scenario set, uncomment below
#i = 1

for(i in seq_along(scenarios)) {
  
  cat("Running scenario", i, "of", length(scenarios), "\n")
  
  # pick the i-th scenario
  scenario_i <- scenarios[[i]]
  
  # pick just one scenario-set for now 
  # scenario_i <- scenarios[[1]]
  # but use all draws of the carbon estimates
  draw <- all_draws
  
  #__________________________________________
  # combine scenarios with carbon values
  #__________________________________________
  
  # function to add carbon estimates to scenario schedule
  add_carbon_fun <- function(x, draw_dt) {
    # Ensure both are data.tables
    setDT(x)
    setDT(draw_dt)
    
    # Separate primary and non-primary to handle joins cleanly
    x_primary <- x[functional_habitat == "primary"]
    x_other   <- x[functional_habitat != "primary"]
    
    # Join scenario primary carbon only - only by functional habitat as functionalhabAge is NA for primary 
    joined_primary <- x_primary[
      draw_dt[functional_habitat == "primary"],
      on = .(functional_habitat),
      nomatch = 0
    ]
    
    # Remove additional column created by join
    joined_primary[, i.functionalhabAge := NULL]
    
    # Join non-primary habitats to carbon values 
    joined_other <- x_other[
      draw_dt[functional_habitat != "primary"],
      on = .(functional_habitat, functionalhabAge),
      nomatch = 0
    ]
    
    # Combine primary and non-primary to get carbon results of each row
    result <- rbindlist(list(joined_primary, joined_other), use.names = TRUE)
    
    return(result)
  }
  
  # Now apply the add carbon function for each posterior draw of the datatable list:
  # Apply add_carbon_fun() to each draw
  scenario_i <- lapply(draw, \(d) add_carbon_fun(scenario_i, d))
  #check <- scenario_i[[1]] %>% scenario_row_order_visualisation_fun()
  
  #________________________________________________________________
  # Calculate ACD per scenario year 
  #________________________________________________________________
  # next - apply the next step in the carbon calculation for the scenario
  
  # --- Function to calculate number of unique delays per transition ---
  calculate_delays_per_transition_fun <- function(x) {
    setDT(x)  # ensure it's a data.table
    
    # Count unique harvest_delay values per transition
    transitions <- x[, .(num_transition_delays = uniqueN(harvest_delay)), 
                     by = .(index, production_target, original_habitat, habitat)]
    
    # Join back to original table
    x[transitions, on = .(index, production_target, original_habitat, habitat), 
      num_transition_delays := i.num_transition_delays]
    
    return(x)
  }
  
  # --- Apply function to each list element ---
  scenario_i <- lapply(scenario_i, calculate_delays_per_transition_fun)
  #check <- scenario_i[[1]] %>% scenario_row_order_visualisation_fun()
  
  #_______________________________________________
  # Function: scenario_ACD_fun
  # Purpose:  Calculate aboveground carbon density (ACD) at multiple scales:
  #           (1) Parcel-level (10 km² units)
  #           (2) Habitat-transition by year
  #           (3) Scenario-level by year
  # Inputs:
  #    x — data.table containing columns:
  #         index, production_target, original_habitat, habitat,
  #         true_year, ACD, num_parcels, num_transition_delays,
  #         scenarioName, scenarioStart
  # Outputs:
  #    data.table with total scenario ACD per (index, production_target, true_year)
  #_______________________________________________
  
  # scenario_test_lowP <- check %>% filter(index == "all_primary_CY_D.csv 10")
  #scenario_test_highP <- check %>% filter(index == "unq2_1055")
  #x <- scenario_test_lowP
  
  scenario_ACD_fun <- function(x) {
    
    # --- Step 0: Ensure input is a data.table ---
    setDT(x)
    
    # --- Step 1: Calculate ACD values at 10 km² parcel scale --------------------
    x[, ACD_10km := ACD * 1000]  # Convert ACD (Mg/ha) to Mg per 10 km²
    x[, ACD_10km2_stag := ACD_10km * (num_parcels / num_transition_delays)]
    
    # --- Step 2: Aggregate to habitat-transition level per year -----------------
    habitat_year_summary <- x[, .(
      hab_ACD_year = sum(ACD_10km2_stag, na.rm = TRUE),
      scenarioName = first(scenarioName),
      scenarioStart = first(scenarioStart)
    ), by = .(index, production_target, original_habitat, habitat, true_year)]
    
    # --- Step 3: Aggregate to full-scenario level per year ----------------------
    scenario_year_summary <- habitat_year_summary[, .(
      scen_ACD_year = sum(hab_ACD_year, na.rm = TRUE),
      scenarioName = first(scenarioName),
      scenarioStart = first(scenarioStart)
    ), by = .(index, production_target, true_year)]
    
    # --- Step 4: Return scenario-level summary ---------------------------------
    return(scenario_year_summary)
  }
  
  # Apply scenario_ACD_fun to all scenarios
  scenario_i <- lapply(scenario_i, scenario_ACD_fun)
  #x <- scenario_i[[1]]
  
  # CALCULATE CARBON STOCK YEARS ####
  carbon_stock_years_fun <- function(x) {
    setDT(x)
    carbon_summary <- x[, .(
      cumaltive_stock_year = sum(scen_ACD_year, na.rm = TRUE),
      scenarioName = first(scenarioName),
      scenarioStart = first(scenarioStart)
    ), by = .(production_target, index)]
    return(carbon_summary)
  }
  
  carbon_yrs <- lapply(scenario_i, carbon_stock_years_fun)
  carbon_yrs <- rbindlist(carbon_yrs, use.names = TRUE, fill = TRUE, idcol = "draw")
  
  # get mean scenario stock years with uncertainty:
  carbon_stock_years_summary <- carbon_yrs[, .(
    mean_cum_stock_year = mean(cumaltive_stock_year, na.rm = TRUE),
    lwr_cum_stock_year_95 = quantile(cumaltive_stock_year, probs = 0.025, na.rm = TRUE),
    upr_cum_stock_year_95 = quantile(cumaltive_stock_year, probs = 0.975, na.rm = TRUE),
    lwr_cum_stock_year_80 = quantile(cumaltive_stock_year, probs = 0.10, na.rm = TRUE),
    upr_cum_stock_year_80 = quantile(cumaltive_stock_year, probs = 0.90, na.rm = TRUE),
    scenarioName = first(scenarioName),
    scenarioStart = first(scenarioStart)
  ), by = .(index, production_target)]
  
  #__________________________________________
  # Get starting landscape carbon trajectories 
  #__________________________________________
  
  all_start_landscape <- all_start_landscape %>%
    uncount(weights = 61) %>%
    group_by(across(everything())) %>%
    mutate(functionalhabAge = 0:60) %>%
    ungroup()
  setDT(all_start_landscape)
  
  all_SL_carbon <- lapply(all_draws, function(draw_dt) {
    setDT(draw_dt)
    setDT(all_start_landscape)
    primary_dt     <- draw_dt[functional_habitat == "primary"]
    non_primary_dt <- draw_dt[functional_habitat != "primary"]
    primary_joined <- all_start_landscape[primary_dt, on = .(functional_habitat), nomatch = 0L]
    non_primary_joined <- all_start_landscape[non_primary_dt, on = .(functional_habitat, functionalhabAge), nomatch = 0L]
    merged_dt <- rbindlist(list(primary_joined, non_primary_joined), use.names = TRUE, fill = TRUE)
    if ("i.functionalhabAge" %in% names(merged_dt)) merged_dt[, `i.functionalhabAge` := NULL]
    return(merged_dt)
  })
  
  SL_ACD_by_age <- map(all_SL_carbon, ~ .x %>%
                         mutate(SL_ACD_year = ACD * 1000 * num_parcels) %>%
                         group_by(functionalhabAge, draw, scenarioStart) %>%
                         summarise(SL_ACD_year = sum(SL_ACD_year, na.rm = TRUE), .groups = "drop") %>%
                         rename(true_year = functionalhabAge)
  )
  
  SL_stockyears_list <- map(all_SL_carbon, ~ .x %>%
                              mutate(yearACD = ACD * 1000 * num_parcels) %>%
                              group_by(scenarioStart) %>%
                              summarise(cumulative_stock_yearSL = sum(yearACD, na.rm = TRUE), .groups = "drop")
  )
  
  SL_stockyears <- bind_rows(SL_stockyears_list, .id = "draw") %>%
    mutate(draw = as.integer(draw))
  setDT(SL_stockyears)
  
  SL_carbon_stock_years_summary <- SL_stockyears[, .(
    mean_cum_stock_year = mean(cumulative_stock_yearSL, na.rm = TRUE),
    lwr_cum_stock_year  = quantile(cumulative_stock_yearSL, probs = 0.025, na.rm = TRUE),
    upr_cum_stock_year  = quantile(cumulative_stock_yearSL, probs = 0.975, na.rm = TRUE),
    scenarioStart = first(scenarioStart)
  ), by = .(scenarioStart )]
  
  #__________________________________________________
  # Calculate fluxes, ACD changes, social cost, summarize, and store for scenario i
  #__________________________________________________
  
  scenario_i <- join_scenario_with_start_ACD(scenario_list = scenario_i, start_ACD_list = SL_ACD_by_age)
  scenario_i <- lapply(scenario_i, ACD_change_function_dt)
  scenario_i <- lapply(scenario_i, flux_conversion_function)
  
  final_carbon2 <- lapply(scenario_i, final_carbon_fun, SDR = socialDR_2) %>% rbindlist(idcol = "draw")
  final_carbon4 <- lapply(scenario_i, final_carbon_fun, SDR = socialDR_4) %>% rbindlist(idcol = "draw")
  final_carbon6 <- lapply(scenario_i, final_carbon_fun, SDR = socialDR_6) %>% rbindlist(idcol = "draw")
  
  posterior_summary2 <- summarize_totcarbon_draws(final_carbon2)
  posterior_summary4 <- summarize_totcarbon_draws(final_carbon4)
  posterior_summary6 <- summarize_totcarbon_draws(final_carbon6)
  
  posterior_summary_comb <- bind_rows(
    posterior_summary2 %>% mutate(discount_rate = "2%"),
    posterior_summary4 %>% mutate(discount_rate = "4%"),
    posterior_summary6 %>% mutate(discount_rate = "6%")
  )
  
  posterior_summary_comb <- posterior_summary_comb %>%  left_join(carbon_stock_years_summary)
  posterior_summary_all[[i]] <- posterior_summary_comb
  
  rm(scenario_i, final_carbon2, final_carbon4, final_carbon6)
  gc()
}
